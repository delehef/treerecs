// Copyright (C) 2018  INRIA
//
// Treerecs is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// Treerecs is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "ReconciledRootedTreeGenerator.h"

#include <treerecs/tools/random.h>
#include <treerecs/tools/IO/IO.h>
#include <treerecs/tools/SpeciesGeneMapper.h>

namespace treerecs {

bool ReconciledRootedTreeGenerator::genetree_consistent_with_compute_distances_option(
    const bpp::PhyloTree &genetree, bool computeDistances
) const {
  auto edges = genetree.getAllEdges();
  for(auto& edge: edges){
    if(not (edge->hasLength() and computeDistances)){
      //std::cerr << "Warning: the tree has branches " << (computeDistances ? "without" : "with") << " length." << std::endl;
      //return false;
    }
  }
  return true;
}

ReconciledRootedTree ReconciledRootedTreeGenerator::operator()(
    const ReconcilableRootedTree &reconcilableRootedTree, DistanceMatrix dmatrix,
    const bool add_losses, const bool computeDistances, const bool printProgression, const bool verbose
) const {
  // Clone the genetree to reconcile.
  std::shared_ptr<bpp::PhyloTree> genetree = PhyloTreeToolBox::cloneTree(*reconcilableRootedTree.genetree(), nullptr, computeDistances);

  // Get polytomy roots to reconcile.
  std::list<Node> internal_nodes;
  internal_nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(*genetree,
                                                                                [&genetree](const Node &node) {
                                                                                  return not genetree->isLeaf(node);
                                                                                }
  );

  // Get the species tree as guide tree.
  auto guidetree = reconcilableRootedTree.guidetree();

  // Get speciesGeneMap
  auto genemap = reconcilableRootedTree.map();

  // Get duplication and loss numbers.
  std::size_t bifurcation_duplication_number = reconcilableRootedTree.bifurcation_duplication_number();
  std::size_t bifurcation_loss_number = reconcilableRootedTree.bifurcation_loss_number();

  std::size_t multifurcation_duplication_number = 0;
  std::size_t multifurcation_loss_number = 0;

  // Then, solve each subtree with a polytomy_root.
  Timer<std::size_t> progression(internal_nodes.size());
  if(printProgression)
    utils::progressionBar(std::cout, progression, true, "Generate reconciled tree");

  for(auto& internal_node: internal_nodes) {
    if(reconcilableRootedTree.has_costtable_at(internal_node)) {
      // If the internal node is a polytomy(= multifurcation)
      // Get polytomy leaves
      std::vector<Node> nodes = genetree->getSons(internal_node);

      // Reduce the map to the polytomy
      SpeciesGeneMap polytomy_map =
          SpeciesGeneMapper::reduceMapAccordingToGeneList(genemap, nodes);

      if (verbose) {
        std::cout << "...multi-furcation map: " << polytomy_map.u_map()
                  << std::endl << std::endl;
      }

      //Then create the guide_tree of species.
      std::shared_ptr<bpp::PhyloTree> polytomy_guidetree =
          PhyloTreeToolBox::cloneTree(
              *guidetree,
              PhyloTreeToolBox::getLastCommonAncestor(
                  polytomy_map.getSpecies(), *guidetree));

      assert(reconcilableRootedTree.has_costtable_at(internal_node));
      auto &cost_table = reconcilableRootedTree.get_costtable_at(internal_node);

      assert(genemap.getAssociatedSpecies(
          PhyloTreeToolBox::getLastCommonAncestor(polytomy_map.getGenes(), *genetree)) ==
             genemap.getAssociatedSpecies(internal_node)); // check if the polytomy root is correct.

      //Deduce the count vector which contains the number of genes from each speciation.
      //This vector will define the history of the genes.
      if(verbose) std::cout << "Cost matrix:" << std::endl << cost_table << std::endl;
      Polytomysolver polytomysolver_;
      std::unordered_map<Node, std::size_t> v_count =
          polytomysolver_.getCountVector(*polytomy_guidetree, cost_table,
                                         polytomy_map);
      if (verbose)
        std::cout << "...count vector: " << v_count << std::endl << std::endl;

      //Finally, apply the ProfileNJ algorithm on the gene tree
      auto events = profileNJ_(*genetree, internal_node, *polytomy_guidetree,
                               dmatrix, v_count, genemap, computeDistances,
                               verbose);

      // Update gene events numbers
      if (events.find(duplication) != events.end()) {
        multifurcation_duplication_number += events.at(duplication).size();
      }

      if (events.find(loss) != events.end()) {
        multifurcation_loss_number += events.at(loss).size();
      }
    } else { // If the internal node is a bifurcation.
      //update the distance matrix
      std::vector<Node> sons = genetree->getSons(internal_node);
      nj_.updateDistances(dmatrix, internal_node, sons.front(), sons.back());
    }

    progression.next();
    if(printProgression)
      utils::progressionBar(std::cout, progression, true, "Generate reconciled tree");
  }

  if(not add_losses) {
    PhyloTreeToolBox::removeArtificialGenes(*genetree);
  }

  double total_execution_time = progression.time_past_since_last_update()
                                + reconcilableRootedTree.execution_time();

  assert(genetree_consistent_with_compute_distances_option(*genetree, computeDistances));

  return ReconciledRootedTree(
      genetree
      , genemap
      , bifurcation_duplication_number + multifurcation_duplication_number
      , bifurcation_loss_number + multifurcation_loss_number
      , reconcilableRootedTree.duplication_cost()
      , reconcilableRootedTree.loss_cost()
      , reconcilableRootedTree.rerooted()
      , total_execution_time
  );
}

std::vector<ReconciledRootedTree> ReconciledRootedTreeGenerator::operator()(
    const std::vector<ReconcilableRootedTree> &reconcilableRootedTrees, const DistanceMatrix &dmatrix
    , const size_t sample_size, const bool add_losses, const bool computeDistances, const bool printProgression
    , const bool verbose
) const {
  std::size_t s_i = 0;
  std::vector<double> occurrences(reconcilableRootedTrees.size());
  std::generate(occurrences.begin(), occurrences.end(), [&s_i, &reconcilableRootedTrees] { return reconcilableRootedTrees.at(s_i++).occurrence(); });
  auto index_sample = WeightedRandomPick(occurrences, sample_size);


  Timer<std::size_t> progression(index_sample.size());
  if(printProgression)
    utils::progressionBar(std::cout, progression,
                          true, "Generate reconciled trees");

  std::vector<ReconciledRootedTree> sample;
  sample.reserve(sample_size);

  #if defined(_OPENMP)
  # pragma omp parallel for
  #endif
  for(std::size_t i = 0; i < index_sample.size(); ++i) {
    const ReconcilableRootedTree& reconcilableRootedTree =
        reconcilableRootedTrees.at(index_sample.at(i));
    ReconciledRootedTree solution =
        operator()(reconcilableRootedTree, dmatrix, add_losses,
                   computeDistances,
                   printProgression and index_sample.size() == 1, verbose);

    // Test if ok
    // the resulting tree has no monofurcation
    assert(not PhyloTreeToolBox::hasMonofurcations(*solution.genetree()));

    // has no polytomy obviously
    assert(PhyloTreeToolBox::findPolytomies(*solution.genetree()).size() == 0);

    #if defined(_OPENMP)
    #pragma omp critical
    #endif
    {
      sample.emplace_back(std::move(solution));
      progression.next();
    }

    if (printProgression)
      utils::progressionBar(std::cout, progression,
                            true, "Generate reconciled trees");
  }

  return sample;
}

} // namespace treerecs

