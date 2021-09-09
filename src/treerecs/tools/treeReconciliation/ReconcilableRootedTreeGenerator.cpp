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
#include "ReconcilableRootedTreeGenerator.h"

// Treerecs tools
#include <treerecs/tools/PhyloTreeToolBox.h>
#include <treerecs/tools/SpeciesGeneMapper.h>

namespace treerecs {

CostTable ReconcilableRootedTreeGenerator::applyPolytomysolver(
    const bpp::PhyloTree &genetree, const Node &polytomy_root, const SpeciesGeneMap &genemap
    , const bpp::PhyloTree &speciestree, const bool verbose
) const {
  /// Solve a polytomy given by the polytomy root and returns the number of events deduced (duplications as losses).
  if(verbose) std::cout << "## Solving multi-furcation at node " << polytomy_root << "..." <<  std::endl;

  //Get nodes of the polytomy (leaves only)
  std::vector<Node> nodes = genetree.getSons(polytomy_root);

  //First, reduce the genemap to the genes of interest.
  SpeciesGeneMap polytomy_map = SpeciesGeneMapper::reduceMapAccordingToGeneList(genemap, nodes);
  if(verbose) std::cout << "...multi-furcation map: " << polytomy_map.u_map() << std::endl << std::endl;

  //Then create the guide_tree of species.
  std::shared_ptr<bpp::PhyloTree> guidetree = PhyloTreeToolBox::cloneTree(speciestree, PhyloTreeToolBox::getLastCommonAncestor(polytomy_map.getSpecies(), speciestree));//phyloTreeToolBox_.cloneTree(speciestree);
  //phyloTreeToolBox_.pruneTree_v0(*guidetree, polytomy_map, false);

  //Compute the cost matrix.
  return solver_.computeCostTable(*guidetree, polytomy_map, duplication_cost_, loss_cost_, verbose);
}

ReconcilableRootedTree ReconcilableRootedTreeGenerator::generate(
    const bpp::PhyloTree &genetree, const bpp::PhyloTree &speciestree, const SpeciesGeneMap &speciesGeneMap
    , const Node &new_root, const bool polytomysolver_on_bifurcations, const bool printProgression, const bool verbose
) {
  // Copy all informations as map, and initial gene tree and reroot.
  SpeciesGeneMap rerooted_genetree_map(speciesGeneMap);

  // Check if we are going to change the root
  bool rerooting = false;//new_root == genetree.getRoot();

  // Create the gene tree and reroot it if asked
  std::shared_ptr<bpp::PhyloTree> rooted_genetree_copy
      = PhyloTreeToolBox::cloneTree(genetree);

  Node original_root = genetree.getRoot();
  if(new_root != nullptr and new_root != rooted_genetree_copy->getRoot()) {
    rooted_genetree_copy->rootAt(new_root);
    rerooting = true;

    // Because of changes in edge direction with rerooting, inner nodes has not
    // the same sons as before. We have to update their species association
    // deduced by the common species ancestor of their son species.
    SpeciesGeneMapper::updateInnerNodeMappingAfterRerooting(
        rerooted_genetree_map, *rooted_genetree_copy, speciestree, original_root);

    // After a change of root, the ancient root can produce a monofurcation (
    // only one son for a given father node).
    if(rooted_genetree_copy->getNumberOfSons(original_root) == 1) {
      PhyloTreeToolBox::removeNode(*rooted_genetree_copy, original_root);
      rerooted_genetree_map.deleteGene(original_root);
    }
  }

  PhyloTreeToolBox::resetNodeIdInPostOrder(*rooted_genetree_copy);

  assert(PhyloTreeToolBox::allEdgesAreCorrect(*rooted_genetree_copy));

  // assertion: the original_genetree has actually no monofurcation.
  assert(not PhyloTreeToolBox::hasMonofurcations(*rooted_genetree_copy));

  // Create the guidetree: which is a pruned tree of the speciestree.
  std::shared_ptr<bpp::PhyloTree> guidetree = PhyloTreeToolBox::cloneTree(
      speciestree,
      PhyloTreeToolBox::getLastCommonAncestor(rerooted_genetree_map.getSpecies(), speciestree));

  // Then, estimate events in bifurcations
  auto roots = PhyloTreeToolBox::getFurcations(*rooted_genetree_copy);
  auto& bifurcation_roots = roots.first;
  auto& multifurcation_roots = roots.second; // if everything is ok, nodes are in post_order.
  if(verbose) {
    std::cout << "...found " << multifurcation_roots.size()
              << " polytomie(s) at: " << std::endl;
    std::cout << multifurcation_roots << std::endl << std::endl;
  }

  // Create progression
  Timer<std::size_t> progression(roots.first.size() + roots.second.size());
  if(printProgression)
    utils::progressionBar(std::cout, progression,
                          true, "Generate reconcilable tree");

  // Create a dictionnary of matrices to solve multi-furcations
  std::unordered_map<Node, CostTable> costMatrices;

  // Estimate events in bifurcations
  std::size_t bifurcation_duplication_number = 0;
  std::size_t bifurcation_loss_number = 0;

  for (auto p_root_it = bifurcation_roots.begin(); p_root_it != bifurcation_roots.end(); p_root_it++) {
      // In reconciliation mode, we are seeking and re-contructing hidden gene losses as duplications.
      // So we solve bifurcation like we do with a multifurcation.
      if(polytomysolver_on_bifurcations) {
        auto p_root_father = rooted_genetree_copy->hasFather(*p_root_it)
                             ? rooted_genetree_copy->getFather(*p_root_it)
                             : nullptr;

        auto polytomy_key = std::pair<Node, Node>(p_root_father, *p_root_it);

        #if defined(_OPENMP)
        # pragma omp critical
        #endif
        {
          bool polytomyNotAlreadySolved =
              polytomies_solved.find(polytomy_key) == polytomies_solved.end();

          if (polytomyNotAlreadySolved) {
            costMatrices[*p_root_it] =
                applyPolytomysolver(*rooted_genetree_copy,
                                    *p_root_it, rerooted_genetree_map,
                                    speciestree);
            polytomies_solved[polytomy_key] = costMatrices.at(*p_root_it);
          } else {
            costMatrices[*p_root_it] = polytomies_solved.at(polytomy_key);
          }
        }
      } else {
        // In non-reconciliation mode, we only needs to guess gene events in a bifurcation.
        std::map<Event, std::size_t> events = PhyloTreeToolBox::getGeneEventsInBifurcation(*p_root_it,
                                                                                           *rooted_genetree_copy,
                                                                                           *guidetree,
                                                                                           rerooted_genetree_map);
        bifurcation_duplication_number += events.at(duplication);
        bifurcation_loss_number += events.at(loss);

        progression.next();
        if (printProgression)
          utils::progressionBar(std::cout, progression,
                                true, "Generate reconcilable tree");
      }

  }

  // Then, create matrices which estimates possible events in multifurcations with polytomysolver
  for(auto& multifurcation_root: multifurcation_roots) {
    auto p_root_father = rooted_genetree_copy->hasFather(multifurcation_root) ? rooted_genetree_copy->getFather(
        multifurcation_root) : nullptr;

    auto polytomy_key = std::pair<Node, Node>(p_root_father, multifurcation_root);

    #if defined(_OPENMP)
    # pragma omp critical
    #endif
    {
      bool polytomyNotAlreadySolved = polytomies_solved.find(polytomy_key) == polytomies_solved.end();
      if (polytomyNotAlreadySolved) {
        costMatrices[multifurcation_root] =
            applyPolytomysolver(
                *rooted_genetree_copy, multifurcation_root, rerooted_genetree_map, speciestree
            );

        polytomies_solved[polytomy_key] = costMatrices.at(multifurcation_root);

      } else {
        costMatrices[multifurcation_root] = polytomies_solved.at(polytomy_key);
      }
    }

    progression.next();
    if (printProgression)
      utils::progressionBar(std::cout, progression,
                            true, "Generate reconcilable tree");
  }

  //Finally we can build the ReconcilableRootedTree

  double execution_time = progression.time_past_since_last_update();
  return ReconcilableRootedTree(
      rooted_genetree_copy // The genetree to solve (with polytomies).
      , guidetree // The species tree used as guide tree.
      , rerooted_genetree_map // The updated map.
      , bifurcation_duplication_number // number of duplications in bifurcations.
      , bifurcation_loss_number // number of losses in bifurcations.
      , duplication_cost_ // parameter: cost for one gene duplication.
      , loss_cost_ // parameter: cost for one gene loss.
      , costMatrices // cost matrices which gives minimal number of gene duplications and losses in polytomies.
      , rerooting // indicates if the tree has been rerooted.
      , execution_time // time past to create this instance.
  );
}

std::vector<ReconcilableRootedTree> ReconcilableRootedTreeGenerator::operator()(
    const bpp::PhyloTree& original_genetree,
    const bpp::PhyloTree& speciestree,
    const SpeciesGeneMap& map,
    const bool reroot,
    const bool polytomysolver_on_bifurcations,
    const bool printProgression,
    const bool verbose) {
  if (verbose) std::cout << "> Current gene tree to reconcile: " << original_genetree << std::endl;

  // Seek roots to reconcile.
  std::vector<Node> roots;
  if (reroot) {
    std::list<Node> internal_nodes = PhyloTreeToolBox::getInternalNodes(original_genetree);
    roots = {internal_nodes.begin(), internal_nodes.end()};
  } else {
    roots.reserve(1);
    roots.push_back(original_genetree.getRoot());
  }

  //Then create ReconcilableRootedTrees
  std::vector<ReconcilableRootedTree> reconcilableRootedTrees;
  reconcilableRootedTrees.reserve(roots.size());

  double minimal_cost;
  bool init_minimal_cost = true;

  Timer<double> progression(roots.size());
  if(printProgression and roots.size() > 1)
    utils::progressionBar(std::cout, progression,
                          true, "Generating reconcilable gene trees");

  #if defined(_OPENMP)
  #pragma omp parallel for
  #endif
  for(std::size_t i = 0; i < roots.size(); ++i) {
    auto& root = roots.at(i);

    auto reconcilableRootedTree =
        generate(
            original_genetree// Genetree to reconcile.
            , speciestree // Speciestree used a guide tree.
            , map // Species to gene map.
            , (reroot) ? root : nullptr // New root of the tree (the node must be created and connected in the genetree).
            , polytomysolver_on_bifurcations
            , printProgression and roots.size() == 1
            , verbose // Print operations.
        );

    double current_cost = reconcilableRootedTree.cost();

    // Check if the tree can be saved -> the associated cost must be equal of
    // inferior than the current minimal solution cost.
    #if defined(_OPENMP)
    #pragma omp critical
    #endif
    {
      if (init_minimal_cost) {
        // Assign a value to the minimal cost because of non default value.
        minimal_cost = current_cost;
        init_minimal_cost = false;
      }

      // If the current cost is minimal, clear the vector and update the minimal
      // cost.
      if (current_cost < minimal_cost) {
        reconcilableRootedTrees.clear(); // delete all solutions.
        minimal_cost = current_cost;
      }

      // Then add the reconcilable rooted tree.
      if (utils::double_equivalence(minimal_cost, current_cost))
        reconcilableRootedTrees.emplace_back(std::move(reconcilableRootedTree));


      progression.next();
      if (printProgression and roots.size() > 1)
        utils::progressionBar(std::cout, progression,
                              true, "Generating reconcilable gene trees");
    }
  }

  // Test if all elements have the same cost:
  for(auto& solvable: reconcilableRootedTrees){
    if(not utils::double_equivalence(minimal_cost, solvable.cost())){
      std::cerr << "Bad selection of solutions" << std::endl;
      assert(utils::double_equivalence(minimal_cost, solvable.cost()));
    }
  }

  polytomies_solved.clear();

  return reconcilableRootedTrees;
}

} // namespace treerecs
