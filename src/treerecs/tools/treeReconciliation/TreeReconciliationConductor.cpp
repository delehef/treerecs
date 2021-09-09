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

#include "TreeReconciliationConductor.h"

#include <treerecs/tools/ALE/ALEevaluation.h>
#include <treerecs/tools/Statistics.h>

namespace treerecs {

std::vector<ReconciledRootedTree> TreeReconciliationConductor::operator()(
    const bpp::PhyloTree &original_genetree,
    const bpp::PhyloTree &speciestree,
    const SpeciesGeneMap &original_map,
    const LibpllAlignmentInfo *alignmentInfo,
    bool estimate_costs,
    const double duplication_cost,
    const double loss_cost,
    const double supportThreshold,
    const bool strict_thresholds,
    const bool reroot,
    const size_t sample_size,
    const bool compute_ale_loglk,
    const bool keep_max_ale_loglk,
    const bool add_loss_events,
    const bool computeDistances,
    const bool printProgression,
    const bool verbose) {
  if (estimate_costs)
    return costEstimations(original_genetree,
                           speciestree,
                           original_map,
                           alignmentInfo,
                           duplication_cost,
                           loss_cost,
                           supportThreshold,
                           strict_thresholds,
                           reroot,
                           sample_size,
                           compute_ale_loglk,
                           keep_max_ale_loglk,
                           add_loss_events,
                           computeDistances,
                           printProgression,
                           verbose);
  else
    return reconcile(original_genetree,
                     speciestree,
                     original_map,
                     alignmentInfo,
                     duplication_cost,
                     loss_cost,
                     supportThreshold,
                     strict_thresholds,
                     reroot,
                     sample_size,
                     compute_ale_loglk,
                     keep_max_ale_loglk,
                     add_loss_events,
                     computeDistances,
                     printProgression,
                     verbose);
}

std::vector<ReconciledRootedTree> TreeReconciliationConductor::costEstimations(
    const bpp::PhyloTree &original_genetree, const bpp::PhyloTree &speciestree
    , const SpeciesGeneMap &original_map
    , const LibpllAlignmentInfo *alignmentInfo
    , const double duplication_cost
    , const double loss_cost, const double supportThreshold
    , const bool contraction_inferior_than_threshold_only, const bool reroot
    , const std::size_t sample_size, const bool compute_ale_loglk
    , const bool keep_max_ale_loglk, const bool add_loss_events
    , const bool computeDistances, const bool printProgression
    , const bool verbose
) {

  double current_dupcost = duplication_cost;
  double current_losscost = loss_cost;
  double current_score = std::numeric_limits<double>::infinity();
  unsigned long while_max_guard = 20;
  unsigned long while_iteration = 0;

  std::list<
      std::tuple<double, double, std::vector<ReconciledRootedTree>>
  > history;

  Timer<unsigned long> progression(while_max_guard);

  if(printProgression) {
    std::string message =
        "Optimizations (D=" + std::to_string(current_dupcost);
    message += ", L=" + std::to_string(current_losscost) + ")";
    utils::progressionBar(std::cout, progression, false, message);
  }

  while(while_iteration < while_max_guard) {

    auto reconciledTrees = reconcile(
      original_genetree, speciestree, original_map, alignmentInfo
      , current_dupcost, current_losscost, supportThreshold
      , contraction_inferior_than_threshold_only, reroot , sample_size
      , compute_ale_loglk, keep_max_ale_loglk, add_loss_events
      , computeDistances, false, verbose

    );

    auto solution = reconciledTrees.front();
    auto this_score = solution.cost();

    auto duplication_number = solution.ndup();
    auto loss_number = solution.nloss();

    double total_events_number = duplication_number + loss_number;

    history.push_back(
        std::make_tuple(current_dupcost, current_losscost, reconciledTrees));

    if(utils::double_equivalence(this_score, current_score)) {
      break;
    } else {
      current_score = this_score;
      current_dupcost = loss_number/total_events_number;
      current_losscost = duplication_number/total_events_number;
    }

    progression.next();
    if(printProgression) {
      std::string message =
          "Optimizations (D=" + std::to_string(current_dupcost);
      message += ", L=" + std::to_string(current_losscost) + ")";
      utils::progressionBar(std::cout, progression, false, message);
    }

    while_iteration++;
  }

  if(printProgression)
    RefreshablePrinter::clean(std::cout);

  estimated_costs_ =
      std::make_unique<std::pair<double, double>>(std::get<0>(history.back()),
                                                  std::get<1>(history.back()));

  return std::get<2>(history.back());
}

std::vector<ReconciledRootedTree> TreeReconciliationConductor::reconcile(
    const bpp::PhyloTree& original_genetree,
    const bpp::PhyloTree& speciestree,
    const SpeciesGeneMap& original_map,
    const LibpllAlignmentInfo *alignmentInfo,
    const double duplication_cost,
    const double loss_cost,
    const double supportThreshold,
    const bool contraction_inferior_than_threshold_only,
    const bool reroot,
    const std::size_t sample_size,
    const bool compute_ale_loglk,
    const bool keep_max_ale_loglk,
    const bool add_loss_events,
    const bool computeDistances,
    const bool printProgression,
    const bool verbose) const {
  // Create generators
  ReconcilableRootedTreeGenerator treegen(duplication_cost, loss_cost);
  ReconciledRootedTreeGenerator recgen;

  // Compute distance matrix
  DistanceMatrix dmatrix = PhyloTreeToolBox::distanceMatrix(
      original_genetree, true, printProgression);

  if(verbose)
    std::cout << "Distance matrix:" << std::endl << dmatrix << std::endl;

  // Clone genetree
  std::shared_ptr<bpp::PhyloTree> genetree =
      PhyloTreeToolBox::cloneTree(original_genetree);

  // Clone speciesGeneMap
  SpeciesGeneMap map(original_map);
  SpeciesGeneMapper::updateInnerNodeMapping(map, *genetree, speciestree);

  // Delete monofurcations in genetree if exists (internal nodes with only one
  // son)
  // Delete in the map too.
  PhyloTreeToolBox::applyInPOT(
      *genetree,
      [&genetree, &map](const Node &node) {
        if (not genetree->isLeaf(node)) {
          if (genetree->getSons(node).size() == 1) {
            if(genetree->getRoot() == node) {
              genetree->rootAt(
                  genetree->getSons(node).at(0));
            }
            PhyloTreeToolBox::removeNode(*genetree, node);
            map.deleteGene(node);
          }
        }
      });

  // Contract tree and update structures
  if (verbose) std::cout << "> Gene tree before contraction: " << *genetree << std::endl;

  auto deletedNodes = PhyloTreeToolBox::mergeWeakBranches(
      *genetree
      , supportThreshold
      , contraction_inferior_than_threshold_only
  );

  for(auto& deletedNode: deletedNodes) {
    map.deleteGene(deletedNode);
  }

  if (verbose) std::cout << "> Contracted tree (support threshold = " << supportThreshold << "): " << *genetree << std::endl;
  assert(original_genetree.getAllLeaves().size() == genetree->getAllLeaves().size());

  // The reconciliation and reconstruction is made in two steps:
  // * first, generate reconcilablesRootedTrees: this step creates matrices and scores to solve the tree (rerooted or not) ;
  // * then, solve each tree by generating reconciledRootedTrees.
  std::vector<ReconcilableRootedTree> reconcilablesTrees = treegen(
      *genetree // Gene tree to solve.
      , speciestree // Species tree used as guidetree for reconciliation.
      , map // Map which indicates relations between genes and species.
      , reroot // Indicates if the tree needs to be rerooted or not (default to false).
      , add_loss_events
      , printProgression // Print progression bar in standard output stream.
      , verbose // Print in standard output stream operations.
  );

  assert(reconcilablesTrees.size() > 0);

  // Then reconcile/ reconstruct gene trees according to data in reconcilable gene trees
  std::vector<ReconciledRootedTree> reconciledTrees = recgen(
      reconcilablesTrees // Reconcilable trees to solve
      , dmatrix // Distance matrix of the gene tree
      , sample_size // Sample size, default to one: only one tree is going to be generated.
      , add_loss_events // add artificial genes in gene trees.
      , computeDistances // compute new branch lengths using Neighbor-Joining.
      , printProgression // print progression bar in standard output stream.
      , verbose // verbose mode.
  );

  if (alignmentInfo) {
    evaluate_with_libpll_(reconciledTrees, *alignmentInfo);
  }

  if(compute_ale_loglk or keep_max_ale_loglk) {
    auto ale_loglks = evaluate_with_ale_(reconciledTrees, speciestree,
                                         printProgression);

    if(verbose) {
      std::cout << "Resulting ALE log likelihoods : ";
      utils::write(std::cout, ale_loglks.begin(), ale_loglks.end(), "", ", ", ".");
      std::cout << std::endl;
    }

    if(verbose and ale_loglks.size() > 1){
      std::cout << "> Distribution:" << std::endl;
      auto ale_loglks_frequencies = Statistics::frequencies(ale_loglks.begin(), ale_loglks.end(), DEFAULT_NB_FREQUENCY_CLASSES);
      Statistics::plotFrequencies(ale_loglks_frequencies, std::cout);
      std::cout << std::endl;
    }

    if(keep_max_ale_loglk) {
      //sort reconciled trees by ale log likelihood.
      std::sort(reconciledTrees.begin(), reconciledTrees.end(),
                [](const ReconciledRootedTree& rrtA, const ReconciledRootedTree& rrtB){
                  return rrtA.ale_loglikelihood() > rrtB.ale_loglikelihood();
                });

      // get maximum
      double maximum_ale_loglk = reconciledTrees.front().ale_loglikelihood();
      std::size_t i = 0;
      while(i < reconciledTrees.size() and utils::double_equivalence(reconciledTrees.at(i).ale_loglikelihood(), maximum_ale_loglk)){
        i++;
      }

      assert(i > 0 and i <= reconciledTrees.size());

      if(verbose) std::cout << "  Keep " << i << " trees upon " << reconciledTrees.size() << "." << std::endl;

      // eliminate worst solutions
      reconciledTrees.erase(reconciledTrees.begin() + i, reconciledTrees.end());
      assert(reconciledTrees.size() == i);
    }
  }

  return reconciledTrees;
}

/**
 * Evaluate the given ReconciledRootedTree with the Phylogeny Likelihood Library
 *
 * \param reconciledTrees Reconstructed gene trees.
 *        Each gene tree will be modified by adding the AE log likelihood.
 * \param info
 * \return the likelihood of each tree in order
 */
std::vector<double> TreeReconciliationConductor::evaluate_with_libpll_(
    std::vector<ReconciledRootedTree>& reconciledTrees,
    const LibpllAlignmentInfo& info) const {
  std::vector<double> res;
  // todobenoit add OpenMP
#if defined(_OPENMP)
  # pragma omp parallel for ordered
#endif
  for (std::size_t i = 0; i < reconciledTrees.size(); ++i) {
    auto& tree = reconciledTrees.at(i);
    std::shared_ptr<LibpllEvaluation> evaluation =
        LibpllEvaluation::buildFromPhylo(tree.genetree(),
                                         info);
    double libpll_loglk = evaluation->optimizeAllParameters();
    res.push_back(libpll_loglk);
    tree.set_libpll_loglikelihood(libpll_loglk);
  }

  return res;
}

std::vector<double> TreeReconciliationConductor::evaluate_with_ale_(
    std::vector<ReconciledRootedTree>& reconciledTrees,
    const bpp::PhyloTree& speciestree, const bool printProgressionBar
) const {
  std::vector<double> res;
  res.reserve(reconciledTrees.size());

  Timer<std::size_t> progression(reconciledTrees.size());
  if(printProgressionBar)
    utils::progressionBar(std::cout, progression, true, "Computing ALE logLk");

  #if defined(_OPENMP)
    # pragma omp parallel for ordered
  #endif
  for(std::size_t i = 0; i < reconciledTrees.size(); ++i){
    auto& recGenetree = reconciledTrees.at(i);
    // Get data as genetree and map. The genetree needs to have no artificial gene.
    auto& genetree = *recGenetree.genetree();
    auto genetree_copy = PhyloTreeToolBox::cloneTree(genetree);
    PhyloTreeToolBox::removeArtificialGenes(*genetree_copy);

    auto genemap = recGenetree.map();

    // Compute ALE log likelihood. Tau (transfer rate) parameter = 0.01
    double ale_loglk = ALEevaluation::evaluate(*genetree_copy, speciestree, genemap, 1, 1, 0.01, 0.00, 0.1);
    recGenetree.set_ale_loglikelihood(ale_loglk);
    res.push_back(ale_loglk);

    progression.next();

    // Print progression bar
    if(printProgressionBar)
      utils::progressionBar(std::cout, progression, true, "Computing ALE logLk");
  }
  return res;
}

} // namespace treerecs

