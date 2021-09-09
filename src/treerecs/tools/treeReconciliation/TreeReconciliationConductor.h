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

#ifndef PHYLASOLVER_TREERECONCILIATIONCONDUCTOR_H
#define PHYLASOLVER_TREERECONCILIATIONCONDUCTOR_H

#include <memory>
#include <stdexcept>

#include <treerecs/tools/treeReconciliation/ReconcilableRootedTreeGenerator.h>
#include <treerecs/tools/treeReconciliation/ReconciledRootedTreeGenerator.h>
#include <treerecs/tools/LibpllEvaluation.h>

namespace treerecs {

/*!
 * @class TreeReconciliationConductor
 * @brief Reconcile a gene tree according to a species tree with Polytomysolver and ProfileNJ.
 */
class TreeReconciliationConductor {
 public:
  /// \brief Reconcile a given gene tree with a species tree used as guide tree.
  /// \return A vector of reconciled/ reconstructed gene trees with some additional data (see ReconciledRootedTree).
  /// \param original_genetree Gene tree to reconcile.
  /// \param speciestree Species tree used as guide tree. Must be fully bifurcated.
  /// \param original_map Map associating genes and species.
  /// \param alignmentInfo Alignment info (null to ignore libpll likelihood)
  /// \param duplication_cost Cost of one gene duplication.
  /// \param loss_cost Cost of one gene loss.
  /// \param supportThreshold Maximal branch support value to define contracted branches.
  /// \param strict_thresholds contract branches with a support strictly lower
  ///        than the threshold.
  /// \param reroot Explore new roots for the tree. The program will chooses existing internal nodes as new possible root.
  /// \param sample_size Size of the sampling, default to one.
  /// \param compute_ale_loglk Compute ale log likelihood.
  /// \param keep_max_ale_loglk Keep gene trees maximizing ale loglikelihood.
  /// \param add_loss_events compute a full reconciliation by adding artificial genes as loss events.
  /// \param computeDistances Compute new branch lengths.
  /// \param printProgression Print in terminal algorithm progressions.
  /// \param verbose Verbose (default = false).
  std::vector<ReconciledRootedTree> operator()(
      const bpp::PhyloTree& original_genetree
      , const bpp::PhyloTree& speciestree
      , const SpeciesGeneMap& original_map
      , const LibpllAlignmentInfo *alignmentInfo
      , bool estimate_costs
      , const double duplication_cost
      , const double loss_cost
      , const double supportThreshold = 0.0
      , const bool strict_thresholds = DEFAULT_STRICT_SUPPORT_THRESHOLDS
      , const bool reroot = false
      , const std::size_t sample_size = 1
      , const bool compute_ale_loglk = false
      , const bool keep_max_ale_loglk = false
      , const bool add_loss_events = false
      , const bool computeDistances = false
      , const bool printProgression = true
      , const bool verbose = false
  );

  /// Get estimated costs
  std::pair<double, double> estimated_costs() {
    if (estimated_costs_) return *estimated_costs_;
    throw std::logic_error("Costs have not been estimated");
  };

 protected:
  std::vector<ReconciledRootedTree> reconcile(
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
      const bool verbose) const;

  /// \brief Reconcile a given gene tree with a species tree used as guide tree by optimizing gene duplication and loss costs.
  /// \return std::tuple of three elements: the resulting duplication and loss costs and a vector of reconciled/ reconstructed gene trees with some additional data (see ReconciledRootedTree).
  /// \param original_genetree Gene tree to reconcile.
  /// \param speciestree Species tree used as guide tree. Must be fully resolved.
  /// \param original_map Map associating genes and species.
  /// \param duplication_cost Cost of one gene duplication.
  /// \param loss_cost Cost of one gene loss.
  /// \param supportThreshold Maximal branch support value to define contracted branches.
  /// \param strict_thresholds contract branches with a support strictly lower
  ///        than the threshold.
  /// \param reroot Explore new roots for the tree. The program will chooses existing internal nodes as new possible root.
  /// \param sample_size Size of the sampling, default to one.
  /// \param compute_ale_loglk Compute ale log likelihood.
  /// \param keep_max_ale_loglk Keep gene trees maximizing ale loglikelihood.
  /// \param add_loss_events compute a full reconciliation by adding artificial genes as loss events.
  /// \param computeDistances Compute new branch lengths.
  /// \param printProgression Print in terminal algorithm progressions.
  /// \param verbose Verbose (default = false).
  std::vector<ReconciledRootedTree> costEstimations(
      const bpp::PhyloTree& original_genetree
      , const bpp::PhyloTree& speciestree
      , const SpeciesGeneMap& original_map
      , const LibpllAlignmentInfo *alignmentInfo
      , const double duplication_cost
      , const double loss_cost
      , const double supportThreshold = 0.0
      , const bool strict_thresholds = DEFAULT_STRICT_SUPPORT_THRESHOLDS
      , const bool reroot = false
      , const std::size_t sample_size = 1
      , const bool compute_ale_loglk = false
      , const bool keep_max_ale_loglk = false
      , const bool add_loss_events = false
      , const bool computeDistances = false
      , const bool printProgression = true
      , const bool verbose = false
  );

  /// Evaluate the given ReconciledRootedTree with the Phylogeny Likelihood
  /// Library
  std::vector<double> evaluate_with_libpll_(
      std::vector<ReconciledRootedTree>& reconciledTrees /// Reconstructed gene trees. Each gene tree will be modified by adding the AE log likelihood.
      , const LibpllAlignmentInfo& info)
  const;

  /// \brief Evaluate and modify each reconstructed gene trees with an ALE log likelihood.
  /// \return A vector of ALE loglikelihood.
  std::vector<double> evaluate_with_ale_(
      std::vector<ReconciledRootedTree>& reconciledTrees /// Reconstructed gene trees. Each gene tree will be modified by adding the AE log likelihood.
      ,
      const bpp::PhyloTree& speciestree /// Species tree used for reconciliation.
      ,
      const bool printProgressionBar = false /// Print a progression bar during computations.
  ) const;

  /// Estimated duplication (first) and loss (second) costs
  std::unique_ptr<std::pair<double, double>> estimated_costs_ = nullptr;
};

} // namespace treerecs

#endif //PHYLASOLVER_TREERECONCILIATIONCONDUCTOR_H
