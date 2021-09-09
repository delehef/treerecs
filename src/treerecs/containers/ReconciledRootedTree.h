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

#ifndef PHYLASOLVER_RECONCILEDROOTEDTREE_H
#define PHYLASOLVER_RECONCILEDROOTEDTREE_H


#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <treerecs/containers/GeneMap.h>

namespace treerecs {

/*!
 * @class ReconciledRootedTree
 * @brief Contains the reconciled tree with duplication and loss events.
 */
class ReconciledRootedTree {
  friend class ReconciledRootedTreeGenerator;

private:
  /// Genetree solved.
  std::shared_ptr<bpp::PhyloTree> genetree_;

  /// GeneMap.
  SpeciesGeneMap map_;

  /// Events in the genetree.
  std::size_t duplication_number_;
  std::size_t loss_number_;

  /// Cost of one duplication.
  double duplication_cost_;

  /// Cost of one loss.
  double loss_cost_;

  /// Gene tree has been rerooted
  bool rerooted_;

  /// Execution_time.
  double execution_time_;

  /// Ale log likelihood
  double ale_loglk_;

  /// libpll log likelihood
  double libpll_loglk_;

  /// Indicated if ale log likelihood has been computed
  bool ale_loglk_computed_;

  /// Indicated if libpll log likelihood has been computed
  bool libpll_loglk_computed_;


protected:

public:
/****************
 * Constructors
 */
  ReconciledRootedTree(
        const std::shared_ptr<bpp::PhyloTree>& genetree /// Reconciled gene tree, needs to be solved (without multifurcation)
      , const SpeciesGeneMap& map /// Map used to reconcile tree.
      , const std::size_t duplication_number /// Number of duplications that occurs in the gene tree.
      , const std::size_t loss_number /// Number of losses that occurs in the gene tree.
      , const double duplication_cost /// Cost of one duplication.
      , const double loss_cost /// Cost of one loss.
      , const bool rerooted /// Gene tree has been rerooted.
      , const double execution_time /// Time used to compute the solution.
  ):    genetree_(genetree)
      , map_(map)
      , duplication_number_(duplication_number)
      , loss_number_(loss_number)
      , duplication_cost_(duplication_cost)
      , loss_cost_(loss_cost)
      , rerooted_(rerooted)
      , execution_time_(execution_time)
      , ale_loglk_(0.0)
      , ale_loglk_computed_(false)
      , libpll_loglk_computed_(false)
  {}


/****************
 * Destructor
 */
  ~ReconciledRootedTree() = default;


/****************
 * Getters
 */
  /// Compute the total cost of the solution.
  double cost() const { return duplication_cost_ * (double)ndup() + loss_cost_ * (double)nloss(); }

  /// Get the number of duplication in the solution.
  std::size_t ndup() const { return duplication_number_; }

  /// Get the number of losses in the solution.
  std::size_t nloss() const { return loss_number_; }

  /// Get the reconciled rooted gene tree.
  std::shared_ptr<bpp::PhyloTree> genetree() const { return genetree_; }

  /// Map used to link species tree to gene tree.
  SpeciesGeneMap map() const { return map_; }

  /// Check if the tree has been rerooted.
  bool rerooted() const { return rerooted_; }

  /// Time used to compute the solution.
  double execution_time() const { return execution_time_; }

  /// Indicates if ALE log likelihood has been computed.
  bool evaluated_with_ale() const { return ale_loglk_computed_; }

  /// Return ALE log likelihood.
  double ale_loglikelihood() const {
    if(not evaluated_with_ale()){
      std::cerr << "Error, tree has not been computed with ALE." << std::endl;
      exit(EXIT_FAILURE);
    }
    return ale_loglk_;
  }

  /// Set ALE evaluation of the gene tree.
  void set_ale_loglikelihood(const double value) {
    ale_loglk_ = value;
    ale_loglk_computed_ = true;
  }

  /// Indicates if libpll log likelihood has been computed.
  bool evaluated_with_libpll() const { return libpll_loglk_computed_; };

  /// Return libpll log likelihood.
  double libpll_loglikelihood() const {
    if(not evaluated_with_libpll()){
      std::cerr << "Error, tree has not been computed with libpll." << std::endl;
      exit(EXIT_FAILURE);
    }
    return libpll_loglk_;
  };

  /// Set libpll evaluation of the gene tree.
  void set_libpll_loglikelihood(const double value) {
    libpll_loglk_ = value;
    libpll_loglk_computed_ = true;
  }

};

} // namespace treerecs

#endif //PHYLASOLVER_RECONCILEDROOTEDTREE_H
