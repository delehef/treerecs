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

#ifndef PHYLASOLVER_RECONCILIABLEROOTEDTREEGENERATOR_H
#define PHYLASOLVER_RECONCILIABLEROOTEDTREEGENERATOR_H

// Treerecs algorithms
#include <treerecs/algorithms/Polytomysolver.h>

// Treerecs containers
#include <treerecs/containers/ReconcilableRootedTree.h>

namespace treerecs {

/*!
 * @class ReconcilableRootedTreeGenerator
 * @brief Generate ReconcilableRootedTrees from a gene tree and a species tree with Polytomysolver.
 * @details A ReconcilableRootedTreeGenerator creates ReconcilableRootedTrees which contains infos to
 *          reconcile/ reconstruct gene trees (expected the species/guide tree).
 */
class ReconcilableRootedTreeGenerator {
private:

  /// gene duplication cost.
  double duplication_cost_;

  /// gene loss cost.
  double loss_cost_;

  /// polytomies already solved.
  std::map<std::pair<Node, Node>, CostTable > polytomies_solved;

  /// algorithm to estimate minimal costs to generate a tree with ProfileNJ
  Polytomysolver solver_;

protected:
  CostTable applyPolytomysolver(const bpp::PhyloTree& genetree,
                           const Node &polytomy_root, const SpeciesGeneMap &genemap,
                           const bpp::PhyloTree &speciestree, const bool verbose = false) const;


  /// Generate a ReconcilableRootedTree. If the optional option new_root is specified, the tree is rerooted.
  ReconcilableRootedTree generate(const bpp::PhyloTree& genetree /// Genetree to reconcile.
      , const bpp::PhyloTree& speciestree /// Speciestree used a guide tree. Must be fully resolved (no multifurcations).
      , const SpeciesGeneMap& speciesGeneMap /// Species to gene map.
      , const Node& new_root = nullptr /// New root of the tree (the node must be created and connected in the genetree).
      , const bool polytomysolver_on_bifurcations = false /// add artificial genes to make a full-reconciliation.
      , const bool printProgression = false /// Print a progression bar.
      , const bool verbose = false /// Print operations.
  );

public:
/****************
 * Constructors
 */
  ReconcilableRootedTreeGenerator(const double duplication_cost = DEFAULT_DUPLICATION_COST, const double loss_cost = DEFAULT_LOSS_COST):
      duplication_cost_(duplication_cost), loss_cost_(loss_cost)
  {
  };


/****************
 * Getters
 */
  /// Returns a vector of ReconcilableRootedTree. Each ReconcilableRootedTree contains elements to solve a the tree.
  /// The returned vector has a size equal to one if there is no rerooting.
  std::vector<ReconcilableRootedTree> operator()(
        const bpp::PhyloTree& original_genetree /// Gene tree to solve.
      , const bpp::PhyloTree& speciestree /// Species tree used as guidetree for reconciliation.
      , const SpeciesGeneMap& map /// Map which indicates relations between genes and species.
      , const bool reroot  = false /// Indicates if the tree needs to be rerooted or not (default to false).
      , const bool polytomysolver_on_bifurcations = false /// Add artificial genes in bifurcations to make a full-reconciliation.
      , const bool printProgression = true /// Print algorithm progression.
      , const bool verbose = false /// Print in standard output operations.
  );

};

} // namespace treerecs

#endif //PHYLASOLVER_RECONCILIABLEROOTEDTREEGENERATOR_H
