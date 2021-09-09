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

#ifndef PHYLASOLVER_RECONCILIABLEROOTEDTREE_H
#define PHYLASOLVER_RECONCILIABLEROOTEDTREE_H

#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <unordered_map>
#include <treerecs/containers/GeneMap.h>
#include <treerecs/algorithms/Polytomysolver.h>

namespace treerecs {

/*!
 * @class ReconcilableRootedTree
 * @brief Class which contains all elements to reconcile a given tree, excepted the species tree.
 */
class ReconcilableRootedTree {
  friend class ReconcilableRootedTreeGenerator;

private:
  /// Rooted multifurcating gene tree.
  std::shared_ptr<bpp::PhyloTree> rooted_multifurcating_genetree_;

  /// Rooted species tree as guide tree.
  std::shared_ptr<bpp::PhyloTree> guidetree_;

  /// Species gene map which describe relation between a gene tree and species tree (each species has genes).
  SpeciesGeneMap speciesGeneMap_;

  /// Duplications counts in bifurcations.
  std::size_t bifurcation_duplication_number_;

  // Losses counts in bifurcations.
  std::size_t bifurcation_loss_number_;

  /// Cost of one duplication.
  double duplication_cost_;

  /// Cost of one loss.
  double loss_cost_;

  /// CostMatrix associated with polytomy roots.
  std::unordered_map<Node, CostTable> costMatrices_;

  /// Sub-tree roots occurences
  std::unordered_map<Node, std::size_t> subtree_root_occurrences_;

  /// Re-rooting
  bool rerooted_;

  /// Execution time of the solution
  double_t execution_time_;

protected:

  /// Propagate subtree root occurrence in father nodes.
  void setSubTreeRootOccurrence_recursive_propagation(const Node& node, const std::size_t& n);

public:
  ReconcilableRootedTree(
        const std::shared_ptr<bpp::PhyloTree>& genetree
      , const std::shared_ptr<bpp::PhyloTree>& guidetree
      , const SpeciesGeneMap& speciesGeneMap
      , const std::size_t& bifurcation_duplication_number
      , const std::size_t& bifurcation_loss_number
      , const double& duplication_cost
      , const double& loss_cost
      , const std::unordered_map<Node, CostTable>& costMatrices
      , const bool rerooted
      , const double& execution_time
  );

  ~ReconcilableRootedTree() = default;

  /// Returns the total cost of a tree according to the cost associated with an event (dup and loss).
  double cost() const;

  /// Check if a CostTable exists/ or has been computed for a given node.
  bool has_costtable_at(const Node &node) const;

  /// Get a CostTable associated with a given Node (std::shared_ptr<bpp::PhyloNode>).
  const CostTable& get_costtable_at(const Node &node) const;

  /// Get the time used to compute the current object.
  double execution_time() const;

  /// Get the genetree to solve.
  std::shared_ptr<const bpp::PhyloTree> genetree() const;

  /// Get the guide tree used to solve the genetree.
  std::shared_ptr<const bpp::PhyloTree> guidetree() const;

  /// Get speciesGenemap to link genes to their species.
  const SpeciesGeneMap& map() const;

  /// Get number of polytomies
  std::size_t n_multifurcations() const;

  /// Set subtree root occurence when this one in a bifurcation
  /// and propagate this count to fathers (while they are not multifurcations).
  void setSubTreeRootOccurence(const Node& node, const std::size_t& n);

  /// Get occurence of the tree.
  std::size_t occurrence(void) const;

  /// Get bifurcation_duplication_number.
  std::size_t bifurcation_duplication_number() const;

  /// Get bifurcation_loss_number.
  std::size_t bifurcation_loss_number() const;

  /// Get duplication_cost used.
  double duplication_cost() const;

  /// Get loss_cost used.
  double loss_cost() const;

  /// Check if the reconciliation is within a rerooting
  bool rerooted() const;

  void print() const;
};

} // namespace treerecs

#endif //PHYLASOLVER_RECONCILIABLEROOTEDTREE_H
