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

#ifndef PHYLASOLVER_RECONCILEROOTEDTREEGENERATOR_H
#define PHYLASOLVER_RECONCILEROOTEDTREEGENERATOR_H

#include <treerecs/containers/ReconciledRootedTree.h>
#include <treerecs/containers/ReconcilableRootedTree.h>
#include <treerecs/algorithms/ProfileNJ.h>

namespace treerecs {

/*!
 * @class ReconciledRootedTreeGenerator
 * @brief Solve ReconcilableRootedTrees with ProfileNJ by generating ReconciledRootedTrees.
 */
class ReconciledRootedTreeGenerator {
private:
  ProfileNJ profileNJ_;
  NeighborJoining nj_;
protected:
  bool genetree_consistent_with_compute_distances_option(
      const bpp::PhyloTree& genetree
      , bool computeDistances) const;

public:

  /// Solve one ReconcilableRootedTree through ProfileNJ.
  ReconciledRootedTree operator()(const ReconcilableRootedTree& reconcilableRootedTree /// Reconcilable rooted tree to solve.
      , DistanceMatrix dmatrix /// Distance matrix of the gene tree.
      , const bool add_losses = false /// Add gene losses in resulting trees.
      , const bool computeDistances = false /// Compute new branch lengths.
      , const bool printProgression = false /// Print a progression bar.
      , const bool verbose = false /// Verbose (default = false).
  ) const;

  /// Solve each given ReconcilableRootedTrees with ProfileNJ.
  std::vector<ReconciledRootedTree> operator()(
        const std::vector<ReconcilableRootedTree>& reconcilableRootedTrees /// Reconcilable trees to solve
      , const DistanceMatrix& dmatrix /// Distance matrix of the gene tree
      , const std::size_t sample_size = 1 /// Sample size, default to one: only one tree is going to be generated.
      , const bool add_losses = false /// Add gene losses in resulting trees.
      , const bool computeDistances = false /// Compute new branch lengths using Neighbor-Joining.
      , const bool printProgression = true /// Print progession.
      , const bool verbose = false /// Verbose mode.
  ) const;
};

/// Print ReconciledRootedTree infos.
inline std::ostream& operator<<(std::ostream& os,
                                const treerecs::ReconciledRootedTree& rrt) {
  os << "Reconciled tree ("
     << "total cost = " << rrt.cost()
     << ", duplications = " << rrt.ndup()
     << ", losses = " << rrt.nloss()
     << ", execution time = " << rrt.execution_time()
     << "s.)" << std::endl;
  os << *rrt.genetree();
  return os;
}

} // namespace treerecs

#endif //PHYLASOLVER_RECONCILEROOTEDTREEGENERATOR_H
