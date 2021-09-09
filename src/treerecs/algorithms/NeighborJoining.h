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

#ifndef PHYLASOLVER_NEIGHBORJOINING_H
#define PHYLASOLVER_NEIGHBORJOINING_H

// Includes std
#include <vector>
#include <utility>
#include <memory>

// Includes Bpp
#include <Bpp/Phyl/Tree/PhyloTree.h>

// Includes Treerecs
#include <treerecs/containers/Table.h>

namespace treerecs {

using Node = std::shared_ptr<bpp::PhyloNode>;
using DistanceMatrix = Table<double, Node, Node>;

/**
 * \class NeighborJoining
 * \brief This class provides methods using Neighbor-Joining algorithm for tree reconciliation.
 */
class NeighborJoining {
 private:
  void updateBifurcation(bpp::PhyloTree& tree
                         , const Node& father
                         , const Node& son_left
                         , const Node& son_right
                         , DistanceMatrix& dmatrix
  ) const;

  std::pair<double, std::size_t> sum_of_distances(const DistanceMatrix& dmatrix
                                                  , const Node& node
  ) const;

 public:
  double DistanceBetweenPairAndOuterNodes(const DistanceMatrix &D
                                          , const Node &x
                                          , const Node &y
                                          , const Node &t) const;

  double DistanceFromPairAncestor(const DistanceMatrix &D
                                  , const Node &x
                                  , const Node &y) const;

  std::pair<Node, Node> findPairMinimizingQ(const DistanceMatrix& D
                                            , const std::vector<Node>& sl_genes
                                            , const std::vector<Node>& sr_genes
                                            , const bool verbose = false) const;

  void updateDistances(DistanceMatrix& D
                       , const Node& father
                       , const Node& son_left
                       , const Node& son_right) const;

  void operator()(bpp::PhyloTree& tree
                  , DistanceMatrix dmatrix) const;
};

} // namespace treerecs

#endif //PHYLASOLVER_NEIGHBORJOINING_H
