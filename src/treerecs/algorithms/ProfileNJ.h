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

#ifndef TREERECS_PROFILENJ_H
#define TREERECS_PROFILENJ_H

//Include BPP
#include <Bpp/Phyl/Tree/PhyloNode.h>

//Include containers
#include <treerecs/containers/Table.h>
#include <treerecs/containers/Grid.h>
#include <treerecs/containers/NodeProperty.h>
#include <treerecs/containers/GeneMap.h>

//Include tools
#include "NeighborJoining.h"

namespace treerecs {

/**
 * \class ProfileNJ
 * \brief This class provides methods using ProfileNJ algorithm for tree reconciliation and reconstruction.
 */
class ProfileNJ {
 private:
  /// Neighbor joining algorithm tools.
  NeighborJoining nj_;

 protected:
  void updateDistances(DistanceMatrix &D, const Node &father, const Node &son_left, const Node &son_right) const;

  Node updateTree_GeneMap_DistanceMatrix(
      bpp::PhyloTree &tree
      , const Node &species
      , const Node &species_left
      , const Node &gene_left
      , const Node &species_right
      , const Node &gene_right
      , SpeciesGeneMap &genemap_temp
      , SpeciesGeneMap &genemap
      , DistanceMatrix &dmatrix
      , const Node &polytomy_root
      , const bool computeBranchLengths = false
      , const bool verbose = false
  ) const;

 public:
  std::map<Event, std::vector<Node>> operator()(
      bpp::PhyloTree &genetree /// Gene tree.
      , const Node &polytomy_root /// Root of the polytomy in the gene tree.
      , const bpp::PhyloTree &speciestree /// Guide/Species tree.
      , DistanceMatrix &distances /// Distance Matrix for the given polytomy leaves
      , const std::unordered_map<std::shared_ptr<bpp::PhyloNode>, std::size_t> &V /// Count vector.
      , SpeciesGeneMap &genemap /// Map which associates gene of the polytomy to species.
      , const bool computeBranchLengths = false /// Indicates if branch lengths has to be computed.
      , const bool verbose = false /// Show operations in terminal.
  ) const;
};

} // namespace treerecs

#endif //TREERECS_PROFILENJ_H
