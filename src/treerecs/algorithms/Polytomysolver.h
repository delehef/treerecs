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

#ifndef TREERECS_POLYTOMYSOLVER_H
#define TREERECS_POLYTOMYSOLVER_H

// Include bpp
#include <Bpp/Phyl/Tree/PhyloTree.h>

// Include Treerecs containers
#include <treerecs/containers/Table.h>
#include <treerecs/containers/Cost.h>
#include <treerecs/containers/GeneMap.h>

namespace treerecs {

using CostTable = Table<Cost, std::shared_ptr<bpp::PhyloNode>, int>;

using Node = std::shared_ptr<bpp::PhyloNode>;

/**
 * \class Polytomysolver
 * \brief Solve polytomies in gene trees.
 * \details Used to get the CostMatrix which gives the path of minimal cost in
 *          gene events in a given polytomy.
 *
 *          A CostMatrix is a dynamic programming table associating species
 *          node with a number of genes which are coming into them.
 *          Each cell of this matrix is associated with a cost to get this
 *          event (to have n genes coming into a species), see Cost.
 */
class Polytomysolver {
 public:
  /*!
   * @brief Choose a node path that conduct to this cost (node, m).
   * @param table Cost table which contains the costs for each species at each count.
   * @param node Species
   * @param m Number of genes
   */
  std::size_t chooseANodePath(const CostTable& table
                            , const Node& node
                            , const std::size_t& m) const;

  /*!
   * @brief Generate the polytomysolver's matrix of costs
   * @param speciestree species tree.
   * @param genemap map associating genes to the species in species tree.
   *        See GeneMap and SpeciesGeneMapper.
   * @param dupCost cost of a duplication.
   * @param lossCost cost of a loss.
   * @param verbose print all operations.
   * @return
   */
  CostTable computeCostTable(
      const bpp::PhyloTree& speciestree,
      const SpeciesGeneMap& genemap,
      double dupCost,
      double lossCost,
      bool verbose = false) const;

  /*!
   * The count vector gives the path to get a reconciliation accroding to the
   * matrix of costs (see computeCostTable)
   * @param speciestree the speciestree
   * @param cost_table the matrix of costs
   * @param map map associating genes to species
   * @return
   */
  std::unordered_map<std::shared_ptr<bpp::PhyloNode>, std::size_t>
  getCountVector(
      const bpp::PhyloTree &speciestree
      , const CostTable &cost_table
      , const SpeciesGeneMap& map) const;
};

} // namespace treerecs

#endif //TREERECS_POLYTOMYSOLVER_H
