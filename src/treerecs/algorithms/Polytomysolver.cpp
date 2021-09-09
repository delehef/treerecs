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

#include "Polytomysolver.h"

// Include std
#include <queue>
#include <limits>
#include <list>
#include <vector>

// Include Treerecs tools
#include <treerecs/tools/random.h>
#include <treerecs/tools/SpeciesGeneMapper.h>
#include <treerecs/tools/utils.h>

namespace treerecs {

std::size_t Polytomysolver::chooseANodePath(
    const CostTable &table, const Node &node, const std::size_t &m
) const {
  auto& occurences = table(node, m).occurence;
  if (occurences.size() == 1) return 0;
  return (WeightedRandomPick(occurences));
}


CostTable Polytomysolver::computeCostTable(
    const bpp::PhyloTree &speciestree, const SpeciesGeneMap &genemap
    , double dupCost, double lossCost, bool verbose
) const {
  /*!
   * Build the cost matrix of a given polytomy according to the species tree.
   * Genemap contains only genes of the polytomy.
   */

  if (verbose)
    std::cout << "Computing cost table with Polytomysolver..." << std::endl;

  auto species_list =
      PhyloTreeToolBox::getNodesInPostOrderTraversalRecursive_list(speciestree);
  std::vector<Node> species{species_list.begin(), species_list.end()};

  /* Looking for the maximum number of genes (max_count) for species */
  int max_count = 0;
  for (auto pair : genemap.u_map()) {
    if (static_cast<int>(pair.second.size()) > max_count)
      max_count = static_cast<int>(pair.second.size());
  }

  /* Build the Table */
  int M_MAX = max_count + 1;  // Column maximum value
  CostTable table(species.size(), static_cast<std::size_t>(M_MAX+1));
  // With Event constructors, values are initialized to 0.0 and path to {none}
  std::vector<int> columns;
  for (int m = 0; m < M_MAX+1; m++) {
    columns.push_back(m);
    // m, the column is the number of genes for one species.
  }
  table.setColIndexes(columns);  // set indexes
  table.setRowIndexes(species);  // set indexes

  /* Fill the Table */
  for (auto& specie : species) {
    // We have zeropos when the number of node from a specie is the same as the
    // column number
    // Fill the table, using the next/previous case cost
    // The node is a leaf, just fill with dupcost and losscost
    int zeropos = static_cast<int>(genemap.ngenes(specie));
    if (speciestree.isLeaf(specie)) {
      // find the column with a cost of zero (zeropos) and fill the table
      // according to this position
      // by default, all the position in the table are 0
      int m = zeropos - 1;  // m is the number og genes
      while (m >= columns[0]) {
        if (m == 0)
          table(specie, m).value =  std::numeric_limits<double>::infinity();
        else
          table(specie, m).value = table(specie, m + 1).value + dupCost;
        table(specie, m).path = {duplication};
        table(specie, m).occurence = {1};
        m -= 1;
      }
      m = zeropos + 1;
      while (m <= M_MAX) {
        table(specie, m).value = table(specie, m - 1).value + lossCost;
        table(specie, m).path = {loss};
        table(specie, m).occurence = {1};
        m += 1;
      }
    } else {
      // Here we have an internal node (not a leaf in the genetree)
      auto sons = speciestree.getSons(specie);
      assert(sons.size() < 3);  // There is no polytomy in the tree.
      auto& son_left = sons[0];
      auto& son_right = sons[1];

      // Fill the table using only the speciation cost (sum of the children's
      // cost of this node specie)
      for (int k = 0; k <= M_MAX; k++) {
        if (k >= zeropos) {
          table(specie, k).value = table(son_left, k - zeropos).value +
                                   table(son_right, k - zeropos).value;
          table(specie, k).path = {speciation};
          table(specie, k).occurence = {
              utils::sum(table(son_left, k - zeropos).occurence) *
              utils::sum(table(son_right, k - zeropos).occurence)};
        } else {
          table(specie, k).value = std::numeric_limits<double>::infinity();
        }
      }

      // Find min of the line and save its position
      std::list<int> minpos = {1};
      auto minvalue = table.getRow(specie)[minpos.front()].value;
      for (auto k : columns) {
        if (table(specie, k).value < minvalue) {
          minvalue = table(specie, k).value;
          minpos.clear();
          minpos.push_back(k);
        } else if (table(specie, k).value == minvalue) {
          minpos.push_back(k);
        }
      }

      // Looking for duplications and losses of genes according to the minimum
      // values of the line.
      for (auto& pos : minpos) {
        auto m = pos - 1;
        while (m > 0) {  // Check duplications
          if (table(specie, m).value == table(specie, m + 1).value + dupCost) {
            // If the current value is equal to the equivalent one with a
            // duplication
            auto &cost = table(specie, m);
            // Add a duplication event if there is none.
            if (std::find(cost.path.begin(), cost.path.end(), duplication) ==
                cost.path.end()) {
              cost.path.push_back(duplication);
              cost.occurence.push_back(
                  utils::sum(table(specie, m + 1).occurence));
            }
          } else if (table(specie, m).value >
                     table(specie, m + 1).value + dupCost) {
            table(specie, m).value = table(specie, m + 1).value + dupCost;
            table(specie, m).path = {duplication};
            table(specie, m).occurence = {
                utils::sum(table(specie, m + 1).occurence)};
          }
          m -= 1;
        }

        m = pos + 1;
        while (m <= M_MAX) {  // Check losses
          if (table(specie, m).value == table(specie, m - 1).value + lossCost) {
            auto &cost = table(specie, m);
            if (std::find(cost.path.begin(), cost.path.end(), loss) ==
                cost.path.end()) {
              cost.path.push_back(loss);
              cost.occurence.push_back(
                  utils::sum(table(specie, m - 1).occurence));
            }
          } else if (table(specie, m).value >
                     table(specie, m - 1).value + lossCost) {
            table(specie, m).value = table(specie, m - 1).value + lossCost;
            table(specie, m).path = {loss};
            table(specie, m).occurence = {
                utils::sum(table(specie, m - 1).occurence)};
          }
          m += 1;
        }
      }
    }
  }
  if (verbose)
    std::cout << "...cost table computed:" << std::endl
              << table << std::endl;

  return table;
}

std::unordered_map<std::shared_ptr<bpp::PhyloNode>, std::size_t>
Polytomysolver::getCountVector(
    const bpp::PhyloTree &speciestree, const CostTable &cost_table
    , const SpeciesGeneMap &map
) const {
  /// Returns the count vector, resulting from the backtracking path of the
  /// cost_table.
  using Node = std::shared_ptr<bpp::PhyloNode>;
  std::unordered_map<Node, std::size_t> counts;
  std::vector<Node> nodes = cost_table.getRowIndexes();

  // Not a need if the costTable is built with
  // Polytomysolver::computeCostTable(...)
  PhyloTreeToolBox::sortNodesByIndex(nodes, speciestree, true);
  Node &node = nodes.back();  // The last element is the root of the tree (with
  // the higher index in speciestree).

  std::queue<Node> nodeQueue;  // Queue of the nodes
  nodeQueue.push(node);  // We start the algorithm by the root...
  counts[node] = 1;  // and its count to 1 gene copy.

  while (not nodeQueue.empty()) {
    node = nodeQueue.front();  // get the next node to work.
    if (not speciestree.isLeaf(node)) {
      std::size_t m = counts.at(node);
      Event event = none;
      unsigned int stop = 0;  // Prevents infinite loop
      while (event != speciation) {
        // while we are not in a speciation event, continue to read the line and
        // accumulate events like duplications or losses.
        std::size_t nodePathIndex = 0;
        auto &nodepaths = cost_table(node, m).path;
        if (nodepaths.size() > 1) {
          nodePathIndex = chooseANodePath(cost_table, node, m);
        }

        event = nodepaths.at(nodePathIndex);

        if (event == duplication) {
          m += 1;
        } else if (event == loss) {
          m -= 1;
        } else if (event != speciation) {
          std::cout << std::endl;
          std::cerr << "Error during node " << node
                    << " backtracking inconsistent event." << std::endl;
          std::cerr << "Cost table:" << std::endl;
          std::cerr << cost_table << std::endl;
          std::cerr << "Guide tree:" << speciestree << std::endl;
          exit(EXIT_FAILURE);
        }

        stop++;

        if (stop > cost_table.ncol() + 1) {
          std::cout << std::endl;
          std::cerr << "Error during node  " << node << " backtracking loop."
                    << std::endl;
          break;
        }
      }

      for (auto son : speciestree.getSons(node)) {
        counts[son] = m - map.getGenes(node).size();
        if (not speciestree.isLeaf(son))
          nodeQueue.push(son);
      }
    }
    nodeQueue.pop();
  }

  return counts;
}

} // namespace treerecs
