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

#include "utils_stream.h"

#include <treerecs/tools/PhyloTreeToolBox.h>

namespace treerecs {

std::ostream &operator<<(std::ostream &os, const bpp::PhyloTree &tree) {
  /// Print a bpp::PhyloTree in newick format.
  bpp::Newick newick = bpp::Newick();
  newick.write(tree, os);
  return os;
}

std::ostream &
operator<<(std::ostream &os, const std::shared_ptr<bpp::PhyloNode> &node_ptr) {
  /// Print the name of a std::shared_ptr<bpp::PhyloNode>.
  if (node_ptr) {
    if (node_ptr->hasName())
      os << node_ptr->getName();
    else
      os << "Noname";
  } else
    os << "nullptr";
  return os;
}

namespace utils {

void printNodeContent(
    const bpp::PhyloTree &tree, const std::shared_ptr<bpp::PhyloNode> &node
    , std::ostream &os
) {
  ///Print Node Content (infos, father, sons).

  // First, print father and itself
  if (tree.getRoot() != node and tree.hasFather(node)) {
    auto father = tree.getFather(node);
    auto edgeFatherNode = tree.getEdgeLinking(father, node);
    if (edgeFatherNode) {
      if (edgeFatherNode->hasLength())
        os << father << " -(" << tree.getEdgeLinking(father, node)->getLength()
           << ")-> " << node << std::endl;
      else
        os << father << " ---> " << node << std::endl;
    } else
      os << father << " -x-> " << node << std::endl;

  } else {
    os << node << " (root)" << std::endl;
  }

  // Then print sons
  if (not tree.isLeaf(node)) {
    auto sons = tree.getSons(node);
    for (auto son : sons) {
      auto edgeNodeSon = tree.getEdgeLinking(node, son);
      if(edgeNodeSon) {
        if (edgeNodeSon->hasLength())
          os << "         '-(" << tree.getEdgeLinking(node, son)->getLength()
             << ")-> " << son << std::endl;
        else
          os << "         '---> " << son << std::endl;
      }
      else
        os << "         '-x-> " << son << std::endl;
    }
  }
}

void print_temp(bpp::PhyloTree &tree, std::ostream& os) {
  /// Print all bifurcations in a given tree.
  PhyloTreeToolBox::applyInPOT(tree,
                               [&tree, &os](
                                   std::shared_ptr<bpp::PhyloNode> &node
                               ) {
                                 printNodeContent(tree, node, os);
                               });
}

std::string time_to_str(const double time) {
  std::string message;
  double clocks_remains = time;

  auto days = floor(clocks_remains / 86400.);
  clocks_remains -= 86400. * days;

  auto hours = floor(clocks_remains / 3600.);
  clocks_remains -= 3600. * hours;

  auto minutes = floor(clocks_remains / 60.);
  clocks_remains -= 60. * minutes;

  if (days > 0)
    message += std::to_string(static_cast<unsigned long>(days))
               + std::string(" d. ");
  if (hours > 0)
    message += std::to_string(static_cast<unsigned long>(hours))
               + std::string(" h. ");
  if (minutes > 0)
    message += std::to_string(static_cast<unsigned long>(minutes))
               + std::string(" m. ");

  message += std::to_string(static_cast<unsigned long>(clocks_remains))
             + std::string(" s.");
  return message;
}

} // namespace utils

} // namespace treerecs
