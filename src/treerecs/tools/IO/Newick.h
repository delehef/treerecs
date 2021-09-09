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

#ifndef TREERECS_NEWICK_H
#define TREERECS_NEWICK_H

#include <Bpp/Phyl/Io/Newick.h>

namespace treerecs {

/*!
 * @class Newick
 * @brief Newick provides functions to load a bpp::PhyloTree or write it into a Newick file (inheritance from bpp::Newick).
 */
class Newick: public bpp::Newick {
  // Prevents hiding valid overrides
  using bpp::Newick::write;
public:
  Newick(const bool allows_comments = true): bpp::Newick(allows_comments) {};

  /// Returns bpp::PhyloTree into parenthesis (newick) format.
  std::string treeToParenthesis(
      const bpp::PhyloTree& tree /// Tree to translate into parenthesis format.
      , bool reconciliation_mode = false /// Add flags specified in the node property "flag" in branches of the parenthesis format.
  ) const {
    if(reconciliation_mode)
      return Newick::treeToParenthesis(tree, false, "flag");
    return bpp::Newick::treeToParenthesis(tree, bpp::Newick::writeId_);
  }

  /// @brief Returns bpp::PhyloTree into parenthesis (newick) format. Used to get Newick::treeToParenthesis public.
  /// @param tree The tree to convert.
  /// @param bootstrap Tell is bootstrap values must be writen.
  ///   If so, the content of the property with name "bootstrap" will be written as bootstrap value.
  ///   The property should be a Number<double> object.
  ///   Otherwise, the content of the property with name 'propertyName' will be written.
  ///   In this later case, the property should be a String object.
  /// @param propertyName The name of the property to use. Only used if bootstrap = false.
  std::string treeToParenthesis(const bpp::PhyloTree& tree, bool bootstrap, const std::string& propertyName) const {
    return bpp::Newick::treeToParenthesis(tree, bootstrap, propertyName);
  }

  /// Writes bpp::PhyloTree into an output stream.
  void write(
      const bpp::PhyloTree& tree /// Tree to translate into parenthesis format.
      , std::ostream& os /// std::ostream where tree is being printed.
      , bool reconciliation_mode = false /// Add flags specified in the node property "flag" in branches of the parenthesis format.
  ) const {
    if (! os) { throw bpp::IOException ("Newick::writeTree: failed to write to stream"); }
    std::string content;

    // If the reconciliation mode is activated, bootstrap are replaced by the name of the internal node.
    content = Newick::treeToParenthesis(tree, reconciliation_mode);

    os << content;
  }
};

} // namespace treerecs

#endif //TREERECS_NEWICK_H
