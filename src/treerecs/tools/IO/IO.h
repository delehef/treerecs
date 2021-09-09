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

#ifndef PHYLASOLVER_IO_H
#define PHYLASOLVER_IO_H

// Include libs
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <memory>
#include <stdexcept>

// Include Bpp
#include <Bpp/Phyl/Tree/PhyloTree.h>

// Include containers
#include <treerecs/containers/ReconciledRootedTree.h>
#include <treerecs/containers/NodeProperty.h>
#include <treerecs/containers/Table.h>

//Include tools
#include "Nhx.h"
#include "Newick.h"
#include "PhyloXML.h"

namespace treerecs {

/// Format of tree's inputs/outputs.
enum class TextFormat {
  newick,
  nhx,
  phyloxml,
  recphyloxml,
  svg,
  unknown
};

inline TextFormat operator++(TextFormat& x) {
  return x = static_cast<TextFormat>((static_cast<int>(x) + 1));
}

/*!
 * @class IO
 * @brief IO provides functions to load a bpp::PhyloTree or write it into a Newick, Nhx or PhyloXML (and RecPhyloXML) file.
 * @details Trees are returned as std::shared_ptr<bpp::PhyloTree> and each node is named.
 */
class IO {
 private:

 protected:
  /// Newick module to write/ read Newick parenthesis format.
  static const Newick newick_;
  /// NHX module to write/ read NHX parenthesis format.
  static const Nhx nhx_;
  /// PhyloXML module to write in PhyloXML and RecPhyloXML format.
  static const PhyloXML phyloXML_;

  /// \brief Set a name for each unnamed internal node, check if the tree has branch lengths in options.
  /// \param tree tree to modify
  /// \param check_branch_length check if tree has correct branch lengths (default = true)
  /// \param force_branch_support set branch support when the branch has none
  static void preProcesses(bpp::PhyloTree& tree
      , const bool check_branch_length = true
      , const bool force_branch_support = false
  );

  /// \brief Rename homonyms in tree by adding a number.
  /// \tparam NodesIterator Iterator on std::shared_ptr<bpp::PhyloNode>
  /// \param nodes_begin
  /// \param nodes_end
  template<typename NodesIterator>
  static void renameNodeHomonyms(const NodesIterator &nodes_begin,
                                 const NodesIterator &nodes_end);

 public:
  static constexpr char preferred_separator =
    #ifdef _WIN32
      '\\';
    #else
      '/';
    #endif
  /// Create directory
  static void MakeDir(const std::string path, uint32_t mode = 0);

  /// Check file existence.
  static bool exists(const std::string& filename);
  static bool IsRegularFile(const char* path);
  static bool IsDirectory(const char* path);
  static bool IsDirEmpty(const char* path);

  /// Compute number of lines in a file. The file is going to be read and reinit.
  static unsigned long nlines(std::ifstream& is);

  /****************
   * LOAD
   */
  /// Translates newick std::string to bpp::PhyloTree.
  static std::shared_ptr<bpp::PhyloTree> nhxToPhyloTree(const std::string& description);

  /// Translates newick std::string to bpp::PhyloTree.
  static std::shared_ptr<bpp::PhyloTree> newickToPhyloTree(const std::string& description, bool bootstrap = false, const std::string& propertyName = "", bool withId = false, bool verbose = false);

  /// Read a tree file in any supported format
  static std::shared_ptr<bpp::PhyloTree> readTreeFile(
      const std::string &filename
      , const bool support = false
      , const bool check_branch_length = true
      , const bool verbose = false);

  /// Read a tree file and return a std::shared_ptr<bpp::PhyloTree>.
  static std::shared_ptr<bpp::PhyloTree> readTreeFile(const std::string &filename /// name of the newick file.
      , const TextFormat format /// file format
      , const bool support = false /// get support values.
      , const bool check_branch_length = true /// check branch lengths. Prints warning when there is no branch length.
      , const bool verbose = false /// verbose mode.
  );

  /// Read a tree file and return a vector of std::shared_ptr<bpp::PhyloTree>.
  static std::vector<std::shared_ptr<bpp::PhyloTree>>
  readTreesFile(const std::string &filename
                , const TextFormat format = TextFormat::newick
                , const int index = -1
                , const bool support = true
                , const bool check_branch_length = true
                , const bool printProgression = true
                , const bool verbose = false
  );

  /// Load a Distance matrix from a file. The tree is important to attribute each node of the matrix to the tree.
  static Table<double, std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>>
  readDistanceMatrix(const std::string& filename /// filename.
      , const bpp::PhyloTree& tree /// tree associated with the distance matrix.
      , const bool verbose = false /// verbose mode.
  );

/****************
 * SAVE
 */
  /// \brief Translate bpp::PhyloTree into newick's parenthesis format.
  /// \return Std::string which contains the tree in the newick parenthesis format.
  static std::string PhyloTreeToNewick(const bpp::PhyloTree &tree);

  /// \brief Translate bpp::PhyloTree into newick's parenthesis format.
  /// @param tree The tree to convert.
  /// @param bootstrap Tell is bootstrap values must be writen.
  ///   If so, the content of the property with name "bootstrap" will be written as bootstrap value.
  ///   The property should be a Number<double> object.
  ///   Otherwise, the content of the property with name 'propertyName' will be written.
  ///   In this later case, the property should be a String object.
  /// @param propertyName The name of the property to use. Only used if bootstrap = false.
  /// \return Std::string which contains the tree in the newick parenthesis format.
  static std::string PhyloTreeToNewick(const bpp::PhyloTree& tree, bool bootstrap, const std::string& propertyName);

  /// \brief Translate bpp::PhyloTree into nhx parenthesis format.
  /// \return Std::string which contains the tree in the NHX parenthesis format.
  static std::string PhyloTreeToNhx(const bpp::PhyloTree &tree);

  /// \brief Translate bpp::PhyloTree into PhyloXML format.
  /// \return Std::string which contains the tree in the PhyloXML format.
  static std::string PhyloTreeToPhyloXML(const bpp::PhyloTree &tree);

  /// \brief Get Reconciliation into PhyloXML format.
  /// \return Std::string which contains the reconciliation in the PhyloXML format.
  static std::string ReconciliationToRecPhyloXML(const bpp::PhyloTree& genetree, const bpp::PhyloTree& speciestree, const SpeciesGeneMap& map);

  /// \brief Print a tree in an ostream in newick (parenthesis) format.
  /// \return Nothing but edit an output stream by writing the tree in parenthesis format (default as Newick).
  static void write(const bpp::PhyloTree& tree /// tree to print.
      , std::ostream& os /// stream object.
      , const TextFormat output = TextFormat::newick /// print tree in format.
      , const std::string& description = "" /// Tree description.
  );

  /// \brief Write Reconciled trees in an std::ostream.
  static void writeReconciledTrees(
      std::ostream &os,
      const std::unordered_map<std::shared_ptr<bpp::PhyloTree>,
       std::map<double, std::vector<ReconciledRootedTree>>>& reconciledTrees,
      const bool& strict_thresholds,
       const std::vector<std::shared_ptr<bpp::PhyloTree>>& trees,
       const bpp::PhyloTree& speciestree,
       const bool printTree,
       const TextFormat output,
       const bool description = true);

  /// \brief Write Reconciled trees in a file.
  static void writeReconciledTrees(
      const std::string& filename,
      const std::unordered_map<std::shared_ptr<bpp::PhyloTree>,
          std::map<double, std::vector<ReconciledRootedTree>>>& reconciledTrees,
      const bool& strict_thresholds,
      const std::vector<std::shared_ptr<bpp::PhyloTree>>& trees,
      const bpp::PhyloTree& speciestree,
      const TextFormat output,
      const bool description = true
  );


  /// \brief Write an orthologous/paralogous file.
  static void writeEventSummary(
      const std::string& filename,
      const std::unordered_map<std::shared_ptr<bpp::PhyloTree>,
      std::map<double, std::vector<ReconciledRootedTree>>> &reconciledTrees,
      const std::vector<std::shared_ptr<bpp::PhyloTree>>& trees,
      const bool full_reconciliation = false
  );

  /// \brief Write an orthologous/paralogous file.
  static void writeEventSummary(std::ostream& os, const ReconciledRootedTree& recgtree, bool full_reconciliation = false);
};

/// Write a TextFormat into an output stream.
std::ostream& operator<<(std::ostream& os,
                         const TextFormat& output_format);

} // namespace treerecs

#endif //PHYLASOLVER_IO_H
