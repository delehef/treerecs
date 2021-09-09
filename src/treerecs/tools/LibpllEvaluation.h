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
#ifndef TREERECS_LIBPLLEVALUATE_H
#define TREERECS_LIBPLLEVALUATE_H

extern "C" {
#include <pll.h>
#include <pllmod_algorithm.h>
#include <pll_binary.h>
#include <pll_msa.h>
#include <pll_optimize.h>
#include <pll_tree.h>
#include <pllmod_util.h>  
}

#include <string>
#include <memory>
#include <exception>
#include <vector>
#include <Bpp/Phyl/Tree/PhyloTree.h>

namespace treerecs {

class PhyloTree;
struct pll_sequence;
using pll_sequence_ptr = std::shared_ptr<pll_sequence>;
using pll_sequences = std::vector<pll_sequence_ptr>;

/*
 *  @brief Exception thrown from LibpllEvaluation
 */
class LibpllException: public std::exception {
public:
  LibpllException(const std::string &s): msg_(s) {}
  LibpllException(const std::string &s1, 
      const std::string s2): msg_(s1 + s2) {}
  virtual const char* what() const throw() { return msg_.c_str(); }
  void append(const std::string &str) {msg_ += str;}

private:
  std::string msg_;
};

struct LibpllAlignmentInfo {
  std::string alignmentFilename;
  std::string model;
};

/*
 * @brief Libpll wraper to compute the phylogenetic likelihood of a tree.
 * 
 * LibpllEvaluation provides methods to compute the phylogenetic
 * likelihood of tree using Felsenstein pruning algorithm.
 */
class LibpllEvaluation {
public:
  /*
   * @brief Build a LibpllEvaluation instance
   * @param newickString the tree in newick format
   * @param alignmentFilename path to the msa file
   * @param modelStr a string representing the model (GTR, DAYOFF...)
   * @return a shared pointer wraping the LibpllEvaluation instance
   */
  static std::shared_ptr<LibpllEvaluation> buildFromString(const std::string &newickString,
      const std::string& alignmentFilename,
      const std::string &modelStr);

  static std::shared_ptr<LibpllEvaluation> buildFromFile(const std::string &newickTree,
      const LibpllAlignmentInfo &info);

  static std::shared_ptr<LibpllEvaluation> buildFromPhylo(std::shared_ptr<bpp::PhyloTree> phyloTree,
      const LibpllAlignmentInfo &info);

  static std::vector<LibpllAlignmentInfo> parseAlignmentInfo(
      const std::string &filename,
      const int tree_index = -1);
  /*
   *  @brief Compute the likelihood of the tree given the alignment
   *  @return the likelihood of the tree
   */
  double computeLikelihood();

  double optimizeAllParameters();
private:
  
  /**
   * Constructors
   */
  LibpllEvaluation():treeinfo_(nullptr) {}
  LibpllEvaluation(const LibpllEvaluation &) = delete;
  
  /**
   * @brief set all the null branch lenghts to length
   */
  static void setMissingBL(pll_utree_t * tree, 
    double length);

  /**
   *  @brief parse sequences and pattern weights from fasta file
   *  @param fasta_file Input file
   *  @param map state map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parseFasta(const char *fasta_file, 
    const pll_state_t *map,
    pll_sequences &sequences,
    unsigned int *&weights);

  /**
   *  @brief parse sequences and pattern weights from phylip file
   *  @param phylip_file Input file
   *  @param map state map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parsePhylip(const char *phylip_file, 
    const pll_state_t *map,
    pll_sequences &sequences,
    unsigned int *&weights);

  static double optimizeAllParametersOnce(pllmod_treeinfo_t *treeinfo);
  
private:
  std::shared_ptr<pllmod_treeinfo_t> treeinfo_;
  std::shared_ptr<pll_utree_t> utree_;
};

} // namespace treerecs

#endif
