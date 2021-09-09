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

#include "IO.h"

#include <sys/stat.h>
#include <dirent.h>
#ifndef _WIN32
  #include <err.h>
#endif

#include <string>

// Include Bpp
#include <Bpp/BppString.h>

// Include Treerecs-code
#include <treerecs/tools/PhyloTreeToolBox.h>
#include <treerecs/tools/ALE/ALEevaluation.h>
#include <treerecs/tools/utils.h>
#include <treerecs/tools/Timer.h>
#include "RecPhyloTreeToSVG.h"
#include "RefreshablePrinter.h"

using std::cerr;
using std::endl;

namespace treerecs {

// Directory separator
constexpr char IO::preferred_separator;

/// bpp::Newick module.
const Newick IO::newick_ = Newick();
/// bpp::NHX module.
const Nhx IO::nhx_ = Nhx(false);
/// bpp::PhyloXML module.
const PhyloXML IO::phyloXML_ = PhyloXML();

std::shared_ptr<bpp::PhyloTree> IO::nhxToPhyloTree(const std::string &description){
  return std::shared_ptr<bpp::PhyloTree>(nhx_.parenthesisToPhyloTree(description));
}

std::shared_ptr<bpp::PhyloTree>
IO::newickToPhyloTree(const std::string &description, bool bootstrap, const std::string &propertyName, bool withId,
                      bool verbose){
  return std::shared_ptr<bpp::PhyloTree>(newick_.parenthesisToPhyloTree(description, bootstrap, propertyName, withId, verbose));
}

/// Read a tree file in any supported format
std::shared_ptr<bpp::PhyloTree> IO::readTreeFile(
    const std::string &filename,
    const bool support,
    const bool check_branch_length,
    const bool verbose) {
  std::shared_ptr<bpp::PhyloTree> tree;

  // Check species tree provided and exists
  if (not(IO::exists(filename))) {
    if (filename.empty())
      std::cerr
          << "Please provide a species tree file with parameter -s or --speciestree"
          << std::endl;
    else
      std::cerr << "Error: species tree file \"" << filename
                << "\" does not exist." << std::endl;

    exit(EXIT_FAILURE);
  }

  // Read species tree
  for (TextFormat format = TextFormat::newick ;
       format != TextFormat::recphyloxml ; ++format) {
    assert(format != TextFormat::unknown);
    try {
      tree = IO::readTreeFile(filename, format, support, check_branch_length, verbose);
      break;
    } catch (std::exception& ee) {}
  }

  if (not tree) {
    std::cerr << "Error: there is no species tree read in "
              << filename << "." << std::endl;
    std::cerr
        << "       Please try with an other format (Newick, NHX or PhyloXML)."
        << std::endl;
    exit(EXIT_FAILURE);
  }

  return tree;
}

std::shared_ptr<bpp::PhyloTree>
IO::readTreeFile(const std::string &filename, const TextFormat format, const bool support, const bool check_branch_length,
                 const bool verbose){
  return IO::readTreesFile(filename, format, 0, support, check_branch_length, 0, verbose).front();
}

std::vector<std::shared_ptr<bpp::PhyloTree>>
IO::readTreesFile(const std::string &filename, const TextFormat format, const int index, const bool support,
                  const bool check_branch_length, const bool printProgression, const bool verbose) {

  if(verbose) std::cout << "Read " << filename << " with " << format << " format." << std::endl;
  std::vector<std::shared_ptr<bpp::PhyloTree>> trees;

  Timer<std::size_t> progression(1);

  if(printProgression)
    utils::progressionBar(std::cout, progression, false, "Reading tree(s)");

  // Open the file
  std::ifstream genefile(filename.c_str(), std::ios::in);
  try {
    if (not genefile.is_open()) {
      throw std::invalid_argument(filename + std::string(" does not exist."));
    }
  } catch(std::exception const& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  // Get file content and strip line feeds.
  std::string file_content((std::istreambuf_iterator<char>(genefile)), std::istreambuf_iterator<char>() ); // content of the file

  /*
  unsigned int corresponding_nchar_with_linefeed = 0;
  unsigned int nchar_linefeed = strlen(LINE_FEED);
  for( auto str_it = file_content.begin() ; str_it != file_content.end() ; str_it++) {
      if(corresponding_nchar_with_linefeed == nchar_linefeed) {
          assert(*(str_it -corresponding_nchar_with_linefeed) == LINE_FEED[0]);
          file_content.erase(str_it - corresponding_nchar_with_linefeed, str_it);
          corresponding_nchar_with_linefeed = 0;
      }

      assert(corresponding_nchar_with_linefeed < nchar_linefeed);

      if((*str_it) == LINE_FEED[corresponding_nchar_with_linefeed]){
          corresponding_nchar_with_linefeed++;
      } else {
          if(corresponding_nchar_with_linefeed > 0) {
              corresponding_nchar_with_linefeed = 0;
          }
      }
  }

  if(corresponding_nchar_with_linefeed == nchar_linefeed) {
      assert(*(file_content.end() -corresponding_nchar_with_linefeed) == LINE_FEED[0]);
      file_content.erase(file_content.end() - corresponding_nchar_with_linefeed, file_content.end());
      corresponding_nchar_with_linefeed = 0;
  }
  */

  file_content.erase(
      std::remove(file_content.begin(), file_content.end(), '\n')
      , file_content.end());

  #if defined _WIN32 || defined __CYGWIN__
  file_content.erase(
      std::remove(file_content.begin(), file_content.end(), '\r')
      , file_content.end());
  #endif

  utils::trim_str(file_content);

  genefile.close();

  std::size_t number_of_trees = 0;
  if(format == TextFormat::newick or format == TextFormat::nhx) {
    number_of_trees = utils::count(file_content, ";");
  } else if(format == TextFormat::phyloxml) {
    number_of_trees = utils::count(file_content, "phylogeny")/2;
  } else {
    std::cerr << "Error: not supported tree format." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(verbose) std::cout << "Number of trees found in " << filename << ": " << number_of_trees << std::endl;

  if(number_of_trees == 0) {
    throw std::runtime_error("IO: no tree read in " + filename + ".");
    return {};
  };

  if(index == -1) {
    trees.reserve(number_of_trees);
    progression.setMaximum(number_of_trees);
  }
  else if (index >= 0) {
    if(static_cast<std::size_t>(index) >= number_of_trees) {
      std::cerr << "Error: index " << index + 1
                << " is greater than the number of trees in " << filename
                << " (" << number_of_trees << ")." << std::endl;
      exit(EXIT_FAILURE);
    }
    trees.reserve(1);
    progression.setMaximum((std::size_t)index);
  }

  if(format == TextFormat::newick) {
    auto contents = utils::splitString(file_content, ";");
    int i = 0;
    for(auto content: contents){
      if(index == -1 or (index >= 0 and i == index)) {
        auto tree = std::shared_ptr<bpp::PhyloTree>(IO::newickToPhyloTree(content + ";", support));
        if(tree) {
          trees.emplace_back(tree);
          if (verbose) std::cout << *trees.back() << std::endl;
          progression.next();
          if (printProgression)
            utils::progressionBar(std::cout, progression, true, "Reading tree(s)");
        } else {
            std::cout << std::endl << "Error: bad format for \"" << content << "\"" << std::endl;
          throw std::runtime_error("IO: bad format in " + filename + ".");
        }
      }
      i++;
    }
  }
  else if(format == TextFormat::nhx) {
    auto contents = utils::splitString(file_content, ";");
    int i = 0;
    for(auto content: contents){
      if(index == -1 or (index >= 0 and i == index)) {
        auto tree = std::shared_ptr<bpp::PhyloTree>(IO::nhxToPhyloTree(content + ";"));
        if(tree) {
          trees.emplace_back(tree);
          progression.next();
          if (printProgression)
            utils::progressionBar(std::cout, progression,
                                  true, "Reading tree(s)");
        }else{
          throw std::runtime_error("IO: bad format in " + filename + ".");
        }
      }
      i++;
    }
  }
  else if(format == TextFormat::phyloxml) {
    std::vector<bpp::PhyloTree*> temp_trees;
    if(printProgression)
      utils::progressionBar(std::cout, progression,
                            false, "Reading "
                                   + std::to_string(number_of_trees)
                                   + " tree(s)");
    temp_trees = phyloXML_.phyloXMLToPhyloTrees(file_content);
    int i = 0;
    for(bpp::PhyloTree* tree: temp_trees) {
      if(index == -1 or (index >= 0 and i == index)) {
        trees.emplace_back(std::shared_ptr<bpp::PhyloTree>(tree));
        progression.next();
        if(printProgression)
          utils::progressionBar(std::cout, progression,
                                false, "Reading "
                                       + std::to_string(number_of_trees)
                                       + " tree(s)");
      } else
        delete tree;
      i++;
    }
  }

  if(trees.size() == 0) {
    std::cerr << "No tree found in " << filename << " with " << format << " format." << std::endl;
    throw std::runtime_error("IO: no tree read in " + filename + ".");
    return {nullptr};
  };

  for(std::size_t i = 0; i < trees.size(); ++i) {
    // Read the content and control

    auto& tree = trees.at(i);

    if(printProgression)
      utils::progressionBar(std::cout, progression, false, "Finish read");
    assert(tree != nullptr);
    if(tree->getName().empty()){
      if(index >= 0)
        tree->setName("Tree" + std::to_string(index + 1));
      else
        tree->setName("Tree" + std::to_string(i + 1));
    }

    // Set a name to each internal node (only leaves has a name yet).
    preProcesses(*tree, check_branch_length, false);
  }

  return trees;
}

Table<double, std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>>
IO::readDistanceMatrix(const std::string &filename, const bpp::PhyloTree &tree, const bool verbose){
  // Open file
  if (verbose) std::cout << "Opening distance matrix file " << filename << "..." << std::endl;
  std::ifstream file(filename);
  try {
    if (not file.is_open()) {
      //std::cerr << "Error in opening " << filename << std::endl;
      throw std::invalid_argument(filename + std::string(" does not exist."));
    }
  } catch(std::exception const& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  // First line = node names
  std::string line;
  getline(file, line);
  auto nodenames = utils::splitString(line, " ");
  std::vector<std::shared_ptr<bpp::PhyloNode>> nodes;
  for(auto& name: nodenames){
    nodes.push_back(PhyloTreeToolBox::getNodeFromName(tree, name, true));
  }

  // Initialize the matrix
  Table<double, std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>> tab(nodenames.size(), nodenames.size(), 0.0);
  tab.setRowIndexes(nodes);
  tab.setColIndexes(nodes);

  // Fill the matrix
  for(std::size_t i = 0; i < nodes.size(); i++){
    for(std::size_t j = 0; j < nodes.size(); j++){
      file >> std::setprecision(5) >> tab(nodes[i], nodes[j]);
    }
  }

  if (verbose) std::cout << "\t...matrix:" << std::endl << tab << std::endl;

  // Close the file and return the matrix
  file.close();

  if(verbose) std::cout << "\t...done." << std::endl;

  return tab;
}

void IO::preProcesses(bpp::PhyloTree &tree
    , const bool check_branch_length /// Check if branch lengths are correct, if there is no branch length it will be replaced with a default value and warning will be printed.
    , const bool force_branch_support /// Add support when there is no support/bootstrap on a branch.
) {
  // First of all: find if a leaf name is a digit or a number, because of each node will be named as the index defined
  // by Bio++...
  int increment = 0;
  for(auto node: tree.getAllLeaves()){
    node->setProperty("isArtificialGene", IsArtificialGene(false));
    if(utils::stringIsNumber(node->getName())) {
      // std::cerr << "The tree can't have an elements named as a number: " << node->getName() << "." << std::endl;
      int num = std::stoi(node->getName());
      if(num > increment) {
        increment = (int) (num + 1);
      }
    }
  }

  // Then, set a name to the internal nodes, the name is going to be an integer which increments as the post order traversal.
  for(auto node: PhyloTreeToolBox::getInternalNodes(tree)){
    node->setProperty("isArtificialGene", IsArtificialGene(false));
    if (!node->hasName()) {
      //auto nodeIndex = tree.getNodeIndex(node);
      //node->setName(std::to_string(nodeIndex + increment));
      node->setName(std::to_string(increment));
      increment++;
    }

    if (force_branch_support and (node != tree.getRoot())) {
      auto edgeToFather = tree.getEdgeToFather(node);
      if (not edgeToFather->hasBootstrapValue()) {
        edgeToFather->setProperty(
            "bootstrap",
            bpp::Number<double>(DEFAULT_BRANCH_SUPPORT_VALUE));
      }
    }
  }

  if(check_branch_length) {
    bool tree_has_edge_without_length = false;
    std::size_t n_edges_without_length = 0;
    double init_edge_length_value = MINIMAL_BRANCH_LENGTH;
    for (auto edge: tree.getAllEdges()) {
      if (not edge->hasLength()) {
        if (not tree_has_edge_without_length) tree_has_edge_without_length = true;
        edge->setLength(init_edge_length_value);
        n_edges_without_length += 1;
      }
    }

    if (tree_has_edge_without_length) {
      std::cerr << "Warning: tree has " << n_edges_without_length << " branch"
                << (n_edges_without_length > 1 ? "es" : "")
                << " without length." << std::endl;
      std::cerr << "Branch length set to " << init_edge_length_value << std::endl;
    }
  }
}

template<typename NodesIterator>
void IO::renameNodeHomonyms(const NodesIterator& nodes_begin,
                            const NodesIterator& nodes_end) {
  std::unordered_map<std::string, std::list<Node>> nodes_by_name;

  // First we classified nodes by their names
  NodesIterator nodes_it;
  for(nodes_it = nodes_begin ; nodes_it != nodes_end ; nodes_it++) {
    auto& node = *nodes_it;
    if(node->hasName()) {
      std::string gene_name = node->getName();
      nodes_by_name[gene_name].push_back(node);
    }
  }

  // Then, if one list has a length > 1, we rename all genes in this list.
  auto nbn_it = nodes_by_name.begin();
  for(; nbn_it != nodes_by_name.end(); nbn_it++) {
    if(nbn_it->second.size() > 1) {
      unsigned int occurrence = 0;
      for (auto node: nbn_it->second) {
        auto node_name = node->getName();
        node_name += ("_" + std::to_string(occurrence));
        node->setName(node_name);
        occurrence++;
      }
    }
  }
}

std::string IO::PhyloTreeToNewick(const bpp::PhyloTree &tree){
  return newick_.treeToParenthesis(tree, false);
}

std::string IO::PhyloTreeToNewick(const bpp::PhyloTree &tree, bool bootstrap, const std::string &propertyName){
  return newick_.treeToParenthesis(tree, bootstrap, propertyName);
}

std::string IO::PhyloTreeToNhx(const bpp::PhyloTree &tree){
  return nhx_.treeToParenthesis(tree);
}

std::string IO::PhyloTreeToPhyloXML(const bpp::PhyloTree &tree) {
  std::ostringstream res;
  phyloXML_.write(tree, res);
  return res.str();
}

void IO::write(
    const bpp::PhyloTree &tree
    , std::ostream &os
    , const TextFormat output
    , const std::string& description
){
  if(not description.empty() and (output != TextFormat::phyloxml and output != TextFormat::recphyloxml)) {
    os << "> " << description << std::endl;
  } else if(not description.empty()) {
    os << "<!-- " << description << " -->" << std::endl;
  }

  if(output == TextFormat::nhx)
    nhx_.write(tree, os);
  else if(output == TextFormat::phyloxml)
    phyloXML_.phylogenyToPhyloXML_(os, tree, "", 1);
  else if(output == TextFormat::recphyloxml)
    phyloXML_.geneTreeToRecPhyloXML(os, tree, -1, 1);
  else
    newick_.write(tree, os, false);
}

std::string IO::ReconciliationToRecPhyloXML(
    const bpp::PhyloTree &genetree
    , const bpp::PhyloTree &speciestree
    , const SpeciesGeneMap &map
) {
  std::ostringstream res;
  res << "<recPhylo>" << std::endl;
  phyloXML_.speciesTreeToPhyloXML(res, speciestree, 1);
  for(auto node: genetree.getAllNodes()) node->setProperty("Species name", bpp::BppString(map.getAssociatedSpecies(node)->getName()));
  phyloXML_.geneTreeToRecPhyloXML(res, genetree, 1, 1);
  res << "</recPhylo>" << std::endl;
  res << std::endl;
  return res.str();
}

void IO::writeReconciledTrees(
    std::ostream &os,
    const std::unordered_map<std::shared_ptr<bpp::PhyloTree>,
    std::map<double, std::vector<ReconciledRootedTree>>>& reconciledTrees,
    const bool& strict_thresholds,
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& trees,
    const bpp::PhyloTree& speciestree,
    const bool printTree,
    const TextFormat output,
    const bool description
) {
  if(output == TextFormat::svg) {
    RecPhyloTreeToSVG::write(os, reconciledTrees, speciestree);
  } else {
    if (output == TextFormat::recphyloxml and printTree) {
      // Used for recphyloxml format.
      os << "<recPhylo>" << std::endl;
      phyloXML_.speciesTreeToPhyloXML(os, speciestree, 1);
    } else if (output == TextFormat::phyloxml and printTree) {
      os << "<phyloxml>" << std::endl;
    }
    std::size_t family_i = 1;

    for (const auto &tree: trees) {
      auto &solutions_for_one_threshold = reconciledTrees.at(
          tree); // gives an std::unordered_map<double, std::vector<ReconciledRootedTree>>
      for (auto &solutions: solutions_for_one_threshold) { // gives std::vector<ReconciledRootedTree> for one Threshold
        std::size_t tree_i = 1;
        auto &contraction_threshold = solutions.first;
        for (auto &solution: solutions.second) {
          std::ostringstream os_description;

          if (description) {
            os_description << "family " << family_i << " tree " << tree_i
                           << (solution.rerooted() ? " re-rooting" : "");
            os_description << " (total cost = " << solution.cost();
            os_description << ", duplications = " << solution.ndup();
            os_description << ", losses = " << solution.nloss();
            if (contraction_threshold >= 0) {
              os_description << ", contraction threshold = "
                             << contraction_threshold;
              if (not strict_thresholds) {
                os_description << "+epsilon";
              }
            }
            else
              os_description << ", contraction threshold = no";
            if (solution.evaluated_with_ale())
              os_description << ", ALE logLk = "
                             << solution.ale_loglikelihood();
            if (solution.evaluated_with_libpll())
              os_description << ", libpll logLk = "
                             << solution.libpll_loglikelihood();
            os_description << ")";
          }

          if (printTree) { //Print tree into Nhx or Newick format
            auto genetree = solution.genetree(); // tree to print
            auto genetree_copy = PhyloTreeToolBox::cloneTree(*genetree);
            auto root = genetree_copy->getRoot(); // root of the tree
            auto map = solution.map(); // map which indicates if two gene nodes
            // belongs to the same species

            if (output == TextFormat::newick
                or output == TextFormat::nhx
                or output == TextFormat::phyloxml) {
              PhyloTreeToolBox::removeArtificialGenes(*genetree_copy);
            } else {
              // First of all eliminate loss duplicates
              auto artificial_nodes =
                  PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(
                      *genetree_copy,
                      [](const Node& node) {
                        return PhyloTreeToolBox::isArtificalGene(node);
                      });

              renameNodeHomonyms(
                  artificial_nodes.begin(), artificial_nodes.end());

              PhyloTreeToolBox::removeBranchLengths(*genetree_copy);
            }

            bool add_tags_in_tree = output == TextFormat::nhx
                                    or output == TextFormat::recphyloxml;

            if (add_tags_in_tree) {
              // If we are in reconciliation mode, each branch needs to be
              // tagged to get duplications in parenthesis format.

              bpp::BppString default_duplication_tag_value("N");

              auto nodes = genetree_copy->getAllNodes();

              for (auto &node: nodes) {
                bpp::BppString duplication_tag_value(
                    default_duplication_tag_value);
                auto sons = genetree_copy->getSons(node);
                if (not PhyloTreeToolBox::isArtificalGene(node) and
                    sons.size() == 2) {
                  auto geneEvents = PhyloTreeToolBox::getGeneEventsInBifurcation(
                      node, *genetree_copy, speciestree, map);
                  if ( geneEvents.at(duplication) > 0 ) {
                    duplication_tag_value = bpp::BppString(
                        "Y"); // means "yes", it is a duplication event.
                  } else {}
                } else {}

                node->setProperty("Species name", bpp::BppString(
                    map.getAssociatedSpecies(node)->getName()));
                assert(node->hasProperty("Species name"));
                node->setProperty(DUPLICATION_STR_FLAG, duplication_tag_value);
              }
            }

            write(*genetree_copy, os, output, os_description.str());

            if (add_tags_in_tree) {
              // Then delete, added properties in nodes.
              auto nodes = genetree_copy->getAllNodes();
              for (auto &node: nodes) {
                assert(node->hasProperty("Species name"));
                node->deleteProperty("Species name");
                assert(node->hasProperty(DUPLICATION_STR_FLAG));
                node->deleteProperty(DUPLICATION_STR_FLAG);
              }
            }
          } else os << "> " << os_description.str() << std::endl;
          tree_i += 1;
        }
      }
      family_i += 1;
    }

    if (output == TextFormat::recphyloxml and printTree)
      os << "</recPhylo>" << std::endl;
    else if (output == TextFormat::phyloxml and printTree)
      os << "</phyloxml>" << std::endl;
  }
}

void IO::writeReconciledTrees(
    const std::string& filename,
    const std::unordered_map<std::shared_ptr<bpp::PhyloTree>,
        std::map<double, std::vector<ReconciledRootedTree>>>& reconciledTrees,
    const bool& strict_thresholds,
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& trees,
    const bpp::PhyloTree& speciestree,
    const TextFormat output,
    const bool description
) {
  std::ofstream genefile(filename.c_str());
  writeReconciledTrees(genefile, reconciledTrees, strict_thresholds,
      trees, speciestree, true, output, description);
  genefile.close();
}


void IO::writeEventSummary(const std::string& filename,
    const std::unordered_map<std::shared_ptr<bpp::PhyloTree>,
        std::map<double, std::vector<ReconciledRootedTree>>> &reconciledTrees,
    const std::vector<std::shared_ptr<bpp::PhyloTree>> &trees,
    const bool full_reconciliation
) {
  std::ofstream file(filename);
  for(auto input_tree : trees) {
    for(auto solutions_per_threshold : reconciledTrees.at(input_tree)) {
      for(auto solution: solutions_per_threshold.second) {
        writeEventSummary(file, solution, full_reconciliation);
      }
    }
  }
  file.close();
}

void IO::writeEventSummary(
    std::ostream &os
    , const ReconciledRootedTree &recgtree
    , bool full_reconciliation
) {
  const bpp::PhyloTree& genetree_topology = *recgtree.genetree();
  const SpeciesGeneMap map = recgtree.map();
  auto nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(*recgtree.genetree());

  os << "DUPLICATIONS : " << recgtree.ndup() << "\tLOSSES : " << recgtree.nloss() << std::endl;

  for(auto node: nodes) {
    if(not genetree_topology.isLeaf(node)) {
      auto sons = genetree_topology.getSons(node);
      assert(sons.size());

      auto left_leaves = genetree_topology.getLeavesUnderNode(sons.front());
      if(not full_reconciliation){
        std::remove_if(left_leaves.begin(), left_leaves.end(), [](const Node& node) {return PhyloTreeToolBox::isArtificalGene(node);});
      }

      auto right_leaves = genetree_topology.getLeavesUnderNode(sons.back());
      if(not full_reconciliation){
        std::remove_if(right_leaves.begin(), right_leaves.end(), [](const Node& node) {return PhyloTreeToolBox::isArtificalGene(node);});
      }

      Event event_type = PhyloTreeToolBox::getEvent(genetree_topology, map, node);

      if (event_type == Event::duplication) {
        os << "    PARALOGY RELATIONSHIP: ";
        utils::write(os, left_leaves.begin(), left_leaves.end(), "", ", ", "");
        os << "    <====>    ";
        utils::write(os, right_leaves.begin(), right_leaves.end(), "", ", ", "");
        os << std::endl;
      } else if(event_type == Event::speciation) {
        os << "    ORTHOLOGY RELATIONSHIP: ";
        utils::write(os, left_leaves.begin(), left_leaves.end(), "", ", ", "");
        os << "    <====>    ";
        utils::write(os, right_leaves.begin(), right_leaves.end(), "", ", ", "");
        os << std::endl;
      } else if(event_type == Event::loss and full_reconciliation) {
        os << "    LOSSES: " << map.getAssociatedSpecies(node);
        os << std::endl;
      }
    }
  }
}

#ifndef _WIN32
void IO::MakeDir(const std::string path, uint32_t mode) {
  if (mkdir(path.c_str(), mode) == -1) {
    err(EXIT_FAILURE, path.c_str(), errno);
  }
}
#else
void IO::MakeDir(const std::string path, uint32_t) {
  if (mkdir(path.c_str()) == -1) {
    cerr << path << ": Permission denied" << endl;
    exit(EXIT_FAILURE);
  }
}
#endif

bool IO::exists(const std::string &filename) {
  return ( access( filename.c_str(), F_OK ) != -1 );
    /*!
  if (FILE *file = fopen(filename.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
     */
}

bool IO::IsRegularFile(const char* path) {
  struct stat buf;
  stat(path, &buf);
  return S_ISREG(buf.st_mode);
}

bool IO::IsDirectory(const char* path) {
  struct stat buf;
  stat(path, &buf);
  return S_ISDIR(buf.st_mode);
}

bool IO::IsDirEmpty(const char* path) {int n = 0;
  struct dirent* d;
  DIR* dir = opendir(path);
  if (dir == NULL) { // Not a directory or doesn't exist
    throw(std::string(path) + ": not a directory");
  }
  while ((d = readdir(dir)) != NULL) {
    if(++n > 2) {
      break;
    }
  }
  closedir(dir);
  return (n <= 2);
}

unsigned long IO::nlines(std::ifstream &is) {
  // Reset the file stream and start the mapping.
  is.clear();
  is.seekg(0, is.beg);

  // Compute the number of lines.
  unsigned long res;
  std::string line_content;
  for (res = 0; std::getline(is, line_content); ++res)
    ;

  // Reset the file stream and start the mapping.
  is.clear();
  is.seekg(0, is.beg);

  return res;
}

/**
 * Write a TextFormat into an output stream.
 *
 * \param os
 * \param output_format
 * \return
 */
std::ostream& operator<<(std::ostream& os,
                                const TextFormat& output_format) {
  switch (output_format) {
    case(TextFormat::nhx): {
      os << "NHX";
      break;
    }
    case(TextFormat::newick): {
      os << "Newick";
      break;
    }
    case(TextFormat::phyloxml): {
      os << "PhyloXML";
      break;
    }
    case(TextFormat::recphyloxml): {
      os << "RecPhyloXML";
      break;
    }
    case(TextFormat::svg): {
      os << "SVG";
      break;
    }
    default:
      os << "Unknown";
  }
  return os;
}

} // namespace treerecs
