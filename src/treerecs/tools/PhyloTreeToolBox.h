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

#ifndef TREERECS_PHYLOTREETOOLBOX_H
#define TREERECS_PHYLOTREETOOLBOX_H

//Includes Bpp
#include <Bpp/Phyl/Tree/PhyloTree.h>

//Includes standards
#include <string>
#include <unordered_set>
#include <stack>
#include <queue>
#include <cmath>

//Includes containers
#include <treerecs/containers/Cost.h>
#include <treerecs/containers/BipartitionList.h>
#include <treerecs/containers/GeneMap.h>
#include <treerecs/containers/Table.h>

//Includes tools
#include "BipartitionTools.h"

namespace treerecs {

//===============================================
/// Typedef Node as a std::shared_ptr<bpp::PhyloNode>.
using Node = std::shared_ptr<bpp::PhyloNode>;
/// Typedef Branch as a std::shared_ptr<std::PhyloBranch>.
using Branch = std::shared_ptr<bpp::PhyloBranch>;
/// Typedef DistanceMatrix as a Table of doubles with Nodes as row
/// and column keys.
using DistanceMatrix = Table<double, Node, Node>;
//===============================================

class PhyloTreeToolBox {
  /*!
   * \class PhyloTreeToolBox
   * \brief PhyloTreeToolBox provides methods to explore and modify a
   *        bpp::PhyloTree as pruning, Post-order traversal edition,
   *        bpp::PhyloTree cloning, etc.
   * \details
   *   Get all nodes according to a post-order traversal:
   *
   *        auto nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(my_PhyloTree);
   *
   *   Get all nodes according to a name pattern:
   *
   *        auto nodes_selected = PhyloTreeToolBox::getNodesFromNamePattern(nodes, "homo_*");
   *
   */

private:
  /************************************
   * Private methods
   */


  /// \brief Used to test good exploration of the visiting matrix in
  ///        computeDistanceMatrix(...) method.
  /// \param m matrix to evaluate
  /// \param check_diag check diagonal (default = true)
  /// \return Return true if matrix contains ones.
  template<typename T, typename Rowkey, typename Colkey>
  static bool isUnary(
      const Table<T, Rowkey, Colkey> &m
      , const bool check_diag = true);

  /// \brief Recursive function for in-order traversal.
  /// \param tree tree to explore.
  /// \param node node which is visited.
  /// \param nodes list of nodes in in-order traversal.
  /// \param conditions condition to push back the node in the resulting list.
  /// \return Nothing, used to fill a list of nodes according to a visiting Node
  ///         node and a Condition condition.
  template<typename Condition>
  static void fillInOrderRecursive(
      const bpp::PhyloTree &tree
      , const Node &node ///
      , std::list<Node> &nodes
      , const Condition &condition
  );

  /// \brief Recursive function for pre-order traversal.
  /// \param tree tree to explore.
  /// \param node node which is visited.
  /// \param nodes list of nodes in pre-order traversal.
  /// \param condition condition to push back the node in the resulting list.
  /// \return Nothing, used to fill a list of nodes according to a visiting Node
  /// node and a Condition condition.
  template<typename Condition>
  static void fillInPreOrderRecursive(
      const bpp::PhyloTree &tree
      , const Node &node
      , std::list<Node> &nodes
      , const Condition &condition
  );

  /// \brief Recursive function for post-order traversal.
  /// \param tree tree to explore.
  /// \param node node which is visited.
  /// \param nodes list of nodes in pre-order traversal.
  /// \param condition condition to push back the node in the resulting list.
  /// \return Nothing, used to fill a list of nodes according to a visiting Node
  ///         node and a Condition condition.
  template<typename Condition>
  static void fillInPostOrderRecursive(
      const bpp::PhyloTree &tree
      , const Node &node
      , std::list<Node> &nodes
      , const Condition &condition
  );

  /// \brief Recursive tree cloning. Clone only tree topology: clone and
  ///        original tree shared the same nodes.
  /// \param tree Tree to clone.
  static void recursiveTreeCloning(
      const bpp::PhyloTree &tree
      , const std::shared_ptr<bpp::PhyloNode> &current_node
      , const std::shared_ptr<bpp::PhyloNode> &father_node
      , bpp::PhyloTree &newtree
      , const bool keep_branch_lengths
      , std::unordered_set<std::shared_ptr<bpp::PhyloNode>> &created
      , const bool verbose
  );

  /// \brief Returns a list of node (in Post-order traversal) according to a
  ///        "push_back" condition and an other to stop exploration.
  ///        Iterative version.
  /// \param tree PhyloTree to visit.
  /// \param get_condition Getter condition. Takes a Node in parameter and
  ///        returns a boolean. If the condition returns true, the node is
  ///        retained.
  /// \param stop_condition Stop condition. Takes a Node in parameter and
  ///        returns a boolean. If the condition return true, the traversal is
  ///        stopped.
  /// \param root Node to start the traversal, default is the root of the tree.
  /// \return A list of nodes according to a post-order traversal, a condition
  /// and a subtree root.
  template<typename Get_condition, typename Stop_condition>
  static std::list<Node> PostOrderTraversalIterative_list(
      const bpp::PhyloTree &tree
      , const Get_condition &get_condition
      , const Stop_condition &stop_condition
      , const Node &root = nullptr);

protected:

public:
  //=========================================================
  // SETTERS, transform a tree, modify, remove and add nodes.
  //=========================================================

  /// \brief Apply a function (from lambda or functor) to nodes in a given tree
  ///        according to a post-order traversal.
  /// \return Nothing but modify nodes in tree.
  template<typename Func>
  static void applyInPOT(
      bpp::PhyloTree &tree, const Func &fun, const Node &node = nullptr);

  /// \brief Add father node and its sons in a tree. Be careful nodes needs to
  ///        be created before.
  /// \param tree Tree to modify.
  /// \param father Father node.
  /// \param son_left Son left of the father.
  /// \param son_right Son right of the father.
  /// \param son_left_branch Branch between father and son_left.
  /// \param son_right_branch Branch between fahter and son_right.
  /// \param verbose Verbosity on std::cout
  static void addNodeWithSonsInPostOrder(
      bpp::PhyloTree &tree
      , const std::shared_ptr<bpp::PhyloNode> &father
      , const std::shared_ptr<bpp::PhyloNode> &son_left
      , const std::shared_ptr<bpp::PhyloNode> &son_right
      , const std::shared_ptr<bpp::PhyloBranch> &son_left_branch = nullptr
      , const std::shared_ptr<bpp::PhyloBranch> &son_right_branch = nullptr
      , const bool verbose = false
  );

  /// \brief Contract tree branches according to a maximal threshold of branch
  ///        support. Branch with a support below the
  ///        threshold are contracted.
  /// \param tree tree to contract
  /// \param threshold minimum value of branch support to keep.
  /// \param strict_thresholds do not contract branches if branch support is
  ///        equal to the threshold.
  /// \return Returns a list of nodes which are removed from the tree.
  static std::list<Node> mergeWeakBranches(
      bpp::PhyloTree &tree
      , const double threshold
      , const bool strict_thresholds = DEFAULT_STRICT_SUPPORT_THRESHOLDS
  );

  /// \brief Prune tree of artificial genes.
  static void removeArtificialGenes(bpp::PhyloTree &tree);

  /// \brief Removes nodes of a tree matching with a specific given condition.
  /// \param tree Tree to modify.
  /// \param root Root of the subtree.
  /// \param condition Condition to remove a node.
  /// \param verbose Print verbosity in std::cout.
  template<typename Func>
  static void removeNodesInPostOrderTraversal(
      bpp::PhyloTree &tree
      , const Node &root
      , const Func &condition
      , const bool verbose = false);

  /// \brief Removes a node from a given tree.
  /// \param tree Tree to modify.
  /// \param node Node to remove.
  /// \return Nothing but modify the tree by removing a node.
  static void removeNode(
      bpp::PhyloTree &tree
      , const Node &node
      , const bool verbose = false);

  /// \brief Reset Node and edges ids in a bpp::Phylotree.
  /// \return Nothing but changes Node Ids in tree.
  static void resetNodeIdInPostOrder(bpp::PhyloTree &tree);

  /// \brief Remove branch length information.
  /// \return Nothing but changes tree.
  static void removeBranchLengths(bpp::PhyloTree &tree);

  // Simplifying trees, pruning
  /// \brief Prune a species tree of leaves which are not specified in the
  ///        SpeciesGeneMap.
  /// \return Nothing but edit the tree.
  static void pruneTree_v0(
      bpp::PhyloTree &speciestree, const SpeciesGeneMap &genemap
      , const bool verbose = false);

  //=========================================================
  // GETTERS
  // Get or find nodes according to conditions or order, get node properties,
  // generate distance matrix, clone tree, etc.
  //=========================================================

  /// \brief Get Nodes (according to a particular given condition) in order,
  ///        returned in a list. Recursive version.
  /// \param tree Tree to retrieve nodes.
  /// \param condition Condition to append node in the resulting list.
  /// \param root Root of the   subtree, default to the root of the entire tree.
  /// \return List of nodes according to an in-order traversal.
  template<typename Condition>
  static std::list<Node> InOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree
      , const Condition &condition
      , const Node &root = nullptr);

  /// \brief Get Nodes (according to a particular given condition) in pre order,
  ///        returned in a list. Recursive version.
  /// \param tree Tree to retrieve nodes.
  /// \param condition Condition to append node in the resulting list.
  /// \param root Root of the subtree, default to the root of the entire tree.
  /// \return List of nodes according to a pre-order traversal.
  template<typename Condition>
  static std::list<Node>
  PreOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree
      , const Condition &condition
      , const Node &root = nullptr);

  /// \brief Get Nodes (according to a particular given condition) in post order
  ///        , returned in a list. Recursive version.
  /// \param tree Tree to retrieve nodes.
  /// \param condition Condition to append node in the resulting list.
  /// \param root Root of the subtree, default to the root of the entire tree.
  /// \return List of nodes according to a post-order traversal.
  template<typename Condition>
  static std::list<Node>
  PostOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree
      , const Condition &condition
      , const Node &root = nullptr);

  /// \brief Get nodes in order traversal. Recursive version.
  /// \param tree Tree to explore.
  /// \param root Root of the subtree, default to the root of the entire tree.
  /// \return List of all nodes of a subtree according to an in-order traversal.
  static std::list<Node>
  getNodesInOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree
      , const Node &root = nullptr);

  /// \brief Get nodes in a pre-order traversal. Recursive version.
  /// \param tree Tree to explore.
  /// \param root Root of the subtree, default to the root of the entire tree.
  /// \return List of all nodes of a subtree according to a pre-order traversal.
  static std::list<Node>
  getNodesInPreOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree
      , const Node &root = nullptr);

  /// \brief Get nodes in a post-order traversal. Recursive version.
  /// \param tree Tree to explore.
  /// \param root Root of the subtree, default to the root of the entire tree.
  /// \return List of all nodes of a subtree according to a post-order traversal
  static std::list<Node>
  getNodesInPostOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree
      , const Node &root = nullptr);

  /// \brief Check if almost one node of the tree returns true according to a
  ///        particular condition.
  /// \param tree Tree to check.
  /// \param condition Lambda condition to check, like
  ///        [](const std::shared_ptr<bpp::PhyloNode>& node){...}.
  /// \param root Root of the subtree to explore, default is the root of the
  ///        entire tree.
  template<typename Condition>
  static bool hasNode(
      const bpp::PhyloTree &tree
      , const Condition &condition
      , const Node &root = nullptr);

  /// \brief Get nodes in a post-order traversal. Iterative version
  /// \param tree Tree to explore.
  /// \param root Root of the subtree, default to the root of the entire tree.
  /// \return List of all nodes of a subtree according to a post-order traversal
  static std::list<Node>
  getNodesInPostOrderTraversalIterative_list(
      const bpp::PhyloTree &tree
      , const Node &root = nullptr);

  /// \brief Get nodes in a post-order traversal according to a specific
  ///        condition.
  /// \param tree Tree to check.
  /// \param condition Condition to append node in the resulting list.
  /// \param root Root of the subtree to explore, default is the root of the
  ///        entire tree.
  /// \return List of all nodes of a subtree according to a post-order traversal
  ///         and which in respect of a given condition.
  template<typename Condition>
  static std::list<Node>
  getNodesInPostOrderTraversalIterative_list(
      const bpp::PhyloTree &tree
      , const Condition &condition
      , const Node &root = nullptr);

  /// \brief Get names from a container of nodes using two iterators.
  /// \return Vector of names of all nodes in a container.
  template<typename Container_iterator>
  static std::vector<std::string>
  getNodeNames(Container_iterator begin, const Container_iterator &end);

  /// \brief Returns all internal nodes under a given node (default = tree's
  ///        root).
  /// \return List of all internal nodes under a node.
  static std::list<Node>
  getInternalNodes(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Returns all internal node names under a give node (default = tree's
  ///        root).
  /// \return Vector of names.
  static std::vector<std::string>
  getInternalNodeNames(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Returns all leaves under a node (default = tree's root).
  /// \return List of std::shared_ptr<bpp::PhyloNode>.
  static std::list<Node>
  getLeaves(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Get all leaves names under a node (default = tree's root).
  /// \return Vector of std::string.
  static std::vector<std::string>
  getLeavesNames(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Get the distance between two nodes.
  static double getDistanceBetweenTwoNodes(
      const Node &nodeA, const Node &nodeB, const bpp::PhyloTree &tree);

  /// \brief Check if there is a distance between two nodes.
  static bool hasDistanceBetweenTwoNodes(
      const Node &nodeA, const Node &nodeB, const bpp::PhyloTree &tree);

  /// \brief Sort Nodes vector by their index in the bpp::PhyloTree.
  ///        Useful when nodes are associated with an id according to a
  ///        post-order traversal.
  static void sortNodesByIndex(
      std::vector<std::shared_ptr<bpp::PhyloNode>> &nodes
      , const bpp::PhyloTree &tree, const bool ascending = true);

  /// \brief Sort nodes according to their name.
  static void sortNodesByNameLength(
      std::vector<std::shared_ptr<bpp::PhyloNode>> &nodes
      , const bool ascending = true);

  /// \brief Create a distance matrix of a tree.
  static DistanceMatrix computeDistanceMatrix(
      const bpp::PhyloTree &tree, const Node &root = nullptr
      , bool printProgressionBar = false);

  /// \brief Generate a DistanceMatrix of a given tree. The matrix contains
  ///        leaves only by default.
  static DistanceMatrix distanceMatrix(
      const bpp::PhyloTree &tree, const bool leavesOnly = true
      , const bool printProgressionBar = false);

  /// \brief Get a DistanceMatrix of leaves given a specific tree.
  static DistanceMatrix distanceMatrixBetweenNodes_unoptimized(
      const bpp::PhyloTree &tree, const std::vector<Node> nodes
      , const bool printProgressionBar = false);

  /// \brief Returns the node which has the same name (first occurence).
  static std::shared_ptr<bpp::PhyloNode>
  getNodeFromName(
      const std::vector<Node> nodes, const std::string &name
      , const bool case_sensitive = false);

  /// \brief Returns the node which has the same name (first occurence).
  static std::shared_ptr<bpp::PhyloNode>
  getNodeFromName(
      const bpp::PhyloTree &tree, const std::string &name
      , const bool case_sensitive = false, const bool leaves_only = true);

  /// \brief Check if a node name has a match with a pattern (supports only "*"
  ///        pattern).
  static bool node_name_match_with_pattern(
      const std::string &node_name, const std::string &pattern
      , const bool case_sensitive = false);

  /// \brief Returns indexes of nodes (in a std::list) with a name following the
  ///        given name pattern.
  static std::list<std::size_t>
  getNodeIndexesFromNamePattern(
      const std::vector<Node> &nodes, const std::string &pattern
      , const bool case_sensitive = false);

  /// \brief Returns indexes of nodes with a name following the given name
  ///        pattern.
  static std::vector<std::shared_ptr<bpp::PhyloNode>>
  getNodesFromNamePattern(
      const std::vector<Node> &nodes, const std::string &pattern
      , const bool case_sensitive = false);

  /// \brief Returns the last common ancestor (lca) of two nodes given nodes.
  static Node getCommonAncestor(
      const bpp::PhyloTree &tree, const Node &nodeA, const Node &nodeB);

  /// \brief Returns the last common ancestor (lca) of a set of nodes.
  static Node
  getLastCommonAncestor(std::vector<Node> nodes, const bpp::PhyloTree &tree);

  /// \brief Get brothers of a given node.
  static std::vector<Node>
  getNodeBrothers(const bpp::PhyloTree &tree, const Node &node);

  /// \brief Clone a tree topology without cloning nodes. Clone and original
  ///        tree share same bpp::PhyloNode.
  /// \param tree Original tree to clone.
  /// \param root Clone a subtree, according to a specific root (default is the
  ///        root of the tree).
  /// \param keep_branch_lengths Keep branch lengths in resulting tree.
  /// \param verbose Verbose.
  static std::shared_ptr<bpp::PhyloTree> cloneTree(
      const bpp::PhyloTree &tree
      , const Node &root = nullptr
      , const bool keep_branch_lengths = true
      , const bool verbose = false);

  /// \brief Get branch support frequencies for a set of trees. The returned map
  /// contains in key: the branch support and the frequency as value.
  static std::map<double, double> getSupportFrequencies(
      const std::vector<std::shared_ptr<bpp::PhyloTree>> &trees, size_t nclasses
  );

  /// \brief Get branch support frequencies. The returned map contains in key:
  ///        the branch support and the frequency as value.
  static std::map<double, double>
  getSupportFrequencies(const bpp::PhyloTree &tree, std::size_t nclasses);

  /// \brief Indicates if the given tree is binary (true) or not (false).
  static bool isBinary(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Compare two tree topologies. Returns true if trees are similar, and
  ///        false if there is a difference.
  ///        Branch lengths are not compared.
  static bool compare_trees_with_days_algorithm(
      const bpp::PhyloTree &treeA, const bpp::PhyloTree &treeB
      , const bool verbose = false);

  /// \brief Returns true if the tree contains at least one artificial gene.
  static bool
  hasArtificialGenes(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Check if a node is an artificial gene. This property is defined in
  ///        bpp::PhyloNode properties.
  static bool isArtificalGene(const Node &node);

  /// \brief Get Node occurrence. Occurrence is a particular integer used in
  ///        Polytomysolver algorithm.
  static std::size_t getNodeOccurence(const Node &node);

  /// \brief Returns a list of nodes matching with a specified outdegree
  ///        (number of sons).
  static std::list<Node> find_nodes_according_to_its_outdegree(
      const bpp::PhyloTree &tree, const size_t &outdegree
      , const Node &root = nullptr);

  /// \brief Returns a list of multifurcation roots.
  ///        Note: polytomy = multifurcation.
  /// \note We consider a polytomy (= multifurcation) as a node which has 3 or
  ///       more sons.
  static std::list<Node>
  findPolytomies(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Returns a list of bifurcation roots.
  /// \note We consider a bifurcation as a node which has exactly two sons.
  static std::list<Node>
  findBifurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Returns a list of monofurcation roots.
  ///        A monofurcation is an internal node with only one son.
  /// \note We consider a monofurcation as a node which has only one son.
  static std::list<Node>
  findMonofurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Check if a tree contains at least one monofurcation.
  ///        A monofurcation is an internal node with only one son.
  static bool
  hasMonofurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Sort internal nodes according to their outdegree. The first element
  ///        of the pair is a list of bifurcations and the
  ///        second a list of multifurcations (= polytomy).
  /// \note We consider bifurcation as a node which contains exactly two sons
  ///       and a multifurcation (= polytomy) as a
  ///       node which has more than two sons.
  static std::pair<std::list<Node>, std::list<Node>>
  getFurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Deduce events in a given bifurcation according to the species tree,
  ///        and a map.
  /// \param node Gene node of interest
  /// \param genetree Gene tree
  /// \param speciestree Species tree
  /// \param map Map
  /// \return A map of Events (duplication, loss) with their number.
  static std::map<Event, std::size_t> getGeneEventsInBifurcation(
      const Node &node
      , const bpp::PhyloTree &genetree
      , const bpp::PhyloTree &speciestree
      , const SpeciesGeneMap &map);

  /// Check if each node of a tree is connected to a father.
  static bool allEdgesAreCorrect(const bpp::PhyloTree &tree);

  /// Deduce event in a gene tree node. GeneTree needs to be reconciled.
  static Event getEvent(
      const bpp::PhyloTree &genetree, const SpeciesGeneMap &map
      , const Node &node);

/**
 * @brief Counts the total number of occurrences of every bipartition from the
 *      input trees
 * @warning Not tested.
 *
 * Returns the list of distinct bipartitions found at least once in the set of
 * input trees,
 * and writes the number of occurrence of each of these bipartitions in vector
 * bipScore.
 *
 * @author Nicolas Galtier, modified by Nicolas Comte
 * @param vecTr Vector of input trees (must share a common set of leaves - not
 *        checked in this function)
 * @param bipScore Output as the numbers of occurrences of the returned distinct
 *        bipartitions
 * @return A BipartitionList object including only distinct bipartitions
 * @note from Bpp-phyl library.
 */
  static BipartitionList *bipartitionOccurrences(
      const std::vector<std::shared_ptr<bpp::PhyloTree>> &vecTr
      , std::vector<std::size_t> &bipScore);

  /**
   * @brief Compute bootstrap values.
   * @warning Not tested.
   *
   * @param bpp::PhyloTree Input tree. the BOOTSTRAP banch property of the tree
   *        will be modified if it already exists.
   * @param vecTr A list of trees to compare to 'tree'.
   * @param verbose Tell if a progress bar should be displayed.
   * @param format The output format of the tree.
   * @note from Bpp-phyl library.
   */
  static void computeBootstrapValues(
      bpp::PhyloTree &tree
      , const std::vector<std::shared_ptr<bpp::PhyloTree>> &vecTr
      , bool verbose = true, int format = 0);

  /**
   * @brief General greedy consensus tree method.
   * @warning Not tested.
   *
   * Calculates the consensus tree of a set of trees defined from the number of
   * occurrences of bipartitions. Bipartitions are considered in decreasing
   * score order.
   * A bipartition is included if it is compatible with all previously included
   * bipartitions, and if its score is higher than a threshold.
   *
   * @author Nicolas Galtier, modified by Nicolas Comte for Treerecs.
   * @param vecTr Vector of input trees (must share a common set of leaves -
   *        checked if checkNames is true)
   * @param threshold Minimal acceptable score =number of occurrence of a
   *        bipartition/number of trees (0.<=threshold<=1.)
   * @param checkNames Tell whether we should check the trees first.
   * @note from Bpp-phyl library.
   */
  static std::shared_ptr<bpp::PhyloTree> thresholdConsensus(
      const std::vector<std::shared_ptr<bpp::PhyloTree>> &vecTr
      , double threshold, bool checkNames = true
  ) noexcept(false);
};

template<typename T, typename Rowkey, typename Colkey>
bool PhyloTreeToolBox::isUnary(
    const Table<T, Rowkey, Colkey> &m
    , const bool check_diag
) {
  auto rowkeys = m.getRowIndexes();
  auto colkeys = m.getColIndexes();

  for (auto &rowkey : rowkeys) {
    for (auto &colkey : colkeys) {
      if (rowkey != colkey or check_diag) {
        if (m(rowkey, colkey) != 1) return false;
      }
    }
  }
  return true;
}

template<typename Condition>
void PhyloTreeToolBox::fillInOrderRecursive(
    const bpp::PhyloTree &tree
    , const Node &node
    , std::list<Node> &nodes
    , const Condition &condition
) {
  bool isLeaf = tree.isLeaf(node);

  std::vector<Node> sons;

  if (not isLeaf) {
    sons = tree.getSons(node);
    fillInOrderRecursive(tree, sons.front(), nodes, condition);
  }

  if (condition(node))
    nodes.push_back(node);

  if (sons.size() > 1) {
    for (std::size_t i = 1; i < sons.size(); ++i)
      fillInOrderRecursive(tree, sons.at(i), nodes, condition);
  }
}

template<typename Condition>
void PhyloTreeToolBox::fillInPreOrderRecursive(
    const bpp::PhyloTree &tree
    , const Node &node
    , std::list<Node> &nodes
    , const Condition &condition
) {
  if (condition(node))
    nodes.push_back(node);

  if (not tree.isLeaf(node))
    for (auto &son : tree.getSons(node))
      fillInPreOrderRecursive(tree, son, nodes, condition);
}

template<typename Condition>
void PhyloTreeToolBox::fillInPostOrderRecursive(
    const bpp::PhyloTree &tree
    , const Node &node
    , std::list<Node> &nodes
    , const Condition &condition
) {
  if (not tree.isLeaf(node))
    for (auto &son : tree.getSons(node))
      fillInPostOrderRecursive(tree, son, nodes, condition);

  if (condition(node))
    nodes.push_back(node);
}

template<typename Get_condition, typename Stop_condition>
std::list<Node> PhyloTreeToolBox::PostOrderTraversalIterative_list(
    const bpp::PhyloTree &tree
    , const Get_condition &get_condition
    , const Stop_condition &stop_condition
    , const Node &root
) {
  std::stack<Node> s1;
  std::list<Node> s2;  // used as a second stack which contains nodes in the
  // right order.
  Node current = (root) ? root : tree.getRoot();
  s1.push(current);
  while (not s1.empty()) {
    current = s1.top();
    s1.pop();
    if (get_condition(current))
      s2.push_front(current);

    if (stop_condition(current))
      return s2;

    if (not tree.isLeaf(current))
      for (auto son : tree.getSons(current))
        s1.push(son);
  }

  return s2;
}

template<typename Func>
void PhyloTreeToolBox::applyInPOT(
    bpp::PhyloTree &tree, const Func &fun, const Node &node
) {
  Node root = node ? node : tree.getRoot();
  auto nodes = getNodesInPostOrderTraversalIterative_list(tree, root);
  std::for_each(nodes.begin(), nodes.end(), fun);
}

template<typename Func>
void PhyloTreeToolBox::removeNodesInPostOrderTraversal(
    bpp::PhyloTree& tree,
    const Node&,
    const Func& condition,
    const bool verbose
) {
  PhyloTreeToolBox::applyInPOT(tree, [&tree, &condition, &verbose](
      const Node &node
  ) {
    if (condition(node, tree))
      PhyloTreeToolBox::removeNode(tree, node, verbose);
  });
}

template<typename Condition>
std::list<Node> PhyloTreeToolBox::InOrderTraversalRecursive_list(
    const bpp::PhyloTree &tree
    , const Condition &condition
    , const Node &root
) {
  std::list<Node> nodes;
  fillInOrderRecursive(tree, (root ? root : tree.getRoot()), nodes,
                       condition);
  return nodes;
}

template<typename Condition>
std::list<Node>
PhyloTreeToolBox::PreOrderTraversalRecursive_list(
    const bpp::PhyloTree &tree
    , const Condition &condition
    , const Node &root
) {
  std::list<Node> nodes;
  fillInPreOrderRecursive(tree, root ? root : tree.getRoot(), nodes,
                          condition);
  return nodes;
}

template<typename Condition>
std::list<Node>
PhyloTreeToolBox::PostOrderTraversalRecursive_list(
    const bpp::PhyloTree &tree
    , const Condition &condition
    , const Node &root
) {
  std::list<Node> nodes;
  fillInPostOrderRecursive(tree, root ? root : tree.getRoot(), nodes,
                           condition);
  return nodes;
}

template<typename Condition>
bool PhyloTreeToolBox::hasNode(
    const bpp::PhyloTree &tree
    , const Condition &condition
    , const Node &root
) {
  return PostOrderTraversalIterative_list(tree, condition, condition,
                                          root).size() > 0;
}

template<typename Condition>
std::list<Node>
PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(
    const bpp::PhyloTree& tree,
    const Condition& condition,
    const Node& root
) {
  return PostOrderTraversalIterative_list(
      tree,
      condition,
      [](const Node&) { return false; },
      root);
}

template<typename Container_iterator>
std::vector<std::string>
PhyloTreeToolBox::getNodeNames(
    Container_iterator begin, const Container_iterator &end) {
  std::vector<std::string> node_names(
      (unsigned long) std::distance(begin, end));

  std::generate(node_names.begin(), node_names.end(), [&begin]() {
    std::string name = "";

    if ((*begin)->hasName())
      name = (*begin)->getName();

    begin++;
    return name;
  });

  return node_names;
}

} // namespace treerecs

#endif //TREERECS_PHYLOTREETOOLBOX_H
