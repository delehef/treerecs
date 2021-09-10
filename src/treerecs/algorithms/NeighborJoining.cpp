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

#include "NeighborJoining.h"

#include <treerecs/tools/PhyloTreeToolBox.h>

namespace treerecs {

void NeighborJoining::updateBifurcation(
    bpp::PhyloTree &tree, const Node &father, const Node &son_left
    , const Node &son_right, DistanceMatrix &dmatrix
) const {
  auto branch_left = tree.getEdgeLinking(father, son_left);
  branch_left->setLength(DistanceFromPairAncestor(dmatrix, son_left, son_right));

  auto branch_right = tree.getEdgeLinking(father, son_right);
  branch_right->setLength(dmatrix(son_left, son_right) - branch_left->getLength());

  updateDistances(dmatrix, father, son_left, son_right);
}

std::pair<double, std::size_t> NeighborJoining::sum_of_distances(
    const DistanceMatrix &dmatrix, const Node &node
) const {
  /// Make a sum of a distance matrix column without artificial genes.
  /// The number of artificial genes in saved as a second value of the returned pair.
  double res = 0.0;
  std::size_t n_artificial_genes = 0;
  for(auto& n: dmatrix.getRowIndexes()){
    if(not PhyloTreeToolBox::isArtificalGene(n)){
      res += dmatrix.get(node, n);
    } else {
      n_artificial_genes++;
    }
  }
  return std::make_pair(res, n_artificial_genes);
}

double NeighborJoining::DistanceBetweenPairAndOuterNodes(
    const DistanceMatrix &D, const Node &x, const Node &y, const Node &t
) const {
  /// Compute new distances for a pair to other nodes.
  if(PhyloTreeToolBox::isArtificalGene(x)){
    return D(y, t);
  } else if(PhyloTreeToolBox::isArtificalGene(y)){
    return D(x, t);
  }
  return 0.5*(D(x, t) + D(y, t) - D(x, y));
}

double NeighborJoining::DistanceFromPairAncestor(
    const DistanceMatrix &D, const Node &x, const Node &y
) const {
  /// Compute distance for a node x inside the pair.
  double sum_d_gene_left;
  double sum_d_gene_right;
  sum_d_gene_left = utils::sum(D.getRow(x));
  sum_d_gene_right = utils::sum(D.getRow(y));
  return (D(x, y) / 2.0) + (sum_d_gene_left - sum_d_gene_right) / (2.0*((double) D.nrow() - 2.0));
}

std::pair<Node, Node> NeighborJoining::findPairMinimizingQ(
    const DistanceMatrix &D, const std::vector<Node> &sl_genes
    , const std::vector<Node> &sr_genes, const bool verbose
) const {
  /// Returns the first pair sl_genes x sr_genes minimizing Q_(...).

  //Init
  std::vector<std::pair<Node, Node>> minimal_pairs;
  minimal_pairs.emplace_back(std::make_pair(nullptr, nullptr));
  double qmin = std::numeric_limits<double>::infinity();

  Table<double, Node, Node> Q;
  if(verbose) { // The verbose mode create a matrix which prints all results of the minimum qvalue search.
    Q = Table<double, Node, Node>(sl_genes.size(), sr_genes.size(), 0.0);
    Q.setRowIndexes(sl_genes);
    Q.setColIndexes(sr_genes);
  }

  for(auto gl: sl_genes){
    //auto d_gl_sum = sum_of_distances(D, gl);
    auto d_gl_sum = 0.0;
    if(not PhyloTreeToolBox::isArtificalGene(gl)) {
      d_gl_sum = utils::sum(D.getRow(gl));
    }
    for(auto gr: sr_genes){
      if(gl != gr) {
        double qvalue;
        if(PhyloTreeToolBox::isArtificalGene(gl) or PhyloTreeToolBox::isArtificalGene(gr))
          qvalue = std::numeric_limits<double>::max();
        else {
          auto d_gr_sum = utils::sum(D.getRow(gr));
          qvalue = (D.ncol() - 2.0) * D(gl, gr) - (d_gl_sum + d_gr_sum);
        }

        if(verbose) Q(gl, gr) = qvalue;

        if (qvalue < qmin) {
          qmin = qvalue;
          minimal_pairs.clear();
          minimal_pairs = { std::pair<Node, Node>(gl, gr) };
        } else if (qvalue == qmin) {
          minimal_pairs.emplace_back(std::pair<Node, Node>(gl, gr));
        }
      }
    }
  }

  if(minimal_pairs.size() > 1) {
    std::random_shuffle(minimal_pairs.begin(), minimal_pairs.end());
  }

  if(verbose) std::cout << "Q matrix:" << std::endl << Q << std::endl;
  return minimal_pairs.front();
}

void NeighborJoining::updateDistances(
    DistanceMatrix &D, const Node &father, const Node &son_left
    , const Node &son_right
) const {
  /// Add the new node father and remove its sons.
  std::vector<double> line_init(D.nrow(), 0.0);
  D.addCol(line_init, father);
  line_init.resize(D.ncol(), 0.0);
  D.addRow(line_init, father);
  for(auto& node: D.getColIndexes()){
    if(node != father and node != son_left and node != son_right) {
      D(father, node) = DistanceBetweenPairAndOuterNodes(D, son_left, son_right, node);
      D(node, father) = D(father, node); //Because of the symetric property of D.
    }
  }
  D.removeRow(son_left);
  D.removeRow(son_right);
  D.removeCol(son_left);
  D.removeCol(son_right);
}

void NeighborJoining::operator()(bpp::PhyloTree &tree, DistanceMatrix dmatrix) const {
  /// Solve a the given tree in parameter. If the tree is bifurcated, all branches are recomputed.
  // Check polytomies
  auto polytomies = PhyloTreeToolBox::findPolytomies(tree);

  if(polytomies.size() == 0){

    // Get all inner nodes of the tree (or subtree if specified)
    auto inner_nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(
        tree, [&tree](const Node& node) { return not tree.isLeaf(node); } );

    // Then, recompute all distances with Neighbor-Joining
    for(auto& inode: inner_nodes){
      auto inode_sons = tree.getSons(inode);
      updateBifurcation(tree, inode, inode_sons.front(), inode_sons.back(), dmatrix);
    }

  } else {
    auto leaves_list = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(
        tree, [&tree](const Node& node) { return tree.isLeaf(node); } );

    std::vector<Node> leaves {leaves_list.begin(), leaves_list.end()};

    auto polytomy_distances = dmatrix.extractRows(leaves);
    polytomy_distances = polytomy_distances.extractCols(leaves);

    while(polytomy_distances.nrow() > 2){
      leaves = polytomy_distances.getColIndexes();

      std::pair<Node, Node> pair = findPairMinimizingQ(polytomy_distances, leaves, leaves, false);

      Node son_left = pair.first;
      Node son_right = pair.second;
      Node father(new bpp::PhyloNode());

      tree.createNode(father);

      PhyloTreeToolBox::addNodeWithSonsInPostOrder(
          tree
          , father, son_left, son_right,
          std::shared_ptr<bpp::PhyloBranch>(new bpp::PhyloBranch()),
          std::shared_ptr<bpp::PhyloBranch>(new bpp::PhyloBranch()), false);

      updateBifurcation(tree, father, son_left, son_right, polytomy_distances);

      //utils::printNodeContent(tree, father, std::cout);
    }

    // Last lap: connect the last two nodes.
    leaves = polytomy_distances.getColIndexes();
    Node& nodeA = leaves.front();
    Node& nodeB = leaves.back();

    Node old_root = tree.getRoot();
    tree.link(old_root, nodeA);
    tree.link(old_root, nodeB);

    tree.rootAt(nodeA);

    tree.deleteNode(old_root);

    tree.link(nodeA, nodeB, std::shared_ptr<bpp::PhyloBranch>(new bpp::PhyloBranch(polytomy_distances(nodeA, nodeB))));

  }
}

} // namespace treerecs
