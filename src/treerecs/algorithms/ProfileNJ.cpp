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

#include "ProfileNJ.h"

//Include libs
#include <limits>
#include <cmath>

//Include tools
#include <treerecs/tools/utils.h>
#include <treerecs/tools/SpeciesGeneMapper.h>
#include <treerecs/algorithms/Polytomysolver.h>

namespace treerecs {

void
ProfileNJ::updateDistances(
    DistanceMatrix &D
    , const Node &father
    , const Node &son_left
    , const Node &son_right)
const {
  /// Add the new node father and remove its sons.
  std::vector<double> line_init(D.nrow(), 0.0);
  //if(D.hasRowIndex(father)) std::cout << "Overwriting in the column for the same node." << std::endl;
  D.addCol(line_init, father);
  line_init.resize(D.ncol(), 0.0);
  D.addRow(line_init, father);
  for (auto &node: D.getColIndexes()) {
    if (node != father and node != son_left and node != son_right) {
      D(father, node) = nj_.DistanceBetweenPairAndOuterNodes(D, son_left, son_right, node);
      D(node, father) = D(father, node); //Because of the symetric property of D.
    }
  }
  if(not PhyloTreeToolBox::isArtificalGene(son_left)) {
    D.removeRow(son_left);
    D.removeCol(son_left);
  }
  if(not PhyloTreeToolBox::isArtificalGene(son_right)) {
    D.removeRow(son_right);
    D.removeCol(son_right);
  }

}

Node ProfileNJ::updateTree_GeneMap_DistanceMatrix(
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
    , const bool computeBranchLengths
    , const bool verbose)
const {
  /// Update the new tree with NJ, distance matrix and SpeciesGeneMap.
  // Initialize the new node father, gene, common ancestor of gene_left and gene_right
  Node gene_father = nullptr;
  if (genemap_temp.getGenes().size() > 2) { // If we are in the end of the polytomy reconciliation (dmatrix.ncol() = 2), set the new-father as the existing gene root polytomy
    gene_father = Node(new bpp::PhyloNode);

    std::string father_name = species->getName() + std::to_string(genemap.ngenes(species) + 1);

    if(species_left == species_right){
      father_name += "_";
      father_name += DUPLICATION_STR_FLAG;
    }
    else {
      father_name += "_";
      father_name += SPECIATION_STR_FLAG;
    }
    gene_father->setName(father_name);
    gene_father->setProperty("isArtificialGene", IsArtificialGene(false));

    tree.createNode(gene_father);
    genemap.addGene(species, gene_father);

    // Change father name according to the event.
  } else {
    gene_father = polytomy_root;
  }

  // Compute distance between father node gene and its sons gene_left and gene_right
  std::shared_ptr<bpp::PhyloBranch> new_branch_left = std::shared_ptr<bpp::PhyloBranch>(new bpp::PhyloBranch());
  std::shared_ptr<bpp::PhyloBranch> new_branch_right = std::shared_ptr<bpp::PhyloBranch>(new bpp::PhyloBranch());

  if (computeBranchLengths) {
    double d_gene_geneleft; // distance from father to the son left.
    double d_gene_generight; // distance from father to the son right.

    if (PhyloTreeToolBox::isArtificalGene(gene_left) or PhyloTreeToolBox::isArtificalGene(gene_right)) {
      // If our nodes are artificial, there is no real branch length.
      d_gene_geneleft = 0.0;
      d_gene_generight = 0.0;
    } else if(gene_father == tree.getRoot()) {
      // NJ cannot compute branch length with only two nodes in the tree and find a root. So, we divide the branch
      // between these two nodes by a half and set the root in the middle.
      d_gene_geneleft = dmatrix(gene_left, gene_right) * 0.50;
      d_gene_generight = dmatrix(gene_left, gene_right) - d_gene_geneleft;
    } else {
      // We use NJ methods to compute new branch lengths.
      d_gene_geneleft = nj_.DistanceFromPairAncestor(dmatrix, gene_left, gene_right);
      d_gene_generight = dmatrix(gene_left, gene_right) - d_gene_geneleft;
    }

    // To avoid any branch length <= 0.0 we can set a value set in eps.
    double eps = MINIMAL_BRANCH_LENGTH;

    //if(utils::double_equal_or_inferior(d_gene_geneleft, 0.0) and utils::double_equal_or_inferior(d_gene_generight, 0.0)){
    if(d_gene_geneleft < eps and d_gene_generight < eps) {
      d_gene_geneleft = eps;
      d_gene_generight = eps;
    } else if(d_gene_geneleft < eps){
      d_gene_generight += fabs(d_gene_geneleft);
      d_gene_geneleft = eps;
      d_gene_generight -= eps;
    } else if (d_gene_generight < eps){
      d_gene_geneleft += fabs(d_gene_generight);
      d_gene_generight = eps;
      d_gene_geneleft -= eps;
    }

    new_branch_left->setLength(d_gene_geneleft);
    new_branch_right->setLength(d_gene_generight);
  }

  // Update resulting tree
  PhyloTreeToolBox::addNodeWithSonsInPostOrder(tree, gene_father, gene_left, gene_right, new_branch_left,
                                               new_branch_right, verbose);

  // Update distance matrix
  updateDistances(dmatrix, gene_father, gene_left, gene_right);

  if(verbose) std::cout << "Distance matrix state:" << std::endl << dmatrix << std::endl;

  // Update genemap
  genemap_temp.deleteGene(species_left, gene_left);
  genemap_temp.deleteGene(species_right, gene_right);
  genemap_temp.addGene(species, gene_father);

  return gene_father;
}

std::map<Event, std::vector<Node>>
ProfileNJ::operator()(
    bpp::PhyloTree &genetree
    , const Node &polytomy_root
    , const bpp::PhyloTree &speciestree
    , DistanceMatrix &distances
    , const std::unordered_map<std::shared_ptr<bpp::PhyloNode>, std::size_t> &V
    , SpeciesGeneMap &genemap
    , const bool computeBranchLengths
    , const bool verbose)
const {
  /// Implements Neighbor Joining algorithm of ProfileNJ.

  // Get all polytomy leaves
  auto polytomy_leaves = genetree.getSons(polytomy_root);

  // This dictionnary, events, store each created node and the event associated.
  std::map<Event, std::vector<Node>> events;

  // Make a copy of the genemap, the speciesgenemap will be modified...
  SpeciesGeneMap polytomy_map = genemap;

  polytomy_map = SpeciesGeneMapper::reduceMapAccordingToGeneList(genemap, polytomy_leaves);

  //Then create the guide_tree of species.
  //std::shared_ptr<bpp::PhyloTree> polytomy_guidetree = phyloTreeToolBox_.cloneTree(speciestree);
  //phyloTreeToolBox_.pruneTree_v0(*polytomy_guidetree, polytomy_map, false);
  std::shared_ptr<bpp::PhyloTree> polytomy_guidetree
      = PhyloTreeToolBox::cloneTree(
          speciestree,
          PhyloTreeToolBox::getLastCommonAncestor(polytomy_map.getSpecies(),
                                                  speciestree));


  // If verbose, print current parameters
  if (verbose) {
    std::cout << "Starting ProfileNJ algorithm on polytomy rooted at " << polytomy_root << " (species "
              << speciestree.getRoot() << ")..." << std::endl;
    utils::printNodeContent(genetree, polytomy_root);
    std::cout << "> Guide tree : " << *polytomy_guidetree;
  }

  // Get all species nodes
  auto nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(
      *polytomy_guidetree); //S, nodes of the species tree
  Node new_node = nullptr;
  Event new_event = none;

  // NJ run
  for (auto s: nodes) {//For each node s of S in a bottom-up traversal of S Do

    if (verbose) std::cout << "  > Working on: " << s << std::endl;
    if (verbose) std::cout << "    Infos: V[s] = " << V.at(s) << " ; M[s] = " << polytomy_map.ngenes(s) << std::endl;

    if (not speciestree.isLeaf(s)) { // If s is an internal node of S with children sl, sr Do

      auto sons = speciestree.getSons(s);
      auto &sl = sons[0];// son left of s
      auto &sr = sons[1];// son right of s

      //{By construction, V[sl] = V[sr] = n}
      assert(V.at(sl) == V.at(sr));
      std::size_t n = V.at(sl);

      if (verbose) std::cout << "\t* Perform " << n << " speciation(s)" << std::endl;

      for (std::size_t i = 1; i <= n; ++i) {
        //Choose in E(sl)*E(sr) the gene pair (gl, gr) minimizing Q_(gl, gr) and create the node g = (gl, gr).

        // Create E(sl)*E(sr)
        auto sl_genes = polytomy_map.getGenes(sl); //list of genes of the sl
        auto sr_genes = polytomy_map.getGenes(sr); //list of genes of the sr
        assert(sl_genes.size() > 0);
        assert(sr_genes.size() > 0);

        // Choose the pair gl_gr minimizing Q_
        std::pair<Node, Node> gl_gr = nj_.findPairMinimizingQ(distances, sl_genes, sr_genes, verbose);
        assert(gl_gr.first != gl_gr.second);

        // Update tree, maps and distance matrix
        new_node = updateTree_GeneMap_DistanceMatrix(genetree, s, sl, gl_gr.first, sr, gl_gr.second, polytomy_map, genemap,
                                                     distances, polytomy_root, computeBranchLengths, verbose);
      }
    }

    // Get the number m of genes in a specific species s.
    std::size_t m = polytomy_map.ngenes(s);

    if (m > V.at(s)) {
      //{Perform m(s) - V[s] duplications}
      if (verbose) std::cout << "\t* Perform " << m - V.at(s) << " duplications" << std::endl;
      for (std::size_t i = 1; i <= m - V.at(s); i++) {
        //Choose in E(s)*E(s) the gene pair (g1, g2) minimizing Q_(...), and create the node g = (g1, g2)
        auto s_genes = polytomy_map.getGenes(s); //list of genes of s

        std::pair<Node, Node> gl_gr = nj_.findPairMinimizingQ(distances, s_genes, s_genes, verbose);
        assert(gl_gr.first != gl_gr.second);

        new_node = updateTree_GeneMap_DistanceMatrix(genetree, s, s, gl_gr.first, s, gl_gr.second, polytomy_map, genemap,
                                                     distances, polytomy_root,
                                                     computeBranchLengths, verbose);

        if (not(PhyloTreeToolBox::isArtificalGene(gl_gr.first) or PhyloTreeToolBox::isArtificalGene(gl_gr.second))) {
          new_event = duplication;
          events[new_event].push_back(new_node);
        }
      }
    } else if (m < V.at(s)) {
      //{Perform V[s] - m(s) losses}
      if (verbose) std::cout << "\t* Perform " << V.at(s) - m << " loss(es) in " << s << std::endl;
      std::size_t i = 0;
      while (i < V.at(s) - m) {
        // Create the artificial gene as node
        new_node = Node(new bpp::PhyloNode);
        new_node->setProperty("isArtificialGene", IsArtificialGene(true));

        // Set a name to the new artificial node.
        std::string new_node_name = LOSS_STR_FLAG;
        new_node->setName(new_node_name);

        // Put the node into the gene tree.
        genetree.createNode(new_node);

        // Update map
        polytomy_map.addGene(s, new_node);
        genemap.addGene(s, new_node);

        // Update events.
        new_event = loss;
        events[new_event].push_back(new_node);
        i++;
      }
    }
  }

  // Reset nodes id of the gene tree according to these changes.
  PhyloTreeToolBox::resetNodeIdInPostOrder(genetree);

  // Test of the bpp::PhyloTree is valid.
  assert(genetree.isValid());

  // End !
  if (verbose) std::cout << "...done." << std::endl;
  if (verbose) std::cout << "Result: " << genetree << std::endl;

  return events;
}

} // namespace treerecs

