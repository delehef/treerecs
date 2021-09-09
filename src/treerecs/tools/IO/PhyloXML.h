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

#ifndef PHYLOXML_H
#define PHYLOXML_H

#include <unordered_map>

#include <Bpp/Phyl/Io/IoTree.h>

#include <treerecs/containers/Cost.h>

namespace treerecs {

using Node = std::shared_ptr<bpp::PhyloNode>;

/*!
 * @class IPhyloXML
 * @brief Read a PhyloXML file.
 */
class IPhyloXML:
  public bpp::AbstractITree,
  public bpp::AbstractIMultiTree
{
private:

protected:

  /*!
   * @brief Edit bpp::PhyloNode with taxonomy infos in PhyloXML description.
   * @details Save each info in Node property. Each taxonomy property has a name like "taxonomy:<property name>".
   *          For example: "taxonomy:scientific_name".
   *
   * @param description PhyloXML description. String which contains the XML taxonomy.
   * @param node Node to edit with taxonomy infos.
   * @param cursor Position of the taxonomy in the phyloxml description.
   */
  void setTaxonomyPhyloXMLTagToPhyloNodeProperty (
      const std::string& description // String which contains the XML taxonomy.
      , bpp::PhyloNode& node // Node to edit with taxonomy infos.
      , std::size_t& cursor // Position to start XML taxonomy read. Change value to position of the end tag character '>'.
  ) const ;

public:

  IPhyloXML(): bpp::AbstractITree(), bpp::AbstractIMultiTree() {};
  virtual ~IPhyloXML() = default;

  const std::string getFormatName() const;

  const std::string getFormatDescription() const;

  /*!
   * @brief Read a tree from a file.
   *
   * @param in The input stream.
   * @return A new tree object.
   * @throw Exception If an error occured.
   */
  bpp::PhyloTree* readP(std::istream& in) const;

  /*!
     * @brief Read trees from a file.
     *
     * @param path The file path.
     * @param trees The output trees container.
     * @throw Exception If an error occured.
     */
  void read(const std::string& path, std::vector<bpp::PhyloTree*>& trees) const;

  /*!
   * @brief Read trees from a stream.
   *
   * @param in The input stream.
   * @param trees The output trees container.
   * @throw Exception If an error occured.
   */
  void read(std::istream& in, std::vector<bpp::PhyloTree*>& trees) const;

  /// PhyloXML to PhyloTrees.
  std::vector<bpp::PhyloTree*> phyloXMLToPhyloTrees(const std::string& description) const ;

  /// PhyloXML to bpp::PhyloTree.
  bpp::PhyloTree* phyloXMLToPhyloTree (
      const std::string& description /// Tree/phylogeny in PhyloXML.
  ) const ;

  /// Read a Phylogeny in a phyloXML, create a bpp::PhyloTree*.
  bpp::PhyloTree* phylogenyToPhyloTree(
      const std::string& description /// Tree/ phylogeny in PhyloXML.
      , std::size_t& cursor /// Position in description of the PhyloXML tree/phylogeny.
  )const ;

  /// Read a PhyloXML clade, create a bpp::PhyloNode, put it into the tree. Returns a smart pointer to bpp:PhyloNode.
  std::shared_ptr<bpp::PhyloNode> phyloXMLCladeToPhyloNode(
      bpp::PhyloTree& tree /// bpp::PhyloTree to edit.
      , const std::shared_ptr<bpp::PhyloNode>& father_node /// bpp::PhyloNode father of the new node.
      , const std::string& description /// Tree in PhyloXML format.
      , std::size_t& nodeCounter /// Node counter. Incremented for each node creation.
      , std::size_t& cursor /// Position of the new Clade info in PhyloXML tree given by description.
      ) const ;
};

/*!
 * @class OPhyloXML
 * @brief Write a PhyloXML and RecPhyloXML file.
 * @details RecPhyloXML is an extension of the PhyloXML format for Reconciled gene trees with species tree.
 */
class OPhyloXML:
    public bpp::AbstractOTree,
    public bpp::AbstractOMultiTree
{
  // Prevents hiding valid overrides
  using bpp::AbstractOTree::write;
  using bpp::AbstractOMultiTree::write;
private:
  /// Indent line.
  std::string line_indentation(const std::size_t indent_level) const;

  /// Print an opening tag.
  std::string opening_tag(const std::string& field, const std::string& param = "", const std::string& value = "") const;

  /// Print a closing tag.
  std::string closing_tag(const std::string& field) const;

  static std::string geneNameField;
  static std::string speciesField;


protected:
  /// Print a gene event in the RecPhyloXML format.
  void eventToPhyloXML(std::ostream& output
                       , const Node& gene_node
                       , std::unordered_map<Node, Event>& associated_events
                       , const std::size_t indent_level) const;

  /// Print a reconciled gene node in the RecPhyloXML format.
  void nodeToRecPhyloXML(std::ostream &output, const Node &gene_node, const bpp::PhyloTree &genetree,
                         std::unordered_map<Node, Event> &associated_events, std::size_t indent_level) const;

  /// Print a node in the PhyloXML format.
  void nodeToPhyloXML(std::ostream &output, const bpp::PhyloTree &tree, Node node,
                      std::size_t indent_level) const;

public:

  OPhyloXML(): bpp::AbstractOTree(), bpp::AbstractOMultiTree() {};

  virtual ~OPhyloXML() = default;

  /// Print in output a gene tree in the RecPhyloXML format.
  void geneTreeToRecPhyloXML(std::ostream &output, const bpp::PhyloTree &tree,
                             int index = -1, std::size_t indent_level = 0) const;

  /// Print in output a species tree in the RecPhyloXML format.
  void speciesTreeToPhyloXML(std::ostream& output
                             , const bpp::PhyloTree& speciestree
                             , std::size_t indent_level = 0) const;

  /// Print in output a tree in the PhyloXML format.
  void treeToPhyloXML(std::ostream& output
                      , const bpp::PhyloTree& tree
                      , const std::string& description = ""
                      , std::size_t indent_level = 0) const;

  /// Print a tree in the PhyloXML format.
  void phylogenyToPhyloXML_(std::ostream &output, const bpp::PhyloTree &tree, const std::string &description = "",
                            std::size_t indent_level = 0) const;

  /// Print tree in stream in PhyloXML format.
  virtual void write(const bpp::PhyloTree& tree, std::ostream& out) const;

  /// Print trees in stream in PhyloXML format.
  void write(const std::vector<const bpp::PhyloTree*>& trees, std::ostream& out) const;
};

/*!
 * @class PhyloXML
 * @brief Provides methods to read or write a PhyloXML file. Can write Reconciled trees in RecPhyloXML format.
 * @details See IPhyloXML and OPhyloXML for more details.
 */
class PhyloXML:
  public IPhyloXML,
  public OPhyloXML {
public:
  PhyloXML(): IPhyloXML(), OPhyloXML() {}
  ~PhyloXML() = default;
};

} // namespace treerecs

#endif //PHYLOXML_H
