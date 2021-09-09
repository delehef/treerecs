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

#include "PhyloXML.h"

// Bpp includes
#include <Bpp/BppString.h>

// Treerecs includes
#include "XMLUtils.h"
#include <assert.h>
#include <treerecs/tools/utils.h>
#include <treerecs/tools/PhyloTreeToolBox.h>

namespace treerecs {

void IPhyloXML::setTaxonomyPhyloXMLTagToPhyloNodeProperty(const std::string &description, bpp::PhyloNode &node,
                                                         std::size_t &cursor) const {
  XMLTag tag_xml = XMLUtils::readNextXMLTag(description, cursor);

  assert(tag_xml.name == "taxonomy" and tag_xml.type == XMLTagType::start);

  if(tag_xml.attributes.find("id_source") != tag_xml.attributes.end()){
    node.setProperty("taxonomy:id_source", bpp::BppString(tag_xml.attributes.at("id_source")));
  }

  tag_xml = XMLUtils::readNextXMLTag(description, cursor);
  while(cursor < description.size() and not(tag_xml.name == "taxonomy" and tag_xml.type == XMLTagType::end)) {
    if(tag_xml.name == "id" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string id = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:id", bpp::BppString(id));
    }

    else if(tag_xml.name == "code" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string code = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:code", bpp::BppString(code));
    }

    else if(tag_xml.name == "scientific_name" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string scientific_name = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:scientific_name", bpp::BppString(scientific_name));
    }

    else if(tag_xml.name == "authority" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string authority = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:authority", bpp::BppString(authority));
    }

    else if(tag_xml.name == "common_name" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string common_name = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:common_name", bpp::BppString(common_name));
    }

    else if(tag_xml.name == "synonym" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string synonym = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:synonym", bpp::BppString(synonym));
    }

    else if(tag_xml.name == "rank" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string rank = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:rank", bpp::BppString(rank));
    }

    else if(tag_xml.name == "uri" and  tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string uri = XMLUtils::readXMLElement(description, cursor);
      node.setProperty("taxonomy:uri", bpp::BppString(uri));
    }

    else {
      //if(tag_xml.type != XMLTagType::end)
      //  std::cerr << "Warning, taxonomy " << tag_xml.name << " tag: bad name or not supported yet." << std::endl;

      if(tag_xml.type == XMLTagType::start)
        XMLUtils::goToEndTag(tag_xml.name, description, cursor);
    }

    tag_xml = XMLUtils::readNextXMLTag(description, cursor);
  }
}

std::shared_ptr<bpp::PhyloNode>
IPhyloXML::phyloXMLCladeToPhyloNode(bpp::PhyloTree &tree, const std::shared_ptr<bpp::PhyloNode> &father_node,
                                   const std::string &description, std::size_t &nodeCounter, std::size_t &cursor) const {
  /*!
   * @details 'Cursor' gives position (in 'description') of the clade start xml tag of the new clade. This will place a
   *          new bpp::PhyloNode in tree as a son of 'father_node'.
   *          If there is no 'father_node' ('= nullptr'), the new node will be inserted without father (like a root).
   */

  XMLTag tag_xml = XMLUtils::readNextXMLTag(description, cursor);

  cursor++;

  assert(tag_xml.type == XMLTagType::start and tag_xml.name == "clade");

  std::shared_ptr<bpp::PhyloNode> current_node (new bpp::PhyloNode());
  std::shared_ptr<bpp::PhyloBranch> edge_to_father (new bpp::PhyloBranch());

  // Check attributes
  // if branch length attribute, set it into the branch/edge to father.
  if(tag_xml.attributes.find("branch_length") != tag_xml.attributes.end()) {
    edge_to_father->setLength(std::atof(tag_xml.attributes.at("branch_length").c_str()));
  }

  // Add node and/ or branch into the tree.
  if(father_node == nullptr){
    tree.createNode(current_node);
  }
  else {
    tree.createNode(father_node, current_node, edge_to_father);
  }

  // Read next tag until the xml closing tag of the current clade.
  tag_xml = XMLUtils::readNextXMLTag(description, cursor);

  cursor++;

  while(not (tag_xml.name == "clade" and tag_xml.type == XMLTagType::end)) {
    // Read inner tags.

    // Add clade with its sons by recursion.
    if (tag_xml.name == "clade" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      phyloXMLCladeToPhyloNode(tree, current_node, description, nodeCounter, cursor);
    }

    // Add clade name.
    else if(tag_xml.name == "name" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string node_name = XMLUtils::readXMLElement(description, cursor);

      current_node->setName(node_name);
    }

    // Add branch to father length.
    else if(tag_xml.name == "branch_length" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string branch_length_str = XMLUtils::readXMLElement(description, cursor);

      double branch_length = atof(branch_length_str.c_str());
      edge_to_father->setLength(branch_length);
    }

    // Add confidence in the branch to father.
    else if(tag_xml.name == "confidence" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string confidence_str = XMLUtils::readXMLElement(description, cursor);

      if(edge_to_father)
        edge_to_father->setProperty("bootstrap", bpp::Number<double>(bpp::TextTools::toDouble(confidence_str)));
    }

    // Add clade description in node property.
    else if(tag_xml.name == "description" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string description_str = XMLUtils::readXMLElement(description, cursor);

      current_node->setProperty("description", bpp::BppString(description_str));
    }

    // Add taxonomy infos in node property.
    else if(tag_xml.name == "taxonomy" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      setTaxonomyPhyloXMLTagToPhyloNodeProperty(description, *current_node, cursor);
    }

    // Default: not supported.
    else {
      //std::cerr << "Unrecognized Clade PhyloXML element called '" << tag_xml.name << "'. Bad name or not supported yet." << std::endl;
      if(tag_xml.type == XMLTagType::start)  cursor = XMLUtils::goToEndTag(tag_xml.name, description, cursor);
    }

    tag_xml = XMLUtils::readNextXMLTag(description, cursor);
    cursor++;
  }

  // Set indexes in tree.
  tree.setNodeIndex(current_node, (unsigned int) nodeCounter);

  if(edge_to_father) tree.setEdgeIndex(edge_to_father, (unsigned int) nodeCounter);

  // If there is no name for the new node, give the id.
  if(not current_node->hasName()) current_node->setName(std::to_string(nodeCounter));

  // Increment nodeCounter.
  nodeCounter++;

  return current_node;
}

bpp::PhyloTree *IPhyloXML::phylogenyToPhyloTree(const std::string &description, std::size_t &cursor) const {
  /*!
   * @details 'Cursor' gives the position in 'description' of the start tag of the new phylogeny. This method will
   *          create a pointer to a new bpp::PhyloTree.
   */

  XMLTag tag_xml;

  tag_xml = XMLUtils::readNextXMLTag(description, cursor);

  assert(tag_xml.name == "phylogeny" and tag_xml.type == XMLTagType::start);
  if(not (tag_xml.name == "phylogeny" and tag_xml.type == XMLTagType::start)){
    throw bpp::IOException ("PhyloXML: failed to read phylogeny.");
    return NULL;
  }

  // Search of the phylogeny tag.
  while(not (tag_xml.name == "phylogeny" and tag_xml.type == XMLTagType::start)) {
    tag_xml = XMLUtils::readNextXMLTag(description, cursor);
    if(cursor >= (description.size() - 1)){
      std::cerr << "No phylogeny found." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  cursor++;

  if(cursor >= description.size()) return NULL;

  //else there is a phylogeny, so create a tree.

  bool rooted = false;

  bpp::PhyloTree* tree = new bpp::PhyloTree(rooted);

  std::shared_ptr<bpp::PhyloNode> first_node = nullptr;

  tag_xml = XMLUtils::readNextXMLTag(description, cursor);
  cursor++;

  while(not (tag_xml.name == "phylogeny" and tag_xml.type == XMLTagType::end)) {

    if(tag_xml.name == "name" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::string tree_name = XMLUtils::readXMLElement(description, cursor);

      tree->setName(tree_name);
    }

    else if(tag_xml.name == "clade" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      std::size_t nodeCounter = 0;
      std::shared_ptr<bpp::PhyloNode> created_node = phyloXMLCladeToPhyloNode(*tree, first_node, description, nodeCounter, cursor);
      if(not first_node) first_node = created_node;
    }

    else {
      if(tag_xml.type == start) {
        //std::cerr << "Unrecognized Phylogeny PhyloXML element called '" << tag_xml.name
        //          << "...'. Bad name or not supported yet." << std::endl;
        cursor = XMLUtils::goToEndTag(tag_xml.name, description, cursor, description.size());
      } else if(tag_xml.type != end) {
        //std::cerr << "Unrecognized Phylogeny PhyloXML element called '" << tag_xml.name
        //          << "...'. Bad name or not supported yet." << std::endl;
      }
    }

    tag_xml = XMLUtils::readNextXMLTag(description, cursor);
    cursor++;
  }

  tree->rootAt(first_node);

  return tree;

}

bpp::PhyloTree *IPhyloXML::phyloXMLToPhyloTree(const std::string &description) const {
  /// @details Read and create only the first phylogeny in PhyloXML 'description'.
  std::size_t cursor = 0;

  XMLTag tag_xml = XMLUtils::readNextXMLTag(description, cursor);

  assert(tag_xml.name == "phyloxml" and tag_xml.type == XMLTagType::start);

  if(not(tag_xml.name == "phyloxml" and tag_xml.type == XMLTagType::start)){
    throw bpp::IOException ("PhyloXML: failed to read phyloxml description.");
    return NULL;
  }

  return phylogenyToPhyloTree(description, cursor);
}

std::vector<bpp::PhyloTree *> IPhyloXML::phyloXMLToPhyloTrees(const std::string &description) const {
  /*!
   * @details Read and create trees in a PhyloXML 'description'. Returns std::vector of pointers to bpp::PhyloTree.
   */
  std::size_t cursor = 0;

  XMLTag tag_xml = XMLUtils::readNextXMLTag(description, cursor);

  while(tag_xml.type != XMLTagType::start and tag_xml.name != "phyloxml" and cursor < description.size()) {
    tag_xml = XMLUtils::readNextXMLTag(description, cursor);
  }

  //assert(tag_xml.name == "phyloxml" and tag_xml.type == XMLTagType::start);
  if(not(tag_xml.name == "phyloxml" and tag_xml.type == XMLTagType::start)) {
    throw bpp::IOException ("PhyloXML: error during phyloxml read.");
    return {};
  }

  std::list<bpp::PhyloTree*> trees;

  tag_xml = XMLUtils::readNextXMLTag(description, cursor);
  while(not (tag_xml.name == "phyloxml" and tag_xml.type == XMLTagType::end) ) {
    if(tag_xml.name == "phylogeny" and tag_xml.type == XMLTagType::start) {
      XMLUtils::reachChar(cursor, description, '<', false);
      trees.emplace_back(phylogenyToPhyloTree(description, cursor));
    } else {
      //std::cerr << "Unknown phyloxml tag \"" << tag_xml.name << "\", bad name or not supported yet." << std::endl;
    }

    tag_xml = XMLUtils::readNextXMLTag(description, cursor);
  }
  return {trees.begin(), trees.end()};
}

const std::string IPhyloXML::getFormatName() const { return "PhyloXML"; }

const std::string IPhyloXML::getFormatDescription() const {
  return getFormatName()
         + ", XML for evolutionary biology and comparative genomics. See http://www.phyloxml.org for more info.";
}

bpp::PhyloTree *IPhyloXML::readP(std::istream &in) const {
  // Checking the existence of specified file
  if (! in) { throw bpp::IOException ("PhyloXML::read: failed to read from stream"); }

  std::string description((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

  if (bpp::TextTools::isEmpty(description))
    throw bpp::IOException("PhyloXML::read: no tree has been found!");
  return phyloXMLToPhyloTree(description);
}

void IPhyloXML::read(const std::string &path, std::vector<bpp::PhyloTree *> &trees) const {
  std::ifstream in(path.c_str(), std::ios::in);
  if (! in) { throw bpp::IOException ("PhyloXML::read: failed to read from stream"); }
  read(in, trees);
}

void IPhyloXML::read(std::istream &in, std::vector<bpp::PhyloTree *> &trees) const {
  if (! in) { throw bpp::IOException ("PhyloXML::read: failed to read from stream"); }

  std::string description((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

  if (bpp::TextTools::isEmpty(description))
    throw bpp::IOException("PhyloXML::read: no tree was found!");

  trees = phyloXMLToPhyloTrees(description);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string OPhyloXML::geneNameField = "geneName=";
std::string OPhyloXML::speciesField = "speciesLocation=";

/// Create an indentation.
std::string OPhyloXML::line_indentation(const std::size_t indent_level) const {
  return std::string(indent_level, ' ');
  //std::string res;
  //while(res.size() < indent_level) res += ' ';
  //return res;
}

/// Create an opening_tag.
std::string OPhyloXML::opening_tag(const std::string &field, const std::string& param, const std::string& value) const {
  return "<" + field
         + (not param.empty() ? (" " + param + (not value.empty() ? ("=" + value) : "")) : "")
         + ">";
}

/// Create a closing tag.
std::string OPhyloXML::closing_tag(const std::string &field) const {
  return "</" + field + ">";
}

/// Print a RecPhyloXML event in output (std::ostream).
void OPhyloXML::eventToPhyloXML(std::ostream &output, const Node &gene_node, std::unordered_map<Node, Event> &associated_events,
                                  const std::size_t indent_level) const {
  Event event = associated_events.at(gene_node);

  if(event == none) return;

  bpp::BppString& species_node = *dynamic_cast<bpp::BppString*>(gene_node->getProperty("Species name"));

  // Begin line
  output << line_indentation(indent_level);

  // First opening field <...>
  output << "<" << event ;

  if(event == extant) {
    if(gene_node->hasName()) {
      output << " " << geneNameField << "\"" << gene_node->getName() << "\"";
    }
  }

  if(event != bifurcationOut) {
    output << " ";
    output << speciesField;
    output << "\"";
    output << species_node;
    output << "\"";
  }

  output << ">";

  // Second ending field </>
  output << closing_tag(event_to_str(event)) << std::endl;

  return;
}

/// Print a RecPhyloXML node in output.
void OPhyloXML::nodeToRecPhyloXML(std::ostream &output, const Node &gene_node, const bpp::PhyloTree &genetree,
                                  std::unordered_map<Node, Event> &associated_events,
                                  std::size_t indent_level) const {

  Event event = associated_events.at(gene_node);
  if(event == none)// ignoring any Null event.
  {
    auto children = genetree.getSons(gene_node);
    nodeToRecPhyloXML(output, children.front(), genetree, associated_events, indent_level);//recursion on the first (and only) son of the null node
  }
  else
  {
    output << line_indentation(indent_level);
    output << opening_tag("clade") << std::endl;
    indent_level++;

    //writing the different informations:

    //1. name:
    output << line_indentation(indent_level);
    output << opening_tag("name");
    if(gene_node->hasName())
    {
      output << gene_node->getName();
    }
    else
    {
      output << "noname";
    }
    output << closing_tag("name") << std::endl;

    //2. the reconciliation event(s)
    output << line_indentation(indent_level);
    output << opening_tag("eventsRec") << std::endl;
    indent_level++;

    auto current_gene_node = gene_node;
    auto next_gene_node = current_gene_node;
    bool donext = true;

    while( donext )
    {
      current_gene_node = next_gene_node; //iteration

      auto children = genetree.getSons(current_gene_node);

      if(children.size() == 0)
        donext = false;

      for(unsigned i = 0 ; i < children.size(); i++)
      {
        auto event = associated_events.at(children[i]);
        if(event != loss) {
          if(next_gene_node  == current_gene_node)//this is the first non-loss child of current gene node -> we accept it as nextid
            next_gene_node = children[i];
          else //this is not the first non-loss child of current id -> we get out of the loop
          {
            donext = false;
            break;
          }
        }
      }
      //we write the event of the current node
      eventToPhyloXML(output, current_gene_node, associated_events, indent_level);
    }

    indent_level--;
    output << line_indentation(indent_level);
    output << closing_tag("eventsRec") << std::endl;


    //3. Sons, if any
    auto children = genetree.getSons(current_gene_node);

    for(unsigned i = 0 ; i < children.size(); i++)
    {
      nodeToRecPhyloXML(output, children.at(i), genetree, associated_events, indent_level);
    }

    indent_level--;
    output << line_indentation(indent_level);
    output << closing_tag("clade") << std::endl;
  }
}

/// Print a tree in PhyloXML format.
void OPhyloXML::nodeToPhyloXML(std::ostream &output, const bpp::PhyloTree &tree, Node node,
                               std::size_t indent_level) const {
  //cout << N->getInfos().postOrder << " " << N->getInfos().timeSlice << " " << N->getInfos().realPostOrder << " " << nodeid << endl;

  output << line_indentation(indent_level);
  output << opening_tag("clade") << std::endl;
  indent_level++;

  double distance = 0.0;
  double bootstrap = -1;

  Node current_node = node;

  if (tree.hasFather(current_node)){
    auto edgeToFather = tree.getEdgeToFather(current_node);
    if(edgeToFather->hasLength()){
      distance += tree.getEdgeToFather(current_node)->getLength();
    }

    if(edgeToFather->hasBootstrapValue()){
      bootstrap = edgeToFather->getBootstrapValue();
    }
  }

  auto children = tree.getSons(current_node);

  /*
  while(children.size() == 1)
  {
    current_node = children.front();
    children = tree.getSons(current_node);
    if (tree.hasFather(current_node) and tree.getEdgeToFather(current_node)->hasLength())
        distance += tree.getEdgeToFather(current_node)->getLength();
  }
   */

  node = current_node;

  //1. name:
  output << line_indentation(indent_level);
  output << opening_tag("name");
  if(tree.isLeaf(node)) //internal nodes actually have names that do not correspond to their RPO but the TS of the speciation
  {
    output << node->getName();
  }
  else
  {
    output << node->getName();
  }

  output << closing_tag("name") << std::endl;

  //2. distance to father and bootstrap
  if(not utils::double_equivalence(distance, 0.0))
  {
    output << line_indentation(indent_level) ;
    output << opening_tag("branch_length") << std::to_string(distance) << closing_tag("branch_length") << std::endl ;
  }

  if(not utils::double_equivalence(bootstrap, -1))
  {
    output << line_indentation(indent_level);
    output << opening_tag("confidence", "type", "\"bootstrap\"") << std::to_string(bootstrap) << closing_tag("confidence") << std::endl;
  }

  //3. id
  //output << line_indentation(indent_level) << std::endl;
  //output << "<node_id>" << N->getInfos().realPostOrder << "</node_id>\n";

  for(unsigned i = 0 ; i < children.size(); i++)
  {
    nodeToPhyloXML(output, tree, children.at(i), indent_level);
  }

  indent_level--;
  output << line_indentation(indent_level);
  output << closing_tag("clade") << std::endl;
}

void OPhyloXML::geneTreeToRecPhyloXML(std::ostream &output, const bpp::PhyloTree &tree, const int index, std::size_t indent_level) const {
  output << line_indentation(indent_level);
  output << opening_tag("recGeneTree") << std::endl;

  indent_level++;

  output << line_indentation(indent_level);
  output << opening_tag("phylogeny", "rooted", "\"true\"") << std::endl;

  indent_level++;

  if( index != -1)
  {
    output << line_indentation(indent_level);
    output << opening_tag("id") << index << closing_tag("id") << std::endl;
  }

  auto root = tree.getRoot();

  // Deduce gene events in the gene tree.
  auto nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(tree);
  std::unordered_map<Node, Event> events;

  for(auto& node: nodes) {
    if(tree.isLeaf(node)){
      if(PhyloTreeToolBox::isArtificalGene(node)) {
        events[node] = loss;
      } else {
        events[node] = extant;
      }
    } else {
      auto sons = tree.getSons(node);
      auto& son_left = sons.front();
      auto& son_right = sons.back();
      assert(son_left->hasProperty("Species name"));
      assert(son_right->hasProperty("Species name"));

      bpp::BppString& son_left_species = *dynamic_cast<bpp::BppString*>(son_left->getProperty("Species name"));
      bpp::BppString& son_right_species = *dynamic_cast<bpp::BppString*>(son_right->getProperty("Species name"));

      if(son_left_species.toSTL() == son_right_species.toSTL()){
        events[node] = duplication;
      } else {
        if(events.at(son_left) == loss or events.at(son_right) == loss){
          events[node] = speciationLoss;
        } else {
          events[node] = speciation;
        }
      }
    }
  }

  nodeToRecPhyloXML(output, root, tree, events, indent_level);


  indent_level--;
  output << line_indentation(indent_level);
  output << closing_tag("phylogeny") << std::endl;

  indent_level--;
  output << line_indentation(indent_level);

  output << closing_tag("recGeneTree") << std::endl;

}

void OPhyloXML::phylogenyToPhyloXML_(std::ostream &output, const bpp::PhyloTree &tree, const std::string &description,
                                     std::size_t indent_level) const {
  output << line_indentation(indent_level);
  if(tree.isRooted())
    output << opening_tag("phylogeny", "rooted", "\"true\"") << std::endl;
  else
    output << opening_tag("phylogeny", "rooted", "\"false\"") << std::endl;

  if(not description.empty())
    output << line_indentation(indent_level) << opening_tag("description") << description << closing_tag("description") << std::endl;

  indent_level++;
  auto root = tree.getRoot();
  nodeToPhyloXML(output, tree, root, indent_level);

  indent_level--;
  output << line_indentation(indent_level);
  output << closing_tag("phylogeny") << std::endl;
}

void OPhyloXML::speciesTreeToPhyloXML(std::ostream& output, const bpp::PhyloTree& speciestree, std::size_t indent_level) const {
  output << line_indentation(indent_level);
  output << opening_tag("spTree") << std::endl;

  indent_level++;
  phylogenyToPhyloXML_(output, speciestree, "", indent_level);

  indent_level--;
  output << line_indentation(indent_level);
  output << closing_tag("spTree") << std::endl;
}

void OPhyloXML::treeToPhyloXML(std::ostream& output, const bpp::PhyloTree& tree, const std::string& description, std::size_t indent_level) const {
  output << line_indentation(indent_level);
  output << opening_tag("phyloxml") << std::endl;

  indent_level++;
  phylogenyToPhyloXML_(output, tree, description, indent_level);

  indent_level--;
  output << line_indentation(indent_level);
  output << closing_tag("phyloxml") << std::endl;
}

void OPhyloXML::write(const bpp::PhyloTree &tree, std::ostream &out) const {
  treeToPhyloXML(out, tree);
}

void OPhyloXML::write(const std::vector<const bpp::PhyloTree *> &trees, std::ostream &out) const {
  std::size_t indent_level = 0;
  out << line_indentation(indent_level);
  out << opening_tag("phyloxml") << std::endl;

  for(auto tree: trees) {
    indent_level++;
    phylogenyToPhyloXML_(out, *tree, "", indent_level);
    indent_level--;
  }

  out << line_indentation(indent_level);
  out << closing_tag("phyloxml") << std::endl;
}

} // namespace treerecs
