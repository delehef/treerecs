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

#ifndef TREERECS_PHYLOXMLTOSVG_H
#define TREERECS_PHYLOXMLTOSVG_H

// Treerecs includes
#include <treerecs/containers/Grid.h>
#include <treerecs/tools/random.h>
#include <treerecs/tools/PhyloTreeToolBox.h>
#include <treerecs/containers/ReconciledRootedTree.h>

// External includes.
#include <simple_svg_1.0.0.hpp>

namespace treerecs {

// Colors chart:
// * leaf = 27, 153, 139
// * dupl = 119, 125, 167
// * loss = 254, 95, 85
// * species = 241, 222, 222
// * gene branches = 161,130,118
// * spec = 187, 172, 193

/*!
 * @class PhyloXMLToSVG
 * @brief Draw a full reconciliation in a svg file.
 * @details Use G. Gence and W. Duchemin algorithm from RecPhyloXML-visu.
 *
 */
class RecPhyloTreeToSVG {
private:
  static constexpr double historySize = 25;
  static constexpr double eventSize = 30;

  class Coordinates {
  public:
    double x;
    double y;

    Coordinates() : x(-1), y(-1) {
    }

    Coordinates(double x_, double y_) : x(x_), y(y_) {
    }

    template<typename T>
    Coordinates(std::initializer_list<T> l) {
      auto current = l.begin();
      x = *current;
      current++;
      y = *current;
    }

    svg::Point toPoint() const {
      return {x, y};
    }

  };

  typedef struct ContainerStruct {
    /*!
   * @struct Container
   * @brief Contains coordinates to draw a species node or gene node.
   * @details A species node is a square delimited by four points. A gene node
     * is given only one point (and will don't use a ContainerUpDown).
   */
    typedef struct ContainerUpDownStruct {
      /*!
       * @struct ContainerUpDown
       * @brief Contains coordinates use by species nodes.
       * @details Up is the entering border of a gene in a species and Down the
       * leaving side of genes from the species.
       */
      Coordinates up;
      Coordinates down;
    } ContainerUpDown;

    ContainerUpDown start; //start.up and start.down
    ContainerUpDown stop;  // stop.up and stop.down

    Coordinates coordinates;
  } Container;

  typedef struct TreeLayoutInfoStruct {
    std::unordered_map<Node, Container> species_containers;
    std::unordered_map<std::shared_ptr<bpp::PhyloTree>,
        std::unordered_map<Node, Container>> genes_containers;
    std::unordered_map<Node,
        std::unordered_map<Node, std::pair<Event, unsigned long>>> events;
    double sizeX = 0;
    double sizeY = 0;
  } TreeLayoutInfo;

protected:

  /*!
   * @brief Draw a cross in svg document.
   * @details The cross is a symbol for a gene loss.
   * @param doc Svg document
   * @param center Position of the center of the cross.
   * @param width
   * @param height
   * @param stroke
   */
  static svg::Document &drawCross(
      svg::Document &doc, const Coordinates &center, const double width
      , const double height, const svg::Stroke &stroke = svg::Stroke()) {
    doc << svg::Line(svg::Point(center.x + width / 2, center.y + height / 2),
                     svg::Point(center.x - width / 2, center.y - height / 2),
                     stroke);
    doc << svg::Line(svg::Point(center.x + width / 2, center.y - height / 2),
                     svg::Point(center.x - width / 2, center.y + height / 2),
                     stroke);
    return doc;
  }

  /*!
   * @brief Draw a star in svg document
   * @param doc
   * @param center
   * @param diameter
   * @param branch_number Number of branches for the star (default = 5)
   * @param fill svg::Fill for the resulting polygon
   * @param stroke svg::Stroke for the resulting polygon
   * @return
   */
  static svg::Document &drawStar(
      svg::Document & doc, const Coordinates &center, const double diameter
      , const unsigned int branch_number = 5, const svg::Fill& fill = svg::Fill(svg::Color::White), const svg::Stroke & stroke = svg::Stroke()) {
    std::vector<double> diameters {1, 0.5};

    std::vector<Coordinates> coordinates;
    coordinates.reserve(branch_number * 2);

    svg::Polygon shape(fill, stroke);

    std::size_t i;
    for (i = 0 ; i < branch_number * 2 ; ++i) {
      double current_radius = diameter * (diameters[i%2]) * 0.5;
      double trigoPortion =  2.0 * PI<double> * i/(branch_number * 2.0);
      double x = (cos(trigoPortion) * current_radius) + center.x;
      double y = (sin(trigoPortion) * current_radius) + center.y;
      Coordinates current(x, y);
      coordinates.emplace_back(current);

      shape << svg::Point(current.x, current.y);
    }

    shape << svg::Point(coordinates.front().x, coordinates.front().y);
    doc << shape;

    return doc;
  }


  /*!
   * @brief Init a layout with the species tree. There is no gene tree yet.
   * @param speciestree species tree.
   * @param nbSpeciesCorridors Number of corridors by species.
   * @param geneIdEvents gene event ids by species.
   * @return TreeLayoutInfo which is all infos to draw species tree and gene
   * trees inside.
   */
  static TreeLayoutInfo speciestree_layout(
      const bpp::PhyloTree &speciestree
      , std::unordered_map<Node, unsigned long> nbSpeciesCorridors
      , std::unordered_map<Node,
      std::unordered_map<Node, std::pair<Event, unsigned long>>> geneIdEvents
  ) {

    TreeLayoutInfo res;
    std::unordered_map<Node, Container> containers;

    double max_depth = 0;
    auto speciesNodesInOrder
        = PhyloTreeToolBox::getNodesInOrderTraversalRecursive_list(
        speciestree);
    auto speciesNodesInPreOrder
        = PhyloTreeToolBox::getNodesInPreOrderTraversalRecursive_list(
        speciestree);
    double y = 10;

    // Compute Y positions
    std::unordered_map<Node, double> speciesY;
    std::unordered_map<Node, double> speciesHeights;
    for (const auto &speciesNode: speciesNodesInOrder) {
      double speciesHeight;
      long nbCorridors = nbSpeciesCorridors.at(speciesNode);
      if (nbCorridors == 0) nbCorridors++;

      speciesHeight = historySize + (nbCorridors * historySize);
      speciesHeights[speciesNode] = speciesHeight;

      y += speciesHeight / 2.0;
      speciesY[speciesNode] = y;
      y += speciesHeight / 2.0;
    }

    // Compute X positions
    double x_root = 50;
    std::unordered_map<Node, double> speciesX;
    std::unordered_map<Node, double> speciesWidths;
    for (const auto &speciesNode: speciesNodesInPreOrder) {
      unsigned long nbGnEvents = 0;
      if (geneIdEvents.find(speciesNode) != geneIdEvents.end()) {
        for (const auto &geneEvents: geneIdEvents.at(speciesNode)) {
          if (nbGnEvents < geneEvents.second.second)
            nbGnEvents = geneEvents.second.second;
        }
      }

      double speciesWidth = eventSize + (nbGnEvents * eventSize);

      speciesWidths[speciesNode] = speciesWidth;

      Node parent = speciestree.hasFather(speciesNode) ? speciestree.getFather(
          speciesNode) : nullptr;

      if (parent) {
        speciesX[speciesNode] =
            speciesX[parent] + speciesWidths[parent] + speciesHeights[parent];
      } else {
        speciesX[speciesNode] = x_root;
      }
    }

    std::unordered_map<Node, Container> speciesContainers;
    for (const auto &speciesNode: speciesNodesInPreOrder) {
      double speciesHeight = speciesHeights[speciesNode];
      double speciesWidth = speciesWidths[speciesNode];
      double speciesx = speciesX[speciesNode];
      double speciesy = speciesY[speciesNode];

      Container container;

      container.start.up.x = speciesx;
      container.start.up.y = speciesy - (speciesHeight / 2.0);

      container.start.down.x = speciesx;
      container.start.down.y = speciesy + (speciesHeight / 2.0);

      container.stop.up.x = speciesx + speciesWidth;
      container.stop.up.y = container.start.up.y;

      container.stop.down.x = container.stop.up.x;
      container.stop.down.y = container.start.down.y;

      if (container.stop.up.x > max_depth) {
        max_depth = container.stop.up.x;
      }
      speciesContainers[speciesNode] = container;
    }

    auto leaves = speciestree.getAllLeaves();
    for (auto leaf: leaves) {
      Container &c = speciesContainers[leaf];
      c.coordinates.x = max_depth + (2 * eventSize);
      c.coordinates.y = y;
    }

    // save tree size.
    res.sizeX = max_depth;
    res.sizeY = y;
    res.species_containers = speciesContainers;
    res.events = geneIdEvents;

    return res;
  }

  /*!
   * @brief Computes all coordinates of genes in a gene tree.
   * @param genetree
   * @param speciestree
   * @param map Map which associates species with genes.
   * @param speciesEvents All gene with their event id by Species (first key is
   *        the species, the gene as second).
   * @param corridorIdGene All gene corridor ids.
   * @param speciesGeometry Species tree nodes coordinates.
   * @return
   */
  static std::unordered_map<Node, Container> genetree_layout(
      const bpp::PhyloTree& genetree,
      const SpeciesGeneMap& map,
      const std::unordered_map<
          Node,
          std::unordered_map<
              Node,
              std::pair<Event, unsigned long>>>& speciesEvents,
      const std::unordered_map<Node, unsigned long>& corridorIdGene,
      const std::unordered_map<Node, Container>& speciesGeometry) {
    std::unordered_map<Node, Container> res;
    for (const auto &gene: genetree.getAllNodes()) {
      Container geneGeo;

      auto species = map.getAssociatedSpecies(gene);
      const Container &speciesGeo = speciesGeometry.at(species);

      // X position
      const std::pair<Event, unsigned long> &event
          = speciesEvents.at(species).at(gene);

      if (event.first == Event::extant) // Evaluate event
        geneGeo.coordinates.x = speciesGeo.stop.up.x;
      else {
        double speciesTopStartX = speciesGeo.start.up.x;
        auto idEvent = event.second;
        geneGeo.coordinates.x =
            speciesTopStartX + ((double) idEvent * eventSize);
      }

      // Y position
      //auto speciesTopStartY = speciesGeo.start.up.y;
      auto speciesBottomY = speciesGeo.start.down.y;
      unsigned long idCorridor = corridorIdGene.at(gene);

      geneGeo.coordinates.y =
          speciesBottomY - ((double) idCorridor * historySize);
      res[gene] = geneGeo;
    }
    return res;
  }

  template<typename GeneTreeIterator, typename MapIterator>
  static void fillLayoutWithGeneTrees(
      TreeLayoutInfo &layout
      , const GeneTreeIterator &genetrees_begin
      , const GeneTreeIterator &genetrees_end
      , const MapIterator &map_begin
      , const std::unordered_map<Node, unsigned long> &idCorridors
      , const std::unordered_map<Node, std::unordered_map<Node,
      std::pair<Event, unsigned long>>> &eventsSummary
  ) {
    MapIterator map_it = map_begin;
    for (auto genetree_it = genetrees_begin;
         genetree_it != genetrees_end; genetree_it++) {
      auto geneRes = genetree_layout(**genetree_it, *map_it,
                                     eventsSummary, idCorridors,
                                     layout.species_containers);
      layout.genes_containers[*genetree_it] = geneRes;
      map_it++;
    }
  }

  /*!
   * @brief Read each gene tree and their nodes to summarize their events
   * (duplications, losses) in a map.
   * @tparam GeneTreeIterator
   * @tparam MapIterator
   * @param genetree_begin
   * @param genetree_end
   * @param speciestree
   * @param map_begin
   * @param map_end
   * @return A map of events. The key is the species. The value is a dict with
   * the gene associated with the nature of the event and an id which is the id
   * of this one.
   */
  template<typename GeneTreeIterator, typename MapIterator>
  static std::unordered_map<Node,
      std::unordered_map<Node, std::pair<Event, unsigned long>>>
  summarizeSpeciesEvents(
      const GeneTreeIterator &genetree_begin
      , const GeneTreeIterator &genetree_end
      , const MapIterator &map_begin
  ) {
    // The key is the species. The value is a dict with the gene associated with
    // the nature of the event and an id which is the id of this one.
    std::unordered_map<Node,
        std::unordered_map<Node, std::pair<Event, unsigned long>>> res;
    MapIterator map_it = map_begin;

    // First we have to determine the number of events per species by an map.
    std::unordered_map<Node, unsigned long> species_events_counts;

    for (auto genetree_it = genetree_begin;
         genetree_it != genetree_end; genetree_it++) {

      // Then work by gene. Each id is given according to a PreOrder traversal.
      const auto &genetree = *genetree_it;
      const auto &map = *map_it;
      auto gene_nodes
          = PhyloTreeToolBox::getNodesInPreOrderTraversalRecursive_list(
              *genetree);
      for (auto &gene_node: gene_nodes) {
        auto associated_species = map.getAssociatedSpecies(gene_node);

        if (species_events_counts.find(associated_species) ==
            species_events_counts.end()) {
          species_events_counts[associated_species] = 1; // if we have not
          // recorded the species, init.
        }

        Event associated_event = PhyloTreeToolBox::getEvent(*genetree, map,
                                                            gene_node);

        res[associated_species][gene_node] = std::pair<Event, unsigned long>(
            associated_event, species_events_counts[associated_species]);

        if (associated_event == none)
          std::cout << "Warning during SVG creation, incorrect event found!"
                    << std::endl;

        if (associated_event == Event::duplication
            or associated_event == Event::speciation
            or associated_event == Event::speciationLoss
            or associated_event == Event::loss
            ) {
          // Fixme: two last conditions are not present in GGence algorithm.
          // Update if we are in duplication event
          species_events_counts[associated_species]++;
        }
      }
      map_it++;
    }
    return res;
  }

  /*!
   * @brief Summarize all paths in a species. One gene is associated with a
   * corridor and gets its Id.
   * @tparam GeneTreeIterator
   * @tparam MapIterator
   * @param genetree_begin
   * @param genetree_end
   * @param speciestree
   * @param map_begin
   * @param map_end
   * @return A pair which summarizes corridors. The first element of the pair is
   * the number of corridors by species nodes. The second is the id of the
   * corridor by gene nodes.
   */
  template<typename GeneTreeIterator, typename MapIterator>
  static std::pair<std::unordered_map<Node, unsigned long>,
      std::unordered_map<Node, unsigned long>>
  summarizeCorridors(
      const GeneTreeIterator &genetree_begin
      , const GeneTreeIterator &genetree_end, const bpp::PhyloTree &speciestree
      , const MapIterator &map_begin
  ) {
    // The first element of the pair is the nb of corridors by species nodes.
    // The second element of the pair is the id of the corridor by gene nodes.

    // For the nb of corridors:
    // The key is the species. The value is the number of corridors inside the
    // species.
    std::unordered_map<Node, unsigned long> nbCorridors;

    // For the corridor ids by genes
    // The key is the gene. The value is the id of the corridor.
    std::unordered_map<Node, unsigned long> corridorIdGenes;

    for (auto species: speciestree.getAllNodes())
      nbCorridors[species] = 0;

    MapIterator map_it = map_begin;
    for (auto genetree_it = genetree_begin;
         genetree_it != genetree_end; genetree_it++) {
      // Then work by gene. Each id is given according to a PostOrder traversal.
      const auto &genetree = *genetree_it;
      const auto &map = *map_it;
      auto gene_nodes
          = PhyloTreeToolBox::getNodesInPostOrderTraversalRecursive_list(
          *genetree);
      for (const auto &gene_node: gene_nodes) {
        auto species = map.getAssociatedSpecies(gene_node);

        if (PhyloTreeToolBox::getEvent(*genetree, map, gene_node) !=
            Event::duplication)
          nbCorridors[species]++;

        corridorIdGenes[gene_node] = nbCorridors[species];
      }
      map_it++;
    }

    return std::make_pair(nbCorridors, corridorIdGenes);
  }

  /*!
   * @brief Creates a layout which contains all coordinates to draw the species
   * tree and its nested gene trees.
   * @tparam GeneTreeIterator
   * @tparam MapIterator
   * @param speciestree
   * @param genetrees_begin
   * @param genetrees_end
   * @param map_begin
   * @param map_end
   * @return TreeLayoutInfo.
   */
  template<typename GeneTreeIterator, typename MapIterator>
  static TreeLayoutInfo createLayout(
      const bpp::PhyloTree &speciestree
      , const GeneTreeIterator &genetrees_begin
      , const GeneTreeIterator &genetrees_end
      , const MapIterator &map_begin
  ) {
    // First, init data
    // Corridors helps to get the column in species to place a gene event.
    auto corridors = summarizeCorridors(genetrees_begin, genetrees_end,
                                        speciestree, map_begin);

    // Events helps to get the floor in species to place a gene event.
    auto events = summarizeSpeciesEvents(genetrees_begin, genetrees_end,
                                         map_begin);

    // Then compute coordinates for each elements
    auto layout = speciestree_layout(speciestree, corridors.first, events);

    alignNodes(layout, speciestree);

    // Finally, the layout needs to include gene names at the top right of the
    // svg file.
    // We need to find the maximum length of a gene name and take care of it in
    // the maximum size of the drawing in x axis.
    std::size_t max_gene_name_size = 0;
    for (auto genetree_it = genetrees_begin;
         genetree_it != genetrees_end; genetree_it++) {
      for (auto leaf_name : (*genetree_it)->getAllLeavesNames()) {
        if (leaf_name.size() > max_gene_name_size)
          max_gene_name_size = leaf_name.size();
      }
    }

    layout.sizeX += max_gene_name_size * 25;

    fillLayoutWithGeneTrees(layout, genetrees_begin, genetrees_end,
                            map_begin, corridors.second, events);

    return layout;
  }

  /*!
   * @brief Draw a line between two points in a matrix of strings.
   * @param res A matrix which is going to be printed. This matrix represents
   *        (with chars or strings) the printed trees.
   * @param p1
   * @param p2
   * @param c Character to draw. This char represents the line.
   */
  static void addEdgeBetweenTwoPointsInGrid(
      Grid<std::string> &res, const Coordinates &p1, const Coordinates &p2
      , const char c = '.'
  ) {
    auto start_x = p1.x;
    auto start_y = p1.y;

    auto end_x = p2.x;
    auto end_y = p2.y;

    int diff_x = (int) end_x - (int) start_x;
    int diff_y = (int) end_y - (int) start_y;

    int n_points = (std::abs(diff_x) > std::abs(diff_y) ? std::abs(diff_x)
                                                        : std::abs(diff_y));

    for (int i = 0; i <= n_points; i++) {
      int current_x = static_cast<int>(
          (end_x - start_x) * ((double) i / n_points) + start_x);
      int current_y = static_cast<int>(
          (end_y - start_y) * ((double) i / n_points) + start_y);
      if (current_x < 0) {}
      else if (current_y < 0) {}
      else if (current_x >= static_cast<int>(res.nrow())) {}
      else if (current_y >= static_cast<int>(res.ncol())) {}
      else
        res(static_cast<const size_t &>(current_x),
            static_cast<const size_t &>(current_y)) = c;
    }
  }

  /*!
   * @brief Translate all leaves to the maximum value in x.
   * @param layout
   * @param speciestree
   */
  static void
  alignNodes(TreeLayoutInfo &layout, const bpp::PhyloTree &speciestree) {
    auto leaves = speciestree.getAllLeaves();
    auto &containers = layout.species_containers;
    for (const auto &leaf: leaves) {
      auto &container = containers.at(leaf);
      container.stop.up.x = layout.sizeX;
      container.stop.down.x = layout.sizeX;
    }
  }

  /*!
   * @brief Draw a species node (but not its name).
   * @tparam SVGObject SVG element to update.
   * @param sgvObject SVG element to update.
   * @param species Species to draw.
   * @param speciestree
   * @param nodes_layout Species tree layout which contains all coordinates.
   */
  template<typename SVGObject>
  static void _drawSpNode(
      SVGObject &sgvObject, const Node &species
      , const bpp::PhyloTree &speciestree
      , const std::unordered_map<Node, Container> &nodes_layout
  ) {
    const Container &species_container = nodes_layout.at(species);
    if (not speciestree.isLeaf(species)) {
      auto sons = speciestree.getSons(species);

      const auto &son_left = sons.front();
      const Container &son_left_container = nodes_layout.at(son_left);

      const auto &son_right = sons.back();
      const Container &son_right_container = nodes_layout.at(son_right);

      sgvObject << species_container.start.up.toPoint()
                << svg::Point(species_container.start.up.x,
                              son_left_container.start.up.y);

      _drawSpNode(sgvObject, son_left, speciestree, nodes_layout);

      sgvObject << svg::Point(species_container.stop.up.x,
                              son_left_container.start.down.y)
                << species_container.stop.up.toPoint()
                << species_container.stop.down.toPoint()
                << svg::Point(species_container.stop.down.x,
                              son_right_container.start.up.y);

      _drawSpNode(sgvObject, son_right, speciestree, nodes_layout);

      sgvObject << svg::Point(species_container.start.down.x,
                              son_right_container.start.down.y)
                << species_container.start.down.toPoint();
    } else {
      sgvObject << species_container.start.up.toPoint()
                << species_container.stop.up.toPoint()
                << species_container.stop.down.toPoint()
                << species_container.start.down.toPoint();
    }
  }

  /*!
   * @brief Draw a gene node (with names)
   * @param svg Document to update.
   * @param gene
   * @param genetree
   * @param nodes_layout All coordinates of the current gene tree.
   * @param map Object to associate a gene with a species.
   * @param color svg::Color of the gene node in SVG.
   */
  static void _drawGnNode(
      svg::Document &svg, const Node &gene, const bpp::PhyloTree &genetree
      , const std::unordered_map<Node, Container> &nodes_layout
      , const SpeciesGeneMap &map, const svg::Color &color
  ) {
    const Coordinates &gene_pos = nodes_layout.at(gene).coordinates;

    double diameter = 25;

    svg::Stroke stroke = svg::Stroke(3, color);

    Event event = PhyloTreeToolBox::getEvent(genetree, map, gene);

    if (not genetree.isLeaf(gene)) {
      auto sons = genetree.getSons(gene);

      const auto &son_left = sons.front();
      const Coordinates &son_left_pos = nodes_layout.at(son_left).coordinates;

      const auto &son_right = sons.back();
      const Coordinates &son_right_pos = nodes_layout.at(son_right).coordinates;

      svg << svg::Line(svg::Point(gene_pos.x, gene_pos.y),
                       svg::Point(gene_pos.x, son_left_pos.y), stroke);
      svg << svg::Line(svg::Point(gene_pos.x, son_left_pos.y),
                       svg::Point(son_left_pos.x, son_left_pos.y), stroke);

      _drawGnNode(svg, son_left, genetree, nodes_layout, map, color);

      svg << svg::Line(svg::Point(gene_pos.x, gene_pos.y),
                       svg::Point(gene_pos.x, son_right_pos.y), stroke);
      svg << svg::Line(svg::Point(gene_pos.x, son_right_pos.y),
                       svg::Point(son_right_pos.x, son_right_pos.y), stroke);

      _drawGnNode(svg, son_right, genetree, nodes_layout, map, color);
    } else {
      if (gene->hasName() and event == Event::extant) // Print gene leaf name.
        svg << svg::Text(svg::Point(gene_pos.x + diameter, gene_pos.y),
                         gene->getName(), svg::Fill(color), svg::Font(18));
    }

    switch (event) {
      case Event::duplication : {
        svg << svg::Rectangle(
            svg::Point(gene_pos.x - diameter / 2, gene_pos.y + diameter / 2),
            diameter, diameter, svg::Fill(color),
            stroke);
        break;
      }
      case Event::loss :
        drawCross(svg, gene_pos, diameter, diameter,
                  svg::Stroke(diameter / 3, color));
        break;
      case Event::speciationLoss :
      case Event::speciation :
        //svg << svg::Circle(gene_pos.toPoint(), diameter,
        // svg::Fill(color), stroke);
        break;
      case Event::extant : {
        //svg << svg::Circle(gene_pos.toPoint(), diameter,
        // svg::Fill(color), stroke);
        break;
      }
      default:
        break;
    }

    if(genetree.getRoot() == gene) {
      // If the node is the root, draw a star at its position.
      drawStar(svg, gene_pos, diameter * 1.25, 5, svg::Fill(svg::Color::White),
               svg::Stroke(diameter / 4, color));
    }
  }

  /*!
   * @brief Draw a species tree.
   * @param doc
   * @param speciestree
   * @param speciestree_layout_info
   */
  static void drawSpTree(
      svg::Document &doc
      , const bpp::PhyloTree &speciestree
      , const TreeLayoutInfo &speciestree_layout_info
  ) {
    const std::unordered_map<Node, Container> &speciestree_layout
        = speciestree_layout_info.species_containers;

    // First: createDocument the shape of the species tree.
    auto polygon = svg::Polygon(
        svg::Color::White, svg::Stroke(10., svg::Color(120, 120, 120)));

    _drawSpNode<svg::Polygon>(polygon, speciestree.getRoot(), speciestree,
                              speciestree_layout);

    doc << polygon;

    // Then add names of the species leaves.
    for (const auto &leaf: speciestree.getAllLeaves()) {
      if (leaf->hasName()) {
        const Container &leaf_container = speciestree_layout.at(leaf);
        auto leaf_name = leaf->getName();
        Coordinates name_pos;
        name_pos.x = leaf_container.stop.down.x - ((leaf_name.size()) * 13);
        name_pos.y = leaf_container.stop.down.y + 18;
        doc << svg::Text(name_pos.toPoint(), leaf_name,
                         svg::Fill(svg::Color::Black), svg::Font(24));
      }
    }
  }

  /*!
   * Draw a gene tree inside a species tree.
   * @param doc svg::Document to update.
   * @param genetree
   * @param layout Coordinates of the gene tree.
   * @param map
   * @param color
   */
  static void drawGnTree(
      svg::Document &doc
      , const std::shared_ptr<bpp::PhyloTree> &genetree
      , const TreeLayoutInfo &layout
      , const SpeciesGeneMap &map
      , const svg::Color &color
  ) {
    _drawGnNode(doc, genetree->getRoot(), *genetree,
                layout.genes_containers.at(genetree), map, color);
  }

  static void drawGlobalLegend(
      svg::Document &doc, const TreeLayoutInfo layout, const double o_x
      , const double o_y
  ) {
    unsigned long loss_number = 0;
    unsigned long dupl_number = 0;

    const auto &events = layout.events;
    for (const auto &species_gene_events: events) {
      for (const auto &gene_event: species_gene_events.second) {
        const Event &event = gene_event.second.first;
        if (event == Event::duplication) {
          dupl_number++;
        } else if (event == Event::loss) {
          loss_number++;
        }
      }
    }

    double x = o_x;
    double y = o_y;
    double line_height = 25;
    double char_width = 15;
    double space_length = 10;

    auto color = svg::Color(120, 120, 120);
    doc << svg::Text(svg::Point(x, y), "Losses", svg::Fill(color),
                     svg::Font(24));
    x += (strlen("Losses") * char_width) + space_length;
    drawCross(doc, {x, y + line_height / 3}, 20, 20, svg::Stroke(8, color));
    x += 20 + space_length;
    doc << svg::Text(svg::Point(x, y), std::to_string(loss_number),
                     svg::Fill(color), svg::Font(24));

    x = o_x;
    y += line_height;
    doc << svg::Text(svg::Point(x, y), "Duplications", svg::Fill(color),
                     svg::Font(24));
    x += (strlen("Duplications") * char_width) + space_length;
    doc << svg::Rectangle(svg::Point(x - 10, y + 10 + line_height / 3), 20, 20,
                          svg::Fill(color), svg::Stroke(3, color));
    x += 20 + space_length;
    doc << svg::Text(svg::Point(x, y), std::to_string(dupl_number),
                     svg::Fill(color), svg::Font(24));

  }

public:

  /*!
   * @brief Print a reconciliation in stream.
   * @tparam GeneTreeIterator
   * @param os Output stream.
   * @param speciestree
   * @param genetrees_begin
   * @param genetrees_end
   * @param layout See TreeLayoutInfo
   * @param width Width of the print
   * @param height Height of the print
   * @return
   */
  template<typename GeneTreeIterator>
  static std::ostream &print(
      std::ostream &os
      , const bpp::PhyloTree &speciestree
      , const GeneTreeIterator &genetrees_begin
      , const GeneTreeIterator &genetrees_end
      , const TreeLayoutInfo &layout
      , const unsigned long width = 50
      , const unsigned long height = 50
  ) {
    const std::unordered_map<Node, Container> &speciestree_layout = layout.species_containers;
    double x_scaler = ((double) width / (double) layout.sizeX);
    double y_scaler = ((double) height / (double) layout.sizeY);


    Grid<std::string> res(width + 5, height + 5,
                          " "); // All prints are going to be saved in a grid where each pos is a character

    //
    for (const auto &species_container: speciestree_layout) {
      Container container = species_container.second;
      Coordinates start_up;
      start_up.x = container.start.up.x * x_scaler;
      start_up.y = container.start.up.y * y_scaler;

      Coordinates start_down;
      start_down.x = container.start.down.x * x_scaler;
      start_down.y = container.start.down.y * y_scaler;

      Coordinates stop_up;
      stop_up.x = container.stop.up.x * x_scaler;
      stop_up.y = container.stop.up.y * y_scaler;

      Coordinates stop_down;
      stop_down.x = container.stop.down.x * x_scaler;
      stop_down.y = container.stop.down.y * y_scaler;

      Coordinates label;
      label.x = (stop_up.x + start_up.x) / 2.0;
      label.y = (start_down.y + start_up.y) / 2.0;
      res((std::size_t) label.x, (std::size_t) label.y) = std::toupper(
          species_container.first->getName().front());


      // Trace boxes which represents species tree
      addEdgeBetweenTwoPointsInGrid(res, stop_down, stop_up, '-');

      if (speciestree.isLeaf(species_container.first)) {
        addEdgeBetweenTwoPointsInGrid(res, start_down, stop_down, '|');
        addEdgeBetweenTwoPointsInGrid(res, stop_up, start_up, '|');
      }

      if (speciestree.hasFather(species_container.first)) { // If son left
        auto father = speciestree.getFather(species_container.first);
        const Container &father_container = speciestree_layout.at(father);
        if (speciestree.getSons(father).front() == species_container.first) {
          Coordinates temp_point(
              {father_container.start.up.x * x_scaler, start_up.y});
          addEdgeBetweenTwoPointsInGrid(res, start_up, temp_point, '|');
          addEdgeBetweenTwoPointsInGrid(res, temp_point,
                                        {father_container.start.up.x * x_scaler,
                                         father_container.start.up.y *
                                         y_scaler}, '-');

          addEdgeBetweenTwoPointsInGrid(res, start_down,
                                        {father_container.stop.up.x * x_scaler,
                                         start_down.y}, '|');
          addEdgeBetweenTwoPointsInGrid(res,
                                        {father_container.stop.up.x * x_scaler,
                                         start_down.y},
                                        {father_container.stop.up.x * x_scaler,
                                         father_container.stop.up.y * y_scaler},
                                        '-');
        } else { // If son right
          addEdgeBetweenTwoPointsInGrid(res, start_up,
                                        {father_container.stop.down.x *
                                         x_scaler,
                                         father_container.stop.down.y *
                                         y_scaler}, '|');
          Coordinates temp_point(
              {father_container.start.down.x * x_scaler, start_down.y});
          addEdgeBetweenTwoPointsInGrid(res, start_down, temp_point, '|');
          addEdgeBetweenTwoPointsInGrid(res, temp_point,
                                        {father_container.start.down.x *
                                         x_scaler,
                                         father_container.start.down.y *
                                         y_scaler}, '-');
        }
      }

    }

    for (auto genetree_it = genetrees_begin;
         genetree_it != genetrees_end; genetree_it++) {
      const bpp::PhyloTree &genetree = *(*genetree_it);

      const auto &genetree_layout = layout.genes_containers.at(*genetree_it);

      auto genes = genetree.getAllNodes();
      for (const auto &gene: genes) {
        const Container &gene_container = genetree_layout.at(gene);
        Coordinates gene_pos;
        gene_pos.x = gene_container.coordinates.x * x_scaler;
        gene_pos.y = gene_container.coordinates.y * y_scaler;

        if (not genetree.isLeaf(gene)) {
          auto sons = genetree.getSons(gene);
          for (auto son: sons) {
            auto edge_x = gene_pos.x;
            auto edge_y = gene_pos.y;

            const Container &son_container = genetree_layout.at(son);

            //auto son_x = son_container.coordinates.x * x_scaler;
            auto son_y = son_container.coordinates.y * y_scaler;

            if (edge_y < son_y) {
              for (edge_y = edge_y + 1; edge_y < son_y; edge_y++) {
                res((std::size_t) edge_x, (std::size_t) edge_y) = "·";
              }
            } else {
              for (edge_y = edge_y - 1; edge_y > son_y; edge_y--) {
                res((std::size_t) edge_x, (std::size_t) edge_y) = "·";
              }
            }
          }
        }

        if (genetree.hasFather(gene)) {
          // Draw a line to the father level
          auto edge_x = gene_pos.x;
          auto edge_y = gene_pos.y;
          auto parent = genetree.getFather(gene);
          const Container &parent_container = genetree_layout.at(parent);
          auto parent_x = parent_container.coordinates.x * x_scaler;

          if (edge_x < parent_x) {
            for (edge_x = edge_x + 1; edge_x < parent_x; edge_x++) {
              res((std::size_t) edge_x, (std::size_t) edge_y) = ":";
            }
          } else {
            for (edge_x = edge_x - 1; edge_x > parent_x; edge_x--) {
              res((std::size_t) edge_x, (std::size_t) edge_y) = ":";
            }
          }
        }
      }

      for (const auto &gene: genes) {
        const Container &gene_container = genetree_layout.at(gene);
        Coordinates gene_pos;
        gene_pos.x = gene_container.coordinates.x * x_scaler;
        gene_pos.y = gene_container.coordinates.y * y_scaler;
        res((std::size_t) gene_pos.x,
            (std::size_t) gene_pos.y) = PhyloTreeToolBox::isArtificalGene(gene)
                                        ? "x" : "o";
      }
    }

    os << res << std::endl;
    return os;
  }

  /*!
   * @brief Draw a reconciliation in svg::Document.
   * @return svg::Document
   * @tparam GeneTreesIterator Iterator to a gene trees container.
   * @tparam MapIterator Iterator to a SpeciesGeneMap container.
   * @param filename Name of the output SVG file.
   * @param genetrees_begin
   * @param genetrees_end
   * @param speciestree
   * @param map_begin
   * @param map_end
   * @param allows_progression_print
   * @param verbose
   */
  template<typename GeneTreesIterator, typename MapIterator>
  static svg::Document createDocument(
      const GeneTreesIterator &genetrees_begin
      , const GeneTreesIterator &genetrees_end
      , const bpp::PhyloTree &speciestree
      , const MapIterator &map_begin
      , const std::string &filename = ""
      , const bool allows_progression_print = true
      , const bool verbose = false
  ) {

    Timer<std::size_t> progression(
        (std::size_t) std::distance(genetrees_begin, genetrees_end));

    if (allows_progression_print)
      utils::progressionBar(std::cout, progression,
                            false, "Computing SVG Layout");

    auto layout = createLayout(speciestree, genetrees_begin, genetrees_end,
                               map_begin);

    if (verbose)
      print(std::cout, speciestree, genetrees_begin, genetrees_end, layout, 40,
            75);

    layout.sizeX += 50;
    layout.sizeY += 50; // Correction according to the last species name which is out of bounds.

    svg::Dimensions dimensions(layout.sizeX, layout.sizeY);
    svg::Document doc(filename,
                      svg::Layout(dimensions, svg::Layout::BottomLeft));

    // Background of the svg.
    //doc << svg::Rectangle(svg::Point(0, (layout.sizeY)), layout.sizeX + 100, layout.sizeY,
    //    svg::Fill(svg::Color(255, 255, 255)));

    if (allows_progression_print)
      utils::progressionBar(std::cout, progression, false, "Drawing speciestree in SVG");

    // Draw species tree
    drawSpTree(doc, speciestree, layout
    );

    // Draw legend
    /*
    const auto species_root = speciestree.getRoot();
    const auto species_son_left = speciestree.getSons(species_root).front();
    auto species_root_start_up_pos = layout.species_containers.at(speciestree.getRoot()).start.up;
    auto species_son_left_up_pos = layout.species_containers.at(species_son_left).start.up;
    drawGlobalLegend(
        doc, layout, species_root_start_up_pos.x + 10, species_son_left_up_pos.y - 50
    );
     */

    // Draw gene trees
    MapIterator map_it = map_begin;
    int original_color_red = 120;
    int original_color_green = 120;
    int original_color_blue = 120;
    int color_variation = 50;
    for (auto genetree_it = genetrees_begin;
         genetree_it != genetrees_end; genetree_it++) {
      auto red = original_color_red + (color_variation * (Random(-1.5, 1.5)));
      auto green = original_color_green +
          (color_variation * (Random(-1.5, 1.5)));
      auto blue = original_color_blue + (color_variation * (Random(-1.5, 1.5)));
      auto color = svg::Color(red, green, blue);
      drawGnTree(doc, *genetree_it, layout, *map_it, color);
      map_it++;
      progression.next();
      if (allows_progression_print)
        utils::progressionBar(std::cout, progression,
                              true, "Drawing genetrees in SVG");
    }
    return doc;
  }

  /*!
   * @brief Draw a reconciliation in SVG file.
   * @tparam GeneTreesIterator Iterator to a gene trees container.
   * @tparam MapIterator Iterator to a SpeciesGeneMap container.
   * @param filename Name of the output SVG file.
   * @param genetrees_begin
   * @param genetrees_end
   * @param speciestree
   * @param map_begin
   * @param map_end
   * @param allows_progression_print
   * @param verbose
   */
  template<typename GeneTreesIterator, typename MapIterator>
  static void write(
      const std::string &filename
      , const GeneTreesIterator &genetrees_begin
      , const GeneTreesIterator &genetrees_end
      , const bpp::PhyloTree &speciestree
      , const MapIterator &map_begin
  ) {
    svg::Document doc = createDocument(genetrees_begin, genetrees_end,
                                       speciestree, map_begin,
                                       filename);
    if (doc.save()) {
      //if(verbose) std::cout << "SVG file created in " << doc.toString() << std::endl;
    } else {
      std::cerr << "Error during SVG creation: file is not created."
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /* @brief Draw a reconciliation in std::ostream.
   * @tparam GeneTreesIterator Iterator to a gene trees container.
   * @tparam MapIterator Iterator to a SpeciesGeneMap container.
   * @param os std::ostream to append svg text.
   * @param genetrees_begin
   * @param genetrees_end
   * @param speciestree
   * @param map_begin
   * @param map_end
   * @param allows_progression_print
   * @param verbose
   */
  template<typename GeneTreesIterator, typename MapIterator>
  static std::ostream &write(
      std::ostream &os
      , const GeneTreesIterator &genetrees_begin
      , const GeneTreesIterator &genetrees_end
      , const bpp::PhyloTree &speciestree
      , const MapIterator &map_begin
  ) {
    svg::Document doc = createDocument(genetrees_begin, genetrees_end,
                                       speciestree, map_begin);
    os << doc.toString();
    return os;
  }

  /*!
   * @brief Draw a ReconciledGeneTrees in SVG file.
   * @tparam ReconciledGeneTreesIterator Iterator to a ReconciledGeneTrees container.
   * @param filename Name of the SVG output file.
   * @param genetrees_begin
   * @param genetrees_end
   * @param speciestree
   * @param allows_progression_print
   * @param verbose
   */
  template<typename ReconciledGeneTreesIterator>
  static void write(
      const std::string &filename
      , const ReconciledGeneTreesIterator &genetrees_begin
      , const ReconciledGeneTreesIterator &genetrees_end
      , const bpp::PhyloTree &speciestree
      , const bool allows_progression_print = true
      , const bool verbose = false
  ) {
    std::vector<std::shared_ptr<bpp::PhyloTree>> genetrees;
    genetrees.reserve(static_cast<unsigned long>(std::distance(genetrees_begin,
                                                               genetrees_end)));

    std::vector<SpeciesGeneMap> maps;
    maps.reserve(static_cast<unsigned long>(std::distance(genetrees_begin,
                                                          genetrees_end)));

    for (auto genetrees_it = genetrees_begin;
         genetrees_it != genetrees_end; genetrees_it++) {
      genetrees.push_back(genetrees_it->genetree());
      maps.push_back(genetrees_it->map());
    }

    RecPhyloTreeToSVG::write(filename, genetrees.begin(), genetrees_end,
                             speciestree,
                             maps.begin(), maps.end(), allows_progression_print,
                             verbose);
  }

  /*!
 * @brief Write a ReconciledGeneTrees in SVG format in std::ostream.
 * @tparam ReconciledGeneTreesIterator Iterator to a ReconciledGeneTrees container.
 * @param os std::ostream to update.
 * @param genetrees_begin
 * @param genetrees_end
 * @param speciestree
 * @param allows_progression_print
 * @param verbose
 */
  template<typename ReconciledGeneTreesIterator>
  static std::ostream &write(
      std::ostream &os
      , const ReconciledGeneTreesIterator &genetrees_begin
      , const ReconciledGeneTreesIterator &genetrees_end
      , const bpp::PhyloTree &speciestree
      , const bool allows_progression_print = true
      , const bool verbose = false
  ) {

    std::vector<std::shared_ptr<bpp::PhyloTree>> genetrees;
    genetrees.reserve(static_cast<unsigned long>(std::distance(genetrees_begin,
                                                               genetrees_end)));

    std::vector<SpeciesGeneMap> maps;
    maps.reserve(static_cast<unsigned long>(std::distance(genetrees_begin,
                                                          genetrees_end)));

    for (auto genetrees_it = genetrees_begin;
         genetrees_it != genetrees_end; genetrees_it++) {
      genetrees.push_back(genetrees_it->genetree());
      maps.push_back(genetrees_it->map());
    }

    return RecPhyloTreeToSVG::write(os, genetrees.begin(), genetrees_end,
                                    speciestree,
                                    maps.begin(), maps.end(),
                                    allows_progression_print, verbose);
  }

  /*!
   * @brief Print Treerecs reconciliations in SVG.
   * @param solutions Solutions generated by Treerecs in main.cpp.
   * @param speciestree
   * @param output_filename Name of the SVG output file.
   * @param allows_progresssion_print
   * @param verbose
   */
  static void write(
      const std::string &output_filename
      , const std::unordered_map<std::shared_ptr<bpp::PhyloTree>, std::map<double, std::vector<ReconciledRootedTree>>> &solutions
      , const bpp::PhyloTree &speciestree
  ) {
    std::list<std::shared_ptr<bpp::PhyloTree>> trees;
    std::list<SpeciesGeneMap> maps;

    for (const auto &solutions_per_genetree: solutions) {
      const auto &first_reconciled_tree = solutions_per_genetree.second.begin()->second.front();
      trees.push_back(first_reconciled_tree.genetree());
      maps.push_back(first_reconciled_tree.map());
    }

    RecPhyloTreeToSVG::write(output_filename, trees.begin(), trees.end(),
                             speciestree,
                             maps.begin());
  }

  /*!
 * @brief Print Treerecs reconciliations in std::ostream.
 * @param os std::ostream to modify.
 * @param solutions Solutions generated by Treerecs in main.cpp.
 * @param speciestree
 * @param allows_progresssion_print
 * @param verbose
 */
  static std::ostream &write(
      std::ostream &os
      , const std::unordered_map<std::shared_ptr<bpp::PhyloTree>, std::map<double, std::vector<ReconciledRootedTree>>> &solutions
      , const bpp::PhyloTree &speciestree
  ) {
    std::list<std::shared_ptr<bpp::PhyloTree>> trees;
    std::list<SpeciesGeneMap> maps;
    for (const auto &solutions_per_genetree: solutions) {
      const auto &first_reconciled_tree = solutions_per_genetree.second.begin()->second.front();
      trees.push_back(first_reconciled_tree.genetree());
      maps.push_back(first_reconciled_tree.map());
    }
    return RecPhyloTreeToSVG::write(os, trees.begin(), trees.end(), speciestree,
                                    maps.begin());
  }

};

} // namespace treerecs

#endif //TREERECS_PHYLOXMLTOSVG_H
