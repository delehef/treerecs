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

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <treerecs/tools/utils.h>
#include <treerecs/tools/IO/RefreshablePrinter.h>
#include <treerecs/tools/Timer.h>
#include <treerecs/tools/IO/IO.h>
#include <treerecs/tools/SpeciesGeneMapper.h>

using namespace treerecs;

static auto& preferred_separator = IO::preferred_separator; // until C++17

void print_header(){
  std::cout << "                                                    " << std::endl;
  std::cout << "              ---- GenetreeEditor ----             " << std::endl;
  std::cout << "    *An editor of gene trees with species trees.*   " << std::endl;
  std::cout << "                                                    " << std::endl;
}

void print_version(char* prog_path) {
  char* prog_name;
  if((prog_name = strrchr(prog_path, preferred_separator))) prog_name++;
  else prog_name = prog_path;
  std::cout << prog_name << " version 0.1" << std::endl;
  std::cout << "                 " << std::endl;
}

void print_usage(char* prog_path) {
  char* prog_name;
  if((prog_name = strrchr(prog_path, preferred_separator))) prog_name++;
  else prog_name = prog_path;
  std::cout << "Usage:   " << prog_name << " -h or --help.                                                       " << std::endl;
  std::cout << "   or:   " << prog_name << " -v or --version.                                                    " << std::endl;
  std::cout << "   or:   " << prog_name << " -g GENETREE_FILE -s SPECIESTREE_FILE [-m MAP_FILE]" << std::endl;
}

void print_utility(char* prog_path) {
  char* prog_name;
  if((prog_name = strrchr(prog_path, preferred_separator))) prog_name++;
  else prog_name = prog_path;
  std::size_t prog_name_size = 0;
  while(prog_name[prog_name_size++] != '\0') {}
  std::cout << prog_name << " provides tools to edit gene trees or create Smap." << std::endl;
  std::cout << std::string(" ", prog_name_size) << " generate genetrees in newick format by inserting associated species name after of before the gene name using a separator." << std::endl;
  std::cout << std::string(" ", prog_name_size) << " generate a map file in standard format." << std::endl;
}

void print_help(char* prog_path) {
  print_header();
  std::cout << "Authors: Nicolas COMTE.                                                                          " << std::endl;
  std::cout << "Team Inria Grenoble Rhône-Alpes Beagle (FRANCE).                                                 " << std::endl;
  std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;
  print_utility(prog_path);
  print_usage(prog_path);
  std::cout << "Options:                                                                                         " << std::endl;
  std::cout << "   -h, --help\n\tprint this help, then exit.\n                                                   " << std::endl;
  std::cout << "   -v, --version\n\tprint version number, then exit.\n                                           " << std::endl;
  std::cout << "   -V, --verbose\n\tverbose.\n                                                                   " << std::endl;
  std::cout << "   -g, --genetree\n\tgene tree filename in a Newick, NHX or PhyloXML format.\n                             " << std::endl;
  std::cout << "   -s, --speciestree\n\tgene tree filename in a Newick, NHX or PhyloXML format.\n                                 " << std::endl;
  std::cout << "   -S, --smap\n\tmap filename. File in the default SMap format.\n                                " << std::endl;
  std::cout << "   -N, --tree-index\n\treduce the selection of gene trees to one defined by its index in the file.\n"<< std::endl;
  std::cout << "   -o, --output\n\toutput filename (default = \"output.txt\").\n                                 " << std::endl;
  std::cout << "   -c, --sep\n\tseparator character in gene names to find the species name (default = '_').\n    " << std::endl;
  std::cout << "   -p, --prefix\n\tposition of the species_name (before the gene_name: yes (y), after no (n), default = yes).\n" << std::endl;
  std::cout << "   -m, --create-map\n\tcreate an smap.\n                                                         " << std::endl;
  std::cout << "   --new-separator\n\tnew separator for gene tree output (default = -c --sep default value)\n    " << std::endl;
  std::cout << "   --new-position\n\tposition of the species_name (prefix/postfix, default = prefix)\n           " << std::endl;
  std::cout << "   -r, --replace-char\n\ttakes a string of size 2. Replace the first character in gene and species names by the second one. Mutiples inputs is possible by separating them with ':'.\n" << std::endl;
  std::cout << "   -x, --multiply-bootstraps\n\tchange bootstrap values by multiplication.\n" << std::endl;
  std::cout << "   -X, --multiply-lengths\n\tchange length values by multiplication.\n" << std::endl;
  std::cout << "   --lengths-precision\n\tnumber of digits after the comma for lengths.\n" << std::endl;
  std::cout << "   --bootstraps-precision\n\tnumber of digits after the comma for bootstraps.\n" << std::endl;
  std::cout << "   --bootstraps-zero-to-one\n\tset bootstraps at zero as value to one.\n" << std::endl;
  std::cout << "   --case-sensitive\n\tuse a case-sensitive mapping.\n                                           " << std::endl;
  std::cout << "                                                                                                 " << std::endl;
  print_version(prog_path);
}

bool isYes(const char* optarg){
  return
      utils::string_comp(optarg, "true", false)
      or utils::string_comp(optarg, "yes", false)
      or utils::string_comp(optarg, "t", false)
      or utils::string_comp(optarg, "y", false);
}

int main(int argc, char * argv[]) {

  // Init objects, defaults
  // Gene tree and species tree.
  std::string genetrees_filename;
  std::vector<std::shared_ptr<bpp::PhyloTree>> genetrees;
  std::string speciestree_filename;
  std::shared_ptr<bpp::PhyloTree> speciestree = nullptr;

  unsigned int genetree_index; // option used, in case of multiple gene trees in genetree file, to select one tree.
  bool genetreeIndexSpecified = false;

  // Gene and species trees mapping.
  std::vector<SpeciesGeneMap> genemaps;
  MappingMethod mappingMethod = MappingMethod::TREES;
  // Mapping with mapping file.
  std::string map_filename; //filename of the map if specified.
  // Mapping with gene names.
  std::string sep = DEFAULT_SEPARATION_CHARACTER_BETWEEN_GENE_AND_SPECIES_NAME;
  bool prefix = false;


  //Outputs parameters
  auto new_separator = sep;
  auto new_position_as_prefix = prefix;
  bool change_separator = false;
  bool change_position = false;
  bool saveMaps = false;
  bool verbose = false;
  std::string output = "output.txt";
  bool allows_progression_print = true;
  TextFormat input_trees_text_format = TextFormat::newick;
  TextFormat output_trees_text_format = TextFormat::unknown;
  bool replace_char_in_names = false;
  std::vector<char> char_old = {'_'};
  std::vector<char> char_new = {'-'};
  double bootstrap_multiplicator = 1.0;
  double length_multiplicator = 1.0;
  int bootstrap_precision = -1;
  int length_precision = -1;
  bool bootstrap_zero_to_one = false;

  // Load options
  static struct option long_options_lists[] = {
      {"help",                  no_argument, NULL,        'h'},
      {"version",               no_argument, NULL,        'v'},
      {"verbose",               no_argument, NULL,        'V'},
      {"genetree",              required_argument, NULL,  'g'},
      {"speciestree",           required_argument, NULL,  's'},
      {"smap",                  required_argument, NULL,  'S'},
      {"sep",                   required_argument, NULL,  'c'},
      {"prefix",                required_argument, NULL,  'p'},
      {"output",                required_argument, NULL,  'o'},
      {"output-format",         required_argument, NULL,  'O'},
      {"tree-index",            required_argument, NULL,  'N'},
      {"create-map",            no_argument, NULL,        'm'},
      {"replace-char",          required_argument, NULL,  'r'},
      {"new-separator",         required_argument, NULL,   0 },
      {"new-position",          required_argument, NULL,   0 },
      {"case-sensitive",        no_argument, NULL,         0 },
      {"translate",             required_argument, NULL,   0 },
      {"multiply-bootstraps",   required_argument, NULL,  'x'},
      {"multiply-lengths",      required_argument, NULL,  'X'},
      {"bootstraps-precision",  required_argument, NULL,   0 },
      {"lengths-precision",     required_argument, NULL,   0 },
      {"bootstraps-zero-to-one",no_argument, NULL,   0 },
      {0,0,0,0}
  };

  int option;
  int option_index = 0;
  const char* options = "hvVmg:s:S:c:p:o:O:N:r:x:X:";

  if(argc == 1){
    // If there is no option command.
    print_usage(argv[0]);
    exit(EXIT_SUCCESS);
  }

  while(((option = getopt_long(argc, argv, options, long_options_lists, &option_index))) != -1){
    switch(option){
      case 0:
      { // We are going to check if the long option given exists
        if (long_options_lists[option_index].flag != 0)
          break;
        if (strcmp(long_options_lists[option_index].name, "translate") == 0){
          if(utils::string_comp("newick", optarg))
            output_trees_text_format = TextFormat::newick;
          else if(utils::string_comp("nhx", optarg))
            output_trees_text_format = TextFormat::nhx;
          else if(utils::string_comp("phyloxml", optarg))
            output_trees_text_format = TextFormat::phyloxml;
          else {
            std::cerr << "Unknown given format given with " << long_options_lists[option_index].name << " option" << std::endl;
            exit(EXIT_FAILURE);
          }
        }
        else if(strcmp(long_options_lists[option_index].name, "case-sensitive") == 0)
          SpeciesGeneMapper::set_case_sensitive_mapping();
        else if(strcmp(long_options_lists[option_index].name, "new-separator") == 0) {
          new_separator = optarg;
          change_separator = true;
        }
        else if(strcmp(long_options_lists[option_index].name, "new-position") == 0) {
          new_position_as_prefix = utils::string_comp(utils::trim(optarg), "prefix", false);
          change_position = true;
        }
        else if(strcmp(long_options_lists[option_index].name, "bootstraps-precision") == 0) {
          bootstrap_precision = atoi(optarg);
        }
        else if(strcmp(long_options_lists[option_index].name, "lengths-precision") == 0) {
          length_precision = atoi(optarg);
        }
        else if(strcmp(long_options_lists[option_index].name, "bootstraps-zero-to-one") == 0) {
          bootstrap_zero_to_one = true;
        }
        break;
      }
      case 'h':
      {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'v':
      {
        print_version(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V':
      {
        verbose = true;
        break;
      }
      case 'g':
      {
        // open genetree.
        genetrees_filename = std::string(optarg);
        break;
      }
      case 's':
      {
        // open speciestree.
        speciestree_filename = std::string(optarg);
        break;
      }
      case 'S':
      {
        // open map.
        map_filename = optarg;
        if (map_filename.find(".emf") != std::string::npos) {
          mappingMethod = MappingMethod::ENSEMBL;
        } else
          mappingMethod = MappingMethod::SMAP;
        break;
      }
      case 'o':
      {
        // save name
        output = optarg;
        break;
      }
      case 'O':
      {
        auto text_format_entries = utils::splitString(utils::trim(optarg), ":");
        for(auto& text_format_entry: text_format_entries) {
          if (utils::string_comp(text_format_entry, "newick", false)) {
            output_trees_text_format = TextFormat::newick;
          } else if (utils::string_comp(text_format_entry, "nhx", false)) {
            output_trees_text_format = TextFormat::nhx;
          } else if (utils::string_comp(text_format_entry, "phyloxml", false)) {
            output_trees_text_format = TextFormat::phyloxml;
          } else if (utils::string_comp(text_format_entry, "recphyloxml",
                                        false)) {
            output_trees_text_format = TextFormat::recphyloxml;
          } else if (utils::string_comp(text_format_entry, "svg",
                                        false)) {
            output_trees_text_format = TextFormat::svg;
          } else {
            std::cerr << "Error: \"" << text_format_entry
                      << "\" is not supported." << std::endl;
            std::cerr << "       please try:" << std::endl;
            for (auto txt_format = TextFormat::newick;
                 txt_format != TextFormat::unknown; ++txt_format) {
              std::cerr << "         * " << txt_format << std::endl;
            }
            exit(EXIT_FAILURE);
          }
        }
        break;
      }
      case 'c':
      {
        // Define separator char.
        sep = optarg;
        mappingMethod = MappingMethod::GENE_NAMES;
        break;
      }
      case 'p':
      {
        // Define name's position.
        prefix = isYes(optarg);
        mappingMethod = MappingMethod::GENE_NAMES;
        break;
      }
      case 'N': {
        genetree_index = atoi(optarg) - 1;
        genetreeIndexSpecified = true;
        break;
      }
      case 'm': {
        saveMaps = true;
        break;
      }
      case 'x': {
        bootstrap_multiplicator = atof(optarg);
        break;
      }
      case 'X': {
        length_multiplicator = atof(optarg);
        break;
      }
      case 'r': {
        replace_char_in_names = true;
        std::string inputs(optarg);
        auto linputs = utils::splitString(inputs, ":");
        for(auto& input: linputs) {
          if(input.size() != 2){
            std::cerr << "Error: two characters minimum "
                      << "needs to be provided with -r --replace-char (\""
                      << input << "\" given)." << std::endl;
          }
          char_old.push_back(input.front());
          char_new.push_back(input.back());
        }
        break;
      }
      default:
      {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
    }
  }

  if(not change_separator) new_separator = sep;
  if(not change_position) new_position_as_prefix = prefix;

  for(TextFormat format = TextFormat::newick; format != TextFormat::recphyloxml; ++format) {
    assert(format != TextFormat::unknown);
    try {
      genetrees = IO::readTreesFile(genetrees_filename, format,
                                    genetreeIndexSpecified ? (int) genetree_index : -1, true, true, allows_progression_print);
      input_trees_text_format = format;
      break;
    } catch (std::exception& ee) {}
  }

  if(output_trees_text_format == TextFormat::unknown)
    output_trees_text_format = input_trees_text_format;

  if(not speciestree_filename.empty()) {
    for (TextFormat format = TextFormat::newick;
         format != TextFormat::recphyloxml; ++format) {
      assert(format != TextFormat::unknown);
      try {
        speciestree = IO::readTreeFile(speciestree_filename, format, false,
                                       false);
        break;
      } catch (std::exception &ee) {}
    }


    // Mapping Genes to Species. There is three methods:
    //   - Reading a file which contains correspondences (SMap format usually);
    //   - Looking in gene names for the species names (with a separator and a given position);
    //   - Idem, but automatically looking the names in species tree and gene tree (default).
    Timer<std::size_t> mappingProgression(1);
    genemaps = SpeciesGeneMapper::map(
        genetrees.begin(), genetrees.end(), *speciestree, mappingMethod,
        map_filename, sep, prefix,
        allows_progression_print
    );
    mappingProgression.next();

    if (verbose and allows_progression_print) std::cout << std::endl;

    if (saveMaps) {
      if (verbose)
        std::cout << "Save map into " << DEFAULT_GLOBAL_MAP_OUTPUT_FILENAME
                  << "." << std::endl;
      SpeciesGeneMapper::save(DEFAULT_GLOBAL_MAP_OUTPUT_FILENAME, genemaps,
                              *speciestree, genetrees,
                              allows_progression_print);
    }
  }

  Timer<std::size_t> editing_progression(genetrees.size());

  for (genetree_index = 0; genetree_index < genetrees.size(); genetree_index++) {
    auto &genetree = genetrees.at(genetree_index);

    // Modify leaves
    for (auto gene_leaf: genetree->getAllLeaves()) {
      std::string gene_leaf_name = gene_leaf->getName();
      std::string associated_species_name;
      if(speciestree) {
         associated_species_name = genemaps.at(
            genetree_index).getAssociatedSpecies(gene_leaf)->getName();
      }

      // Create to copies of strings to find subtr in case of non-sensitive.
      std::string gene_leaf_name_copy = gene_leaf_name;
      std::string associated_species_name_copy = associated_species_name;

      std::transform(gene_leaf_name_copy.begin(), gene_leaf_name_copy.end(), gene_leaf_name_copy.begin(),
                       ::tolower);
      std::transform(associated_species_name_copy.begin(), associated_species_name_copy.end(),
                       associated_species_name_copy.begin(), ::tolower);

      unsigned long species_name_pos = std::string::npos;

      if(speciestree) {
        // Find species name in gene name.
        if (mappingMethod == MappingMethod::GENE_NAMES) {
          if (not prefix) {
            species_name_pos = gene_leaf_name_copy.find(sep) + 1;
          } else {
            //If prefix
            species_name_pos = 0;
          }
        } else {
          species_name_pos = gene_leaf_name_copy.find(
              associated_species_name_copy);
        }

        // Erase species name from gene name.
        if (species_name_pos != std::string::npos) {
          std::string old_separator_candidate = "";

          if (species_name_pos == 0) {
            // Species_name as prefix case.
            old_separator_candidate = gene_leaf_name.at(
                associated_species_name.size());
            bool old_separator_candidate_is_correct = (
                old_separator_candidate == "=" or
                old_separator_candidate == "_" or
                old_separator_candidate == "." or
                old_separator_candidate == "-" or
                old_separator_candidate == "|" or
                old_separator_candidate == "~");
            gene_leaf_name.erase(species_name_pos, species_name_pos +
                                                   associated_species_name.size() +
                                                   old_separator_candidate_is_correct);
          } else {
            // Species_name as postfix case.
            old_separator_candidate = gene_leaf_name.at(species_name_pos - 1);
            bool old_separator_candidate_is_correct = (
                old_separator_candidate == "=" or
                old_separator_candidate == "_" or
                old_separator_candidate == "." or
                old_separator_candidate == "-" or
                old_separator_candidate == "|" or
                old_separator_candidate == "~");
            gene_leaf_name.erase(
                species_name_pos - old_separator_candidate_is_correct,
                species_name_pos + associated_species_name.size());
          }
        }
      }

      // Edit gene name.
      if (replace_char_in_names) {
        for(std::size_t i = 0 ; i < char_old.size() ; i++) {
          associated_species_name =
              utils::replace(associated_species_name,
                             char_old.at(i), char_new.at(i));
          gene_leaf_name =
              utils::replace(gene_leaf_name,
                             char_old.at(i), char_new.at(i));
        }
      }

      std::string resulting_gene_leaf_name;

      if (new_position_as_prefix) {
        resulting_gene_leaf_name = (speciestree ?
                                   associated_species_name
                                   + new_separator : "")
                                   + gene_leaf_name;
      } else {
        resulting_gene_leaf_name = gene_leaf_name
                                    + (speciestree ?
                                   new_separator
                                   + associated_species_name : "");
      }

      gene_leaf->setName(resulting_gene_leaf_name);
    }

    // Modify edges
    for(auto edge: genetree->getAllEdges()) {
      if(edge->hasBootstrapValue()) {
        double bootstrap_value = edge->getBootstrapValue();
        if(not utils::double_equivalence(1.0, bootstrap_multiplicator)) {
          // Then change bootstrap by multiplying with value given
          // by the user.
          bootstrap_value *= bootstrap_multiplicator;
        }

        if(bootstrap_precision != -1) {
          double scale = pow(10, -1.0 * bootstrap_precision);
          bootstrap_value = (int)round(bootstrap_value / scale) * scale;
        }

        if(bootstrap_zero_to_one and (bootstrap_value == 0)) {
          bootstrap_value = 1;
        }

        edge->setProperty("bootstrap", bpp::Number<double>(bootstrap_value));
      }

      if(edge->hasLength()) {
        double length_value = edge->getLength();
        if(not utils::double_equivalence(1.0, bootstrap_multiplicator)) {
          // Then change length by multiplying with value given
          // by the user.
          length_value *= length_multiplicator;
        }

        if(length_precision != -1) {
          double scale = pow(10, -1.0 * length_precision);
          length_value = (int)round(length_value / scale) * scale;
        }

        edge->setLength(length_value);
      }
    }

    editing_progression.next();
    if (allows_progression_print)
      utils::progressionBar(std::cout, editing_progression,
                            true, "Editing gene trees");
  }

  if (allows_progression_print)
    RefreshablePrinter::clean(std::cout);

  // Write resulting genetrees in one file.
  std::ofstream output_genefile(output.c_str());

  if(output_trees_text_format == TextFormat::phyloxml) {
    output_genefile << "<phyloxml>" << std::endl;
  }

  for(auto genetree: genetrees) {
    IO::write(*genetree, output_genefile, output_trees_text_format, "");
  }

  if(output_trees_text_format == TextFormat::phyloxml) {
    output_genefile << "</phyloxml>" << std::endl;
  }

  output_genefile.close();

  return EXIT_SUCCESS;
}
