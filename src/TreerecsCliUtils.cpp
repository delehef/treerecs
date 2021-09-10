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
//

#include "TreerecsCliUtils.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef TREERECS_VERSION_NUMBER
#warning TREERECS_VERSION_NUMBER not defined, falling back to 0.0
#define TREERECS_VERSION_NUMBER "0.0"
#endif

using std::cerr;
using std::cout;
using std::endl;

static auto& preferred_separator = IO::preferred_separator; // until C++17

/**
 * Parse command line
 * \param argc
 * \param argv
 * \return
 */
std::unique_ptr<TreerecsParameters> ParseCommandLine(int argc, char** argv) {
  // Make sure we are reading argv from the beginning
  optind = 1;  // index of the next element to be processed in argv

  std::unique_ptr<TreerecsParameters> params(new TreerecsParameters());

  // Init option list (short and long)
  char options[50] = "hvVYCMqrg:s:a:D:S:d:l:c:p:o:O:ft:n:N:";
  #if defined(_OPENMP)
  strcat(options, "P");
  #endif

  static struct option options_long[] = {
      {"help",              no_argument,       NULL,  'h'},
      {"verbose",           no_argument,       NULL,  'v'},
      {"version",           no_argument,       NULL,  'V'},
      {"superverbose",      no_argument,       NULL,  'Y'},
      {"genetree",          required_argument, NULL,  'g'},
      {"speciestree",       required_argument, NULL,  's'},
      {"alignments",        required_argument, NULL,  'a'},
      {"dmatrix",           required_argument, NULL,  'D'},
      {"smap",              required_argument, NULL,  'S'},
      {"dupcost",           required_argument, NULL,  'd'},
      {"losscost",          required_argument, NULL,  'l'},
      {"threshold",         required_argument, NULL,  't'},
      {"sample-size",       required_argument, NULL,  'n'},
      {"sep",               required_argument, NULL,  'c'},
      {"prefix",            required_argument, NULL,  'p'},
      {"outdir",            required_argument, NULL,  'o'},
      {"output-format",     required_argument, NULL,  'O'},
      {"force",             required_argument, NULL,  'f'},
      {"reroot",            no_argument,       NULL,  'r'},
      {"costs-estimation",  no_argument,       NULL,  'C'},
      {"tree-index",        required_argument, NULL,  'N'},
  #if defined(_OPENMP)
      {"parallelize",       no_argument,       NULL,  'P'},
  #endif
      {"save-map",          no_argument,       NULL,  'M'},
      {"quiet",             no_argument,       NULL,  'q'},
      {"usage",             no_argument,       NULL,   0 },
      {"info",              no_argument,       NULL,   0 },
      {"ale-evaluation",    no_argument,       NULL,   0 },
      {"ale-selection",     no_argument,       NULL,   0 },
      {"case-sensitive",    no_argument,       NULL,   0 },
      {"stats",             no_argument,       NULL,   0 },
      {"fevent",            no_argument,       NULL,   0 },
      {"seed",              required_argument, NULL,   0 },
      {"output-without-description", no_argument, NULL,0 },
      {"number-of-threads", required_argument, NULL,   0 },
      {0, 0, 0, 0}
  };

  int option;
  int option_index = 0;

  if (argc == 1) {
    // If no option given
    throw malformed_cmd("mandatory options missing");
  }

  while (((option = getopt_long(argc, argv, options, options_long,
                                &option_index))) != -1) {
    switch (option) {
      case 'h': {
        params->functionalities.emplace(Functionalities::PRINT_HELP);
        return params;
      }
      case 'v': {
        if (params->verbosity != Verbosity::NORMAL)
          throw malformed_cmd("please use a single verbosity option");
        params->verbosity = Verbosity::VERBOSE;
        break;
      }
      case 'V': {
        params->functionalities.emplace(Functionalities::PRINT_VERSION);
        return params;
      }
      case 'g': {
        // open genetree.
        params->genetrees_filename = std::string(optarg);
        break;
      }
      case 's': {
        // open speciestree.
        params->speciestree_filename = std::string(optarg);
        break;
      }
      case 'a': {
        #ifndef _WIN32
        params->functionalities.emplace(Functionalities::PLL_LIKELIHOOD);
        params->alignments_filename = std::string(optarg);
        #else
        cerr << "Error: the alignments option is not supported on windows"
             << endl;
        exit(EXIT_FAILURE);
        #endif
        break;
      }
      case 'D': {
        params->functionalities.emplace(Functionalities::DMATRIX);
        params->dmatrix_filename = std::string(optarg);
        break;
      }
      case 'S': {
        // open map.
        params->map_filename = optarg;
        if (params->map_filename.find(".emf") != std::string::npos) {
          params->mappingMethod = treerecs::MappingMethod::ENSEMBL;
        } else
          params->mappingMethod = treerecs::MappingMethod::SMAP;
        break;
      }
      case 'd': {
        // set dupcost.
        params->dupcost = atof(optarg);
        break;
      }
      case 'l': {
        // set losscost.
        params->losscost = atof(optarg);
        break;
      }
      case 'C': {
        // change for duplication and loss cost estimation
        params->functionalities.emplace(Functionalities::COST_ESTIMATION);
        params->costs_estimation = true;
        break;
      }
      case 't': {
        // set threshold.
        std::string trim_optarg = treerecs::utils::trim(optarg);
        auto split_trim_optarg = treerecs::utils::splitString(trim_optarg, ":");
        params->thresholds_inputs.insert(
            params->thresholds_inputs.end(),
            split_trim_optarg.begin(), split_trim_optarg.end());
        break;
      }
      case 'o': {
        params->outdir = optarg;
        break;
      }
      case 'O': {
        auto text_format_entries = treerecs::utils::splitString(treerecs::utils::trim(optarg), ":");
        for (auto& text_format_entry: text_format_entries) {
          if (treerecs::utils::string_comp(text_format_entry, "newick", false)) {
            params->trees_output_formats.push_back(treerecs::TextFormat::newick);
          } else if (treerecs::utils::string_comp(text_format_entry, "nhx", false)) {
            params->trees_output_formats.push_back(treerecs::TextFormat::nhx);
          } else if (treerecs::utils::string_comp(text_format_entry, "phyloxml", false)) {
            params->trees_output_formats.push_back(treerecs::TextFormat::phyloxml);
          } else if (treerecs::utils::string_comp(text_format_entry, "recphyloxml",
                                        false)) {
            params->trees_output_formats.push_back(treerecs::TextFormat::recphyloxml);
          } else if (treerecs::utils::string_comp(text_format_entry, "svg",
                                        false)) {
            params->trees_output_formats.push_back(treerecs::TextFormat::svg);
          } else {
            cerr << "Error: \"" << text_format_entry
                 << "\" is not supported." << endl
                 << "       please try:" << endl;
            for (auto txt_format = treerecs::TextFormat::newick ;
                 txt_format != treerecs::TextFormat::unknown ; ++txt_format) {
              cerr << "         * " << txt_format << endl;
            }
            exit(EXIT_FAILURE);
          }
        }
        break;
      }
      case 'f': {
        params->force = true;
        break;
      }
      case 'c': {
        // Define separator char.
        params->sep = optarg;
        params->mappingMethod = treerecs::MappingMethod::GENE_NAMES;
        break;
      }
      case 'p': {
        // Define name position.
        params->prefix = treerecs::utils::isYes(optarg);
        params->mappingMethod = treerecs::MappingMethod::GENE_NAMES;
        break;
      }
      case 'Y': {
        if (params->verbosity != Verbosity::NORMAL)
          throw malformed_cmd("please use a single verbosity option");
        params->verbosity = Verbosity::CHATTY;
        break;
      }
      case 'n': {
        // Keep artificial genes created by profileNJ
        params->sample_size = atoi(optarg);
        break;
      }
      case 'r': {
        // Define name position.
        params->reroot = true;
        params->functionalities.emplace(Functionalities::REROOTING);
        break;
      }
      case 'N': {
        auto given_input = atoi(optarg);
        if (given_input < 1) {
          cerr << "Error: minimal gene tree index accepted is 1 ("
               << given_input << " is given)." << endl;
          exit(EXIT_FAILURE);
        }
        params->genetree_index = given_input - 1;
        break;
      }
      case 'P': {
        // Parallelize.
        #if defined(_OPENMP)
        if (params->nthreads == 1)
          params->nthreads = omp_get_num_procs();
        #endif
        break;
      }
      case 'M': {
        // Save mapping in file.
        params->saveMaps = true;
        break;
      }
      case 'q': {
        if (params->verbosity != Verbosity::NORMAL)
          throw malformed_cmd("please use a single verbosity option");
        params->verbosity = Verbosity::QUIET;
        break;
      }
      case 0: {  // Handle long-only options
        assert(options_long[option_index].flag == 0);
        if (strcmp(options_long[option_index].name, "usage") == 0) {
          params->functionalities.emplace(Functionalities::PRINT_USAGE);
          return params;
        } else if (strcmp(options_long[option_index].name, "info") == 0) {
          if (params->verbosity == Verbosity::QUIET) {
            throw malformed_cmd();
          } else {
            params->functionalities.emplace(
                Functionalities::PRINT_SUPPORT_DIAGRAM);
          }
        } else if (strcmp(options_long[option_index].name,
                          "case-sensitive") == 0) {
          treerecs::SpeciesGeneMapper::set_case_sensitive_mapping();
        } else if (strcmp(options_long[option_index].name,
                          "ale-evaluation") == 0) {
          params->add_ale_evaluation = true;
          params->functionalities.emplace(Functionalities::ALE_LIKELIHOOD);
        } else if (strcmp(options_long[option_index].name,
                          "ale-selection") == 0) {
          params->add_ale_evaluation = true;
          params->ale_selection = true;
          params->functionalities.emplace(Functionalities::ALE_SELECTION);
        } else if (strcmp(options_long[option_index].name, "stats") == 0) {
          params->statistics = true;
        } else if (strcmp(options_long[option_index].name, "seed") == 0) {
          params->seed = (long) atoi(optarg);
          params->seed_set = true;
        } else if (strcmp(options_long[option_index].name,
                          "output-without-description") == 0) {
          params->output_without_description = true;
        } else if (strcmp(options_long[option_index].name, "fevent") == 0) {
          params->write_gene_relationships_summary = true;
        } else if (strcmp(options_long[option_index].name,
                          "number-of-threads") == 0) {
          #if defined(_OPENMP)
          params->nthreads = atoi(optarg);
          #endif
        }
        break;
      }
      case '?':
      case ':': {
        throw malformed_cmd(""); // Let getopt issue the error message
      }
      default: {
        throw malformed_cmd();
      }
    }
  }

  if((not params->genetrees_filename.empty()) and
     (not params->speciestree_filename.empty())) {
    if (not params->thresholds_inputs.empty()) {
      params->functionalities.emplace(Functionalities::CORRECTION);
    } else {
      params->functionalities.emplace(Functionalities::RESOLVING);
    }
  }

  return params;
}

/**
 * Print Minimal header
 */
void PrintMinimalHeader() {
  cout << "Treerecs (" << TREERECS_VERSION_NUMBER << "), Inria - Beagle"
       << endl;
}

/**
 * Print Treerecs current version
 */
void PrintVersion() {
  cout << "Treerecs " << TREERECS_VERSION_NUMBER << endl;
}

/**
 * Print Treerecs usage
 */
void PrintUsage(char* prog_path) {
  char* prog_name;
  if ((prog_name = strrchr(prog_path, preferred_separator))) prog_name++;
  else prog_name = prog_path;
  cout << "Usage:   " << prog_name
       << " -h | --help" << endl
       << "   or:   " << prog_name
       << " -V | --version" << endl
       << "   or:   " << prog_name
       << " --usage" << endl
       << "   or:   " << prog_name
       << " -g GENETREE_FILE -s SPECIESTREE_FILE" << endl
       << "                  [-S MAP_FILE] [-t BRANCH_SUPPORT_THRESHOLD]"
       << " [...]" << endl
       << "   or:   " << prog_name
       << " -g GENETREE_FILE" << " --info" << endl;
}

/**
 * Print Treerecs help
 * \param prog_path
 */
void PrintHelp(char* prog_path) {
  cout
      << "--------------------------------------------------"
      << endl;

  PrintUsage(prog_path);

  cout << endl
       << "Options:" << endl;

  cout << "   -h, --help"
    << endl <<
    "\tprint this help, then exit."
    << endl << endl;

  cout << "   --usage"
    << endl <<
    "\tprint usage, then exit."
    << endl << endl;

  cout << "   -V, --version"
    << endl <<
    "\tprint version number, then exit."
    << endl << endl;

  cout << "   -v, --verbose"
    << endl <<
    "\tverbose mode. Causes Treerecs to print messages about its progress."
    << endl << endl;

  cout << "   -Y, --superverbose"
    << endl <<
    "\tsuper-verbose mode. Print even more messages than in verbose mode."
    << endl << endl;

  cout << "   -g, --genetree GENETREE_FILE"
    << endl <<
    "\tinput gene tree(s) (supported formats: Newick, NHX or PhyloXML)."
    << endl << endl;

  cout << "   -s, --speciestree SPECIESTREE_FILE"
    << endl <<
    "\tinput species tree (supported formats: Newick, NHX or PhyloXML)."
    << endl << endl;

  cout << "   -a, --alignments ALIGNMENTS_FILE"
    << endl <<
    "\tinput alignment file. Must contain:"
    << endl <<
    "\t  * the pll substitution model to use"
    << endl <<
    "\t  * the paths to the multiple alignments (one per gene-tree)"
    << endl << endl;

  cout << "   -D, --dmatrix MATRIX_FILE"
    << endl <<
    "\tinput distance matrix file. Will override brnach lengths as they are specified in the gene tree."
    << endl << endl;

  cout << "   -S, --smap SMAP_FILE"
    << endl <<
    "\tinput gene-to-species mapping file."
    << endl << endl;

  cout << "   -r, --reroot"
    << endl <<
    "\tfind the best root according to the reconciliation cost."
    << endl << endl;

  cout << "   -d, --dupcost VALUE"
    << endl <<
    "\tspecify gene duplication cost (default value = " <<
    DEFAULT_DUPLICATION_COST << ")."
    << endl << endl;

  cout << "   -l, --losscost VALUE"
    << endl <<
    "\tspecify gene loss cost (default value = " << DEFAULT_LOSS_COST << ")."
    << endl << endl;

  cout << "   -t, --threshold EXPRESSION | quantiles[N]"
    << endl <<
    "\tspecify branch support thresholds to use while contracting gene trees."
    << endl << endl <<
    "\t* EXPRESSION can be any colon-separated combination of the following:"
    << endl <<
    "\t  none: no contraction"
    << endl <<
    "\t  all: contract all branches. The tree collapses into a single polytomy"
    << endl <<
    "\t  VALUE: contract branches with support strictly lower than VALUE"
    << endl <<
    "\t  VALUE+epsilon (short VALUE+e): contract branches with support equal to"
    << endl <<
    "\t  or lower than VALUE"
    << endl << endl <<
    "\t* quantiles[N]: use several threshold values: none, all, and the"
    << endl <<
    "\t  quantiles dividing the branch supports into N groups"
    << endl << endl;

  cout << "   -n, --sample-size VALUE"
    << endl <<
    "\tsize of the tree sampling (default value = 1)."
    << endl << endl;

  cout << "   -N, --tree-index VALUE"
    << endl <<
    "\tonly consider the VALUE-th gene tree in the gene tree file."
    << endl << endl;

  cout << "   -o, --outdir OUTPUT_DIR"
    << endl <<
    "\toutput directory (default: treerecs_output)."
    << endl << endl;

  cout << "   -O, --output-format FORMAT"
    << endl <<
    "\toutput format(s): newick(default), nhx, phyloxml, recphyloxml or svg."
    << endl <<
    "\trepeat option or use a colon-separated list of formats to get multiple"
    << endl <<
    "\toutput"
    << endl << endl;

  cout << "   -f, --force"
    << endl <<
    "\tforce possible overwrite of existing files."
    << endl << endl;

  cout << "   -c, --sep CHARACTER"
    << endl <<
    "\tspecify separator character for species names embedded in gene names"
    << endl <<
    "\t(default = '"
    << DEFAULT_SEPARATION_CHARACTER_BETWEEN_GENE_AND_SPECIES_NAME << "')."
      << endl << endl;

  cout << "   -p, --prefix Y/N"
    << endl <<
    "\tspecify whether the species_name is a prefix of gene_name"
    << endl <<
    "\tdefault = "
    << (DEFAULT_SPECIES_NAME_IS_PREFIX_IN_GENE_NAME ? "Y" : "N")
    << ")."
    << endl << endl;

#if defined(_OPENMP)
  cout << "   -P, --parallelize"
    << endl <<
    "\trun in parallel if possible."
    << endl << endl;
#endif

  cout << "   -M, --save-map"
    << endl <<
    "\tsave map(s) used during execution."
    << endl << endl;

  cout << "   -q, --quiet"
    << endl <<
    "\tsilent mode (no print, no progression bar)."
    << endl << endl;

  cout << "   -C, --costs-estimation"
    << endl <<
    "\testimate duplication and loss costs."
    << endl << endl;

  cout << "   --info"
    << endl <<
    "\tprint informations about genetree(s) with a branch support diagram."
    << endl << endl;

  cout << "   --case-sensitive"
    << endl <<
    "\tuse case sensitive mapping."
    << endl << endl;

  cout << "   --ale-evaluation"
  << endl <<
  "\tcompute ALE log likelihood for each solution."
  << endl << endl;

  cout << "   --output-without-description"
    << endl <<
    "\tstrip output from gene tree descriptions."
    << endl << endl;

//  cout << "   --ale-selection"
//    << endl <<
//    "\tkeep solutions maximizing the ALE log likelihood."
//    << endl << endl;

  cout << "   --fevent"
    << endl <<
    "\tcreate a file that summarizes orthology/paralogy relationships."
    << endl << endl;

//  cout << "   --stats"
//    << endl <<
//    "\tprint in file statistics of resulting trees."
//    << endl << endl;

  cout << endl;
  PrintVersion();
}

/**
 * @brief Summarize all inputs and outputs and so, the behaviour of Treerecs for a run.
 * @param os Output stream to change.
 * @param params treerecs parameters
 * @param genetrees gene trees as a vector of ptr to bpp::PhyloTrees.
 * @param sorted_supports Branch supports, sorted.
 * @param supportThresholdsToUse all thresholds used.
 * @param output_map_filename output map filename if asked.
 * @param output_statistics_filename statistics output filename.
 */
void SummarizeParameters(
    std::ostream& os,
    const std::unique_ptr<TreerecsParameters>& params,
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& genetrees,
    const std::vector<double>& sorted_supports,
    const std::vector<double>& supportThresholdsToUse,
    const std::string& output_map_filename,
    const std::string& gene_relationships_summary_filename) {
  os  << "-Inputs--------------------------------------------------------------------------------------------"
      << endl;

  std::string treerecs_mode;
  // Treerecs modes :
  // * Resolve: solve each polytomy without gene tree(s) contraction.
  // * Rearrange: solve each gene tree with contractions in gene tree(s).
  // * Rooting: find best root in un-rooted gene tree(s).
  // * Re-rooting: find best root in rooted gene tree(s).
  // * Reconcile: contructs tree(s) with all gene events (shows gene losses and duplications).
  if (sorted_supports.size() > 0 and
      supportThresholdsToUse.back() < sorted_supports.front())
    treerecs_mode = "resolve";
  else
    treerecs_mode = "rearrange";

  if (params->reroot) treerecs_mode += " and root";

  // Trees inputs
  //   gene tree(s)
  if (genetrees.size() == 1) {
    os << "> 1 tree to " << treerecs_mode << " with "
       << genetrees.front()->getNumberOfNodes() << " nodes." << endl;
  } else
    os << "> " << genetrees.size() << " trees to " << treerecs_mode
       << " and correct." << endl;

  //   species tree
  os << "> Speciestree : " << params->speciestree_filename << endl;

  // Alignments
  if (params->alignments_filename.size()) {
    os << "> Alignments file: " << params->alignments_filename << endl;
  }

  // Distance matrix
  if (params->dmatrix_filename.size()) {
    os << "> Distance matrix file: " << params->dmatrix_filename << endl;
  }

  // Mapping
  os << "> Mapping according to ";
  if (params->mappingMethod == treerecs::MappingMethod::ENSEMBL) {
    os << "an ensembl file map: " << params->map_filename << ".";
  } else if (params->mappingMethod == treerecs::MappingMethod::SMAP) {
    os << "smap file map: " << params->map_filename << ".";
  } else if (params->mappingMethod == treerecs::MappingMethod::GENE_NAMES) {
    os << "gene names with the separator " << params->sep
       << ((params->prefix) ? "after" : "before") << " the species name.";
  } else if (params->mappingMethod == treerecs::MappingMethod::NHX_TAGS) {
    os << "NHX tags.";
  } else {
    assert(params->mappingMethod == treerecs::MappingMethod::TREES);
    os << "gene names with species names likeness.";
  }
  os << endl;

  // Threshold choice
  if (sorted_supports.size() == 0)
    os << "> No branch support/ bootstrap found in gene tree(s)." << endl;
  else if (supportThresholdsToUse.size() == 1 and
           treerecs::utils::double_equal_or_superior(supportThresholdsToUse.back(),
                                           sorted_supports.front()))
    os << "> Contraction of branches with a support less than or equal to "
       << supportThresholdsToUse.front() << "." << endl;
  else if (supportThresholdsToUse.size() == 1 and
           supportThresholdsToUse.back() < sorted_supports.front()) {
    os << "> No contraction with a threshold equal to: "
       << supportThresholdsToUse.back();
    os << ". Minimal branch support: " << sorted_supports.front();
    os << endl;
  } else if (supportThresholdsToUse.size() > 1) {
    os << "> Working with " << supportThresholdsToUse.size()
       << " thresholds of branch supports:" << endl;
    os << "  ";
    treerecs::utils::write(os, supportThresholdsToUse.begin(),
                 supportThresholdsToUse.end(), "", ", ", "");
    os << endl;
  }

  //Other options
  if (params->seed_set)
    os << "> Seed set to " << params->seed << "." << endl;

  if (params->ale_selection)
    os << "> Selection of solutions with ALE log likelihood." << endl;

  if (params->reroot)
    os << "> Re-rooting mode activated." << endl;

  if (params->costs_estimation)
    os << "> Duplication and loss costs estimations activated." << endl;

  if (params->sample_size > 0)
    os << "> Sampling (of size " << params->sample_size << " per tree)." << endl;

  if (params->nthreads > 1)
    os << "> Parallelization activated with " << params->nthreads << " threads."
       << endl;

  os << "> Output(s):" << endl;

  if (params->saveMaps)
    os << "  * An Smap in " << output_map_filename << "." << endl;

  os << "  * Reconstructed tree(s) in " << params->output << endl;

  os << "    in " << params->trees_output_formats << " format." << endl;

  if (params->output_without_description)
    os << "    without description." << endl;

  if (params->write_gene_relationships_summary)
    os << "  * Genes relationships summary in "
       << gene_relationships_summary_filename << endl;


  os  << "---------------------------------------------------------------------------------------------------"
      << endl;
  os << endl;
}
