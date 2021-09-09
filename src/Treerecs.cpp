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

#include <sys/stat.h>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
#include <algorithm>

#if defined(_OPENMP)
  #include <omp.h>
#endif

#include <TreerecsCliUtils.h>
#include <treerecs/tools/SpeciesGeneMapper.h>
#include <treerecs/tools/SupportTools.h>
#include <treerecs/tools/treeReconciliation/TreeReconciliationConductor.h>
#include <treerecs/tools/Statistics.h>
#include <treerecs/tools/IO/IO.h>
#include <treerecs/tools/random.h>
#include <treerecs/tools/Timer.h>

using std::cerr;
using std::cout;
using std::endl;
using std::list;
using std::string;
using namespace treerecs;

void CheckParametersConsistency(
    const std::unique_ptr<TreerecsParameters>& params);
void CheckOutputDir(const std::unique_ptr<TreerecsParameters>& params);
void CheckBranchSupports(
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& genetrees,
    std::vector<std::vector<double>> branch_supports_vector,
    const std::unique_ptr<TreerecsParameters>& params
);

int main(int argc, char* argv[]) {
  // Initialize timer used to compute the total execution time of Treerecs.
  Timer<std::size_t> execution_progression(1);

  // Get parameters from command line arguments
  std::unique_ptr<TreerecsParameters> params;
  try {
    params = ParseCommandLine(argc, argv);
  } catch (malformed_cmd& e) {
    if (strcmp(e.what(), "") != 0) {
      std::cerr << argv[0] << ": " << e.what() << std::endl;
    }
    PrintUsage(argv[0]);
    exit(EXIT_FAILURE);
  }

  // Print treerecs header
  if (params->verbosity >= Verbosity::NORMAL)
    PrintMinimalHeader();

  // Handle trivial requests
  if (params->functionalities.find(Functionalities::PRINT_HELP) !=
      params->functionalities.end()) {
    PrintHelp(argv[0]);
    exit(EXIT_SUCCESS);
  }
  if (params->functionalities.find(Functionalities::PRINT_USAGE) !=
      params->functionalities.end()) {
    PrintUsage(argv[0]);
    exit(EXIT_SUCCESS);
  }
  if (params->functionalities.find(Functionalities::PRINT_VERSION) !=
      params->functionalities.end()) {
    PrintVersion();
    exit(EXIT_SUCCESS);
  }

  CheckParametersConsistency(params);
  CheckOutputDir(params);

  // Random seed
  std::srand(params->seed);
  default_random_generator.seed(params->seed);

  // Check that the provided gene tree is accessible
  if (not IO::exists(params->genetrees_filename)) {
    std::cerr << "Error: gene tree(s) file \"" << params->genetrees_filename
              << "\" does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Check that the provided species tree file exists
  if (not params->speciestree_filename.empty() and
      not IO::exists(params->speciestree_filename)) {
    std::cerr << "Error: species tree file \"" << params->speciestree_filename
              << "\" does not exist." << std::endl;

    exit(EXIT_FAILURE);
  }

  // Check that the provided mapping file exists
  if ((params->mappingMethod == MappingMethod::SMAP or
       params->mappingMethod == MappingMethod::ENSEMBL) and
      not(IO::exists(params->map_filename))) {
    std::cerr << "Error: mapping file \"" << params->map_filename
              << "\" does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Check that the provided alignment file exists
  if (params->alignments_filename.size() and
      not IO::exists(params->alignments_filename)) {
    std::cerr << "Error: " << params->alignments_filename
              << " does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }

#if defined(_OPENMP) && !defined(__clang__)
  // Define the number of threads to use
  // With no parallelization, nthreads is set to 1.
  omp_set_num_threads(params->nthreads);
#endif

  // Set output format to default if the user didn't ask for any
  if (params->trees_output_formats.empty()) {
    params->trees_output_formats.push_back(TextFormat::newick);
  }

  // Open gene trees and try different formats
  std::vector<std::shared_ptr<bpp::PhyloTree>> genetrees;
  TextFormat genetrees_input_format = TextFormat::unknown;
  for (TextFormat format = TextFormat::newick ;
       format != TextFormat::recphyloxml ; ++format) {
    assert(format != TextFormat::unknown);
    try {
      genetrees = IO::readTreesFile(params->genetrees_filename, format,
                                    params->genetree_index, true, true,
                                    params->verbosity >= Verbosity::NORMAL,
                                    params->verbosity >= Verbosity::CHATTY);
      genetrees_input_format = format;
      break;
    } catch (std::exception& ee) {}
  }

  if (genetrees.size() <= 0) {
    std::cerr <<
        "Error: failed to read gene tree from " <<
        params->genetrees_filename << std::endl <<
        "       Please check conformity to one of the supported formats " <<
        std::endl <<
        "       (Newick, NHX or PhyloXML)" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Compute total number of gene nodes.
  std::size_t total_number_of_nodes = 0;
  for (auto& tree: genetrees) {
    total_number_of_nodes += tree->getNumberOfNodes();
  }

  // Opening alignments
  std::vector<LibpllAlignmentInfo> alignmentsInfo;
  if (params->alignments_filename.size()) {
    try {
      alignmentsInfo = LibpllEvaluation::parseAlignmentInfo(
          params->alignments_filename,
          params->genetree_index);
    }
    catch (LibpllException& e) {
      std::cerr << argv[0] << ": " << e.what()
                << ", please check your alignment file" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (alignmentsInfo.size() != genetrees.size()) {
      std::cerr << argv[0] << ": " << "the number of alignments ("
                << alignmentsInfo.size()
                << ") should be equal to the number of input gene trees ("
                << genetrees.size() << ")" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Edit output file names.
  if (params->output == "") {
    // Extract the input file name (without the path) and add postfix.
    params->output = params->outdir + IO::preferred_separator +
                     utils::extractFilename(params->genetrees_filename) +
                     "_" + OUTPUT_POSTFIX;
  }

  std::string output_map_filename =
      params->output +
      (strcmp(DEFAULT_GLOBAL_MAP_OUTPUT_FILENAME, "") == 0 ? "" : ".") +
      DEFAULT_GLOBAL_MAP_OUTPUT_FILENAME + "." + SMAP_EXTENTION;

  std::string output_statistics_filename = params->output
                                           + (strcmp(
      DEFAULT_STATISTICS_OUTPUT_FILENAME, "") == 0 ? "" : ".")
                                           + DEFAULT_STATISTICS_OUTPUT_FILENAME
                                           + "."
                                           + STATS_EXTENTION;

  std::string output_gene_relationships_summary_filename = params->output +
      (strcmp(DEFAULT_FEVENT_OUTPUT_FILENAME, "") == 0 ? "" : ".") +
      DEFAULT_FEVENT_OUTPUT_FILENAME + "." + FEVENT_EXTENTION;

  bool add_all_gene_losses_during_computations = false;
  for (auto& txt_format : params->trees_output_formats) {
    if (txt_format == TextFormat::recphyloxml
        or txt_format == TextFormat::svg) {
      add_all_gene_losses_during_computations = true;
      break;
    }
  }

  // Get the branch support for each of the provided gene trees
  auto branch_supports = GetBranchSupports(genetrees);
  CheckBranchSupports(genetrees, branch_supports, params);
  // Create a global container containing all the support values from all the
  // trees sorted by ascending value
  decltype(branch_supports)::value_type sorted_supports;
  for (const auto& supports : branch_supports) {
    sorted_supports.insert(sorted_supports.end(),
                           supports.begin(),
                           supports.end());
  }
  std::sort(sorted_supports.begin(), sorted_supports.end());

  // Print support diagram if requested
  if (params->functionalities.find(Functionalities::PRINT_SUPPORT_DIAGRAM) !=
      params->functionalities.end()) {
    PrintSupportDiagram(sorted_supports);
  }
  // If no species tree was provided, we are done => Exit
  if (params->speciestree_filename.empty()) {
    // If no species tree was provided, the info flag must have been set or
    // else we should have detected the error earlier
    assert(params->functionalities.find(
        Functionalities::PRINT_SUPPORT_DIAGRAM) !=
           params->functionalities.end());
    exit(EXIT_SUCCESS);
  }

  // Interpret command line threshold option argument
  decltype(ConvertThresholdInputsToDoubles(params, sorted_supports))
      supportThresholdsToUse;
  try {
    supportThresholdsToUse =
        ConvertThresholdInputsToDoubles(params, sorted_supports);
  }
  catch (malformed_cmd& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  // Read species tree
  auto speciestree =
      IO::readTreeFile(params->speciestree_filename, false, false,
                       params->verbosity >= Verbosity::CHATTY);

  // Check that the input species tree is fully resolved
  auto speciestree_multifurcations
      = PhyloTreeToolBox::findPolytomies(*speciestree);

  if (not speciestree_multifurcations.empty()) {
    std::cerr << "Error: there is " << speciestree_multifurcations.size()
              << " multifurcations found in the species tree."
              << std::endl
              << "       Please give a fully resolved species tree."
              << std::endl
              << "       Treerecs currently does not handle species trees "
                 "with multifurcations."
              << std::endl;

    exit(EXIT_FAILURE);
  }

  // Gene <-> species mapping
  Timer<std::size_t> mappingProgression(1);
  auto genemaps =
      SpeciesGeneMapper::GenerateMapping(genetrees_input_format,
                                         params->mappingMethod,
                                         params->map_filename,
                                         params->sep,
                                         params->prefix,
                                         genetrees,
                                         speciestree,
                                         params->verbosity >= Verbosity::NORMAL,
                                         mappingProgression);
  if (params->saveMaps)
    SpeciesGeneMapper::save(output_map_filename, genemaps, *speciestree,
                            genetrees, params->verbosity >= Verbosity::NORMAL);

  // Summarize all infos for users in at least Verbosity::VERBOSE mode
  // before reconciliations and reconstructions.
  if (params->verbosity >= Verbosity::VERBOSE) {
    SummarizeParameters(std::cout,
                        params,
                        genetrees,
                        sorted_supports,
                        supportThresholdsToUse,
                        output_map_filename,
                        output_gene_relationships_summary_filename);
  }

  // Container for all the solutions generated by Treerecs.
  std::unordered_map<
      std::shared_ptr<bpp::PhyloTree>,
      std::map<double, std::vector<ReconciledRootedTree>>> trees_solutions;
  // two keys: the first is a ptr to the original bpp::PhyloTree and
  // the second the threshold used to contract the tree.

  Timer<unsigned int> genetrees_progression(
      (unsigned int) (total_number_of_nodes * supportThresholdsToUse.size()));
  if (params->verbosity >= Verbosity::NORMAL)
    utils::progressionBar(std::cout, genetrees_progression, true,
                          "Reconciling genetrees");

  // Define the number of threads used to solve trees.
  int gtrees_nthreads =
      (params->nthreads > 1 and genetrees.size() > 1) ? params->nthreads : 1;

  // Define the number of threads in the case of multiple thresholds.
  int thresholds_exploration_nthreads =
      (supportThresholdsToUse.size() == 1) ? 1 : params->nthreads;

  // Finally, run rearrangements/reconciliations/rerootings.
  #if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic) num_threads(gtrees_nthreads)
  #endif
  for (std::size_t gtree_index = 0 ;
       gtree_index < genetrees.size() ; ++gtree_index) {
  #if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic) num_threads(thresholds_exploration_nthreads)
  #endif
    for (std::size_t supportThreshold_index = 0 ;
         supportThreshold_index < supportThresholdsToUse.size() ;
         ++supportThreshold_index) {
      auto supportThreshold = supportThresholdsToUse.at(supportThreshold_index);
      if (genetrees.size() > 1 and params->verbosity >= Verbosity::NORMAL) {
        utils::progressionBar(std::cout, genetrees_progression,
                              true, std::string("Reconciling genetree ")
                                    + std::to_string(gtree_index + 1));
      }

      if (params->verbosity >= Verbosity::CHATTY)
        std::cout << "---------------------------------------------------------------------------------------------------"
                  << std::endl;

      bool printProgression =
          genetrees.size() == 1 and supportThresholdsToUse.size() == 1;

      // Get the gene tree
      auto& genetree = genetrees.at(gtree_index);

      // Check if the tree is binary and un-rooted (compute the number of multifurcations)
      auto polytomies = PhyloTreeToolBox::findPolytomies(*genetree);

      bool geneTreeIsBinary =
          polytomies.size() == 0; // A binary gene tree has no polytomy.

      auto* alignment = alignmentsInfo.size() ? &alignmentsInfo[gtree_index]
                                              : nullptr;

      // But it could be unrooted...
      if (polytomies.size() == 1) {
        // If there is only one polytomy at the root with 3 sons the tree is unrooted and binary
        if (polytomies.front() == genetree->getRoot() and
            genetree->getSons(polytomies.front()).size() == 3)
          geneTreeIsBinary = true;
      }

      bool root_tree = false;

      // Apply tree reconciliation.
      if (not params->reroot and
          (geneTreeIsBinary and
           genetree->getSons(genetree->getRoot()).size() == 3)) {
        RefreshablePrinter::clean(std::cout);
        std::cout << "> Note: gene tree " << gtree_index + 1 << " is un-rooted."
                  << std::endl;
        std::cout << "  ...the tree is going to be rooted." << std::endl;

        root_tree = true;
      }

      TreeReconciliationConductor treesReconciliationConductor; // conducts all operations of rearrangments/reconciliation/rerooting
      auto solutions_for_one_tree_and_one_threshold =
          treesReconciliationConductor(
              *genetrees.at(
                  gtree_index) // original gene tree to reconcile/reconstuct.
              , *speciestree // speciestree as guide tree.
              , genemaps.at(
                  gtree_index) // map to link species nodes with gene nodes.
              , alignment // alignment to compute the libpll likelihood
              , params->costs_estimation,
              params->dupcost // cost of one duplication event.
              , params->losscost // cost of one loss event.
              ,
              supportThreshold // branch supports threshold to contract gene tree.
              ,
              params->strict_support_thresholds // contraction branch behaviour.
              , params->reroot or root_tree // re-root tree or not.
              , params->sample_size // number of outputs.
              ,
              params->add_ale_evaluation // add ALE log-likelihood evaluation in output.
              ,
              params->ale_selection // select best solution according to profileNJ cost and ALE log-likelihood.
              ,
              add_all_gene_losses_during_computations // add all artificial genes in gene trees.
              , params->compute_distances // compute new branch lengths.
              , printProgression and params->verbosity >=
                                     Verbosity::NORMAL // print progression bar in standard output.
              , params->verbosity >=
                Verbosity::CHATTY // print all operations in standard output.
          );

      if (params->costs_estimation && params->verbosity >= Verbosity::NORMAL) {
        RefreshablePrinter::clean(std::cout);

        std::cout << "> Estimated costs for gene tree "
                  << gtree_index + 1;
        if (utils::double_equal_or_superior(supportThreshold,
                                            sorted_supports.front())) {
          std::cout << " with contraction threshold at " << supportThreshold;
        }
        std::cout << ":" << std::endl;
        std::cout << "\t* Duplication cost = "
                  << treesReconciliationConductor.estimated_costs().first
                  << std::endl;
        std::cout << "\t* Loss cost = "
                  << treesReconciliationConductor.estimated_costs().second
                  << std::endl;
        std::cout << std::endl;
      }

      assert(solutions_for_one_tree_and_one_threshold.size() == params->sample_size);

#if defined(_OPENMP)
  # pragma omp critical
#endif
      {
        // Save each solution in a map with two keys: the first is the gene tree index, the second
        // the branch support threshold used to solve the gene tree.
        trees_solutions[genetrees.at(
            gtree_index)][supportThreshold] = std::move(
            solutions_for_one_tree_and_one_threshold);

        // Update progression bar.
        genetrees_progression.next(
            (unsigned int) genetrees.at(gtree_index)->getNumberOfNodes());
        if ((genetrees.size() > 1 or supportThresholdsToUse.size() > 1) and
            params->verbosity >= Verbosity::NORMAL)
          utils::progressionBar(std::cout, genetrees_progression,
                                true, std::string("Reconciling genetree ")
                                      + std::to_string(gtree_index + 1));
      }
    } //end for by support threshold
  } //end for by genetree

  // Manage final outputs
  if (params->statistics) { // If user ask for statistics, create statistics.
    auto statistics_table = Statistics::average_scores_by_thresholds(
        trees_solutions);
    statistics_table.csv(
        output_statistics_filename); // Save statistics in file.
    if (params->verbosity >= Verbosity::VERBOSE) {
      std::cout
          << "-Results-statistics--------------------------------------------------------------------------------"
          << std::endl;
      std::cout << statistics_table;
      std::cout
          << "---------------------------------------------------------------------------------------------------"
          << std::endl;
    }
  }

  RefreshablePrinter::clean(std::cout);

  // Print results in file
  if (params->write_gene_relationships_summary)
    IO::writeEventSummary(output_gene_relationships_summary_filename,
                          trees_solutions, genetrees);

  RefreshablePrinter::clean(std::cout);

  std::vector<std::string> written_files;
  written_files.reserve(params->trees_output_formats.size());
  for (const auto& text_format: params->trees_output_formats) {
    std::string output_filename(params->output);

    // Add output extension
    output_filename += ".";
    if (text_format == TextFormat::newick) {
      output_filename += NEWICK_EXTENSION;
    } else if (text_format == TextFormat::nhx) {
      output_filename += NHX_EXTENSION;
    } else if (text_format == TextFormat::phyloxml) {
      output_filename += PHYLOXML_EXTENTION;
    } else if (text_format == TextFormat::recphyloxml) {
      output_filename += RECPHYLOXML_EXTENTION;
    } else if (text_format == TextFormat::svg) {
      output_filename += SVG_EXTENTION;
    }

    IO::writeReconciledTrees(output_filename,
        trees_solutions,
        params->strict_support_thresholds,
        genetrees,
        *speciestree,
        text_format,
        not params->output_without_description);

    written_files.push_back(output_filename);
  }

  // Print infos like results and computations times.
  if (params->verbosity >= Verbosity::VERBOSE) {
    std::cout
        << "-Solutions-----------------------------------------------------------------------------------------"
        << std::endl;

    IO::writeReconciledTrees(std::cout,
        trees_solutions,
        params->strict_support_thresholds,
        genetrees,
        *speciestree,
        params->verbosity >= Verbosity::CHATTY,
        TextFormat::newick,
        true);

    std::cout
        << "---------------------------------------------------------------------------------------------------"
        << std::endl;

    std::cout << std::endl;

    std::cout
        << "-Execution time------------------------------------------------------------------------------------"
        << std::endl;

    std::cout << "> Mapping time: "
              << (mappingProgression.time_past_since_last_update())
              << " seconds."
              << std::endl;

    std::cout << "> Total solver execution time: "
              << (genetrees_progression.time_past_since_last_update())
              << " seconds."
              << std::endl;

    std::cout << "> Average solver execution time for one given tree: "
              << (genetrees_progression.time_past_since_last_update() /
                  ((double_t) genetrees.size() * (double_t) params->sample_size))
              << " seconds."
              << std::endl;
    std::cout
        << "---------------------------------------------------------------------------------------------------"
        << std::endl;
  }

  RefreshablePrinter::clean(std::cout);

  execution_progression.next();

  if (params->verbosity >= Verbosity::NORMAL) {
    std::cout << "Solution" << ((genetrees.size() == 1) ? (
        trees_solutions.begin()->second.size() > 1 ? "s" : "") : "s")
              << " ";
    utils::write(std::cout, written_files.begin(), written_files.end(),
                 "saved in ", ", ", "");
    std::cout << std::endl;
    std::cout << "Total elapsed time " << execution_progression.time_past()
              << " seconds." << std::endl;
  }

  return EXIT_SUCCESS;
}

/**
 * Check that the provided set of parameters is valid
 *
 * if not, issue a message on stderr and exit with failure code
 *
 * \param params
 */
void CheckParametersConsistency(
    const std::unique_ptr<TreerecsParameters>& params) {
  // A gene tree must be provided
  if (params->genetrees_filename.empty()) {
    std::cerr << "Please provide a gene tree file with parameter -g"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // A species tree must be provided unless the --info flag is set
  if (params->speciestree_filename.empty() and
      params->functionalities.find(Functionalities::PRINT_SUPPORT_DIAGRAM) ==
      params->functionalities.end()) {
    cerr << "Please provide a species tree with option -s or --speciestree"
         << endl;
    cerr << "Alternatively, you can ask for branch support information on a "
            "gene tree with the --info switch" << endl;
    exit(EXIT_FAILURE);
  }
}

/**
 * Check that output directory parameters comply with the following rules:
 *
 * outdir does not exist
 * OR
 * outdir exists and is writable and option -f | --force was provided
 *
 * if not, issue a message on stderr and exit with failure code
 *
 * \param params
 */
void CheckOutputDir(const std::unique_ptr<TreerecsParameters>& params) {
  // Create shorthands
  const auto& outdir = params->outdir;

  if (not IO::exists(outdir)) {
    IO::MakeDir(params->outdir, 0755);
  } else if (not IO::IsDirectory(params->outdir.c_str())) {
    cerr << "'" << params->outdir << "' exists and is not a directory" << endl;
    exit(EXIT_FAILURE);
  } else if (not params->force and
             not IO::IsDirEmpty(params->outdir.c_str())) {
    cerr << "Cowardly refusing to write in '" << outdir <<
            "' which is not empty." << endl <<
            "Choose another output directory using -o or use -f to force" << endl;
    exit(EXIT_FAILURE);
  }
}

void CheckBranchSupports(
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& genetrees,
    std::vector<std::vector<double>> branch_supports_vector,
    const std::unique_ptr<TreerecsParameters>& params
) {
  // For each tree,
  // if there is no bootstrap found in the tree, and a contraction was
  // requested (there is at least one positive threshold), issue a
  // warning/note to the user
  for (std::size_t i_tree = 0 ; i_tree < genetrees.size() ; ++i_tree) {
    auto& tree = genetrees[i_tree];
    auto& branch_supports = branch_supports_vector[i_tree];

    if (branch_supports.empty() and
        (params->thresholds_inputs.empty() or
         (params->thresholds_inputs.size() == 1 and
          params->thresholds_inputs.front() != "none" and
          params->thresholds_inputs.front() != "-1")) and
        (not params->reroot)) {
      // If the tree is not a single polytomy
      if (tree->getNumberOfLeaves() < tree->getNumberOfEdges()) {
        std::cerr << "Note: tree " << (i_tree + 1) << "/" << genetrees.size()
                  << " has no branch support, its topology won't be changed";
        auto polytomy_number = PhyloTreeToolBox::findPolytomies(*tree).size();
        if (polytomy_number > 0) {
          std::cerr << " but Treerecs will solve "
                    << polytomy_number
                    << (polytomy_number > 1 ? " polytomies" : " polytomy");
        }
        std::cerr << "." << std::endl;
      }
    }
  }
}
