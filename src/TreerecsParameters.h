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

#ifndef TREERECS_TREERECSPARAMETERS_H
#define TREERECS_TREERECSPARAMETERS_H

#include <string>
#include <list>

#include "treerecs/Constants.h"
#include "treerecs/tools/SpeciesGeneMapper.h"
#include "treerecs/tools/IO/IO.h"

class TreerecsCliUtils;

// Verbosity levels
enum class Verbosity {
  QUIET,
  NORMAL,
  VERBOSE,
  CHATTY
};

// Functionalities
enum class Functionalities {
  PRINT_HELP,
  PRINT_USAGE,
  PRINT_VERSION,
  PRINT_SUPPORT_DIAGRAM,
  RESOLVING,
  REROOTING,
  CORRECTION,
  DMATRIX,
  COST_ESTIMATION,
  ALE_LIKELIHOOD,
  ALE_SELECTION,
  PLL_LIKELIHOOD,
};

using namespace treerecs;

class TreerecsParameters {
  friend std::unique_ptr<TreerecsParameters> ParseCommandLine(int argc,
                                                              char** argv);

 protected:
  TreerecsParameters() = default;
 public:
  TreerecsParameters(const TreerecsParameters&) = delete;
  TreerecsParameters(TreerecsParameters&&) = default;
  TreerecsParameters& operator=(const TreerecsParameters&) = delete;
  TreerecsParameters& operator=(TreerecsParameters&&) = delete;
  virtual ~TreerecsParameters() = default;

 public:
  // Input file names
  /// Gene trees filename
  std::string genetrees_filename;
  /// Species tree filename
  std::string speciestree_filename;
  /// Filename with per-gene-family alignements paths
  std::string alignments_filename;
  /// Filename with distance matrix overriding tree annotations
  std::string dmatrix_filename;
  /// Gene<->species mapping filename
  std::string map_filename;

  // Gene and Species Trees parameters
  int genetree_index = -1;  // To select one gene tree in a file containing
  // multiple gene trees
  /// Mapping method to be used
  treerecs::MappingMethod mappingMethod = treerecs::MappingMethod::TREES;
  /// Character that separates species name and gene name in gene names
  std::string sep = DEFAULT_SEPARATION_CHARACTER_BETWEEN_GENE_AND_SPECIES_NAME;
  /// Indicates the position of the species name in gene names
  /// (prefix as species before the separator)
  bool prefix = DEFAULT_SPECIES_NAME_IS_PREFIX_IN_GENE_NAME;

  // Parameters of the model with duplication cost and loss cost.
  double dupcost = DEFAULT_DUPLICATION_COST;
  double losscost = DEFAULT_LOSS_COST;
  std::list<std::string> thresholds_inputs; // User-defined thresholds
  bool strict_support_thresholds =
      DEFAULT_STRICT_SUPPORT_THRESHOLDS;

  // Functionalities to perform
  std::set<Functionalities> functionalities;
  /// Whether to try to reroot gene tree(s)
  bool reroot = false;
  /// Number of results to compute for each gene tree
  std::size_t sample_size = DEFAULT_SAMPLE_SIZE;
  /// Whether to select the best gene tree in ex-aequos according to ale loglk
  bool add_ale_evaluation = false;
  bool ale_selection = false;
  /// Whether to estimate duplication and loss costs from the inputs
  bool costs_estimation = false;
  /// Whether to compute branch distances
  const bool compute_distances = true;

  // Output parameters
  /// Output directory
  std::string outdir = "treerecs_output";
  /// Output filename
  std::string output;
  /// Force possible overwrite
  bool force = false;
  /// Whether to omit descriptions in the results
  bool output_without_description = false;
  /// Formats in which to write the outputs (newick, nhx, ...)
  std::list<treerecs::TextFormat> trees_output_formats;
  /// Whether to save the gene<->species map
  bool saveMaps = false;
  bool write_gene_relationships_summary = false;
  /// Whether to print statistics
  bool statistics = false;

  // Verbosity
  Verbosity verbosity = Verbosity::NORMAL;

  // PRNG
  /// Seed used to initialize random generator
  long int seed = time(nullptr);
  /// Whether a specific seed is provided to initialize random generator
  bool seed_set = false;

  // Multithreading
  /// Number of threads to use
  int nthreads = 1;
};


#endif //TREERECS_TREERECSPARAMETERS_H
