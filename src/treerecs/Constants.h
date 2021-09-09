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

#ifndef TREERECS_CONSTANTS_H
#define TREERECS_CONSTANTS_H

#include <iostream>
#include <random>

namespace treerecs {

/// Define Pi
template<class T>
constexpr T PI = T(3.1415926535897932385L);

/// Default value of a cost for one duplication event.
constexpr double DEFAULT_DUPLICATION_COST = 2.0;

/// Default value of a cost for one loss event.
constexpr double DEFAULT_LOSS_COST = 1.0;

/// Default output filename.
constexpr const char* OUTPUT_POSTFIX = "recs";

constexpr const char* NEWICK_EXTENSION = "nwk";

constexpr const char* NHX_EXTENSION = "nhx";

constexpr const char* PHYLOXML_EXTENTION = "phylo.xml";

constexpr const char* RECPHYLOXML_EXTENTION = "recphylo.xml";

constexpr const char* SVG_EXTENTION = "svg";

constexpr const char* SMAP_EXTENTION = "smap";

constexpr const char* STATS_EXTENTION = "csv";

constexpr const char* FEVENT_EXTENTION = "txt";

/// Default global map output filename.
constexpr const char* DEFAULT_GLOBAL_MAP_OUTPUT_FILENAME = "global_map";

/// Default statistics output filename.
constexpr const char* DEFAULT_STATISTICS_OUTPUT_FILENAME = "statistics";

/// Default svg output filename.
constexpr const char* DEFAULT_SVG_OUTPUT_FILENAME = "";

/// Default gene relationships summary filename
constexpr const char* DEFAULT_FEVENT_OUTPUT_FILENAME = "relationships_summary";

/// Default character which separates species name from gene name.
constexpr const char* DEFAULT_SEPARATION_CHARACTER_BETWEEN_GENE_AND_SPECIES_NAME = "_";

/// Default position on species name in gene name according to a separator
constexpr bool DEFAULT_SPECIES_NAME_IS_PREFIX_IN_GENE_NAME = false;

/// Initial value for the support thresold.
constexpr double INIT_SUPPORT_THRESHOLD = -1.0;

/// Contraction behaviour when there is no support on a branch
constexpr bool CONTRACT_BRANCHES_WITHOUT_SUPPORT = true;

/// Default branch support value if not provided
constexpr double DEFAULT_BRANCH_SUPPORT_VALUE = 0.0;

/// Default size of the sample.
constexpr std::size_t DEFAULT_SAMPLE_SIZE = 1;

/// Default number of classes in a frequency classification.
constexpr std::size_t DEFAULT_NB_FREQUENCY_CLASSES = 4;

/// Default std::string flag of a gene duplication.
constexpr const char* DUPLICATION_STR_FLAG = "duplication";

/// Default std::string flag of a gene loss.
constexpr const char* LOSS_STR_FLAG = "loss";

/// Default std::string flag of a speciation event.
constexpr const char* SPECIATION_STR_FLAG = "speciation";

/// Default std::string flag of an extant event.
constexpr const char* EXTANT_STR_FLAG = "leaf";

/// Default std::string flag of a bifurcation out.
constexpr const char* BIFURCATION_OUT_STR_FLAG = "bifurcationOut";

/// Default std::string flag of a speciation loss event.
constexpr const char* SPECIATION_LOSS_STR_FLAG = "speciationLoss";

/// Precision of doubles for comparisons.
constexpr double DEFAULT_DOUBLE_EQUIVALENCE_PRECISION = 1.0e-10;

/// Minimal branch length.
constexpr double MINIMAL_BRANCH_LENGTH = 1.0e-6;

/// Case sensitive mapping.
constexpr bool DEFAULT_CASE_SENSITIVE_MAPPING = false;

/// Threshold comparisons inferior only
constexpr bool DEFAULT_STRICT_SUPPORT_THRESHOLDS = true;

/// Maximal plot width.
constexpr int MAXIMAL_PLOT_WIDTH = 20;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
  constexpr const char* LINE_FEED = "\r\n";
  constexpr const char* CARRIAGE_RETURN = "\r";
#else
  constexpr const char* LINE_FEED = "\n";
  constexpr const char* CARRIAGE_RETURN = "\r";
#endif

} // namespace treerecs

#endif //TREERECS_CONSTANTS_H
