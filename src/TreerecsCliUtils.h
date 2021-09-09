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

#ifndef TREERECS_TREERECSCLIUTILS_H
#define TREERECS_TREERECSCLIUTILS_H

#include "TreerecsParameters.h"

#include <getopt.h>
#include <exception>

/// Malformed command exception
class malformed_cmd : public std::exception {
 public:
  malformed_cmd() : what_("unknown error") {}
  malformed_cmd(std::string msg) : what_(msg) {}
  const char* what() const noexcept override {
    return what_.c_str();
  }
 protected:
  std::string what_;
};

/// Parse command line
std::unique_ptr<TreerecsParameters> ParseCommandLine(int argc, char** argv);

/// Print Minimal header
void PrintMinimalHeader();

/// Print Treerecs current version
void PrintVersion();

/// Print Treerecs usage
void PrintUsage(char* prog_path);

/// Print Treerecs help
void PrintHelp(char* prog_path);

/// Summarize all inputs and outputs, thus the behaviour of Treerecs for a run
void SummarizeParameters(
    std::ostream& os,
    const std::unique_ptr<TreerecsParameters>& params,
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& genetrees,
    const std::vector<double>& sorted_supports,
    const std::vector<double>& supportThresholdsToUse,
    const std::string& output_map_filename,
    const std::string& gene_relationships_summary_filename);


#endif //TREERECS_TREERECSCLIUTILS_H
