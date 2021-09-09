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

#include "Statistics.h"

namespace treerecs {

Table<double, double, std::string> Statistics::average_scores_by_thresholds(
    const std::unordered_map<std::shared_ptr<bpp::PhyloTree>, std::map<double, std::vector<ReconciledRootedTree>>> &input) {

  std::map<double, double> libpllLoglk_by_threshold;
  std::map<double, double> aleLoglk_by_threshold;
  std::map<double, double> totalCost_by_threshold;

  for(auto one_tree_results: input) {
    // first : the original tree
    // second : the map<thresholds, solutions>
    for(auto solutions_per_threshold: one_tree_results.second) {
      // first : threshold
      // second : reconciled trees
      double average_totalCost = 0.0;
      double average_aleLoglk = 0.0;
      double average_libpllLoglk = 0.0;
      for(auto& solution: solutions_per_threshold.second) {
        average_totalCost += solution.cost();
        if(solution.evaluated_with_ale())
          average_aleLoglk += solution.ale_loglikelihood();
        if(solution.evaluated_with_libpll())
          average_libpllLoglk += solution.libpll_loglikelihood();
      }
      libpllLoglk_by_threshold[solutions_per_threshold.first] = average_libpllLoglk/solutions_per_threshold.second.size();
      aleLoglk_by_threshold[solutions_per_threshold.first] = average_aleLoglk/solutions_per_threshold.second.size();
      totalCost_by_threshold[solutions_per_threshold.first] = average_totalCost/solutions_per_threshold.second.size();
    }
  }

  // Get all computed data.
  std::vector<double> thresholds;
  thresholds.reserve(aleLoglk_by_threshold.size());

  std::vector<double> aleLoglks;
  aleLoglks.reserve(aleLoglk_by_threshold.size());

  std::vector<double> libpllLoglks;
  aleLoglks.reserve(libpllLoglk_by_threshold.size());

  std::vector<double> costs;
  costs.reserve(aleLoglk_by_threshold.size());

  for(auto& pair: aleLoglk_by_threshold) {
    pair.second/= input.size();
    thresholds.push_back(pair.first);
    aleLoglks.push_back(pair.second);
  }

  for(auto& pair: libpllLoglk_by_threshold) {
    pair.second/= input.size();
    libpllLoglks.push_back(pair.second);
  }

  for(auto& pair: totalCost_by_threshold) {
    pair.second/= input.size();
    costs.push_back(pair.second);
  }

  // Compile all of data into a Table.
  Table<double, double, std::string> results;
  results.addCol(costs, "total_cost");
  results.addCol(aleLoglks, "ale_loglk");
  results.addCol(libpllLoglks, "libpll_loglk");
  results.setRowIndexes(thresholds);

  return results;
}

} // namespace treerecs
