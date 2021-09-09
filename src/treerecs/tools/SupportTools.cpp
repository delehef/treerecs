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

#include "SupportTools.h"

#include <iostream>

#include <TreerecsCliUtils.h>
#include "Statistics.h"

bool AtLeastOneBinaryTree(
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& trees) {
  for (auto& tree : trees) {
    // TODO(dpa): This block should be as simple as the following:
    //if (tree->isBinary()) {
    //  return true;
    //}
    // Such a function exists: PhyloTreeToolBox::isBinary
    // However, its behaviour differs from what is done here
    auto polytomies = PhyloTreeToolBox::findPolytomies(*tree);

    // A tree with no polytomy is binary.
    if (polytomies.size() == 0) {
      return true;
    }

    // A tree that has a single polytomy with 3 children at its "root" is
    // actually an unrooted binary tree
    if (polytomies.size() == 1 and
        polytomies.front() == tree->getRoot() and
        tree->getSons(polytomies.front()).size() == 3) {
      return true;
    }
  }

  // All trees have been considered, none of them are binary
  return false;
}

std::vector<double> GetBranchSupports(
    const std::shared_ptr<bpp::PhyloTree>& tree) {
  std::vector<double> branch_supports;
  for (const auto& edge: tree->getAllEdges()) {
    if (edge->hasBootstrapValue()) {
      branch_supports.push_back(edge->getBootstrapValue());
    }
  }
  return branch_supports;
}

std::vector<std::vector<double>> GetBranchSupports(
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& genetrees) {
  // Init container
  std::vector<std::vector<double>> branch_supports_vector;

  for (const auto& tree : genetrees) {
    branch_supports_vector.push_back(GetBranchSupports(tree));
  }

  return branch_supports_vector;
}

std::vector<double> ConvertThresholdInputsToDoubles(
    std::unique_ptr<TreerecsParameters>& params,
    const std::vector<double>& sorted_supports) {
  // Default (user has not specified a thing
  if (params->thresholds_inputs.empty()) {
    return {INIT_SUPPORT_THRESHOLD};
  }

  // Init resulting vector
  std::vector<double> supportThresholdsToUse;

  // Get thresholds asked by the user
  auto estimated_number_of_thresholds = params->thresholds_inputs.size();
  supportThresholdsToUse.reserve(estimated_number_of_thresholds);
  for (auto& thresholds_input: params->thresholds_inputs) {
    if ("none" == thresholds_input and (not sorted_supports.empty())) {
      // No contraction (default behaviour).
      supportThresholdsToUse.push_back(-1);
    } else if ("all" == thresholds_input and (not sorted_supports.empty())) {
      // Contract all branches
      supportThresholdsToUse.push_back(std::numeric_limits<double>::infinity());
    } else if (thresholds_input.compare(0, 9, "quantiles") == 0
               and (not sorted_supports.empty())) {
      // Set number of "quantile-groups"
      std::size_t number_of_groups = thresholds_input.size() == 9 ?
          DEFAULT_NB_FREQUENCY_CLASSES :
          static_cast<size_t>(
              stoi(thresholds_input.substr(9, thresholds_input.size())));

      if(number_of_groups < 1) {
        std::cerr << "the number of \"quantile-groups\" needs to be >= 1"
                  << " (" << std::to_string(number_of_groups) << " given)"
                  << std::endl;
        exit(EXIT_FAILURE);
      }

      // In the quantiles mode we aggregate different thresholds:
      // * no contraction
      // * the frontiers between branch support quantiles
      // * contract all
      // Add no contraction
      supportThresholdsToUse.push_back(-1);

      // Add quantiles
      auto supports_quantiles = utils::quantile(sorted_supports.begin(),
                                                sorted_supports.end(),
                                                number_of_groups);
      supportThresholdsToUse.insert(supportThresholdsToUse.end(),
                                    supports_quantiles.begin(),
                                    supports_quantiles.end());// Add quantiles

      // Add contract all
      supportThresholdsToUse.push_back(
          std::numeric_limits<double>::infinity());
    } else {
      // No special value found, value should be numeric
      // todo: Check that string is a number.
      supportThresholdsToUse.push_back(std::stod(thresholds_input));
      // If the user has postfixed the value with +epsilon (or a shorthand
      // thereof), set the corresponding flag accordingly
      if (thresholds_input.rfind("+e") != std::string::npos) {
        // Check that either this is the first threshold provided or that the
        // "epsilon" flag has already been set
        if (supportThresholdsToUse.size() > 1 and
            params->strict_support_thresholds) {
          throw malformed_cmd("the epsilon postfix must either be present on "
                              "all threshold values or none");
        }
        params->strict_support_thresholds = false;
      } else if (not params->strict_support_thresholds) {
        throw malformed_cmd("the epsilon postfix must either be present on "
                            "all threshold values or none");
      }
    }
  }

  return supportThresholdsToUse;
}

/**
 * @brief Print frequency diagram for the provided support values
 */
void PrintSupportDiagram(const std::vector<double>& sorted_supports) {
  if (sorted_supports.size() > 0) {
    // First: get frequencies of thresholds which are separated in several classes.
    auto supportFrequency = Statistics::frequencies(
        sorted_supports.begin(), sorted_supports.end(),
        DEFAULT_NB_FREQUENCY_CLASSES);
    std::cout
        << "> Branch supports distribution (support value, frequency):                           "
        << std::endl;
    Statistics::plotFrequencies(supportFrequency,
                                std::cout); // Print in terminal a pseudo-diagram of the frequencies.
    std::cout << std::endl;

    // Compute median support value
    //const auto median_it1 = sorted_supports.begin() + sorted_supports.size() / 2 - 1;
    //const auto median_it2 = sorted_supports.begin() + sorted_supports.size() / 2;
    //auto median_support_value = (sorted_supports.size() % 2 == 0) ? (*median_it1 + *median_it2) / 2 : *median_it2;
    auto median_support_value = utils::quantile(sorted_supports.begin(),
                                                sorted_supports.end(),
                                                2).front();
    //assert(utils::quantile(sorted_supports.begin(), sorted_supports.end(), 2).front() == median_support_value);

    std::cout << "Median value of supports: ";
    std::cout << median_support_value;
    std::cout << " (min = " << sorted_supports.front();
    std::cout << " and max = " << sorted_supports.back() << ")" << std::endl;
    //std::cout << "Please run the program with the option -t (--threshold)." << std::endl;
  } else {
    // User has asked for a support frequency diagram but has only provided
    // gene trees with no support values
    // => issue a message on stderr and exit with failure code
    std::cerr << "It looks like you have asked for information about branch "
        "support but have provided only gene trees with no branch "
        "support values." << std::endl;
    exit(EXIT_FAILURE);
  }
}
