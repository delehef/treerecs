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

#ifndef TREERECS_SUPPORTTOOLS_H
#define TREERECS_SUPPORTTOOLS_H

#include <vector>
#include <memory>

#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <TreerecsParameters.h>


/**
 * @brief Check whether at least one of the provided trees is binary
 * @param trees trees to check
 * @return true if at least one the provided trees is binary, false otherwise
 */
bool AtLeastOneBinaryTree(
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& trees);

std::vector<std::vector<double>> GetBranchSupports(
    const std::vector<std::shared_ptr<bpp::PhyloTree>>& genetrees);

std::vector<double> ConvertThresholdInputsToDoubles(
    std::unique_ptr<TreerecsParameters>& params,
    const std::vector<double>& sorted_supports);

void PrintSupportDiagram(const std::vector<double>& sorted_supports);

#endif //TREERECS_SUPPORTTOOLS_H
