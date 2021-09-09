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

#ifndef TREERECS_RANDOM_H
#define TREERECS_RANDOM_H

#include <cassert>

#include <numeric>
#include <random>

#include "treerecs/tools/utils.h"

namespace treerecs {

/// Random generator engine (standard mersenne_twister_engine)
#if defined(__clang__) // TODO<dpa> Problem only on mac OS => check and correct (check specifically the thread_local support)
static std::mt19937 default_random_generator;
#else
static thread_local std::mt19937 default_random_generator;
#endif

double Random(const double min = 0.0, const double max = 1.0);

/*!
 * @brief random pick of an index according to its weight.
 * @tparam Weight
 * @param weights Weight at an index.
 * @param n sample size
 * @return picked indexes.
 */
template <typename Weight>
std::vector<std::size_t>
WeightedRandomPick(const std::vector<Weight>& weights, const std::size_t n);

/*!
 * @brief random pick of an index according to its weight.
 * @tparam Weight
 * @param weights Weight at an index.
 * @return picked indexe.
 */
template <typename Weight>
std::size_t
WeightedRandomPick(const std::vector<Weight>& weights);

} // namespace treerecs

#include "random.hpp"

#endif //TREERECS_RANDOM_H
