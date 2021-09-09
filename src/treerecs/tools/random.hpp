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

#ifndef TREERECS_RANDOM_HPP
#define TREERECS_RANDOM_HPP

namespace treerecs {

template <typename Weight>
std::vector<std::size_t>
WeightedRandomPick(const std::vector<Weight>& weights, const std::size_t n) {
  /// Uniform distribution in [0, 1)
  static std::uniform_real_distribution<double> uniform_0_1(0, 1);

  /// Random index pick according to the value in weights.
  /// For example, with [0.9, 0.1, 0.0], index 0 will be picked frequently,
  /// index 2 never.
  // Init the vector which contains picked indexes
  std::vector<std::size_t> picked;
  picked.reserve(n);

  Weight total = utils::sum(weights);

  // Get the descending order of each element of weight
  auto weights_indexes_in_order = utils::getOrder(weights, true);

  for (std::size_t i = 0 ; i < n ; i++) {
    double rand_value = uniform_0_1(default_random_generator);
    double cumul_prob_value = 0.0;
    for (std::size_t wi : weights_indexes_in_order) {
      cumul_prob_value += ((double) weights.at(wi) / total);
      if (cumul_prob_value > rand_value) {
        picked.push_back(wi);
        break;
      }
    }
  }

  return picked;
}

template <typename Weight>
std::size_t WeightedRandomPick(const std::vector<Weight>& weights) {
  /// Random index pick according to the value in weights.
  /// For example, with [0.9, 0.1, 0.0], index 0 will be picked frequently,
  /// index 2 never.
  auto values = WeightedRandomPick(weights, 1);
  assert(values.size() == 1);
  return values.at(0);
}

} // namespace treerecs

#endif //TREERECS_RANDOM_HPP
