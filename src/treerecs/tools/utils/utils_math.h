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

#ifndef TREERECS_UTILSMATHS_H
#define TREERECS_UTILSMATHS_H

#include <cassert>
#include <algorithm>

#include <treerecs/Constants.h>
#include "utils_containers.h"

namespace treerecs {

namespace utils {

  template<typename T>
  T sum(const std::vector <T> &x);

  template<typename T>
  inline T linearInterpolation(const T& v0, const T& v1, const T& t)
  {
    return (1.0 - t) * v0 + t * v1;
  }

  /*!
 * @brief Get quantiles of a sample of numbers.
 * @tparam Iterator
 * @param begin First element (as an iterator) of the sample to evaluate.
 * @param end Last element (as an iterator) of the sample to evaluate.
 * @param n Type of quantile (2 = median, 3 = tertile, 4 = quartile, ...)
 * @return A vector containing quantiles (as std::vector<double>).
 */
  template<typename Iterator>
  std::vector<double> quantile(const Iterator &begin, const Iterator &end,
      const std::size_t &n);

  /// Compare two doubles according to a specific precision, defined by delta
  /// (equivalent to operator==).
  inline bool double_equivalence(
      const double a, const double b
      , const double delta = DEFAULT_DOUBLE_EQUIVALENCE_PRECISION
  )
  { return fabs(a - b) < delta; }

  /// Compare two doubles (operator <=).
  inline bool double_equal_or_inferior(const double a, const double b) {
    return double_equivalence(a, b) or (a < b);
  }

  /// Compare two doubles (operator >=).
  inline bool double_equal_or_superior(const double a, const double b) {
    return double_equivalence(a, b) or (a > b);
  }

  /// Factorial(n).
  template<typename T>
  T factorial(T n){
    T res = 1.0;

    for (T i = 1.0 ; i <= n ; ++i) { // warning: double precision cannot be fine with equality
      res *= i;
    }

    return res;
  }

  /// Double factorial function.
  template<typename T>
  T double_factorial(T n) {
    T res = 1.0;

    for(T i = n ; i >= 1.0 ; i-= 2.0){
      // warning: double precision cannot be fine with equality
      res *=i;
    }

    return res;
  }
}

template<typename T>
T utils::sum(const std::vector <T> &x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

template<typename Iterator>
std::vector<double> utils::quantile(
    const Iterator &begin,
    const Iterator &end,
    const std::size_t &n)
{

  // Change name of template value type for "T".
  typedef typename std::iterator_traits<Iterator>::value_type T;

  std::vector<T> data {begin, end};

  if (data.size() == 0 or data.size() == 1 or n == 1)
  {
    return {data.begin(), data.end()};
  } else if(n < 1) {
    return {};
  }

  std::vector<double> probs(n - 1);

  probs[0] = (1.0/(double)n);

  std::size_t i = 1;
  std::generate(probs.begin() + 1, probs.end(),
                [&i, p = probs.front()](){ return (p * (++i)); });

  std::sort(data.begin(), data.end());
  std::vector<double> quantiles;
  quantiles.reserve(probs.size());

  for (i = 0; i < probs.size(); ++i)
  {
    double index_li = linearInterpolation<double>(0, data.size() - 1, probs[i]);

    auto left_index = std::max(int64_t(std::floor(index_li)), int64_t(0));
    auto right_index = std::min(int64_t(std::ceil(index_li)),
                                int64_t(data.size() - 1));

    auto left_element = static_cast<double>(data.at(left_index));
    auto right_element = static_cast<double>(data.at(right_index));

    auto quantile = linearInterpolation<double>(left_element, right_element,
                                                index_li - left_index);

    quantiles.push_back(quantile);
  }

  return quantiles;

} // namespace utils

} // namespace treerecs

#endif //TREERECS_UTILSMATHS_H
