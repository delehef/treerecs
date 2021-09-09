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

#ifndef TREERECS_UTILSCONTAINERS_H
#define TREERECS_UTILSCONTAINERS_H

#include <vector>
#include <list>
#include <algorithm>

namespace treerecs {

namespace utils {
  template<typename T>
  std::vector<std::size_t>
  getOrder(const std::vector<T> &x, const bool ascending = true);

  template<typename IterA, typename IterB>
  bool comp_all(
      IterA a_begin, const IterA &a_end, IterB b_begin, const IterB &b_end
  );

  template<typename IterA, typename IterB, typename Comp>
  bool comp_all(
      IterA a_begin, const IterA &a_end, IterB b_begin, const IterB &b_end
      , const Comp &comparator
  );

  /// Returns a list of value indexes which have a match with a template
  /// condition in a container.
  /// For example:
  ///
  ///   std::vector<int> my_vect {1, 42, 3, 300, 28};
  ///   std::cout << "Values > 30 in my_vect:" ;
  ///   auto indexes = matchesValuesIndexes(my_vect.begin(), my_vect.end(),
  ///                                       [](int x){ return x > 30; });
  ///   for(auto i: indexes) std::cout << " " << my_vect.at(i);
  ///   std::cout << "." << std::endl;
  ///
  /// \tparam Iterator
  /// \tparam Condition
  /// \param begin
  /// \param end
  /// \param condition
  /// \return std::list of indexes.
  template<typename Iterator, typename Condition>
  std::list<std::size_t> matchesValueIndexes(
      Iterator begin, const Iterator& end, const Condition& condition);

  /// Check if a container has an element which match with a condition
  template<typename Iteration, typename Lambda>
  bool
  contains(const Iteration &begin, const Iteration &end, const Lambda &lambda);
}

template<typename T>
std::vector<std::size_t>
utils::getOrder(const std::vector<T> &x, const bool ascending) {
  /// Returns a vector of indexes in the order of their corresponding value.
  /// Ex: std::vector<int> x = {3, 1, 2};
  ///     getOrder(x [, true]) will return {1, 2, 0} if ascending = true. 1 is
  ///     the index of the lowest value in the vector and 0 the index of the
  ///     greatest element.
  ///     getOrder(x, false) will return {0, 2, 1} if ascending = false.

  std::vector<std::size_t> order(x.size());
  std::size_t n(0);
  std::generate(std::begin(order), std::end(order), [&n] { return n++; });
  std::sort(std::begin(order), std::end(order),
            [&x, &ascending](const std::size_t &i, const std::size_t &j) {
              return (ascending) == (x.at(i) < x.at(j));
            }
  );
  return order;
}

template<typename IterA, typename IterB>
bool utils::comp_all(
    IterA a_begin, const IterA &a_end, IterB b_begin, const IterB &b_end
) {
  /// Compare two sorted containers.
  if(std::distance(a_begin, a_end) != std::distance(b_begin, b_end))
    return false;

  while(a_begin != a_end){
    if(*a_begin != *b_begin) return false;
    a_begin++;
    b_begin++;
  }
  return true;
}

template<typename IterA, typename IterB, typename Comp>
bool utils::comp_all(
    IterA a_begin, const IterA &a_end, IterB b_begin, const IterB &b_end
    , const Comp &comparator
) {
  /// Compare two sorted containers using a template comparator.
  if(std::distance(a_begin, a_end) != std::distance(b_begin, b_end))
    return false;

  while(a_begin != a_end){
    if(not comparator(*a_begin, *b_begin)) return false;
    a_begin++;
    b_begin++;
  }
  return true;
}

template<typename Iterator, typename Condition>
std::list<std::size_t> utils::matchesValueIndexes(
    Iterator begin, const Iterator& end, const Condition& condition) {
  std::list<std::size_t> res;
  std::size_t i = 0;
  while(begin != end){
    if(condition(*begin)) {
      res.push_back(i);
    }
    i++;
    begin++;
  }
  return res;
}

/// Check if a container has an element which match with a condition
template<typename Iteration, typename Lambda>
bool utils::contains(
    const Iteration &begin, const Iteration &end, const Lambda &lambda
) {
  auto current = begin;
  while(current != end){
    if(lambda(*current)) return true;
    current++;
  }
  return false;

} // namespace utils

} // namespace treerecs

#endif //TREERECS_UTILSCONTAINERS_H
