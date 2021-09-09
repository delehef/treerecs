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

#ifndef TREERECS_UTILSSTREAM_H
#define TREERECS_UTILSSTREAM_H

//Include std
#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <cstring>

//Include Bpp
#include <Bpp/Phyl/Io/Newick.h>

#include <treerecs/tools/Timer.h>
#include <treerecs/tools/IO/RefreshablePrinter.h>

namespace treerecs {

// Print in std::ostream stl elements.
template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &pair);

template<class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v);

template<class T>
std::ostream &operator<<(std::ostream &os, const std::list<T> &l);

template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::map<T1, T2> &map);

template<class T1, class T2>
std::ostream &
operator<<(std::ostream &os, const std::unordered_map<T1, T2> &map);

std::ostream &operator<<(std::ostream &os, const bpp::PhyloTree &tree);

std::ostream &
operator<<(std::ostream &os, const std::shared_ptr<bpp::PhyloNode> &node_ptr);

namespace utils {
  template<typename T>
  void write_array(std::ostream& os, const T* elements, const std::size_t n);

  template<typename Iterator>
  std::ostream& write(std::ostream &os
                      , Iterator it
                      , const Iterator& itend
                      , const std::string& opening_str = "{"
                      , const std::string& elements_separator = ", "
                      , const std::string& closing_str = "}");

  void printNodeContent(
      const bpp::PhyloTree &tree, const std::shared_ptr<bpp::PhyloNode> &node
      , std::ostream &os = std::cout);

  void print_temp(bpp::PhyloTree &tree, std::ostream &os);

  /// Returns length in stream of a given element.
  template<class T>
  std::size_t getStreamObjectSize(const T &o, const std::ostream &os);

  /// \brief Change time in seconds to a string. The time will be decomposed in
  ///        hours, minutes and seconds.
  /// \param time number of seconds.
  std::string time_to_str(const double time);

  /// \brief Prints a progression bar in a given ostream.
  /// \tparam T
  /// \param os std::ostream to print the progression bar.
  /// \param timer Timer to print.
  /// \param remaining_time print an estimation of remaining time
  /// \param text change text printed with the bar.
  /// \param force_print progression bar is actually printed according to a
  ///                    frequency display. This option set to true allows a
  ///                    print for each call.
  template<typename T>
  void progressionBar(
      std::ostream &os
      , const Timer<T> &timer
      , const bool remaining_time = false
      , const std::string &text = "Progression"
      , const bool force_print = false);

  template<typename T>
  void write_array(std::ostream &os, const T *elements, const std::size_t n) {
    os << "{";
    for (std::size_t i = 0; i < n; i++) {
      os << elements[i];
      if (i < (n - 1)) os << ", ";
    }
    os << "}";
  }

template<typename Iterator>
std::ostream& write(std::ostream &os
    , Iterator it
    , const Iterator& itend
    , const std::string& opening_str
    , const std::string& elements_separator
    , const std::string& closing_str) {
  os << opening_str;
  if (it != itend) {
    os << *it;
    it++;
    while (it != itend) {
      os << elements_separator << *it;
      it++;
    }
  }
  os << closing_str;

  return os;
}

template<class T>
std::size_t getStreamObjectSize(const T &o, const std::ostream &os) {
  std::ostringstream temp;
  auto old_precision = os.precision();
  temp << std::setprecision(old_precision) << o;
  return temp.str().size();
}

template<typename T>
void progressionBar(
    std::ostream &os
    , const Timer<T> &timer
    , const bool remaining_time
    , const std::string &text
    , const bool force_print
){
  if (RefreshablePrinter::printable() or force_print) {
    int progression_bar_size = 20;
    std::string message;
    if(not text.empty()) {
      message += text + " ";
    }
    message += "[";
    double progression = timer.progress();
    if (progression > 1.0) progression = 1.0;
    int pos = (int) (progression_bar_size * progression);
    for (int cursor = 0; cursor < progression_bar_size; ++cursor) {
      if (cursor < pos) message += std::string("=");
      else if (cursor == pos) message += std::string(">");
      else message += std::string(" ");
    }
    message +=
        std::string("]") + std::to_string((int) round(100.0 * progression)) +
        "%";

    // Print remaining time.
    if (remaining_time) {
      message += " (remaining time: ";
      if (timer.current() == timer.min())
        message += "unknown";
      else
        message += time_to_str(timer.remaining_time());
      message += ")";
    }

    RefreshablePrinter::print(os, message, force_print);
  }
}

} // namespace utils

/// Print a std::pair.
template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &pair) {
  os << "(" << pair.first << ", " << pair.second << ")";
  return os;
}

/*!
 * @brief Prints each element of a given std::vector.
 * @param os
 * @param v The vector that contains each element to print
 * @return std::ostream&
 */
template<class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  ///Streams all elements of a vector (std).
  return treerecs::utils::write(os, v.begin(), v.end(), "[", ", ", "]");
}

template<class T>
std::ostream &operator<<(std::ostream &os, const std::list<T> &l) {
  /// Streams all elements of a list (std).
  return treerecs::utils::write(os, l.begin(), l.end(), "[", ", ", "]");
}

template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::map<T1, T2> &map) {
  /// Streams all pairs of a map (std).
  return treerecs::utils::write(os, map.begin(), map.end(), "{", ", ", "}");
}

template<class T1, class T2>
std::ostream &
operator<<(std::ostream &os, const std::unordered_map<T1, T2> &map) {
  /// Streams all pairs of an unordered_map (std).
  return treerecs::utils::write(os, map.begin(), map.end(), "{", ", ", "}");
}

} // namespace treerecs

#endif //TREERECS_UTILSSTREAM_H
