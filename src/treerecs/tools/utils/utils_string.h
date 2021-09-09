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

#ifndef TREERECS_UTILSSTRING_H
#define TREERECS_UTILSSTRING_H

#include <string>
#include <treerecs/Constants.h>
#include "utils_containers.h"

namespace treerecs {

namespace utils {
/// Check if a given const char* means "yes" or "true".
bool isYes(const char* optarg);


inline bool
case_insensitive_char_comp(const unsigned char& a, const unsigned char& b) {
  return std::tolower(a) == std::tolower(b);
}

inline bool
case_sensitive_char_comp(const unsigned char& a, const unsigned char& b) {
  return a == b;
}

bool string_comp(
    const std::string& a, const std::string& b, const bool case_sensitive = true
);

/// Check if a given string matches with a regex pattern.
bool strmatch_regex(const std::string& str, const std::string& pattern);

/// Returns count of non-overlapping occurrences of 'sub' in 'str'.
unsigned int count(const std::string& str, const std::string& sub);

std::vector<std::string> splitString(
    const std::string& str, const char* c,
    const bool concatenate_delimiters = true
);

/// Return a file name, without path and file extension.
std::string extractFilename(const std::string& str);

inline std::string
trunc_string_number(const std::string& s, const std::size_t n_dec = 2) {
  auto dot_pos = std::find(s.begin(), s.end(), '.');
  for (std::size_t i = 0 ; i <= n_dec and dot_pos != s.end() ; ++i)
    dot_pos++;
  return std::string(s.begin(), dot_pos);
}

/// Trim at left.
inline void ltrim_str(std::string& str) {
  str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](int ch) {
    return !std::isspace(ch);
  }));
}

/// Trim at right.
inline void rtrim_str(std::string& str) {
  str.erase(std::find_if(str.rbegin(), str.rend(), [](int ch) {
    return !std::isspace(ch);
  }).base(), str.end());
}

/// Trim at left and right.
inline void trim_str(std::string& str) {
  utils::ltrim_str(str);
  utils::rtrim_str(str);
}

/// Trim at right.
inline std::string ltrim(std::string str) {
  utils::ltrim_str(str);
  return str;
}

/// Trim at left and right.
inline std::string rtrim(std::string str) {
  utils::rtrim_str(str);
  return str;
}

/// Trim at left and right.
inline std::string trim(std::string str) {
  utils::trim_str(str);
  return str;
}

/// Returns a string with a specific substring replaced by an other
/// \param str
/// \param old_substr Substring to replace, can be a regex pattern.
/// \param new_substr New element to put in place of the previous substring/
///        pattern.
/// \return std::string with substrings replaced (all occurrences).
std::string replace(const std::string& str, const std::string& old_substr,
                    const std::string& new_substr);

/// Returns a string with a specific char replaced by an other
/// \param str
/// \param old_char character to replace.
/// \param new_char New element to put in place of the previous character.
/// \return std::string with substrings replaced (all occurrences).
std::string
replace(const std::string& str, const char old_char, const char new_char);

/// Check if a string has a match with a number.
bool stringIsNumber(const std::string& s);

} // namespace utils

} // namespace treerecs

#endif //TREERECS_UTILSSTRING_H
