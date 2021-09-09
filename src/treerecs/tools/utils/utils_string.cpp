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

#include "utils_string.h"

#include <regex>

#include <treerecs/tools/IO/IO.h>

namespace treerecs {

static auto& preferred_separator = IO::preferred_separator; // until C++17

bool utils::stringIsNumber(const std::string &s) {
  /// Check if a string is an int
  if (s.empty())
    return false;
  char* p;
  std::strtol(s.c_str(), &p, 10);
  return (*p == 0);
}

bool utils::string_comp(
    const std::string &a, const std::string &b, const bool case_sensitive
) {
  if (a.length() == b.length()) {
    return std::equal(b.begin(), b.end(),
                      a.begin(),
                      case_sensitive ? utils::case_sensitive_char_comp
                                     : utils::case_insensitive_char_comp);
  } else {
    return false;
  }
}

std::vector<std::string> utils::splitString(
    const std::string &str, const char *c, const bool concatenate_tokens
) {
  /// Split a string according to delimiters defined in c.
  std::string buff{""};
  std::vector<std::string> ll;

  for (auto n : str) {
    // find if n is in the delimiters list in const char* c
    bool n_is_delimiter = false;
    std::size_t i = 0;

    while (!n_is_delimiter && c[i] != '\0') {
      n_is_delimiter = (n == c[i++]);
    }

    if (!n_is_delimiter) {
      buff += n;
    } else if (buff != "") {
      ll.push_back(buff);
      buff = "";
    } else if (not concatenate_tokens) {
      ll.push_back(buff);
    }
  }
  if (buff != "") ll.push_back(buff);

  return ll;
}

bool utils::isYes(const char *optarg) {
  /// Check if a given const char* means "yes" or "true".
  return
      string_comp(optarg, "true", false)
      or string_comp(optarg, "yes", false)
      or string_comp(optarg, "t", false)
      or string_comp(optarg, "y", false);
}

bool utils::strmatch_regex(const std::string &str, const std::string &pattern) {
  auto regex = std::regex(pattern);
  return std::regex_match(str, regex);
}

unsigned int utils::count(const std::string &str, const std::string &sub) {
  if (sub.length() == 0 or sub.length() > str.length()) return 0;
  unsigned int count = 0;
  for (std::size_t offset = str.find(sub); offset != std::string::npos;
       offset = str.find(sub, offset + sub.length()))
  {
    ++count;
  }
  return count;
}

std::string utils::extractFilename(const std::string &str) {
  // First, extract the correct file name (without the path). So, find the last
  // path delimiter (which is the beginning of the filename) and keep only the
  // string after.
  std::size_t last_path_delimiter = str.find_last_of(preferred_separator);

  if (last_path_delimiter != std::string::npos)
    return str.substr(last_path_delimiter + 1);
  else
    return str;
}

std::string utils::replace(const std::string& str
                           , const std::string& old_substr
                           , const std::string& new_substr
) {
  return std::regex_replace(str, std::regex(old_substr), new_substr);
}

std::string utils::replace(const std::string& str
                           , const char old_char
                           , const char new_char
) {
  std::string res(str);
  for (char& c : res) {
    if (c == old_char)
      c = new_char;
  }
  return res;
}

} // namespace treerecs
