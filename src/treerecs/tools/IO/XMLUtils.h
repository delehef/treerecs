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

#ifndef TREERECS_XMLUTILS_H
#define TREERECS_XMLUTILS_H

#include <iostream>       // std::cout
#include <string>         // std::string
#include <tuple>          // std::tuple
#include <map>            // std::map
#include <cassert>        // assert()

namespace treerecs {

/// Defines types of XML tags.
enum XMLTagType {
  start /// opening XML tag.
  , end /// closing XML tag.
  , oneline /// XML tag in one block.
  , comment /// XML tag comment.
  , declaration /// XML tag declaration.
  , undefined /// Undefined/ not supported tag.
};

typedef struct {
  XMLTagType type;
  std::string name;
  std::map<std::string, std::string> attributes;
} XMLTag;

inline std::ostream& operator<<(std::ostream& os, XMLTag tag){
  os << "<Name=" << tag.name << " type=";
  switch (tag.type){
    case (XMLTagType::start):
      os << "start";
      break;
    case (XMLTagType::end):
      os << "end";
      break;
    case (XMLTagType::oneline):
      os << "oneline";
      break;
    case (XMLTagType::declaration):
      os << "declaration";
      break;
    default:
      os << "undefined";
  }
  os << " with " << tag.attributes.size() << " attributes.>" << std::endl;
  return os;
}

class XMLUtils {
public:

  ///Read content from XML element (<...>content</...>). Update i in the last '>' position.
  static std::string readXMLElement(const std::string& str, std::size_t& i);

  /// Reach next char element.
  /// Changes value of 'i' to the index of the character equal to 'c'.
  static void reachChar(std::size_t& i
                        , const std::string& str
                        , const char c
                        , const bool forward = true);

  /// Reach char element which matches with a given template condition from a
  /// value given by the first element (i).
  /// Changes value of the Iterator 'i'.
  ///
  /// \tparam Iterator Value to increment to reach value which is
  ///         validating condition.
  /// \tparam Condition Condition to validate to stop increment of i.
  /// \param i Value to change to reach the next value validating the condition.
  /// \param str
  /// \param stop Condition to validate to stop incrementation.
  /// \param forward Boolean, to look forward (++) or backward (--).
  template<typename Iterator, typename Condition>
  static void reach(
      Iterator& i /// Value to start search of character.
      , const std::string& str /// String where looking for.
      , const Condition& stop /// Condition to stop. Example: '[](const char c) { return c == '>'; }'.
      , const bool forward = true /// Looking forward (true), backward (false).
  ) {
    if(forward) {
      while (i < str.size() and not stop(str.at(i))) { i++; }
    } else {
      while (i >= 0 and not stop(str.at(i))) { i--; }
    }
  }

  /// Read the next tag and put i at the end of the tag (on '>').
  ///
  /// \param str String to explore.
  /// \param i
  /// \return Index used to explore the string.
  ///         If everything is ok, str[i] == '>' in the next tag.
  static XMLTag readNextXMLTag(
      const std::string& str
      , std::size_t& i
  );

  /// Get xml tag from string.
  static XMLTag readXMLTag(const std::string &strtag);

  /// Go to the closing tag according to the name given in first parameter.
  /// Returns the position of the '>' character in the second parameter.
  static std::size_t goToEndTag(
      const std::string& name /// Name of the XML tag.
      , const std::string& content /// String content where closing tag is.
      , const std::size_t& from = 0 /// Start of the search.
      , std::size_t to = 0 /// End of the search.
  );

};

} // namespace treerecs

#endif //TREERECS_XMLUTILS_H
