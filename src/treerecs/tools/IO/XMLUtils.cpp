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

#include "XMLUtils.h"

namespace treerecs {

std::string XMLUtils::readXMLElement(const std::string &str, std::size_t &i) {
  std::size_t begin = i;

  XMLTag tag = readNextXMLTag(str, i);
  if(tag.type == XMLTagType::start){
    i++;
    begin = i;
  } else {
    std::cerr << "XMLUtils::readXMLElement: bad XML Element given with "
              << str.substr(begin, 10)
              << "..." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string tag_name = tag.name;

  i = goToEndTag(tag_name, str, i, str.size());

  reachChar(i, str, '<', false);
  std::string res = str.substr(begin, i - begin);
  reachChar(i, str, '>', true);

  return res;
}

void XMLUtils::reachChar(
    std::size_t &i, const std::string &str, const char c, const bool forward
) {
  return reach(i, str, [c](const char& e) { return e == c; }, forward);
}

XMLTag XMLUtils::readNextXMLTag(const std::string &str, std::size_t &i) {
  XMLUtils::reach(i, str, [](const char& c) { return c == '<'; });
  std::size_t start = i;
  XMLUtils::reach(i, str, [](const char& c) { return c == '>'; });
  std::string substr = str.substr(start, i - start + 1);
  return readXMLTag(substr);
}

XMLTag XMLUtils::readXMLTag(const std::string &strtag) {
  // results
  XMLTag res;
  res.type = undefined;
  res.name = "";
  res.attributes.clear();

  std::size_t i = strtag.find_first_of("<");

  if(i != std::string::npos) {
    i++;

    // Check next character
    reach(i, strtag, [](const char& c) { return c != ' '; });

    // The first character can give a clue on the type of the XML tag
    if (strtag.at(i) == '/') {
      res.type = end;
      ++i;
    } else if(strtag.at(i) == '!') {
      res.type = comment;
      i+=2; //Because of the two following '-'.
    } else if(strtag.at(i) == '?') {
      res.type = declaration;
    } else if(strtag.at(i) == '>'){
      res.type = undefined;
      return res;
    } else {
      res.type = start;
    }

    // get the name of the tag.
    while (i < strtag.size() and strtag.at(i) != '>' and strtag.at(i) != ' ') {
      res.name += strtag.at(i);
      ++i;
    }

    // get all parameters
    bool is_value = false;
    std::string key = "";
    std::string value = "";

    while (i < strtag.size() and strtag.at(i) != '>') {
      if (strtag.at(i) == '\"') {
        is_value = !is_value;
        if (not is_value) {
          // If the read of the value is finished, add key and value and create
          // new keys and values.
          res.attributes[key] = value;
          key = "";
          value = "";
        }
      } else if (is_value) {
        value += strtag.at(i);
      } else {
        if(strtag.at(i) != ' ' and strtag.at(i) != '=' and strtag.at(i) != '/')
          key += strtag.at(i);
        else if(strtag.at(i) == '/') {
          break;
        }
      }
      i++;
    }

    // end line
    if(strtag.at(i) != '>') {
      while (i < strtag.size() and strtag.at(i) != '>') {
        if (strtag.at(i) == '/') {
          res.type = oneline;
          break;
        }
        i++;
      }
    }
  } else {
    res.type = undefined;
  }

  return res;
}

std::size_t XMLUtils::goToEndTag(
    const std::string &name, const std::string &content, const size_t &from
    , size_t to
) {
  if(to <= 0) to = content.size();

  assert(content.size() >= to);

  std::size_t pos = from;

  int open_tags = 0; // Number of tags open tags

  while(pos < to and open_tags >= 0) {
    XMLTag tag_xml = readNextXMLTag(content, pos);
    if(tag_xml.type == XMLTagType::start) {
      open_tags++;
    } else if(tag_xml.type == end) {
      if(open_tags == 0 and tag_xml.name == name) return pos;
      else open_tags--;
    }
    if(open_tags < 0)
      std::cerr << "Error while seeking closing tag for "
                << name << "." << std::endl;
  }

  return to;
}

} // namespace treerecs
