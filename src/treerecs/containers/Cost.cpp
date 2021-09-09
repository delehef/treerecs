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

#include "Cost.h"

namespace treerecs {

std::ostream& operator<< (std::ostream& os, const treerecs::Cost& cost) {
  /*!
   * Streams all pairs of an unordered_map (std).
   */
  os << "{";
  os << cost.value;
  os << " : ";
  if(cost.path.empty())
    os << " ";
  else
    for(std::size_t i = 0; i < cost.path.size(); i++){
      auto& type = cost.path[i];
      auto& n = cost.occurence[i];
      os << type << ":" + std::to_string(n);
      if(i != cost.path.size() - 1)
        os << ", ";
    }
  os << "}";
  return os;
}

} // namespace treerecs
