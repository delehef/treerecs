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

#ifndef PHYLASOLVER_COST_H
#define PHYLASOLVER_COST_H

#include <ostream>
#include <vector>
#include <treerecs/Constants.h>

namespace treerecs {

/// Types of events in genes history.
enum Event { duplication /// Duplication.
              , loss /// Loss.
              , speciation /// Speciation.
              , speciationLoss /// Speciation resulting in one loss son.
              , extant /// Extant.
              , bifurcationOut /// Bifurcation out.
              , none /// Nothing/ Unknown.
              };

/// Get event in std::string.
inline std::string event_to_str(const Event& event) {
  switch (event) {
    case duplication:
      return DUPLICATION_STR_FLAG;
    case loss:
      return LOSS_STR_FLAG;
    case speciation:
      return SPECIATION_STR_FLAG;
    case extant:
      return EXTANT_STR_FLAG;
    case bifurcationOut:
      return BIFURCATION_OUT_STR_FLAG;
    case speciationLoss:
      return SPECIATION_LOSS_STR_FLAG;
    default:
      return "none";
  }
}

class Cost {
  /*!
   * \class CostStruct
   * \brief contains cost informations for a Cost Table
   *
   * contains the value of the cost and what kind of event it corresponds to
   */
 public:
  double value = 0.0;
  std::vector<Event> path = {none};
  std::vector<std::size_t> occurence = {1};
  friend std::ostream& operator<< (std::ostream& os, const Cost& cost);
};

inline std::ostream& operator<< (std::ostream& os, const treerecs::Event& cost_type) {
  os << event_to_str(cost_type);
  return os;
}

std::ostream& operator<< (std::ostream& os, const Cost& cost);

} // namespace treerecs

#endif //PHYLASOLVER_COST_H
