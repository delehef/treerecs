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

#include "Nhx.h"

#include <treerecs/Constants.h>

namespace treerecs {

Nhx::Nhx(const bool useTagsAsPptNames) : bpp::Nhx(useTagsAsPptNames) {
  // Create bootstrap property
  //Property bootstrap_property("bootstrap", "B", true, 2);

  // Create species property
  //Property species_property("species", "S", false, 0);

  // Create duplication property
  Property duplication_property(DUPLICATION_STR_FLAG, "D", false, 0);

  // Register properties.
  //registerProperty(bootstrap_property);
  //registerProperty(species_property);
  registerProperty(duplication_property);
}

} // namespace treerecs
