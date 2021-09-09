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

#ifndef TREERECS_NHX_H
#define TREERECS_NHX_H

#include <Bpp/Phyl/Io/Nhx.h>

namespace treerecs {

/*!
 * @class Nhx
 * @brief Nhx provides functions to load a bpp::PhyloTree or write it into a Nhx file (inheritance from bpp::Nhx).
 */
class Nhx: public bpp::Nhx {
public:
  explicit Nhx(const bool useTagsAsPptNames = false);
};

} // namespace treerecs

#endif //TREERECS_NHX_H
