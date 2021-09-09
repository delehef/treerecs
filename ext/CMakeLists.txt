# Copyright (C) 2018  INRIA
#
# Treerecs is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Treerecs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Add _USE_MATH_DEFINES for cmath.
add_definitions(-D_USE_MATH_DEFINES)

# Add Bio++ and Boost sources and create a libraries
add_subdirectory(boost)
add_subdirectory(bpp)
add_subdirectory(pll)
add_subdirectory(simple_svg)
