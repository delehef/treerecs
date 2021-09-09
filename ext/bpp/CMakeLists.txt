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

cmake_minimum_required(VERSION 3.0.1)

project(bpp)

# Specify C++11
include(CheckCXXCompilerFlag)
if(CYGWIN)
  CHECK_CXX_COMPILER_FLAG(-std=gnu++11 CXX_STANDARD_11_SUPPORT)
  if(CXX_STANDARD_11_SUPPORT)
    add_definitions(-std=gnu++11)
  else(CXX_STANDARD_11_SUPPORT)
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} does not support the -std=gnu++11 flag. Please use a different C++ compiler.")
  endif(CXX_STANDARD_11_SUPPORT)
else(CYGWIN)
  CHECK_CXX_COMPILER_FLAG(-std=c++11 CXX_STANDARD_11_SUPPORT)
  if(CXX_STANDARD_11_SUPPORT)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
  else(CXX_STANDARD_11_SUPPORT)
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif(CXX_STANDARD_11_SUPPORT)
endif(CYGWIN)


# Add Bpp sources
add_subdirectory(bpp-core)
add_subdirectory(bpp-seq)
add_subdirectory(bpp-phyl)
