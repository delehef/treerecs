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

# Specify C++17
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++17 CXX_STANDARD_17_SUPPORT)
if(CXX_STANDARD_17_SUPPORT)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
else()
  message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
endif()

# Export variables as CACHE entries
set(TREERECS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL "treerecs: Include dirs" FORCE)

# Add treerecs-core library
add_subdirectory(treerecs)

# Add executables misc
add_subdirectory(misc)

# Add executable Treerecs
add_executable(treerecs Treerecs.cpp TreerecsCliUtils.h TreerecsCliUtils.cpp TreerecsParameters.h)
target_link_libraries(treerecs
  PRIVATE treerecs-core)
set_target_properties(treerecs
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  )

# Add install target
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
install(TARGETS treerecs RUNTIME DESTINATION bin)

# Add uninstall target
add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake)

# libpll
set(LIBPLL_INCLUDE_DIR ${LIBPLL_INCLUDE_DIR} PARENT_SCOPE)
set(LIBPLL_LIBRARY ${LIBPLL_LIBRARY} PARENT_SCOPE)
