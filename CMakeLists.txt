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

cmake_minimum_required(VERSION 3.4.3)

# Project name
project(Treerecs
        VERSION 1.2)

# Add compile definition TREERECS_VERSION_NUMBER
if(${CMAKE_VERSION} VERSION_LESS "3.12.4")
  add_definitions(-DTREERECS_VERSION_NUMBER="${PROJECT_VERSION}")
else()
  add_compile_definitions(TREERECS_VERSION_NUMBER="${PROJECT_VERSION}")
endif()

# Tell CMake where to look for custom modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Set a default build type if none was specified
include(DefaultBuildType)

# Add ad-hoc build type for testing
set(CMAKE_CXX_FLAGS_TESTING "-Os" CACHE STRING
        "Flags used by the C++ compiler during test builds."
        FORCE)
set(CMAKE_C_FLAGS_TESTING "-Os" CACHE STRING
        "Flags used by the C compiler during test builds."
        FORCE)
set(CMAKE_EXE_LINKER_FLAGS_TESTING
        "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
        "Flags used for linking binaries during test builds."
        FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_TESTING
        "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
        "Flags used by the shared libraries linker during test builds."
        FORCE)
MARK_AS_ADVANCED(
        CMAKE_CXX_FLAGS_TESTING
        CMAKE_C_FLAGS_TESTING
        CMAKE_EXE_LINKER_FLAGS_TESTING
        CMAKE_SHARED_LINKER_FLAGS_TESTING)

# Build-type-specific flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -O0")

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Testing."
        FORCE)

# Add uninstall target (uses a custom script)
configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

# On Windows, we use static linking to avoid dependencies to non-standard dlls
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

#if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
#  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.0")
#    message(FATAL_ERROR "Insufficient AppleClang version")
#  endif()
#elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0")
#    message(FATAL_ERROR "Insufficient gcc version")
#  endif()
#endif()

# Add OpenMP
set(USE_OMP_IF_FOUND TRUE CACHE BOOL
  "Use OpenMP if found on the system.")
if (USE_OMP_IF_FOUND)
  find_package(OpenMP)
  if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    message(WARNING "Compiler does not support OpenMP.")
  endif()
endif()

# We will use CTest
enable_testing()

# Add subdirectories
add_subdirectory(ext)
add_subdirectory(src)
