cmake_minimum_required (VERSION 3.0.1 FATAL_ERROR)

project (Libpll)

add_subdirectory(src)

set(PLL_EXPORTED_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src PARENT_SCOPE)
