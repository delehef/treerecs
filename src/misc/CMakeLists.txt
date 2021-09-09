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

#Create executable Ale.
add_executable(ALEevaluate_undated ALEevaluate_undated.cpp)
target_link_libraries(ALEevaluate_undated
  treerecs-core)
target_include_directories(ALEevaluate_undated
  PRIVATE ${TREERECS_INCLUDE_DIRS})
set_target_properties(ALEevaluate_undated
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/misc"
  )

#Create executable genetreeEditor.
add_executable(genetreeEditor genetreeEditor.cpp)
target_link_libraries(genetreeEditor
  treerecs-core)
target_include_directories(genetreeEditor
  PRIVATE ${TREERECS_INCLUDE_DIRS})
set_target_properties(genetreeEditor
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/misc"
  )

#Move tools README.md to binaries.
file(COPY "README.md" DESTINATION "${CMAKE_BINARY_DIR}/bin/misc")
