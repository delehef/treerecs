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

set(treerecs-core_HEADER_FILES
  # Containers
  containers/AbstractGrid.h
  containers/AbstractGrid.hpp
  containers/Cost.h
  containers/GeneMap.h
  containers/Grid.h
  containers/Table.h
  containers/Singleton.h
  containers/ReconcilableRootedTree.h
  containers/ReconciledRootedTree.h
  containers/NodeProperty.h
  containers/BipartitionList.h
  # libpll
  tools/LibpllEvaluation.h
  # Tools
  tools/SpeciesGeneMapper.h
  tools/utils.h
  tools/utils/utils_stream.h
  tools/utils/utils_math.h
  tools/utils/utils_containers.h
  tools/utils/utils_string.h
  tools/PhyloTreeToolBox.h
  tools/Timer.h
  tools/BipartitionTools.h
  tools/Statistics.h
  # Tree reconciliation and recontruction
  tools/treeReconciliation/ReconcilableRootedTreeGenerator.h
  tools/treeReconciliation/ReconciledRootedTreeGenerator.h
  tools/treeReconciliation/TreeReconciliationConductor.h
  # Input/Output
  tools/IO/IO.h
  tools/IO/Newick.h
  tools/IO/Nhx.h
  tools/IO/XMLUtils.h
  tools/IO/PhyloXML.h
  tools/IO/RecPhyloTreeToSVG.h
  tools/IO/RefreshablePrinter.h
  # Algorithms
  algorithms/NeighborJoining.h
  algorithms/ProfileNJ.h
  algorithms/Polytomysolver.h
  # Constants
  Constants.h
  # ALE tools
  tools/ALE/ALE.h
  tools/ALE/exODT.h
  tools/ALE/fractionMissing.h
  tools/ALE/pairHasher.h
  tools/ALE/ALEevaluation.h
  tools/random.h
  tools/random.hpp
  )

set(treerecs-core_SOURCE_FILES
  # Containers
  containers/Cost.cpp
  containers/GeneMap.cpp
  containers/Grid.cpp
  containers/Table.cpp
  containers/Singleton.cpp
  containers/ReconcilableRootedTree.cpp
  containers/ReconciledRootedTree.cpp
  containers/NodeProperty.cpp
  containers/BipartitionList.cpp
  # libpll
  tools/LibpllEvaluation.cpp
  # Tools
  tools/SpeciesGeneMapper.cpp
  tools/utils.cpp
  tools/utils/utils_stream.cpp
  tools/utils/utils_math.cpp
  tools/utils/utils_containers.cpp
  tools/utils/utils_string.cpp
  tools/PhyloTreeToolBox.cpp
  tools/Timer.cpp
  tools/BipartitionTools.cpp
  tools/Statistics.cpp
  tools/SupportTools.cpp
  # Tree reconciliation and recontruction
  tools/treeReconciliation/ReconcilableRootedTreeGenerator.cpp
  tools/treeReconciliation/ReconciledRootedTreeGenerator.cpp
  tools/treeReconciliation/TreeReconciliationConductor.cpp
  # Input/Output
  tools/IO/IO.cpp
  tools/IO/Newick.cpp
  tools/IO/Nhx.cpp
  tools/IO/XMLUtils.cpp
  tools/IO/PhyloXML.cpp
  tools/IO/RecPhyloTreeToSVG.cpp
  tools/IO/RefreshablePrinter.cpp
  # Algorithms
  algorithms/NeighborJoining.cpp
  algorithms/ProfileNJ.cpp
  algorithms/Polytomysolver.cpp
  # ALE tools
  tools/ALE/ALE.cpp
  tools/ALE/exODT.cpp
  tools/ALE/fractionMissing.cpp
  tools/ALE/ALEevaluation.cpp
  tools/random.cpp)

add_library(treerecs-core
  ${treerecs-core_SOURCE_FILES}
  )
target_include_directories(treerecs-core
  PUBLIC ${BOOST_INCLUDE_DIRS}
  PUBLIC ${BPP-PHYL_INCLUDE_DIRS}
  PUBLIC ${TREERECS_INCLUDE_DIRS}
  PUBLIC ${SIMPLE-SVG_INCLUDE_DIRS}
  PUBLIC ${PLL_INCLUDE_DIRS}
  PUBLIC ${PLLMODULES_INCLUDE_DIRS}
  )
target_link_libraries(treerecs-core
  PUBLIC bpp-phyl
  PUBLIC ${PLL_LIBRARIES}
  PUBLIC ${PLLMODULES_LIBRARIES})

if (APPLE AND OPENMP_FOUND)
  target_link_libraries(treerecs-core
    PUBLIC ${OpenMP_libomp_LIBRARY})
endif(APPLE AND OPENMP_FOUND)
