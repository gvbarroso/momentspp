# CMake script for Moments++
# Authors:
#   Gustavo V. Barroso
# Created: 06/09/2022

SET(momentspp_CPP
  Log.cpp
  Graph.cpp
  Demes.cpp
  Data.cpp
  SumStatsLibrary.cpp
  AbstractOperator.cpp
  Drift.cpp
  Migration.cpp
  Mutation.cpp
  Admixture.cpp
  Selection.cpp
  Recombination.cpp
  Epoch.cpp
  Model.cpp
  OptimizationWrapper.cpp
  main.cpp
)

ADD_EXECUTABLE (momentspp ${momentspp_CPP})

SET(momentspp-targets momentspp)

FOREACH (target ${momentspp-targets})
	TARGET_LINK_LIBRARIES(${target} ${BPP_LIBS_SHARED})
	TARGET_LINK_LIBRARIES(${target} ${YAML_CPP_LIBRARIES})
        #TARGET_LINK_LIBRARIES (${target} Threads::Threads)
ENDFOREACH (target)

INSTALL(TARGETS ${momentspp-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})

