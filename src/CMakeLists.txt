# CMake script for Moments++
# Authors:
#   Gustavo V. Barroso
# Created: 06/09/2022

SET(momentspp_CPP
  OptimizationWrapper.cpp
  Model.cpp
  Epoch.cpp
  SumStatsLibrary.cpp
  Drift.cpp
  Migration.cpp
  Mutation.cpp
  Recombination.cpp
  main.cpp
)

ADD_EXECUTABLE (momentspp ${momentspp_CPP})

SET(momentspp-targets momentspp)

FOREACH (target ${momentspp-targets})
   TARGET_LINK_LIBRARIES(${target} ${BPP_LIBS_SHARED} ${BOOST_LIBS_SHARED} ${EIGEN3_LIBS_SHARED})
   TARGET_LINK_LIBRARIES (${target} ${LIBS})
   #TARGET_LINK_LIBRARIES (${target} Threads::Threads)
ENDFOREACH (target)

INSTALL(TARGETS ${momentspp-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
