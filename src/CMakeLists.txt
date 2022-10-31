# CMake script for Moments++
# Authors:
#   Gustavo V. Barroso
# Created: 06/09/2022

SET(momentspp_CPP
  Log.cpp
  PolymorphismData.cpp
  SumStatsLibrary.cpp
  Drift.cpp
  Migration.cpp
  Mutation.cpp
  Recombination.cpp
  Epoch.cpp
  Model.cpp
  OptimizationWrapper.cpp
  Demes.cpp
  main.cpp
)

ADD_EXECUTABLE (momentspp ${momentspp_CPP})
add_executable(yaml_cpp main.cpp)
target_link_libraries(yaml_cpp ${yaml-cpp_LIBRARIES})

SET(momentspp-targets momentspp)

FOREACH (target ${momentspp-targets})
   TARGET_LINK_LIBRARIES(${target} ${BPP_LIBS_SHARED} ${BOOST_LIBS_SHARED} ${EIGEN3_LIBS_SHARED})
   TARGET_LINK_LIBRARIES (${target} ${LIBS})
   #TARGET_LINK_LIBRARIES (${target} Threads::Threads)
ENDFOREACH (target)

INSTALL(TARGETS ${momentspp-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
