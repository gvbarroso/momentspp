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
  #Migration.cpp
  Mutation.cpp
  #Admixture.cpp
  Selection.cpp
  Recombination.cpp
  Epoch.cpp
  Model.cpp
  OptimizationWrapper.cpp
  main.cpp
)

SET(twoLocusSim_CPP
  twoLocusSim.cpp
)

ADD_EXECUTABLE (momentspp ${momentspp_CPP})
ADD_EXECUTABLE (twoLocusSim ${twoLocusSim_CPP})

SET(momentspp-targets momentspp twoLocusSim)

target_include_directories(momentspp PRIVATE ${mpreal_SOURCE_DIR})

FOREACH (target ${momentspp-targets})
	TARGET_LINK_LIBRARIES(${target} ${BPP_LIBS_SHARED})
	TARGET_LINK_LIBRARIES(${target} ${YAML_CPP_LIBRARIES})
	TARGET_LINK_LIBRARIES(${target} GSL::gsl GSL::gslcblas)
	TARGET_LINK_LIBRARIES(${target} Boost::iostreams)
	TARGET_LINK_LIBRARIES(${target} ${MPFR_LIB} ${GMP_LIB})
        TARGET_LINK_LIBRARIES (${target} ${MPREAL_LIB})
ENDFOREACH (target)

INSTALL(TARGETS ${momentspp-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})

