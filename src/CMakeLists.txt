# CMake script for Moments++
# Authors:
#   Gustavo V. Barroso
# Created: 06/09/2022

# Source files for main executables
set(momentspp_CPP
  Log.cpp
  Graph.cpp
  Demes.cpp
  Data.cpp
  SumStatsLibrary.cpp
  AbstractOperator.cpp
  Drift.cpp
  Mutation.cpp
  Selection.cpp
  Recombination.cpp
  Epoch.cpp
  Model.cpp
  OptimizationWrapper.cpp
  main.cpp
)

set(twoLocusSim_CPP
  twoLocusSim.cpp
)

# Create executables
add_executable(momentspp ${momentspp_CPP})
add_executable(twoLocusSim ${twoLocusSim_CPP})

set(momentspp_targets momentspp twoLocusSim)

# Apply common settings to all targets
foreach(target ${momentspp_targets})
  target_compile_features(${target} PRIVATE cxx_std_20)

  target_include_directories(${target}
    PRIVATE
      ${mpreal_SOURCE_DIR}
      ${Boost_INCLUDE_DIRS}
      ${Eigen3_INCLUDE_DIRS}
      ${YAML_CPP_INCLUDE_DIR}
  )

  target_link_libraries(${target}
    PRIVATE
      ${BPP_LIBS_SHARED}
      ${YAML_CPP_LIBRARIES}
      GSL::gsl
      GSL::gslcblas
      Boost::iostreams
      ${MPFR_LIB}
      ${GMP_LIB}
      ${MPREAL_LIB}
  )
endforeach()

# Install binaries
install(TARGETS ${momentspp_targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
