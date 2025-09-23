# CMake script for Moments++
# Authors:
#   Gustavo V. Barroso
# Created: 06/09/2022
# src/CMakeLists.txt

# 1) Explicit source lists
set(MOMENTS_SOURCES
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

set(TWOLOCUS_SOURCES
  twoLocusSim.cpp
)

# 2) Build momentspp
add_executable(momentspp ${MOMENTS_SOURCES})
target_compile_features(momentspp PRIVATE cxx_std_23)
target_include_directories(momentspp PRIVATE ${EXTRA_INCLUDES})
target_link_libraries(momentspp PRIVATE ${EXTRA_LIBS})

# 3) Build twoLocusSim
add_executable(twoLocusSim ${TWOLOCUS_SOURCES})
target_compile_features(twoLocusSim PRIVATE cxx_std_23)
target_include_directories(twoLocusSim PRIVATE ${EXTRA_INCLUDES})
target_link_libraries(twoLocusSim PRIVATE ${EXTRA_LIBS})

# 4) Install rules
install(
  TARGETS momentspp twoLocusSim
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
