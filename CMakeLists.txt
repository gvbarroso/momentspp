# CMake script for moments++
# Authors:
#   Gustavo V. Barroso
# Created: 28/07/2022

cmake_minimum_required(VERSION 3.23)
project(momentspp VERSION 0.0.1 LANGUAGES CXX)

#------------------------------------------------------------------------------
# 0) Bio++ install location
#------------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/.local")

#------------------------------------------------------------------------------
# 1) Suppress Boost CONFIG warning
#------------------------------------------------------------------------------
if (POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

#------------------------------------------------------------------------------
# 2) User options & C++ standard
#------------------------------------------------------------------------------
option(NativeBuild "Enable CPU-specific optimization flags" OFF)
option(VERBOSE     "Enable verbose mode"                OFF)
option(NAKED_D     "Enable naked D mode"                OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#------------------------------------------------------------------------------
# 3) Compiler warnings & optimization flags
#------------------------------------------------------------------------------
set(BASE_FLAGS
  -Wall
  -Wextra
  -pedantic
  -Wshadow
  -Weffc++
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND BASE_FLAGS -g)
else()
  list(APPEND BASE_FLAGS -O3)
endif()

if (NativeBuild)
  message(STATUS "Native build: adding -march=native")
  list(APPEND BASE_FLAGS -march=native)
endif()

add_compile_options(${BASE_FLAGS})

if (VERBOSE)
  add_compile_definitions(VERBOSE)
endif()

if (NAKED_D)
  add_compile_definitions(NAKED_D)
endif()

#------------------------------------------------------------------------------
# 4) Boost (CONFIG mode → Boost::iostreams)
#------------------------------------------------------------------------------
set(Boost_FIND_CMAKE CONFIG)
find_package(Boost 1.67 REQUIRED COMPONENTS iostreams)

#------------------------------------------------------------------------------
# 5) OpenMP (module first, fallback for macOS)
#------------------------------------------------------------------------------
find_package(OpenMP)

if (OpenMP_CXX_FOUND)
  add_compile_options(${OpenMP_CXX_FLAGS})
  set(OPENMP_LIBS OpenMP::OpenMP_CXX)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  find_path(OMP_INCLUDE_DIR omp.h
    HINTS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include
  )
  find_library(OMP_LIBRARY NAMES omp libomp
    HINTS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib
  )
  if (OMP_INCLUDE_DIR AND OMP_LIBRARY)
    add_compile_options(-Xpreprocessor -fopenmp -I${OMP_INCLUDE_DIR})
    set(OPENMP_LIBS ${OMP_LIBRARY})
  else()
    message(FATAL_ERROR "OpenMP not found. Install with: brew install libomp")
  endif()
endif()

#------------------------------------------------------------------------------
# 6) GSL, Eigen3, yaml-cpp, MPFR, GMP
#------------------------------------------------------------------------------
find_package(GSL       REQUIRED)
find_package(Eigen3    3.4 REQUIRED NO_MODULE)
find_package(yaml-cpp  REQUIRED)

find_library(MPFR_LIB mpfr REQUIRED PATHS "$ENV{HOME}/.local/lib")
find_library(GMP_LIB  gmp  REQUIRED PATHS "$ENV{HOME}/.local/lib")

#------------------------------------------------------------------------------
# 7) Bio++ core3 (includes numeric3)
#------------------------------------------------------------------------------
find_library(BPP_CORE3_LIB bpp-core3 REQUIRED PATHS "$ENV{HOME}/.local/lib")
message(STATUS "Resolved bpp-core3 library: ${BPP_CORE3_LIB}")

#------------------------------------------------------------------------------
# 8) mpreal (header-only)
#------------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
  mpreal
  GIT_REPOSITORY https://github.com/advanpix/mpreal.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(mpreal)

#------------------------------------------------------------------------------
# 9) Collect include paths & libraries
#------------------------------------------------------------------------------
set(EXTRA_INCLUDES
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${YAML_CPP_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${mpreal_SOURCE_DIR}
)

set(EXTRA_LIBS
  Boost::iostreams
  yaml-cpp
  GSL::gsl
  GSL::gslcblas
  ${MPFR_LIB}
  ${GMP_LIB}
  ${OPENMP_LIBS}
  ${BPP_CORE3_LIB}
)

#------------------------------------------------------------------------------
# 10) Build & install src/
#------------------------------------------------------------------------------
add_subdirectory(src)

include(GNUInstallDirs)
install(
  TARGETS momentspp twoLocusSim
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(CPack)

