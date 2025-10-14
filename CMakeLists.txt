# CMake script for Moments++
# Authors:
#   Gustavo V. Barroso
# Created: 06/09/2022

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
option(DEBUG       "Enable debug mode"                      OFF)
option(NAKED_D     "Enable naked D mode"                    OFF)

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

if (DEBUG)
  message(STATUS "DEBUG build: adding -DDEBUG=ON")
  add_compile_definitions(DEBUG)
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
# 5) OpenMP detection with macOS-specific fallback
#------------------------------------------------------------------------------
# Goal: provide OpenMP support on Linux/other platforms via find_package(OpenMP),
# and on macOS prefer Homebrew libomp. Avoid calling find_package(OpenMP) on macOS
# so CMake won't list OpenMP as "not found" in its summary.

set(OPENMP_COMPILE_FLAGS "")
set(OPENMP_LINK_FLAGS "")
set(OPENMP_LIB_TARGET "")

if (APPLE)
  # macOS: prefer Homebrew libomp (keg-only)
  set(HOMEBREW_LIBOMP_PREFIX "/opt/homebrew/opt/libomp")
  if (EXISTS "${HOMEBREW_LIBOMP_PREFIX}/include/omp.h" AND EXISTS "${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib")
    message(STATUS "Using Homebrew libomp at ${HOMEBREW_LIBOMP_PREFIX}")

    # Compile flags must be separate list entries
    set(OPENMP_COMPILE_FLAGS -Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include)
    set(OPENMP_LINK_FLAGS -L${HOMEBREW_LIBOMP_PREFIX}/lib)
    set(OPENMP_LIB "${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib")

    # Create an imported target so other targets can link uniformly
    if (NOT TARGET OpenMP::HomebrewLibomp)
      add_library(OpenMP::HomebrewLibomp UNKNOWN IMPORTED)
      set_target_properties(OpenMP::HomebrewLibomp PROPERTIES
        IMPORTED_LOCATION "${OPENMP_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${HOMEBREW_LIBOMP_PREFIX}/include"
      )
    endif()
    set(OPENMP_LIB_TARGET OpenMP::HomebrewLibomp)

    # Create helper interface that carries flags and link target
    add_library(momentspp_openmp_flags INTERFACE)
    target_compile_options(momentspp_openmp_flags INTERFACE ${OPENMP_COMPILE_FLAGS})
    target_link_options(momentspp_openmp_flags INTERFACE ${OPENMP_LINK_FLAGS})
    target_link_libraries(momentspp_openmp_flags INTERFACE OpenMP::HomebrewLibomp)

    # Mark a variable so downstream code can detect OpenMP support
    set(MOMENTSPP_HAS_OPENMP TRUE CACHE INTERNAL "Moments++: OpenMP support available")
  else()
    message(STATUS "Homebrew libomp not found at ${HOMEBREW_LIBOMP_PREFIX}; continuing without OpenMP on macOS")
    set(MOMENTSPP_HAS_OPENMP FALSE CACHE INTERNAL "Moments++: OpenMP support available")
  endif()
else()
  # Non-macOS: use standard find_package(OpenMP QUIET). If found, create helper interface.
  find_package(OpenMP QUIET)
  if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found via find_package: ${OpenMP_CXX_FLAGS}")
    set(OPENMP_COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    set(OPENMP_LIB_TARGET OpenMP::OpenMP_CXX)

    add_library(momentspp_openmp_flags INTERFACE)
    # CMake's OpenMP flags are often a single string; ensure they are passed correctly.
    foreach(_flag IN LISTS OPENMP_COMPILE_FLAGS)
      target_compile_options(momentspp_openmp_flags INTERFACE ${_flag})
    endforeach()
    target_link_libraries(momentspp_openmp_flags INTERFACE OpenMP::OpenMP_CXX)

    set(MOMENTSPP_HAS_OPENMP TRUE CACHE INTERNAL "Moments++: OpenMP support available")
  else()
    message(STATUS "OpenMP not found on this platform; continuing without OpenMP")
    set(MOMENTSPP_HAS_OPENMP FALSE CACHE INTERNAL "Moments++: OpenMP support available")
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
# 9) Spectra (header-only)
#------------------------------------------------------------------------------
FetchContent_Declare(
  Spectra
  GIT_REPOSITORY https://github.com/yixuan/spectra.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(Spectra)

#------------------------------------------------------------------------------
# 10) Collect include paths & libraries (used when creating targets)
#------------------------------------------------------------------------------
set(EXTRA_INCLUDES
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${YAML_CPP_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${mpreal_SOURCE_DIR}
  ${spectra_SOURCE_DIR}/include
)

set(EXTRA_LIBS
  Boost::iostreams
  yaml-cpp
  GSL::gsl
  GSL::gslcblas
  ${MPFR_LIB}
  ${GMP_LIB}
  ${BPP_CORE3_LIB}
)

#------------------------------------------------------------------------------
# 11) macOS rpath suppression for relocatable installs (optional)
#------------------------------------------------------------------------------
if(APPLE)
  set(CMAKE_INSTALL_RPATH "")
  set(CMAKE_SKIP_RPATH TRUE)
endif()

#------------------------------------------------------------------------------
# 12) Build & install src/
#------------------------------------------------------------------------------
add_subdirectory(src)

include(GNUInstallDirs)
install(
  TARGETS momentspp twoLocusSim
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(CPack)

#------------------------------------------------------------------------------
# NOTE: how to consume OpenMP in src/CMakeLists.txt
#   if (TARGET momentspp_openmp_flags)
#     target_link_libraries(momentspp PRIVATE momentspp_openmp_flags)
#   endif()
#
# Use the cached variable MOMENTSPP_HAS_OPENMP to check for availability in downstream logic.

