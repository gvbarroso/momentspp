# CMake script for moments++
# Authors:
#   Gustavo V. Barroso and Copilot
# Created: 28/07/2022

cmake_minimum_required(VERSION 3.14)
project(momentspp LANGUAGES CXX)

# Options
option(NativeBuild "Enable CPU-specific optimization flags" OFF)
option(NAKED_D "Enable naked D mode" OFF)
option(VERBOSE "Enable verbose mode" OFF)

# Base compiler flags
set(BASE_CXX_FLAGS "-std=c++20 -Weffc++ -Wshadow -Wall -Wextra -g -O3")

# Platform-specific optimization
if(NativeBuild)
  if(APPLE)
    message(STATUS "Native build on macOS: using -mcpu=apple-m1")
    set(NATIVE_OPT_FLAGS "-mcpu=apple-m1")
  else()
    message(STATUS "Native build on Linux: using -march=native")
    set(NATIVE_OPT_FLAGS "-march=native")
  endif()
else()
  message(STATUS "Native build disabled")
  set(NATIVE_OPT_FLAGS "")
endif()

# OpenMP setup for macOS
if(APPLE)
  set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
  set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/include")
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/include")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY "/opt/homebrew/lib/libomp.dylib")
endif()

find_package(OpenMP REQUIRED)

# Final compiler flags
set(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${NATIVE_OPT_FLAGS}")

# Macro definitions
if(VERBOSE)
  add_compile_definitions(VERBOSE)
endif()

if(NAKED_D)
  add_compile_definitions(NAKED_D)
endif()

# Dependencies
include(GNUInstallDirs)

if(CMAKE_INSTALL_PREFIX)
  list(APPEND CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}")
endif()

find_package(GSL REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(yaml-cpp REQUIRED)
find_package(Boost 1.67 REQUIRED CONFIG COMPONENTS iostreams)
find_library(MPFR_LIB mpfr REQUIRED)
find_library(GMP_LIB gmp REQUIRED)
find_package(bpp-core3 1.0.0 REQUIRED)

# mpreal via FetchContent
include(FetchContent)
FetchContent_Declare(
  mpreal
  GIT_REPOSITORY https://github.com/advanpix/mpreal.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(mpreal)

# Include directories
include_directories(
  ${Boost_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS}
  ${YAML_CPP_INCLUDE_DIR}
)

# Link libraries
set(LIBS
  ${Boost_LIBRARIES}
  ${Eigen3_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${GSL_LIBRARIES}
  ${MPFR_LIB}
  ${GMP_LIB}
)

# Source
add_subdirectory(src)

# Packaging
set(CPACK_PACKAGE_NAME "moments++")
set(CPACK_PACKAGE_VENDOR "Gustavo V. Barroso")
set(CPACK_PACKAGE_VERSION "0.0.1")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Fast and accurate prediction of two-locus statistics")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
  "/build/"
  "/\\.git/"
  "/\\.gitignore"
)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}")
include(CPack)

# Custom target (guarded to avoid duplication)
if(NOT TARGET dist)
  add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
endif()
