# CMake script for moments++
# Authors:
#   Gustavo V. Barroso
# Created: 28/07/2022

cmake_minimum_required (VERSION 3.1.0)
project (momentspp CXX)

# Telling the compiler to optimize based on CPU architecture leads to considerable speed-up.
# On Linux, add "-march=native" after -O3 flag
# On MacOS, clang may not support this option (depending on the version); add instead e.g. "-mcpu=apple-M1" after -O3 flag.
# clang --print-supported-cpus may give additional insight

SET(CMAKE_CXX_FLAGS "-std=c++20 -Weffc++ -Wshadow -Wall -Wextra -ffast-math -O3 -march=native")

# Find dependencies (add install directory to search)
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

INCLUDE(GNUInstallDirs)
SET(THREADS_PREFER_PTHREAD_FLAG ON)
FIND_PACKAGE(Threads REQUIRED)

# Boost:
SET(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
SET(Boost_USE_MULTITHREADED OFF)
FIND_PACKAGE( Boost 1.67.0 REQUIRED COMPONENTS "iostreams")
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${Boost_LIBRARIES})
  MESSAGE("-- Boost libraries found here:")
  MESSAGE("   includes: ${Boost_INCLUDE_DIRS}")
  MESSAGE("   dynamic libraries: ${Boost_LIBRARIES}")
ENDIF()

# BLAS:
FIND_PACKAGE(BLAS REQUIRED)
IF(BLAS_FOUND)
  INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${BLAS_LIBRARIES})
  MESSAGE("-- BLAS libraries found here:")
  MESSAGE("   includes: ${BLAS_INCLUDE_DIRS}")
  MESSAGE("   dynamic libraries: ${BLAS_LIBRARIES}")
ENDIF()

# here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  # start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})

  SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
  ELSE()
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF()
  # add the dependency:
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
  
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

# Subdirectories
ADD_SUBDIRECTORY(src)

# Packager
SET(CPACK_PACKAGE_NAME "moments++")
SET(CPACK_PACKAGE_VENDOR "Gustavo V. Barroso")
SET(CPACK_PACKAGE_VERSION "0.0.1")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "PopGen inference using two-locus statistics")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
#SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
#SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES 
       "/build/" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

# This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)



