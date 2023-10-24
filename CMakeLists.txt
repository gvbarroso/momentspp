# CMake script for moments++
# Authors:
#   Gustavo V. Barroso
# Created: 28/07/2022

cmake_minimum_required (VERSION 3.5.0)
project (momentspp CXX)

# Telling the compiler to optimize based on CPU architecture leads to considerable speed-up.
# On Linux, add "-march=native" after -O3 flag
# On MacOS, clang may not support this option (depending on the version); add instead e.g. "-mcpu=apple-M1" after -O3 flag.
# clang --print-supported-cpus may give additional insight

SET(CMAKE_CXX_FLAGS "-std=c++20 -Weffc++ -Wshadow -Wall -Wextra -lgsl -O3 -fopenmp -march=native")

# Find dependencies (add install directory to search)
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

IF(VERBOSE)
  ADD_DEFINITIONS(-DVERBOSE)
ENDIF()

#SET(THREADS_PREFER_PTHREAD_FLAG ON)
#FIND_PACKAGE(Threads REQUIRED)

INCLUDE(GNUInstallDirs)
FIND_PACKAGE(bpp-core3 1.0.0 REQUIRED)

FIND_PACKAGE(GSL REQUIRED)  


SET(Boost_USE_MULTITHREADED OFF)
FIND_PACKAGE(Boost 1.67.0 REQUIRED COMPONENTS "iostreams")
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${Boost_LIBRARIES})
  MESSAGE("-- Boost libraries found here:")
  MESSAGE("   includes: ${Boost_INCLUDE_DIRS}")
  MESSAGE("   dynamic libraries: ${Boost_LIBRARIES}")
ENDIF()

FIND_PACKAGE(Eigen3 3.4 REQUIRED NO_MODULE)
IF(Eigen3_FOUND)
  INCLUDE_DIRECTORIES(${Eigen3_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${Eigen3_LIBRARIES})
  MESSAGE("-- Eigen3 libraries found.")
ENDIF()

FIND_PACKAGE(yaml-cpp REQUIRED)
IF(yaml-cpp_FOUND)
	INCLUDE_DIRECTORIES(${YAML_CPP_INCLUDE_DIR})
  SET(LIBS ${LIBS} ${YAML_CPP_LIBRARY_DIR})
  MESSAGE("-- yaml-cpp libraries found.")
  MESSAGE("   includes: ${YAML_CPP_INCLUDE_DIR}")
  MESSAGE("   dynamic libraries: ${YAML_CPP_LIBRARIES}")
ENDIF()

IF(NAKED_D)
  ADD_DEFINITIONS(-DNAKED_D)
ENDIF()

# Set the CMAKE_PREFIX_PATH for the find_library function when using non
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

