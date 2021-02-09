# This is the main CMake build file for NCEPLIBS-ncio, a NOAA library
# to for reading and writing netCDF data.
#
# 2/9/21

cmake_minimum_required(VERSION 3.15)

file(STRINGS "VERSION" pVersion)

project(ncio
        VERSION ${pVersion}
        LANGUAGES Fortran)

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED Fortran)

add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Determine whether or not to generate documentation.
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  add_subdirectory(docs)  
endif()
