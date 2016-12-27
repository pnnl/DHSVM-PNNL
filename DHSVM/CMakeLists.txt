# SUMMARY:      CMakeLists.txt
# USAGE:        Part of DHSVM

# AUTHOR:       William A. Perkins
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       william.perkins@pnl.gov
# ORIG-DATE:    Dec-2016
# DESCRIPTION:  Configruation rules for CMake
# DESCRIP-END.
# COMMENTS:
#
# Last Change: 2016-12-27 10:41:11 d3g096

# -------------------------------------------------------------
# General set up
# -------------------------------------------------------------

cmake_minimum_required (VERSION 2.8.12)
project (DHSVM C)

# this is where extra CMake functions can be found
list (APPEND CMAKE_MODULE_PATH "${DHSVM_SOURCE_DIR}/cmake")

# -------------------------------------------------------------
# Configuration options
# -------------------------------------------------------------

# Require NetCDF File I/O
option (DHSVM_USE_NETCDF "Look for NetCDF library and require its use" ON)

# Require X11 use on UNIX systems
if (UNIX)
  option (DHSVM_USE_X11 "Look for X11 libraries and require their use" ON)
endif(UNIX)

# Build test programs
option (DHSVM_BUILD_TESTS "Build several module test programs in addition to DHSVM" OFF)

# Limit calculations to snow pack only
option (DHSVM_SNOW_ONLY "Only simulate snow pack (no ET or infiltration)" OFF)
if (DHSVM_SNOW_ONLY) 
  add_definitions(-DSNOW_ONLY)
endif (DHSVM_SNOW_ONLY)

# -------------------------------------------------------------
# system-specific options
# -------------------------------------------------------------
include ( CheckIncludeFile )
include (CheckFunctionExists)

# On most UNIX-like platforms, the math library needs to be explicitly linked
if (UNIX)
  find_library(MATH_LIBRARY m)
endif(UNIX)

# -------------------------------------------------------------
# NetCDF is optional
# -------------------------------------------------------------
if (DHSVM_USE_NETCDF)
  find_package(NetCDF REQUIRED)
  add_definitions(-DHAVE_NETCDF)
  include_directories(AFTER ${NETCDF_INCLUDES})
endif (DHSVM_USE_NETCDF)

# -------------------------------------------------------------
# X11 is optional
# -------------------------------------------------------------
if (DHSVM_USE_X11)
  find_package(X11 REQUIRED)
  add_definitions(-DHAVE_X11)
  include_directories(AFTER ${X11_INCLUDE_DIR})
endif (DHSVM_USE_X11)

# -------------------------------------------------------------
# Use FLEX if it is available
# -------------------------------------------------------------
find_package(FLEX)

# -------------------------------------------------------------
# sub directories
# -------------------------------------------------------------
add_subdirectory(sourcecode)
add_subdirectory(program)
