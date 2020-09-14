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
# Last Change: 2020-05-12 09:47:01 d3g096

# -------------------------------------------------------------
# General set up
# -------------------------------------------------------------

cmake_minimum_required (VERSION 2.8.12)
project (DHSVM)

# this is where extra CMake functions can be found
list (APPEND CMAKE_MODULE_PATH "${DHSVM_SOURCE_DIR}/cmake")

# -------------------------------------------------------------
# Configuration options
# -------------------------------------------------------------

enable_language(C)

# Option to dump topography maps, for debugging
option (DHSVM_DUMP_TOPO
  "Make DHSVM dump input topography maps (in a non-standard way) for debugging"
)

# Build RBM if desired
option (DHSVM_USE_RBM "Build RBM and related programs" OFF)
if (DHSVM_USE_RBM)
  enable_language(Fortran)
endif (DHSVM_USE_RBM)

# Require NetCDF File I/O
option (DHSVM_USE_NETCDF "Look for NetCDF library and require its use" ON)

# Require X11 use on UNIX systems
if (UNIX)
  option (DHSVM_USE_X11 "Look for X11 libraries and require their use" ON)
endif(UNIX)

# Build test programs
option (DHSVM_BUILD_TESTS "Build several module test programs in addition to DHSVM" OFF)

# Limit calculations to snow pack only
option (DHSVM_SNOW_ONLY "Builds an addition executable, DHSVM_SNOW, that simulates snow only" OFF)

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
add_subdirectory(DHSVM/sourcecode)
add_subdirectory(DHSVM/program)
if (DHSVM_USE_RBM) 
  add_subdirectory(RBM)
  add_subdirectory(Create)
endif (DHSVM_USE_RBM)

