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
# Last Change: 2018-10-18 12:15:01 d3g096

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

enable_language(CXX)

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
  option (DHSVM_USE_X11 "Look for X11 libraries and require their use" OFF)
endif(UNIX)

# Build test programs
option (DHSVM_BUILD_TESTS "Build several module test programs in addition to DHSVM" OFF)

# Limit calculations to snow pack only
option (DHSVM_SNOW_ONLY 
  "Builds an addition executable, DHSVM_SNOW, that simulates snow only" OFF)

# Use the General Purpose Timing Library (GPTL) for performance timing
option (DHSVM_USE_GPTL, 
  "Use the General Purpose Timing Library (GPTL) for performance timing" OFF)

set(DHSVM_TIMING_LEVEL "1" CACHE STRING "The level of timing to use (1-4)")
if (NOT (DHSVM_TIMING_LEVEL MATCHES "^[1234]")) 
  message(ERROR "DHSVM_TIMING_LEVEL=\"${DHSVM_TIMING_LEVEL}\" not understood")
endif()

# -------------------------------------------------------------
# MPI is required
# -------------------------------------------------------------
message(STATUS "Checking MPI ...")
find_package (MPI REQUIRED)

message(STATUS "MPI_C_INCLUDES: ${MPI_C_INCLUDES}")
message(STATUS "MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")

# -------------------------------------------------------------
# Global Arrays is required
# -------------------------------------------------------------
message(STATUS "Checking GA ...")
find_package(GA REQUIRED)

# -------------------------------------------------------------
# GPTL is optional, but if asked for, it must be found
# -------------------------------------------------------------
if (DHSVM_USE_GPTL)
  find_package(GPTL REQUIRED)
endif (DHSVM_USE_GPTL)
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

