# - Try to find GPTL libraries. This assumes that GPTL was built with
#   MPI enabled, but OpenMP disabled.
#
# Once done this will define
#
#  GPTL_FOUND        - system has GPTL
#  GPTL_INCLUDE_DIRS - include directories for GPTL
#  GPTL_LIBRARIES    - libraries for GPTL
#
# Variables used by this module. They can change the default behaviour and
# need to be set before calling find_package:
#
#  GTPL_DIR          - Prefix directory of the GTPL installation
#  GTPL_INCLUDE_DIR  - Include directory of the GTPL installation
#                          (set only if different from ${GTPL_DIR}/include)
#  GTPL_LIB_DIR      - Library directory of the GTPL installation
#                          (set only if different from ${GTPL_DIR}/lib)
#  GTPL_TEST_RUNS    - Skip tests building and running a test
#                          executable linked against GTPL libraries
#  GTPL_LIB_SUFFIX   - Also search for non-standard library names with the
#                          given suffix appended

#=============================================================================

find_path(GPTL_INCLUDE_DIR gptl.h
  HINTS ${GPTL_INCLUDE_DIR} ENV GPTL_INCLUDE_DIR ${GPTL_DIR} ENV GPTL_DIR
  PATH_SUFFIXES include
  DOC "Directory where the GPTL header files are located"
)
message(STATUS "GPTL_INCLUDE_DIR: ${GPTL_INCLUDE_DIR}")

find_library(GPTL_LIBRARY
  NAMES gptl GPTL${GPTL_LIB_SUFFIX}
  HINTS ${GPTL_LIB_DIR} ENV GPTL_LIB_DIR ${GPTL_DIR} ENV GPTL_DIR
  PATH_SUFFIXES lib
  DOC "Directory where the GPTL library is located"
)
message(STATUS "GPTL_LIBRARY: ${GPTL_LIBRARY}")

# Try compiling and running test program
if (GPTL_INCLUDE_DIR AND GPTL_LIBRARY)

  # Test requires MPI
  find_package(MPI QUIET REQUIRED)

  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${GPTL_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
  if (NOT MPI_LIBRARY OR NOT MPI_EXTRA_LIBRARY)
    set(CMAKE_REQUIRED_LIBRARIES 
      ${GPTL_LIBRARY} ${GPTL_EXTRA_LIBS}
    )
  else()
    set(CMAKE_REQUIRED_LIBRARIES 
      ${GPTL_LIBRARY} ${GPTL_EXTRA_LIBS} ${MPI_C_LIBRARIES}
    )
  endif()
endif()

# Build and run test program, maybe

set(gptl_test_src "
#include <gptl.h>

int main ()
{
  int ret;

  ret = GPTLinitialize ();
  ret = GPTLstart (\"total\");

  ret = GPTLstop (\"total\");
  ret = GPTLpr (0);
  return ret;
}
")

include(CheckCSourceRuns)
include(CheckCSourceCompiles)
check_c_source_runs("${gptl_test_src}" GPTL_TEST_RUNS)

unset(gptl_test_src)

# Standard package handling
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.2)
  find_package_handle_standard_args(GPTL
    REQUIRED_VARS GPTL_LIBRARY GPTL_INCLUDE_DIR GPTL_TEST_RUNS
    VERSION_VAR GPTL_VERSION_STRING)
else()
  find_package_handle_standard_args(GPTL
    REQUIRED_VARS GPTL_LIBRARY GPTL_INCLUDE_DIR GPTL_TEST_RUNS)
endif()

if(GPTL_FOUND)
  set(GPTL_LIBRARIES ${GPTL_LIBRARY})
  set(GPTL_INCLUDE_DIRS ${GPTL_INCLUDE_DIR})
endif()

mark_as_advanced(GPTL_INCLUDE_DIR GPTL_LIBRARY)
