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
# Last Change: 2019-10-02 08:15:25 d3g096


# -------------------------------------------------------------
# Set up for this subdirectory
# -------------------------------------------------------------

# Look in both source and binary trees for headers
include_directories(BEFORE 
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# -------------------------------------------------------------
# RBM
# -------------------------------------------------------------
add_executable(RBM RBM.f)

# -------------------------------------------------------------
# Installation
# -------------------------------------------------------------
install(TARGETS RBM DESTINATION bin)
