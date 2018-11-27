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
# Last Change: 2018-11-08 14:26:11 d3g096

# -------------------------------------------------------------
# Set up for this subdirectory
# -------------------------------------------------------------

# Look in the both source and binary DHSVM trees for headers
include_directories(BEFORE 
  ${CMAKE_CURRENT_SOURCE_DIR}/../sourcecode
  ${CMAKE_CURRENT_BINARY_DIR}/../sourcecode
)

if (DHSVM_USE_NETCDF)
  set(NetCDF_LIBS ${NetCDF_C_LIBRARIES} ${MATH_LIBRARY})
  if (HDF5_IS_PARALLEL)
    set(NetCDF_LIBS ${NetCDF_LIBS} ${MPI_C_LIBRARIES})
  endif(HDF5_IS_PARALLEL)
else(DHSVM_USE_NETCDF)
  set(NetCDF_LIBS "")
endif(DHSVM_USE_NETCDF)

# -------------------------------------------------------------
# libraries of common IO source
# -------------------------------------------------------------
add_library(locBinIO STATIC
  FileIOBin.c 
  Files.c 
  InitArray.c
  ReportError.c
  SizeOfNT.c
  Calendar.c
  )
if (DHSVM_USE_NETCDF)
  add_library(NetCDFIO STATIC
    FileIONetCDF.c 
    Files.c 
    InitArray.c
    ReportError.c
    SizeOfNetCDF.c
    Calendar.c
    )
endif (DHSVM_USE_NETCDF)    

# -------------------------------------------------------------
# myconvert
# -------------------------------------------------------------
add_executable(myconvert
  myconvert.c
  )

# -------------------------------------------------------------
# MakeModelState
# -------------------------------------------------------------
add_executable(MakeModelStateBin
  MakeModelStateBin.c
)

target_link_libraries(MakeModelStateBin
  locBinIO
  ${MATH_LIBRARY}
)

if (DHSVM_USE_NETCDF)
  add_executable(MakeModelStateNetCDF
    MakeModelStateNetCDF.c
    )
  target_link_libraries(MakeModelStateNetCDF
    NetCDFIO
    ${NetCDF_LIBS}
    )
endif (DHSVM_USE_NETCDF)

# -------------------------------------------------------------
# skyview
# -------------------------------------------------------------
add_executable(skyviewBin
  skyviewBin.c 
  )
target_link_libraries(skyviewBin
  locBinIO
  ${MATH_LIBRARY}
  )

if (DHSVM_USE_NETCDF)
  add_executable(skyviewNetCDF
    skyviewNetCDF.c 
    )
  target_link_libraries(skyviewNetCDF
    NetCDFIO
    ${NetCDF_LIBS}
    )
endif (DHSVM_USE_NETCDF)

# -------------------------------------------------------------
# average_shadow
# -------------------------------------------------------------
add_executable(average_shadow_bin
  average_shadow_bin.c 
  )
target_link_libraries(average_shadow_bin
  locBinIO
  ${MATH_LIBRARY}
  )

if (DHSVM_USE_NETCDF)
  add_executable(average_shadow_netcdf
    average_shadow_netcdf.c
    )
  target_link_libraries(average_shadow_netcdf
    NetCDFIO
    ${NetCDF_LIBS}
    )
endif (DHSVM_USE_NETCDF)




# -------------------------------------------------------------
# find_nearest_channel
# -------------------------------------------------------------
add_executable(find_nearest_channel_bin
  find_nearest_channel_bin.c 
  )
  target_link_libraries(find_nearest_channel_bin
    locBinIO
    ${MATH_LIBRARY}
    )

if (DHSVM_USE_NETCDF)
    add_executable(find_nearest_channel_netcdf
    find_nearest_channel_netcdf.c 
    )
  target_link_libraries(find_nearest_channel_netcdf
    NetCDFIO
    ${NetCDF_LIBS}
    )
endif (DHSVM_USE_NETCDF)

# -------------------------------------------------------------
# make_shade_maps
# -------------------------------------------------------------
add_executable(make_shade_maps_bin
  make_shade_maps_bin.c
  )
  target_link_libraries(make_shade_maps_bin
    locBinIO
    ${MATH_LIBRARY}
    )

if (DHSVM_USE_NETCDF)
  add_executable(make_shade_maps_netcdf
    make_shade_maps_netcdf.c
    )
  target_link_libraries(make_shade_maps_netcdf
    NetCDFIO
    ${NetCDF_LIBS}
    )
endif (DHSVM_USE_NETCDF)

# -------------------------------------------------------------
# fixroads
# -------------------------------------------------------------
add_executable(fixroads
  fixroads.c
)
target_link_libraries(make_shade_maps_bin
  ${MATH_LIBRARY}
)

# -------------------------------------------------------------
# fill_sinks_dhsvm
# -------------------------------------------------------------
add_executable(fill_sinks_dhsvm
  FILL_SINKS_DHSVM.c
)

# -------------------------------------------------------------
# WriteConstantMapBin
# -------------------------------------------------------------
add_executable(WriteConstantMapBin
  WriteConstantMap.c
)  
target_link_libraries(WriteConstantMapBin
  locBinIO
  ${MATH_LIBRARY}
  )

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(TARGETS
  myconvert
  MakeModelStateBin
  skyviewBin
  average_shadow_bin
  find_nearest_channel_bin
  make_shade_maps_bin 
  fixroads 
  fill_sinks_dhsvm
  DESTINATION bin
)

if (DHSVM_USE_NETCDF)
  install(TARGETS
    MakeModelStateNetCDF
    skyviewNetCDF
    average_shadow_bin
    make_shade_maps_netcdf 
    find_nearest_channel_netcdf
    DESTINATION bin
  )
endif (DHSVM_USE_NETCDF)
