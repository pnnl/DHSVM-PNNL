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
# Last Change: 2020-05-12 09:45:58 d3g096

# -------------------------------------------------------------
# Set up for this subdirectory
# -------------------------------------------------------------

# Look in both source and binary trees for headers
include_directories(BEFORE 
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

# Flags for build options
if (DHSVM_DUMP_TOPO)
  add_definitions(-DTOPO_DUMP)
endif (DHSVM_DUMP_TOPO)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("debug mode: turning DEBUG messages on")
    add_definitions(-DDEBUG=1)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

# tablio.c generation
if (FLEX_FOUND)
  FLEX_TARGET(tableio 
    tableio.lex  
    ${CMAKE_CURRENT_BINARY_DIR}/tableio.c
    COMPILE_FLAGS -Ptable_yy
    )
  set_directory_properties(
    PROPERTIES 
    ADDITIONAL_MAKE_CLEAN_FILES "${FLEX_tableio_OUTPUTS}"
    )
else(FLEX_FOUND)
  set(FLEX_tableio_OUTPUTS tableio.c)
endif(FLEX_FOUND)

# -------------------------------------------------------------
# Binary I/O library
# -------------------------------------------------------------
add_library(BinIO STATIC
  Calendar.c 
  InitFileIO.c
  FileIOBin.c 
  FileIONetCDF.c
  Files.c 
  InitArray.c 
  ReportError.c
  SizeOfNT.c 
  VarID.c
)

# -------------------------------------------------------------
# DHSVM Target
# -------------------------------------------------------------

set(DHSVM_SRC
  AdjustStorage.c
  Aggregate.c
  AggregateRadiation.c
  Avalanche.c
  CalcAerodynamic.c
  CalcAvailableWater.c
  CalcDistance.c
  CalcEffectiveKh.c
  CalcKhDry.c
  CalcKinViscosity.c
  CalcSatDensity.c
  CalcSnowAlbedo.c
  CalcSolar.c
  CalcTotalWater.c
  CalcTransmissivity.c
  CalcWeights.c
  CanopyGapEnergyBalance.c
  CanopyGapRadiation.c
  CanopyResistance.c
  ChannelState.c
  CheckOut.c
  CutBankGeometry.c
  DHSVMChannel.c
  Desorption.c
  DistributeSatflow.c
  Draw.c
  EvalExponentIntegral.c
  EvapoTranspiration.c
  ExecDump.c
  FinalMassBalance.c
  GetInit.c
  GetMetData.c
  InArea.c
  InitAggregated.c
  InitConstants.c
  InitDump.c
  InitInterpolationWeights.c
  InitMetMaps.c
  InitMetSources.c
  InitModelState.c
  InitNetwork.c
  InitNewMonth.c
  InitParameterMaps.c
  InitSnowMap.c
  InitTables.c
  InitTerrainMaps.c
  InitUnitHydrograph.c
  InitXGraphics.c
  InterceptionStorage.c
  IsStationLocation.c
  LapseT.c
  LookupTable.c
  MainDHSVM.c
  MakeLocalMetData.c
  MassBalance.c
  MassEnergyBalance.c
  MassRelease.c
  MaxRoadInfiltration.c
  NoEvap.c
  RadiationBalance.c
  ReadMetRecord.c
  ReadRadarMap.c
  ResetAggregate.c
  RiparianShading.c
  RootBrent.c
  Round.c
  RouteSubSurface.c
  RouteSurface.c
  SatVaporPressure.c
  SensibleHeatFlux.c
  SeparateRadiation.c
  SlopeAspect.c
  SnowInterception.c
  SnowMelt.c
  SnowPackEnergyBalance.c
  SnowStats.c
  SoilEvaporation.c
  StabilityCorrection.c
  StoreModelState.c
  SurfaceEnergyBalance.c
  UnsaturatedFlow.c
  WaterTableDepth.c
  channel.c
  channel_grid.c
  channel_complt.c
  deg2utm.c
  equal.c
  errorhandler.c
  globals.c
  ${FLEX_tableio_OUTPUTS}
)

add_executable(DHSVM
  ${DHSVM_SRC}
)

target_link_libraries(DHSVM
  BinIO
  ${NETCDF_LIBRARIES}
  ${X11_LIBRARIES}
  ${MATH_LIBRARY}
)

if(DHSVM_SNOW_ONLY)

  add_executable(DHSVM_SNOW
    ${DHSVM_SRC}
    )

  target_compile_definitions(DHSVM_SNOW
    PUBLIC SNOW_ONLY=1
    )

  target_link_libraries(DHSVM_SNOW
    BinIO
    ${NETCDF_LIBRARIES}
    ${X11_LIBRARIES}
    ${MATH_LIBRARY}
    )
endif(DHSVM_SNOW_ONLY)

# -------------------------------------------------------------
# channel_test
# -------------------------------------------------------------
# if (DHSVM_BUILD_TESTS)
#   add_executable(channel_test
#     channel.c
#     errorhandler.c
#     ${FLEX_tableio_OUTPUTS}
#     )
#   if (MATH_LIBRARY)
#     target_link_libraries(channel_test
#       ${MATH_LIBRARY}
#       )
#   endif (MATH_LIBRARY)
#   set_target_properties(channel_test
#     PROPERTIES
#     COMPILE_DEFINITIONS "TEST_CHANNEL=1"
#     )
# endif (DHSVM_BUILD_TESTS)

# -------------------------------------------------------------
# channel_grid_test
# -------------------------------------------------------------
# if (DHSVM_BUILD_TESTS)
#   add_executable(channel_grid_test
#     channel_grid.c
#     tableio.c
#     )
#   set_target_properties(channel_grid_test
#     PROPERTIES
#     COMPILE_DEFINITIONS "TEST_CHANNEL_GRID=1"
#     )
# endif (DHSVM_BUILD_TESTS)

# -------------------------------------------------------------
# report_error_test
# -------------------------------------------------------------
if (DHSVM_BUILD_TESTS)
  add_executable(report_error_test
    ReportError.c
    )
  set_target_properties(report_error_test
    PROPERTIES
    COMPILE_DEFINITIONS "TEST_REPORTERROR=1"
    )
endif(DHSVM_BUILD_TESTS)

# -------------------------------------------------------------
# calendar_test
# -------------------------------------------------------------
if (DHSVM_BUILD_TESTS)
  add_executable(calendar_test 
    Calendar.c
    ReportError.c
    Round.c
    equal.c
    )
  target_link_libraries(calendar_test
    ${MATH_LIBRARY}
    )
  set_target_properties(calendar_test
    PROPERTIES
    COMPILE_DEFINITIONS "TEST_CALENDAR=1"
    )
endif (DHSVM_BUILD_TESTS)

# -------------------------------------------------------------
# error_handler_test
# -------------------------------------------------------------
if (DHSVM_BUILD_TESTS)
  add_executable(error_handler_test
    errorhandler.c
    )
  set_target_properties(error_handler_test
    PROPERTIES
    COMPILE_DEFINITIONS "TEST_EH_MAIN=1"
    )
endif(DHSVM_BUILD_TESTS)

# -------------------------------------------------------------
# varid_test
# -------------------------------------------------------------
if (DHSVM_BUILD_TESTS)
  add_executable(varid_test
    VarID.c
    ReportError.c
    )
  set_target_properties(varid_test
    PROPERTIES
    COMPILE_DEFINITIONS "TEST_VARID=1"
    )  
endif(DHSVM_BUILD_TESTS)

# -------------------------------------------------------------
# NetCDF test
# -------------------------------------------------------------
if (DHSVM_BUILD_TESTS AND NETCDF_FOUND) 
  add_executable(netcdf_test
    FileIONetCDF.c
    ReportError.c
    VarID.c
    )
  target_link_libraries(netcdf_test
    ${NETCDF_LIBRARIES}
    ${MATH_LIBRARY}
    )
  set_target_properties(netcdf_test
    PROPERTIES
    COMPILE_DEFINITIONS "TEST_FILEIONETCDF=1"
    )
endif()

# -------------------------------------------------------------
# table_test
# -------------------------------------------------------------
if (DHSVM_BUILD_TESTS) 
  add_executable(table_test
    ${FLEX_tableio_OUTPUTS}
    errorhandler.c
    )
  set_target_properties(table_test
    PROPERTIES
    COMPILE_DEFINITIONS "TEST_TABLE=1"
    )
endif (DHSVM_BUILD_TESTS)

# -------------------------------------------------------------
# equal_test
# -------------------------------------------------------------
if (DHSVM_BUILD_TESTS) 
  add_executable(equal_test
    equal.c
    )
  set_target_properties(equal_test
    PROPERTIES
    COMPILE_DEFINITIONS "TEST_EQUAL=1"
    )
endif (DHSVM_BUILD_TESTS)
