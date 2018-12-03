# -*- mode: cmake -*-

#
# Amanzi NetCDF Find Module
#
# Usage:
#    Control the search through NetCDF_DIR or setting environment variable
#    NetCDF_ROOT to the NetCDF installation prefix.
#    By default only searches for C library. To search for the C++ library
#    set 'CXX' in the COMPONENTS option in find_package(NetCDF) 
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    NetCDF_FOUND            (BOOL)       Flag indicating if NetCDF was found
#    NetCDF_INCLUDE_DIR      (PATH)       Path to the NetCDF include file
#    NetCDF_INCLUDE_DIRS     (LIST)       List of all required include files
#    NetCDF_C_LIBRARIES      (LIST)       List of all required libraries to link to 
#                                          the NetCDF C library
#    NetCDF_CXX_LIBRARIES    (LIST)       List of all required libraries to link to 
#                                          the NetCDF C++ library (If CXX is in COMPONENTS)
#
#    Additional variables set
#    NetCDF_C_LIBRARY        (FILE)       NetCDF C library
#    NetCDF_CXX_LIBRARY      (FILE)       NetCDF C++ library (If CXX is in the COMPONENTS)
#    NetCDF_LARGE_DIMS       (BOOL)       Checks the header files for size of 
#                                          NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VARS_DIMS
#                                          Returns TRUE if
#                                          NC_MAX_DIMS >= 655363
#                                          NC_MAX_VARS >= 524288
#                                          NC_MAX_VAR_DIMS >= 8
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)
include(AddImportedLibrary)
# include(AddPackageDependency)

# Macro to handle print out
macro(_netcdf_status mess)
  if (NOT NetCDF_FIND_QUIETLY)
    message(STATUS "${mess}")
  endif()
endmacro()

# Search for the nc-config binary
find_program(NetCDF_CONFIG_BINARY nc-config
             HINTS ${NetCDF_DIR} ${NetCDF_BIN_DIR}
	     PATH_SUFFIXES bin Bin
	     DOC "NetCDF configuration script")


if (NetCDF_C_LIBRARIES AND NetCDF_INCLUDE_DIRS )

  # Do nothing. Variables are set. No need to search again

else(NetCDF_C_LIBRARIES AND NetCDF_INCLUDE_DIRS)

  # Cache variables
  if (NetCDF_DIR)
    set(NetCDF_DIR "${NetCDF_DIR}" CACHE PATH "Path to search for NetCDF include and library files")
  endif()

  if (NetCDF_INCLUDE_DIR)
    set(NetCDF_INCLUDE_DIR "${NetCDF_INCLUDE_DIR}" CACHE PATH "Path to search for NetCDF include files")
  endif()

  if (NetCDF_LIBRARY_DIR)
    set(NetCDF_LIBRARY_DIR "${NetCDF_LIBRARY_DIR}" CACHE PATH "Path to search for NetCDF library files")
  endif()

  if (NOT NetCDF_INCLUDE_DIR )    
    # Search for include files
    # Search order preference:
    #  (1) NetCDF_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) NetCDF_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(netcdf_inc_names "netcdf.h")
    if (NetCDF_INCLUDE_DIR)
      if (EXISTS "${NetCDF_INCLUDE_DIR}")
  
        find_path(cdf_test_include_path
                  NAMES ${netcdf_inc_names}
                  HINTS ${NetCDF_INCLUDE_DIR}
                  NO_DEFAULT_PATH)

        if (NOT cdf_test_include_path)
          message(WARNING "Can not locate ${netcdf_inc_names} in ${NetCDF_INCLUDE_DIR}")
        endif()
        set(NetCDF_INCLUDE_DIR "${cdf_test_include_path}")
  
      else()
        message(WARNING "NetCDF_INCLUDE_DIR=${NetCDF_INCLUDE_DIR} does not exist")
        set(NetCDF_INCLUDE_DIR "NetCDF_INCLUDE_DIR-NOTFOUND")
      endif()
  
    else() 
  
      set(netcdf_inc_suffixes "include")
      if (NetCDF_DIR)
        if (EXISTS "${NetCDF_DIR}" )
  
          find_path(NetCDF_INCLUDE_DIR
                    NAMES ${netcdf_inc_names}
                    HINTS ${NetCDF_DIR}
                    PATH_SUFFIXES ${netcdf_inc_suffixes}
                    NO_DEFAULT_PATH)
  
        else()
          message(WARNING "NetCDF_DIR=${NetCDF_DIR} does not exist")
          set(NetCDF_INCLUDE_DIR "NetCDF_INCLUDE_DIR-NOTFOUND")
        endif()    
  
      else()
  
        find_path(NetCDF_INCLUDE_DIR
                  NAMES ${netcdf_inc_names}
                  PATH_SUFFIXES ${netcdf_inc_suffixes})
  
      endif()
    endif()
  

    if (NOT NetCDF_INCLUDE_DIR)
      message(WARNING "Can not locate NetCDF include directory")
    endif()
  endif(NOT NetCDF_INCLUDE_DIR)  
  
 
  if (NOT NetCDF_C_LIBRARIES)
    # Search for libraries 
    # Search order preference:
    #  (1) NetCDF_LIBRARY_DIR - check existence of path AND if the include files exist
    #  (2) NetCDF_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    if (NetCDF_LIBRARY_DIR)
      if (EXISTS "${NetCDF_LIBRARY_DIR}")
  
        find_library(_NetCDF_C_LIBRARY
                     NAMES netcdf
                     HINTS ${NetCDF_LIBRARY_DIR}
                     NO_DEFAULT_PATH)
  
      else()
        message(WARNING "NetCDF_LIBRARY_DIR=${NetCDF_LIBRARY_DIR} does not exist")
        set(NetCDF_C_LIBRARY NetCDF_C_LIBRARY-NOTFOUND)
      endif()
  
    else() 
      if (NetCDF_DIR)
        if (EXISTS "${NetCDF_DIR}")
  
          find_library(_NetCDF_C_LIBRARY
                       NAMES netcdf
                       HINTS ${NetCDF_DIR}
                       PATH_SUFFIXES "lib" "Lib"
                       NO_DEFAULT_PATH)
  
        else()
          message(WARNING "NetCDF_DIR=${NetCDF_DIR} does not exist")
          set(NetCDF_C_LIBRARY "NetCDF_C_LIBRARY-NOTFOUND")
        endif()    
  
      else()
  
        find_library(_NetCDF_C_LIBRARY
                     NAMES netcdf
                     PATH_SUFFIXES ${netcdf_lib_suffixes})
              
      endif()
    endif()
  
    # Define the NetCDF library targets
    if (_NetCDF_C_LIBRARY) 
      if (NetCDF_CONFIG_BINARY)
        execute_process(COMMAND "${NetCDF_CONFIG_BINARY}" "--has-hdf5"
                        RESULT_VARIABLE _ret_code
                        OUTPUT_VARIABLE _stdout
                        ERROR_VARIABLE  _stderr)
                          
        string(REGEX REPLACE "[\n\r ]" "" _hdf5_answer ${_stdout})
        _netcdf_status("${NetCDF_CONFIG_BINARY} --has-hdf5 returned '${_hdf5_answer}'")
        string(COMPARE EQUAL "${_hdf5_answer}" "yes" _has_hdf5)
        if (${_has_hdf5} ) 
          set(NetCDF_NEEDS_HDF5 TRUE)
          _netcdf_status( "NetCDF requires HDF5")
        else()
          set(NetCDF_NEEDS_HDF5 FALSE)
          _netcdf_status( "NetCDF does not require HDF5")
        endif()    
      endif()
  
      # If NetCDF was built with HDF5 then add that to the target properties
      # NetCDF calls HDF5 HL routines and not all HDF5 installations will have 
      # this library. Warn the user if HL is not detected. Use the 
      # HDF5_C_LIBRARIES to define link needs since it will contain hdf5 and
      # hdf5_hl.
      set(_NetCDF_C_LIBRARY_LINK_LIBS)
      if (NetCDF_NEEDS_HDF5) 

        find_package(HDF5 QUIET REQUIRED COMPONENTS C HL)
        if (HDF5_FOUND)
          if (NOT HDF5_HL_FOUND)
            message(WARNING "NetCDF calls the HDF5 HL library but this HDF5 does not have one")
          endif() 
          set(_NetCDF_C_LIBRARY_LINK_LIBS ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES})
          list(APPEND NetCDF_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})		      
        endif()    
      endif()
  
      # Now create the target 'netcdf'
      add_imported_library(netcdf LOCATION "${_NetCDF_C_LIBRARY}" LINK_LANGUAGES "C")
      if (_NetCDF_C_LIBRARY_LINK_LIBS)
        set_target_properties(netcdf PROPERTIES
                              INTERFACE_LINK_LIBRARIES "${_NetCDF_C_LIBRARY_LINK_LIBS}")
      endif()
      set(NetCDF_C_LIBRARY netcdf)
      set(NetCDF_C_LIBRARIES ${NetCDF_C_LIBRARY} ${_NetCDF_C_LIBRARY_LINK_LIBS})
      set(NetCDF_C_FOUND TRUE)
  
    endif(_NetCDF_C_LIBRARY)  

  endif(NOT NetCDF_C_LIBRARIES)      

endif(NetCDF_C_LIBRARIES AND NetCDF_INCLUDE_DIRS)    


# Large dimension check here
if (NetCDF_INCLUDE_DIR) 
  set(netcdf_h "${NetCDF_INCLUDE_DIR}/netcdf.h" )
  _netcdf_status( "NetCDF include file ${netcdf_h} will be searched for define values")

  file(STRINGS "${netcdf_h}" netcdf_max_dims_string REGEX "^#define NC_MAX_DIMS")
  string(REGEX REPLACE "[^0-9]" "" netcdf_max_dims "${netcdf_max_dims_string}")

  file(STRINGS "${netcdf_h}" netcdf_max_vars_string REGEX "^#define NC_MAX_VARS")
  string(REGEX REPLACE "[^0-9]" "" netcdf_max_vars "${netcdf_max_vars_string}")

  file(STRINGS "${netcdf_h}" netcdf_max_var_dims_string REGEX "^#define NC_MAX_VAR_DIMS")
  string(REGEX REPLACE "[^0-9]" "" netcdf_max_var_dims "${netcdf_max_var_dims_string}")

  if (
       ( (netcdf_max_dims EQUAL 65536)  OR (netcdf_max_dims GREATER 65536)  ) AND
       ( (netcdf_max_vars EQUAL 524288) OR (netcdf_max_vars GREATER 524288) ) AND
       ( (netcdf_max_var_dims EQUAL 8)  OR  (netcdf_max_var_dims GREATER 8) )
     )
    set(NetCDF_LARGE_DIMS TRUE)
  else()
    # this is not important to DHSVM
    # message(WARNING "The NetCDF found in ${NetCDF_DIR} does not have the correct NC_MAX_DIMS, NC_MAX_VARS and NC_MAX_VAR_DIMS\n"
    #                 "It may not be compatible with other TPL libraries such MOAB and ExodusII\n" )
    set(NetCDF_LARGE_DIMS FALSE)
  endif()
endif()    


# Determine the version
if (NOT NetCDF_VERSION)
  if (NetCDF_CONFIG_BINARY)
    execute_process(COMMAND "${NetCDF_CONFIG_BINARY}" "--version"
                    RESULT_VARIABLE _ret_code
                    OUTPUT_VARIABLE _stdout
                    ERROR_VARIABLE  _stderr)
                          
     string(REGEX REPLACE "[\n\r ]" "" _version_answer "${_stdout}")
     string(REGEX REPLACE "netCDF" "" NetCDF_VERSION "${_version_answer}")
  endif()   
endif()


# Search for the C++ libraries if requested
if (NetCDF_FIND_COMPONENTS)
  list(FIND NetCDF_FIND_COMPONENTS "CXX" _search_cxx_library)

  if (( NOT NetCDF_CXX_LIBRARIES ) AND 
       ( NOT ("${_search_cxx_library}" LESS "0") ) )

    if (NetCDF_LIBRARY_DIR)
      if (EXISTS "${NetCDF_LIBRARY_DIR}")
  
        find_library(_NetCDF_CXX_LIBRARY
                     NAMES netcdf_c++
                     HINTS ${NetCDF_LIBRARY_DIR}
                     NO_DEFAULT_PATH)
  
      else()
        message(WARNING "NetCDF_LIBRARY_DIR=${NetCDF_LIBRARY_DIR} does not exist")
        set(NetCDF_CXX_LIBRARY NetCDF_CXX_LIBRARY-NOTFOUND)
      endif()
  
    else() 
      if (NetCDF_DIR)
        if (EXISTS "${NetCDF_DIR}")
  
          find_library(_NetCDF_CXX_LIBRARY
                       NAMES netcdf_c++
                       HINTS ${NetCDF_DIR}
                       PATH_SUFFIXES "lib" "Lib"
                       NO_DEFAULT_PATH)
  
        else()
          message(WARNING "NetCDF_DIR=${NetCDF_DIR} does not exist")
          set(NetCDF_CXX_LIBRARY "NetCDF_CXX_LIBRARY-NOTFOUND")
        endif()    
  
      else()
  
        find_library(_NetCDF_CXX_LIBRARY
                     NAMES netcdf_c++
                     PATH_SUFFIXES ${netcdf_lib_suffixes})
              
      endif()
    endif()

    if (_NetCDF_CXX_LIBRARY)
      add_imported_library(netcdf_c++ LOCATION "${_NetCDF_C_LIBRARY}" LINK_LANGUAGES "CXX")
      set_target_properties(netcdf_c++ PROPERTIES
                           INTERFACE_LINK_LIBRARIES "${NetCDF_C_LIBRARY}")
      set(NetCDF_CXX_LIBRARY netcdf_c++)			 
      set(NetCDF_CXX_FOUND TRUE)
      set(NetCDF_CXX_LIBRARIES ${NetCDF_CXX_LIBRARY} ${NetCDF_C_LIBRARIES})
    endif()
  endif()

endif(NetCDF_FIND_COMPONENTS)    


# Send useful message if everything is found
find_package_handle_standard_args(NetCDF DEFAULT_MSG
                                  NetCDF_VERSION
                                  NetCDF_C_LIBRARIES
                                  NetCDF_INCLUDE_DIRS)

mark_as_advanced(
  NetCDF_INCLUDE_DIR
  NetCDF_INCLUDE_DIRS
  NetCDF_C_LIBRARY
  NetCDF_CXX_LIBRARY
  NetCDF_C_LIBRARIES
)
