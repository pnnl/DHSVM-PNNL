#!/bin/sh
# SUMMARY:      CMakeLists.txt
# USAGE:        Part of DHSVM

# AUTHOR:       William A. Perkins
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       william.perkins@pnl.gov
# ORIG-DATE:    Dec-2016
# DESCRIPTION:  Example DHSVM CMake configuration for some systems
# DESCRIP-END.
# COMMENTS:
#
# Last Change: 2020-05-12 10:09:21 d3g096

set -xue

# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r] [-t] [name]"

opts=`getopt dr $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi
set -- $opts

build="RelWithDebInfo"
for o in $*; do
    case $o in
        -d)
            build="Debug"
            shift
            ;;
        -r)
            build="Release"
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "$0: error: $o: unknown option" >&2
            echo $usage >&2
            exit 2
    esac
done

if [ $# -gt 0 ]; then
    host="$1"
else
    host=`uname -n`
fi

rm -rf CMakeCache.txt CMakeFiles

options="-Wdev --debug-trycompile"

# useful build types: Debug, Release, RelWithDebInfo
common_flags="\
        -D CMAKE_BUILD_TYPE:STRING=$build \
        -D DHSVM_SNOW_ONLY:BOOL=ON \
        -D DHSVM_BUILD_TESTS:BOOL=OFF \
        -D DHSVM_USE_RBM:BOOL=ON \
        -D DHSVM_DUMP_TOPO:BOOL=ON \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
"

if [ $host == "flophouse" ]; then

    cmake $options \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        $common_flags \
        ..

elif [ $host == "WE32673" ]; then

    # this is a Mac system with NetCDF installed using MacPorts
    # using the GNU compiler installed via MacPorts

    CC=gcc-mp-6
    CXX=g++-mp-6
    export CC CXX

    cmake $options \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -D NETCDF_DIR:PATH=/opt/local/include \
        -D DHSVM_USE_X11:BOOL=ON \
        -D DHSVM_USE_NETCDF:BOOL=ON \
        -D DHSVM_USE_RBM:BOOL=OFF \
        $common_flags \
        ..

elif [ $host == "WE32673-clang" ]; then

    # this is a Mac system with NetCDF installed using MacPorts
    # using the system (XCode) compiler

    CC=/usr/bin/clang
    CXX=/usr/bin/clang++
    export CC CXX

    cmake $options \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -D NETCDF_DIR:PATH=/opt/local/include \
        -D DHSVM_USE_X11:BOOL=ON \
        -D DHSVM_USE_NETCDF:BOOL=ON \
        -D DHSVM_USE_RBM:BOOL=OFF \
        $common_flags \
        ..

elif [ $host == "pe10900" ]; then
    
    # this is an older Mac system with Intel compilers and NetCDF
    # installed via MacPorts. This is how you use non-default compilers. 
    CC=icc
    FC=ifort
    export CC FC
    cmake $options \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -D DHSVM_USE_X11:BOOL=ON \
        -D DHSVM_USE_NETCDF:BOOL=ON \
        -D NETCDF_DIR:PATH=/opt/local/include \
        -D DHSVM_USE_RBM:BOOL=ON \
        $common_flags \
        ..

else

    # For an unknown system, turn most options off
    cmake $options \
        -D CMAKE_BUILD_TYPE:STRING=$build \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=FALSE \
        -D DHSVM_SNOW_ONLY:BOOL=OFF \
        -D DHSVM_USE_X11:BOOL=OFF \
        -D DHSVM_USE_NETCDF:BOOL=OFF \
        -D DHSVM_USE_RBM:BOOL=OFF \
        -D DHSVM_BUILD_TESTS:BOOL=OFF \
        ..

    echo "Unknown host: $host"
    exit 2
    
fi
