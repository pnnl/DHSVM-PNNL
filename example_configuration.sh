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
# Last Change: 2017-02-16 07:11:28 d3g096

set -xue

# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r] [name]"

set -- `getopt d $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi

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
        -D DHSVM_SNOW_ONLY:BOOL=OFF \
        -D DHSVM_BUILD_TESTS:BOOL=ON \
"

if [ $host == "flophouse" ]; then

    prefix="/net/flophouse/files0/perksoft/linux64"
    cmake $options \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec" \
        -D GA_DIR:STRING="$prefix" \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        $common_flags \
        ..

elif [ $host == "flophouse48" ]; then

    prefix="/net/flophouse/files0/perksoft/linux64/openmpi48"
    CC="$prefix/bin/gcc"
    CFLAGS="-pthread -Wall"
    export CC CFLAGS
    
    PATH="${prefix}/bin:${PATH}"
    export PATH
    cmake $options \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec" \
        -D GA_DIR:STRING="$prefix/ga-5-4" \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        $common_flags \
        ..

elif [ $host == "WE32673" ]; then

    # this is a Mac system with MPI and NetCDF installed using
    # MacPorts.  You cannot use the Apple CLang because Global Arrays
    # does not work with it.  It's better to use a GNU compiler.
    # Here, gcc-6 and mpich, installed via MacPorts, are used.

    prefix="/opt/local"
    CC="$prefix/bin/clang-mp-3.8"
    export CC

    cmake $options \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc-openmpi-clang38" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec-openmpi-clang38" \
        -D NETCDF_DIR:PATH="$prefix/include" \
        -D DHSVM_USE_X11:BOOL=OFF \
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


    echo "Unknown host: $host"
    exit 2
    
fi
