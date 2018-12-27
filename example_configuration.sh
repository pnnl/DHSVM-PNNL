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
# Last Change: 2018-12-12 11:01:36 d3g096

set -xue

# -------------------------------------------------------------
# handle command line options
# -------------------------------------------------------------
usage="$0 [-d|-r|-t] [name]"

set -- `getopt drt $*`
if [ $? != 0 ]; then
    echo $usage >&2
    exit 2
fi

timing="OFF"
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
	-t)
	    timing="ON"
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
        -D DHSVM_USE_RBM:BOOL=OFF \
        -D DHSVM_DUMP_TOPO:BOOL=ON \
	-D DHSVM_USE_GPTL:BOOL=$timing \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
"

if [ $host == "flophouse" ]; then

    CC=/usr/bin/gcc
    CXX=/usr/bin/g++
    export CC

    # For GPTL (not working)
    # CFLAGS="-finstrument-functions"
    # export CFLAGS

    prefix="/net/flophouse/files0/perksoft/linux64"
    cmake3 $options \
        -D MPI_C_COMPILER:STRING="/usr/lib64/openmpi/bin/mpicc" \
        -D MPI_CXX_COMPILER:STRING="/usr/lib64/openmpi/bin/mpicxx" \
        -D MPIEXEC:STRING="/usr/lib64/openmpi/bin/mpiexec" \
        -D GA_DIR:PATH="${prefix}/ga-c++" \
	-D GA_EXTRA_LIBS:STRING="-lm" \
        -D DHSVM_TIMING_LEVEL:STRING="1" \
        -D GPTL_DIR:PATH="$prefix/gptl-v5.5.3-2-gbb58395" \
        -D DHSVM_USE_NETCDF:BOOL=ON \
	-D NetCDF_BIN_DIR:PATH="/usr/lib64/openmpi/bin" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/dhsvm" \
        $common_flags \
        ..

elif [ $host == "tlaloc" ]; then

    CC=/usr/bin/gcc
    CXX=/usr/bin/g++
    export CC

    prefix="/file0/perksoft"
    cmake $options \
        -D MPI_C_COMPILER:STRING="mpicc" \
        -D MPI_CXX_COMPILER:STRING="mpicxx" \
        -D MPIEXEC:STRING="mpiexec" \
        -D GA_DIR:STRING="$prefix/ga-c++" \
	-D GA_EXTRA_LIBS:STRING="-lm" \
        -D DHSVM_TIMING_LEVEL:STRING="1" \
        -D GPTL_DIR:PATH="$prefix/gptl-v5.5.3-2-gbb58395" \
        -D DHSVM_USE_NETCDF:BOOL=ON \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/dhsvm" \
        $common_flags \
        ..

elif [ $host == "WE32673" ]; then

    # this is a Mac system with MPI and NetCDF installed using
    # MacPorts.  You cannot use the Apple CLang because Global Arrays
    # does not work with it.
    
    prefix="/opt/local"
    CC="$prefix/bin/clang-mp-6.0"
    CXX="$prefix/bin/clang++-mp-6.0"
    export CC CXX

    cmake $options \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc-mpich-clang60" \
        -D MPIEXEC:STRING="$prefix/bin/mpiexec-mpich-clang60" \
        -D GA_DIR:PATH="$HOME/Projects/GridPACK" \
	-D NetCDF_BIN_DIR:PATH="/opt/local/bin" \
        -D HDF5_DIR:PATH="/opt/local" \
        -D HDF5_INCLUDE_DIRS:PATH="/opt/local/include" \
        -D HDF5_PREFER_PARALLEL:BOOL=TRUE \
        -D DHSVM_USE_X11:BOOL=OFF \
        -D DHSVM_USE_NETCDF:BOOL=ON \
        -D NetCDF_BIN_DIR:PATH=/opt/local/bin \
        -D CMAKE_INSTALL_PREFIX:PATH="$HOME/Projects/DHSVM" \
        $common_flags \
        ..

elif [ $host == "WE32673-gnu" ]; then

    # this is a Mac system with MPI and NetCDF installed using
    # MacPorts.  

    prefix="/opt/local"
    CC="$prefix/bin/gcc-mp-6"
    CXX="$prefix/bin/g++-mp-6"
    export CC CXX

    cmake $options \
        -D MPI_C_COMPILER:STRING="$prefix/bin/mpicc-openmpi-gcc6" \
        -D MPIEXEC:STRING="$prefix/bin/mpicxx-openmpi-gcc6" \
        -D NETCDF_DIR:PATH="$prefix/include" \
        -D DHSVM_USE_X11:BOOL=OFF \
        -D DHSVM_USE_NETCDF:BOOL=ON \
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
        -D DHSVM_USE_X11:BOOL=ON \
        -D DHSVM_USE_NETCDF:BOOL=ON \
        -D NETCDF_DIR:PATH=/opt/local/include \
        $common_flags \
        ..

elif [ $host = "briareus" ]; then

    # with these modules (default compilers:

    # module load intel/15.0.1
    # module load intelmpi/2017.4.056

    prefix=/files0/dhsvm

    CC="icc"
    export CC

    cmake $options \
        -D DHSVM_USE_NETCDF:BOOL=ON \
        -D MPI_C_COMPILER:STRING="mpicc" \
        -D GA_DIR:STRING="$prefix" \
        -D MPI_CXX_COMPILER:STRING="mpicxx" \
        -D GA_EXTRA_LIBS:STRING="-libverbs -lm" \
        -D DHSVM_TIMING_LEVEL:STRING="1" \
        -D GPTL_DIR:PATH="$prefix" \
	-D NetCDF_DIR:PATH="$prefix" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix" \
        $common_flags \
        ..

elif [ $host = "briareus-gnu" ]; then

    # with these modules (default compilers:

    # module load gcc
    # module load mpi/openmpi/1.4.3/gnu

    prefix=/files0/dhsvm

    CC="/share/apps/gcc/4.5.0/bin/gcc"
    export CC

    cmake \
    -D DHSVM_USE_NETCDF:BOOL=OFF \
    -D MPI_C_COMPILER:STRING="/share/apps/openmpi/1.4.3/gnu/bin/mpicc" \
    -D GA_DIR:STRING="/files0/dhsvm" \
    -D CMAKE_INSTALL_PREFIX:PATH="/files0/dhsvm" \
    ..

elif [ $host = "constance" ]; then

    # with these modules (default compilers:

    # module load precision/i4
    # module load intel/15.0.1
    # module load intelmpi/2017.1.132
    # module load netcdf/4.3.2
    # module load cmake/2.8.12
    
    # GA installed here:

    prefix=/pic/projects/informed_hydro
    PATH="$prefix/netcdf-intel:$PATH"
    CC=icc
    CXX=icpc
    export CC CXX PATH

    cmake $options \
        -D MPI_C_COMPILER:STRING="mpicc" \
        -D MPIEXEC:STRING="mpiexec" \
        -D GA_DIR:STRING="$prefix/dhsvm-intel" \
	-D GA_EXTRA_LIBS:STRING="-libverbs -lm" \
        -D GPTL_DIR:PATH="$prefix" \
        -D DHSVM_USE_X11:BOOL=OFF \
        -D DHSVM_USE_NETCDF:BOOL=ON \
	-D NetCDF_DIR:PATH="$prefix/netcdf-intel" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/dhsvm-intel" \
        $common_flags \
        ..

elif [ $host = "constance-gnu" ]; then

    # with these modules (default compilers:

    # module load precision/i4
    # module load gcc/4.8.2
    # module load openmpi/1.8.3
    # module load netcdf/4.3.2
    
    # GA installed here:

    prefix=/pic/projects/informed_hydro
    PATH="$prefix/netcdf-gnu:$PATH"
    CC=/share/apps/gcc/4.8.2/bin/gcc
    CXX=/share/apps/gcc/4.8.2/bin/g++
    export CC CXX PATH

    cmake $options \
        -D MPI_C_COMPILER:STRING="/share/apps/openmpi/1.8.3/gcc/4.8.2/bin/mpicc" \
        -D MPIEXEC:STRING="/share/apps/openmpi/1.8.3/gcc/4.8.2/bin/mpiexec" \
        -D GA_DIR:STRING="$prefix/dhsvm-gnu" \
	-D GA_EXTRA_LIBS:STRING="-libverbs -lm -lpthread" \
        -D GPTL_DIR:PATH="$prefix/dhsvm-gnu" \
        -D DHSVM_USE_X11:BOOL=OFF \
        -D DHSVM_USE_NETCDF:BOOL=ON \
	-D NetCDF_DIR:PATH="$prefix/netcdf-gnu" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/dhsvm-gnu" \
        $common_flags \
        ..

elif [ $host = "constance-gnu-pr" ]; then

    # with these modules (default compilers:

    # module load precision/i4
    # module load gcc/4.8.2
    # module load openmpi/1.8.3
    # module load netcdf/4.3.2
    
    # GA installed here:

    prefix=/pic/projects/informed_hydro
    PATH="$prefix/netcdf-gnu:$PATH"
    CC=/share/apps/gcc/4.8.2/bin/gcc
    CXX=/share/apps/gcc/4.8.2/bin/g++
    export CC CXX PATH
    CFLAGS="-pthread"
    CXXFLAGS="-pthread"
    export CFLAGS CXXFLAGS

    cmake $options \
        -D MPI_C_COMPILER:STRING="/share/apps/openmpi/1.8.3/gcc/4.8.2/bin/mpicc" \
        -D MPIEXEC:STRING="/share/apps/openmpi/1.8.3/gcc/4.8.2/bin/mpiexec" \
        -D GA_DIR:STRING="$prefix/dhsvm-gnu-pr" \
	-D GA_EXTRA_LIBS:STRING="-lrt -lm" \
        -D GPTL_DIR:PATH="$prefix/dhsvm-gnu-pr" \
        -D DHSVM_USE_X11:BOOL=OFF \
        -D DHSVM_USE_NETCDF:BOOL=ON \
	-D USE_PROGRESS_RANKS:BOOL=ON \
	-D NetCDF_DIR:PATH="$prefix/netcdf-gnu" \
        -D CMAKE_INSTALL_PREFIX:PATH="$prefix/dhsvm-gnu-pr" \
        $common_flags \
        ..

else

    # For an unknown system, turn most options off
    cmake $options \
        -D DHSVM_SNOW_ONLY:BOOL=OFF \
        -D DHSVM_USE_X11:BOOL=OFF \
        -D DHSVM_USE_NETCDF:BOOL=OFF \
        -D DHSVM_BUILD_TESTS:BOOL=OFF \
        ..

    echo "Unknown host: $host"
    exit 2
    
fi
