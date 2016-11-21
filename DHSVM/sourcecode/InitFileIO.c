/*
 * SUMMARY:      InitFileIO.c - Initialize the file IO functions
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize the file IO functions depending on the file
 *               format to be used (at this time either binary or HDF v3.3
 * DESCRIP-END.
 * FUNCTIONS:    InitFileIO()
 * COMMENTS:     In order to use the NetCDF, you have to define HAVE_NETCDF 
 *               during the build
 * $Id: InitFileIO.c,v 3.1 2013/02/06 19:12 ning Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fileio.h"
#include "fifobin.h"
#include "fifoNetCDF.h"
#include "DHSVMerror.h"

/*******************************************************************************
  Function name: InitFileIO()

  Purpose      : Initialize function pointers for file IO

  Required     :
    int FileFormat - identifier for the file format to be used

  Returns      : void

  Modifies     : function pointers for file IO

  Comments     :

  This function sets the file pointers for file I/O functions to the
  functions that implement the desired file format.  By using file pointers,
  the main routines do not need to be changed if a new file format is to be
  supported in the future (at least that is what I am trying to accomplish).
  The only thing that will need to be done is write the necessary I/O
  functions for the new file format, and add the additional options to this 
  initialization routine.

  For example: Currently three different file formats are supported, i.e.
  (plain vanilla) binary, swapped binary, and NetCDF v 3.4. (HDF v. 3.3 support
  has been discontinued as of Mon Jan 25 1999).  The user can specify which file
  format is to be used, and then the function pointer for the function
  Read2DMatrix will be set to Read2DMatrixBin, Read2DMatrixByteSwapBin, or
  Read2DMatrixNetCDF.  In the remaining part of the program this change is
  transparent, and the function call used to read a matrix is Read2DMatrix for
  either case.

  Information is stored in all files in the following way:
  fastest varying dimension: X (West to East)
  next fastest dimension   : Y (North to South)
          .                : Variable (if more than one)
  slowest varying dimension: Time (if more than one timestep)
  All the information is written out or read in for one timestep at a time.
            
  The following terminology is used (this is not based on anything
  "official", and may (will) not correspond to either dictionary or
  scientific definitions, but there is only so much one can do with variable
  and function names).  Some of the distinctions being made may not always
  seem to make sense for certain file formats, and partly result from the
  fact that the first version of the model used the HDF file format:
   - 2DMatrix :
       a map layer with X and Y dimension.  
   - 2DImage :
       a map layer with X and Y dimension in which the data are stored as 
       unsigned char, i.e. values in the interval [a, b] are mapped to numbers 
       in the range [0, 255].

   In order to use the NetCDF functions the NetCDF library needs to be 
   installed on your system.  If this is the case, HAVE_NETCDF needs to be
   defined at compile time.  If it is not defined the NetCDF functions cannot be
   used, and DHSVM will not try to access the NetCDF libraries.
*******************************************************************************/
void InitFileIO(int FileFormat)
{
  const char *Routine = "InitFileIO";

  printf("Initializing file IO\n");

  /************************* Binary format **********************/
  if (FileFormat == BIN) {
    strcpy(fileext, ".bin");
    CreateMapFile = CreateMapFileBin;
    Read2DMatrix = Read2DMatrixBin;
    Write2DMatrix = Write2DMatrixBin;
  }
  else if (FileFormat == BYTESWAP) {
    strcpy(fileext, ".bin");
    CreateMapFile = CreateMapFileBin;
    Read2DMatrix = Read2DMatrixByteSwapBin;
    Write2DMatrix = Write2DMatrixByteSwapBin;
  }
  /************* NetCDF File Format (version 3.4) ****************/
  else if (FileFormat == NETCDF) {
#ifdef HAVE_NETCDF
    strcpy(fileext, ".nc");
    CreateMapFile = CreateMapFileNetCDF;
    Read2DMatrix = Read2DMatrixNetCDF;
    Write2DMatrix = Write2DMatrixNetCDF;
#else
    ReportError((char *) Routine, 56);
#endif
  }
  else
    ReportError((char *) Routine, 38);
}
