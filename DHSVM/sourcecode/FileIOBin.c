/*
 * SUMMARY:      FileIOBin.c - Functions for binary IO
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Functions for binary IO
 * DESCRIP-END.
 * FUNCTIONS:    CreateMapFileBin()
 *               Read2DMatrixBin()
 *               Read2DMatrixByteSwapBin()
 *               Write2DMatrixBin()
 *		 Write2DMatrixByteSwapBin()
 *               SizeOfNumberType()
 *               byte_swap_long()
 *               byte_swap_short()
 * COMMENTS:
 * $Id: FileIOBin.c,v 1.4 2003/07/01 21:26:14 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "fifobin.h"
#include "fileio.h"
#include "sizeofnt.h"
#include "settings.h"
#include "DHSVMerror.h"

/*****************************************************************************
  Function name: CreateMapFileBin()

  Purpose      : Open and close a new file.  If the file already exists it 
                 will be overwritten.

  Required     : 
    FileName  - Name of the new file

  Returns      : void

  Modifies     : 

  Comments     :
*****************************************************************************/
void CreateMapFileBin(char *FileName, ...)
{
  FILE *NewFile;

  OpenFile(&NewFile, FileName, "w", TRUE);
}

/*****************************************************************************
  Function name: Read2DMatrixBin()

  Purpose      : Function to read a 2D array from a file.

  Required     :
    FileName   - name of input file
    Matrix     - address of array data into
    NumberType - code for number type (taken from HDF, see comments at the
                 beginning of InitFileIO.c for more detail)
    NY         - Number of rows
    NX         - Number of columns
    NDataSet   - number of the dataset to read, i.e. the first matrix in a 
                 file is number 0, etc.
    Any remaining arguments are not used in straight binary

  Returns      : Number of elements read

  Modifies     : Matrix

  Comments     :
*****************************************************************************/
int Read2DMatrixBin(char *FileName, void *Matrix, int NumberType, int NY,
		    int NX, int NDataSet, ...)
{
  FILE *InFile;
  int NElements = 0;		/* number of elements read */
  size_t ElemSize;
  unsigned long OffSet;		/* number of bytes to OffSet (is non-zero when
				   reading matrices other than the first one
				   in the file */

  OpenFile(&InFile, FileName, "rb", FALSE);

  ElemSize = SizeOfNumberType(NumberType);

  OffSet = NY * NX * ElemSize * NDataSet;

  if (fseek(InFile, OffSet, SEEK_SET))
    ReportError(FileName, 39);
  NElements = fread(Matrix, ElemSize, NY * NX, InFile);
  if (NElements != NY * NX)
    ReportError(FileName, 2);

  fclose(InFile);

  return NElements;
}

/******************************************************************************/
int Read2DMatrixByteSwapBin(char *FileName, void *Matrix, int NumberType,
			    int NY, int NX, int NDataSet, ...)
{
  FILE *InFile;
  int NElements = 0;		/* number of elements read */
  size_t ElemSize;
  unsigned long OffSet;		/* number of bytes to OffSet (is non-zero when
				   reading matrices other than the first one
				   in the file */

  OpenFile(&InFile, FileName, "rb", FALSE);

  ElemSize = SizeOfNumberType(NumberType);

  OffSet = NY * NX * ElemSize * NDataSet;

  if (fseek(InFile, OffSet, SEEK_SET)) {
    ReportError(FileName, 39);
  }
  NElements = fread(Matrix, ElemSize, NY * NX, InFile);
  if (NElements != NY * NX) {
    ReportError(FileName, 2);
  }
  fclose(InFile);

  if (ElemSize == 4) {
    byte_swap_long(Matrix, NElements);
  }
  else if (ElemSize == 2) {
    byte_swap_short(Matrix, NElements);
  }
  else if (ElemSize != 1) {
    ReportError(FileName, 61);
  }

  return NElements;
}

/*****************************************************************************
  Function name: Write2DMatrixBin()

  Purpose      : Function to write a 2D array to a file.  Data is appended to
                 the end of the file.

  Required     :   
    FileName   - name of output file
    Matrix     - address of array containing matrix elements
    NumberType - code for number type (see comments at the
                 beginning of InitFileIO.c for more detail)
    NY         - Number of rows
    NX         - Number of columns

  Returns      : Number of elements written 

  Modifies     :

  Comments     :
*****************************************************************************/
int Write2DMatrixBin(char *FileName, void *Matrix, int NumberType, int NY,
		     int NX, ...)
{
  FILE *OutFile;		/* output file */
  size_t ElemSize = 0;		/* size of number type in bytes */

  OpenFile(&OutFile, FileName, "ab", FALSE);
  ElemSize = SizeOfNumberType(NumberType);

  if (!(fwrite(Matrix, ElemSize, NY * NX, OutFile)))
    ReportError(FileName, 41);

  fclose(OutFile);

  return NY * NX;
}

/******************************************************************************/
int Write2DMatrixByteSwapBin(char *FileName, void *Matrix, int NumberType,
			     int NY, int NX, ...)
{
  FILE *OutFile;		/* output file */
  size_t ElemSize = 0;		/* size of number type in bytes */
  int NElements;
  NElements = NX * NY;

  OpenFile(&OutFile, FileName, "ab", FALSE);
  ElemSize = SizeOfNumberType(NumberType);

  if (ElemSize == 4) {
    byte_swap_long(Matrix, NElements);
  }
  else if (ElemSize == 2) {
    byte_swap_short(Matrix, NElements);
  }
  else if (ElemSize != 1) {
    ReportError(FileName, 61);
  }

  if (!(fwrite(Matrix, ElemSize, NY * NX, OutFile))) {
    ReportError(FileName, 41);
  }

  fclose(OutFile);

  return NY * NX;
}

/******************************************************************************/
void byte_swap_short(short *buffer, int number_of_swaps)
{
  short *temp;
  int swap_loop;

  for (swap_loop = 0, temp = buffer; swap_loop < number_of_swaps;
       swap_loop++, temp++) {
    *temp = ((*temp & 0x00ff) << 8) | ((*temp & 0xff00) >> 8);
  }
}

/******************************************************************************/
void byte_swap_long(long *buffer, int number_of_swaps)
{
  long *temp;
  int swap_loop;

  for (swap_loop = 0, temp = buffer; swap_loop < number_of_swaps;
       swap_loop++, temp++) {
    *temp = ((*temp & 0x000000ff) << 24) | ((*temp & 0x0000ff00) << 8) |
      ((*temp & 0x00ff0000) >> 8) | ((*temp & 0xff000000) >> 24);
  }
}
