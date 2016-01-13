/*
 * SUMMARY:      FileIOBin.c - Functions for binary IO
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD: Fri Sep 27 14:30:07 1996 by  William A Perkins <perk@euterpe.muse.pnl.gov>
 * DESCRIPTION:  Functions for binary IO
 * DESCRIP-END.
 * FUNCTIONS:    CreateFileBin()
 *               MakeFileNameBin()
 *               Read2DMatrixBin()
 *               Write2DImageBin()
 *               Write2DMatrixBin()
 *               SizeOfNumberType()
 * COMMENTS:     
 */

#ifndef lint
static char vcid[] = "$Id: FileIOBin.c,v 1.5 1996/10/18 16:52:13 nijssen Exp $";
#endif /* lint */

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
  Function name: CreateFileBin()

  Purpose      : Open and close a new file.  If the file already exists it 
                 will be overwritten.

  Required     : 
    FileName  - Name of the new file
    FileLabel - String describing file contents (currently not used for
                binary IO

  Returns      : void

  Modifies     : 

  Comments     :
*****************************************************************************/
void CreateFileBin(char *FileName, char *FileLabel)
{
  const char *Routine = "CreateFileBin";
  FILE *NewFile;

  OpenFile(&NewFile, FileName, "w", TRUE);
}

/*****************************************************************************
  Function name: MakeFileNameBin()

  Purpose      : Create a new file name ending in extension ".bin" to
                 indicate that the file is binary

  Required     : 
    Path     - directory for files
    Str1     - first part of filename
    Str2     - second part of filename
    FileName - filename

  Returns      : void

  Modifies     : FileName

  Comments     : 
*****************************************************************************/
void MakeFileNameBin(char *Path, char *Str1, char *Str2, char *FileName)
{
  const char *Routine = "MakeFileNameBin";
  
  MakeFileNameGen(Path, Str1, Str2, ".bin", FileName);
}


/*****************************************************************************
  Function name: Write2DMatrixBin()

  Purpose      : Function to write a 2D array to a file.  Data is appended to
                 the end of the file.

  Required     :   
    NY         - Number of rows
    NX         - Number of columns
    NumberType - code for number type (taken from HDF, see comments at the
                 beginning of InitFileIO.c for more detail)
    DataLabel  - string describing the data (currently not used for binary 
                 data)
    Units      - string describing the units of the data (currently not 
                 used for binary data) 
    Matrix     - address of array containing matrix elements
    FileName   - name of output file

  Returns      : Number of elements written 

  Modifies     :

  Comments     :
*****************************************************************************/
int Write2DMatrixBin(int NY, int NX, int NumberType, char *DataLabel, 
                     char *Units, void *Matrix, char *FileName)
{
  const char *Routine = "Write2DMatrixBin";
  FILE *OutFile;                /* output file */
  size_t ElemSize = 0;          /* size of number type in bytes */

  OpenFile(&OutFile, FileName, "ab", FALSE);
  ElemSize = SizeOfNumberType(NumberType);

  if (!(fwrite(Matrix, ElemSize, NY*NX, OutFile)))
    ReportError(FileName, 41);
  
  fclose(OutFile);
  
  return NY*NX;
}

