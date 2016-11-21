/*
 * SUMMARY:      Files.c - File functions
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD:     29-Oct-1996 at 16:37:50 by Bart Nijssen
 * DESCRIPTION:  File and I/O functions
 * DESCRIP-END.
 * FUNCTIONS:    MakeFileNameGen() 
 *               OpenFile() 
 *               ScanInts() 
 *               ScanFloats() 
 *               SkipLines()             
 *               SkipHeader() 
 *               ScanDoubles() 
 *               ScanUChars() 
 * COMMENTS:     
 */

#ifndef lint
static char vcid[] = "nil";
#endif /* lint */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "init.h"
#include "fileio.h"

/*****************************************************************************
  Function name: MakeFileNameGen()

  Purpose      : Generic file name generator, called by format specific file 
                 name generating functions

  Required     :
    Path     - directory where file will be created
    Str1     - first part of file name
    Str2     - second part of file name
    Str3     - third part of file name
    FileName - new file name

  Returns      : void

  Modifies     : FileName

  Comments     : Str3 is intended to be used as an extension indicating the 
                 file format , e.g. ".bin" or ".hdf"
*****************************************************************************/
void MakeFileNameGen(char *Path, char *Str1, char *Str2, char *Str3,
                     char *FileName)
{
  char TempStr[MAXSTRING+1];
  const char *Routine = "MakeFileNameGen";

  sprintf(TempStr, "%s%s%s%s", Path, Str1, Str2, Str3);

  if (strlen(TempStr) > NAMESIZE) 
    TempStr[NAMESIZE] = '\0';

  InitCharArray(FileName, NAMESIZE+1);
  strcpy(FileName, TempStr);
}

/*****************************************************************************
  OpenFile()
*****************************************************************************/
void OpenFile(FILE **FilePtr, char *FileName, char *Mode, 
	      uchar OverWrite)
{
  struct stat FileInfo;

  /* if OverWrite is FALSE and the Mode is not read, check whether the file 
     already exists */ 

  if (!OverWrite && strstr(Mode, "w")) {
    if (!(stat(FileName, &FileInfo))) 
      ReportError(FileName, 4);
  }

  if (!(*FilePtr = fopen(FileName, Mode)))
    ReportError(FileName, 3);
}

