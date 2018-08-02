/*
 * SUMMARY:      fileio.h - header file for file I/O
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for file I/O
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: fileio.h,v3.1.2 2013/07/01 ning Exp $     
 */

#ifndef FILEIO_H
#define FILEIO_H

#include "data.h"

/* define identifiers for different file formats */

#define BIN 1			/* binary IO */
#define NETCDF 2		/* NetCDF format */
#define BYTESWAP 3		/* binary IO but byteswap reads */
void InitFileIO(int FileFormat);

/* global file extension string */
extern char fileext[];

/* function pointers for 2D file IO */

void CreateMapFile(char *FileName, char *FileLabel, MAPSIZE *Map);

int Read2DMatrix(char *FileName, void *Matrix, int NumberType, 
                 MAPSIZE *Map, int NDataSet, char *VarName, int index);

int Read2DMatrixAll(char *FileName, void *Matrix, int NumberType, 
                    MAPSIZE *Map, int NDataSet, char *VarName, int index);

int Write2DMatrix(char *FileName, void *Matrix, int NumberType, 
                  MAPSIZE *Map, MAPDUMP *DMap, int index);


/* generic file functions */
void OpenFile(FILE **FilePtr, char *FileName, char *Mode,
	      unsigned char OverWrite);

#endif
