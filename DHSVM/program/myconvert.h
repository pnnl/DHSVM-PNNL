/*
 * SUMMARY:      convert.h - header for convert.c
 * USAGE:        
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:    18-Oct-1996 at 11:33:25
 * LAST-MOD:     18-Oct-1996 at 17:20:23 by Bart Nijssen
 * DESCRIPTION:  
 * DESCRIP-END.
 * COMMENTS:     
 */

/* $Id: convert.h,v 1.1 1996/10/19 21:37:13 nijssen Exp $ */

#ifndef _CONVERT_H
#define _CONVERT_H

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define N_FORMATS 10

enum ValidFormat {character = 0, ucharacter, shortint, ushortint, integer,
                  uinteger, longint, ulongint, floatp, doublep, ascii};

void Cast(int nCols, int readFormat, void *readArray, 
          int writeFormat, void *writeArray);

void Convert(int nRows, int nCols, int readFormat, FILE *inFile, 
             int writeFormat, FILE *outFile);

int GetFormat(char *formatStr);

int GetNumber(char *numberStr);

void InitSize(void);

void OpenFile(FILE **FilePtr, char *FileName, char *Mode, 
	      unsigned char OverWrite);

int ReadAscii(FILE *inFile, int nCols, int format, void *array);

int ReadBin(FILE *inFile, int nCols, int format, void *array);

void ReportError(char *errorStr1, char *errorStr2);

int WriteAscii(FILE *outFile, int nCols, int format, void *array);

int WriteBin(FILE *outFile, int nCols, int format, void *array);

#endif
