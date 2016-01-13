/*
 * SUMMARY:      convert.c - convert one matrix into another
 * USAGE:        convert <from format> <to format> <infile> <nrows> <ncols>
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:    18-Oct-1996 at 10:56:31
 * LAST-MOD: Wed Jun 24 21:24:16 1998 by Bart Nijssen <nijssen@u.washington.edu>
 * DESCRIPTION:  Converts a matrix from one format to another, e.g.
 *                 short to float 
 *                 float to asc
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 */

#define _POSIX_SOURCE 1

#ifndef lint
static char vcid[] = "$Id: convert.c,v 1.4 1998/06/25 04:25:10 nijssen Exp $";
#endif /* lint */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "myconvert.h"

const char *usage = 
"myconvert source_format target_format source_file target_file\n        number_of_rows number_of_column\n";

int sizeArray[N_FORMATS];     /* Array with sizeof() */

int main(int argc, char **argv)
{
  FILE *inFile;                 /* Input file */
  FILE *outFile;                /* Output file */
  char inFilename[BUFSIZ+1];    /* Name of input file */
  char outFilename[BUFSIZ+1];   /* Name of output file */
  int readFormat;               /* Formatspecifier for input file */
  int nRows;                    /* Number of rows */
  int nCols;                    /* Number of columns */
  int writeFormat;              /* Format specifier for output file */

  InitSize();

  if (argc < 7) {
    fprintf(stderr, usage);
    exit(EXIT_FAILURE);
  }

  readFormat  = GetFormat(argv[1]);
  writeFormat = GetFormat(argv[2]);
  
  strcpy(inFilename, argv[3]);
  strcpy(outFilename, argv[4]);

  nRows = GetNumber(argv[5]);
  nCols = GetNumber(argv[6]);
  
  if (readFormat == ascii)
    OpenFile(&inFile, inFilename, "r", FALSE);
  else
    OpenFile(&inFile, inFilename, "rb", FALSE);

  if (writeFormat == ascii)
    OpenFile(&outFile, outFilename, "w", TRUE);
  else
    OpenFile(&outFile, outFilename, "wb", TRUE);
  
  Convert(nRows, nCols, readFormat, inFile, writeFormat, outFile); 
  
  fclose(inFile);
  fclose(outFile);
  
  return EXIT_SUCCESS;
}

/*****************************************************************************
  Cast()
*****************************************************************************/
void Cast(int nCols, int readFormat, void *readArray, 
          int writeFormat, void *writeArray)
{
  long *longValue = NULL;
  unsigned long *ulongValue = NULL;
  double *doubleValue = NULL;
  int promotion = -1;
  int i;
  
  if (readFormat == writeFormat)
    return;
  
  switch (readFormat) {
  case character:
  case shortint:
  case integer:
  case longint:
    promotion = longint;
    longValue = calloc(nCols, sizeArray[promotion]);
    if (longValue == NULL)
      ReportError("", "Failure to allocate memory");
    break;
  case ucharacter:
  case ushortint:
  case uinteger:
  case ulongint:
    promotion = ulongint;
    ulongValue = calloc(nCols, sizeArray[promotion]);
    if (ulongValue == NULL)
      ReportError("", "Failure to allocate memory");
    break;
  case floatp:
  case doublep:
    promotion = doublep;
    doubleValue = calloc(nCols, sizeArray[promotion]);
    if (doubleValue == NULL)
      ReportError("", "Failure to allocate memory");
    break;
  default:
    ReportError("Unrecognized format specifier:", "Fatal error");
    break;
  }

  switch (readFormat) {
  case character:
    for (i = 0; i < nCols; i++)
      longValue[i] = (long) ((char *)readArray)[i];
    break;  
  case shortint:
    for (i = 0; i < nCols; i++)
      longValue[i] = (long) ((short *)readArray)[i];
    break;  
  case integer:
    for (i = 0; i < nCols; i++)
      longValue[i] = (long) ((int *)readArray)[i];
    break;  
  case longint:
    for (i = 0; i < nCols; i++)
      longValue[i] = (long) ((long *)readArray)[i];
    break;  
  case ucharacter:
    for (i = 0; i < nCols; i++)
      ulongValue[i] = (unsigned long) ((unsigned char *)readArray)[i];
    break;  
  case ushortint:
    for (i = 0; i < nCols; i++)
      ulongValue[i] = (unsigned long) ((unsigned short *)readArray)[i];
    break;  
  case uinteger:
    for (i = 0; i < nCols; i++)
      ulongValue[i] = (unsigned long) ((unsigned int *)readArray)[i];
    break;  
  case ulongint:
    for (i = 0; i < nCols; i++)
      ulongValue[i] = (unsigned long) ((unsigned long *)readArray)[i];
    break;  
  case floatp:
    for (i = 0; i < nCols; i++)
      doubleValue[i] = (double) ((float *)readArray)[i];
    break;  
  case doublep:
    for (i = 0; i < nCols; i++)
      doubleValue[i] = (double) ((double *)readArray)[i];
    break;  
  default:
    ReportError("Unrecognized format specifier:", "Fatal error");
    break;
  }

  switch (promotion) {
  case longint:
    switch (writeFormat) {
    case character:
      for (i = 0; i < nCols; i++)
        ((char *)writeArray)[i] = (char)longValue[i];
      break;
    case ucharacter:
      for (i = 0; i < nCols; i++)
        ((unsigned char *)writeArray)[i] = (unsigned char)longValue[i];
      break;
    case shortint:
      for (i = 0; i < nCols; i++)
        ((short *)writeArray)[i] = (short)longValue[i];
      break;
    case ushortint:
      for (i = 0; i < nCols; i++)
        ((unsigned short *)writeArray)[i] = (unsigned short)longValue[i];
      break;
    case integer:
      for (i = 0; i < nCols; i++)
        ((int *)writeArray)[i] = (int)longValue[i];
      break;
    case uinteger:
      for (i = 0; i < nCols; i++)
        ((unsigned int *)writeArray)[i] = (unsigned int)longValue[i];
      break;
    case longint:
      for (i = 0; i < nCols; i++)
        ((long *)writeArray)[i] = (long)longValue[i];
      break;
    case ulongint:
      for (i = 0; i < nCols; i++)
        ((unsigned long *)writeArray)[i] = (unsigned long)longValue[i];
      break;
    case floatp:
      for (i = 0; i < nCols; i++)
        ((float *)writeArray)[i] = (float)longValue[i];
      break;
    case doublep:
      for (i = 0; i < nCols; i++)
        ((double *)writeArray)[i] = (double)longValue[i];
      break;
    default:
      ReportError("Unrecognized format specifier:", "Fatal error");
      break;
    }
    break;
  case ulongint:
    switch(writeFormat) {
    case character:
      for (i = 0; i < nCols; i++)
        ((char *)writeArray)[i] = (char)ulongValue[i];
      break;
    case ucharacter:
      for (i = 0; i < nCols; i++)
        ((unsigned char *)writeArray)[i] = (unsigned char)ulongValue[i];
      break;
    case shortint:
      for (i = 0; i < nCols; i++)
        ((short *)writeArray)[i] = (short)ulongValue[i];
      break;
    case ushortint:
      for (i = 0; i < nCols; i++)
        ((unsigned short *)writeArray)[i] = (unsigned short)ulongValue[i];
      break;
    case integer:
      for (i = 0; i < nCols; i++)
        ((int *)writeArray)[i] = (int)ulongValue[i];
      break;
    case uinteger:
      for (i = 0; i < nCols; i++)
        ((unsigned int *)writeArray)[i] = (unsigned int)ulongValue[i];
      break;
    case longint:
      for (i = 0; i < nCols; i++)
        ((long *)writeArray)[i] = (long)ulongValue[i];
      break;
    case ulongint:
      for (i = 0; i < nCols; i++)
        ((unsigned long *)writeArray)[i] = (unsigned long)ulongValue[i];
      break;
    case floatp:
      for (i = 0; i < nCols; i++)
        ((float *)writeArray)[i] = (float)ulongValue[i];
      break;
    case doublep:
      for (i = 0; i < nCols; i++)
        ((double *)writeArray)[i] = (double)ulongValue[i];
      break;
    default:
      ReportError("Unrecognized format specifier:", "Fatal error");
      break;
    }
    break;
  case doublep:
    switch (writeFormat) {
    case character:
      for (i = 0; i < nCols; i++)
        ((char *)writeArray)[i] = (char)doubleValue[i];
      break;
    case ucharacter:
      for (i = 0; i < nCols; i++)
        ((unsigned char *)writeArray)[i] = (unsigned char)doubleValue[i];
      break;
    case shortint:
      for (i = 0; i < nCols; i++)
        ((short *)writeArray)[i] = (short)doubleValue[i];
      break;
    case ushortint:
      for (i = 0; i < nCols; i++)
        ((unsigned short *)writeArray)[i] = (unsigned short)doubleValue[i];
      break;
    case integer:
      for (i = 0; i < nCols; i++)
        ((int *)writeArray)[i] = (int)doubleValue[i];
      break;
    case uinteger:
      for (i = 0; i < nCols; i++)
        ((unsigned int *)writeArray)[i] = (unsigned int)doubleValue[i];
      break;
    case longint:
      for (i = 0; i < nCols; i++)
        ((long *)writeArray)[i] = (long)doubleValue[i];
      break;
    case ulongint:
      for (i = 0; i < nCols; i++)
        ((unsigned long *)writeArray)[i] = (unsigned long)doubleValue[i];
      break;
    case floatp:
      for (i = 0; i < nCols; i++)
        ((float *)writeArray)[i] = (float)doubleValue[i];
      break;
    case doublep:
      for (i = 0; i < nCols; i++)
        ((double *)writeArray)[i] = (double)doubleValue[i];
      break;
    default:
      ReportError("Unrecognized format specifier:", "Fatal error");
      break;
    }
    break;
  default:
    ReportError("Unrecognized promotion specifier:", "Fatal error");
    break;
  }

  switch (promotion) {
  case longint: 
    free(longValue);
    break;
  case ulongint:
    free(ulongValue);
    break;
  case doublep:
    free(doubleValue);
    break;
  default:
    ReportError("Unrecognized promotion specifier:", "Fatal error");
    break;
  }
}

/*****************************************************************************
  Convert()
*****************************************************************************/
void Convert(int nRows, int nCols, int readFormat, FILE *inFile, 
             int writeFormat, FILE *outFile)
{
  char errorStr[BUFSIZ+1];
  int nElements;
  int row;
  int (*Read)(FILE *inFile, int nCols, int format, void *Array);
  int (*Write)(FILE *outFile, int nCols, int format, void *Array);
  void *readArray;
  void *writeArray;

  if (readFormat == ascii) {
    Read = ReadAscii;
    readFormat = writeFormat;
  }
  else 
    Read = ReadBin;

  if (writeFormat == ascii) {
    Write = WriteAscii;
    writeFormat = readFormat;
  }
  else 
    Write = WriteBin;
  
  readArray = calloc(nCols, sizeArray[readFormat]);
  if (readArray == NULL)
    ReportError("", "Failure to allocate memory");
  
  if (readFormat != writeFormat) {
    writeArray = calloc(nCols, sizeArray[writeFormat]);
    if (writeArray == NULL)
      ReportError("", "Failure to allocate memory");
  }
  else
    writeArray = readArray;

  for (row = 0; row < nRows; row++) {
    nElements = Read(inFile, nCols, readFormat, readArray);
    if (nElements != nCols) {
      sprintf(errorStr, "Row: %d\tColumn: %d", row, nElements);
      ReportError(errorStr, "Error reading input:");
    }
    Cast(nCols, readFormat, readArray, writeFormat, writeArray);
    nElements = Write(outFile, nCols, writeFormat, writeArray);
    if (nElements != nCols) {
      sprintf(errorStr, "Row: %d\tColumn: %d", row, nElements);
      ReportError(errorStr, "Error writing output:");
    }
  }
  free(readArray);
  if (readFormat != writeFormat)
    free(writeArray);
} 
  
/*****************************************************************************
  GetFormat()
*****************************************************************************/
int GetFormat(char *formatStr) 
{
  char *pStr;
  int format = -1;
  unsigned char notSigned = FALSE;

  for (pStr = formatStr; *pStr; pStr++)
    *pStr = tolower(*pStr);

  pStr = formatStr;

  if (*pStr == 'u') {
    notSigned = TRUE;
    pStr++;
  }

  switch (*pStr) {
  case 'c': 
    format = notSigned ? ucharacter : character;
    break;
  case 's':
    format = notSigned ? ushortint : shortint;
    break;
  case 'i':
    format = notSigned ? uinteger : integer;
    break;
  case 'l':
    format = notSigned ? ulongint : longint;
    break;
  case 'f':
    if (notSigned) 
      ReportError(formatStr, "Unrecognized format specifier:");
    else
      format = floatp;
    break;
  case 'd':
    if (notSigned) 
      ReportError(formatStr, "Unrecognized format specifier:");
    else
      format = doublep;
    break;
  case 'a':
    if (notSigned) 
      ReportError(formatStr, "Unrecognized format specifier:");
    else
      format = ascii;
    break;
  default:
    ReportError(formatStr, "Unrecognized format specifier:");
    break;
  }

  return format;
}

/*****************************************************************************
  GetNumber()
*****************************************************************************/
int GetNumber(char *numberStr) 
{
  char *endPtr;
  int number = 0;

  number = (int) strtol(numberStr, &endPtr, 0);
  if (*endPtr != '\0')
    ReportError(numberStr, "Invalid integer:");

  return number;
}

/*****************************************************************************
  InitSize()
*****************************************************************************/
void InitSize(void) 
{
  sizeArray[character]  = sizeof(char);
  sizeArray[ucharacter] = sizeof(unsigned char);
  sizeArray[shortint]   = sizeof(short);
  sizeArray[ushortint]  = sizeof(unsigned short);
  sizeArray[integer]    = sizeof(int);
  sizeArray[uinteger]  = sizeof(unsigned int);
  sizeArray[longint]    = sizeof(long);
  sizeArray[ulongint]   = sizeof(unsigned long);
  sizeArray[floatp]     = sizeof(float);
  sizeArray[doublep]    = sizeof(double);
}  
    
/*****************************************************************************
  OpenFile()
*****************************************************************************/
void OpenFile(FILE **FilePtr, char *FileName, char *Mode, 
	      unsigned char OverWrite)
{
  struct stat FileInfo;

  /* if OverWrite is FALSE and the Mode is not read, check whether the file 
     already exists */ 

  if (!OverWrite && strstr(Mode, "w")) {
    if (!(stat(FileName, &FileInfo))) 
      ReportError(FileName, 
                  "File already exists, and should not be overwritten:");
  }

  if (!(*FilePtr = fopen(FileName, Mode)))
    ReportError(FileName, "Cannot open file:");
}

/*****************************************************************************
  ReadAscii()
*****************************************************************************/
int ReadAscii(FILE *inFile, int nCols, int format, void *array) 
{
  int i = 0;
  int tempint;
  unsigned int tempuint;

  switch(format) {
  case character: 
    for (i = 0; i < nCols; i++) {
      fscanf(inFile, "%d", &tempint);
      ((char *)array)[i] = (char)tempint;
    }
    break;
  case ucharacter: 
    for (i = 0; i < nCols; i++) {
      fscanf(inFile, "%u", &tempuint);
      ((unsigned char *)array)[i] = (unsigned char)tempuint;
    }
    break;
  case shortint: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%hd", &((short *)array)[i]);
    break;
  case ushortint: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%hu", &((unsigned short *)array)[i]);
    break;
  case integer: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%d", &((int *)array)[i]);
    break;
  case uinteger: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%u",  &((unsigned int *)array)[i]);
    break;
  case longint: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%ld", &((long *)array)[i]);
    break;
  case ulongint: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%lu", &((unsigned long *)array)[i]);
    break;
  case floatp: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%f", &((float *)array)[i]);
    break;
  case doublep: 
    for (i = 0; i < nCols; i++) 
      fscanf(inFile, "%lf", &((double *)array)[i]);
    break;
  default:
    ReportError("Unknown format specifier in switch statement", 
                "Fatal error:");
    break;
  }
  return i;
}
              
/*****************************************************************************
  ReadBin()
*****************************************************************************/
int ReadBin(FILE *inFile, int nCols, int format, void *array) 
{
  int elementSize = sizeArray[format];

  return fread(array, elementSize, nCols, inFile);    
}

/*****************************************************************************
  ReportError()
*****************************************************************************/
void ReportError(char *errorStr1, char *errorStr2)
{
  fprintf(stderr, "%s %s\n", errorStr2, errorStr1);

  exit(EXIT_FAILURE);
}

/*****************************************************************************
  WriteAscii()
*****************************************************************************/
int WriteAscii(FILE *outFile, int nCols, int format, void *array) 
{
  int i = 0;
  
  switch(format) {
  case character: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%d\t",(int)((char *)array)[i]);
    break;
  case ucharacter: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%u\t",(unsigned int)((unsigned char *)array)[i]);
    break;
  case shortint: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%hd\t",((short *)array)[i]);
    break;
  case ushortint: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%hu\t",((unsigned short *)array)[i]);
    break;
  case integer: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%d\t",((int *)array)[i]);
    break;
  case uinteger: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%u\t", ((unsigned int *)array)[i]);
    break;
  case longint: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%ld\t",((long *)array)[i]);
    break;
  case ulongint: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%lu\t",((unsigned long *)array)[i]);
    break;
  case floatp: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%f\t", ((float *)array)[i]);
    break;
  case doublep: 
    for (i = 0; i < nCols; i++) 
      fprintf(outFile, "%f\t", ((double *)array)[i]);
    break;
  default:
    ReportError("Unknown format specifier in switch statement", 
                "Fatal error:");
    break;
  }

  fprintf(outFile, "\n");
  return i;
}

/*****************************************************************************
  WriteBin()
*****************************************************************************/
int WriteBin(FILE *outFile, int nCols, int format, void *array)
{
  int elementSize = sizeArray[format];

  return fwrite(array, elementSize, nCols, outFile);
}
