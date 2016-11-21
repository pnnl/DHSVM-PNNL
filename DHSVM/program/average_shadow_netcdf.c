/*
 * SUMMARY:      average shadow images from fine to coarse times
 * USAGE:        
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    March-2000
 * Last Change:  Feb-2013
 *               
 * $Id:          average_shadow_netcdf.c, v3.1  2013/2/5  Ning Exp $
 * COMMENTS:
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fifoNetCDF.h"
#include "sizeofNetCDF.h"
#include "data.h"

int CopyDouble(double *Value, char *Str, const int NValues);
unsigned char IsLeapYear(int Year);
int GetNumber(char *numberStr);
void ReverseUCHARMatrix(unsigned char *Array1, unsigned char *Array2, 
						int nrow, int ncol);
int DayOfYear(int Year, int Month, int Day);

int main(int argc, char **argv)
{
  char   infilename[255],outfilename[255];
  char   VarName[255];
  int nRows;                    /* Number of rows */
  int nCols;                    /* Number of columns */
  int nIn,nOut;
  unsigned char *input, *input2;
  unsigned char  *output;
  float *temp;               /* to trap overflows of unsigned char *output */
  int i,j,k, q;
  float dx;
  int flag;
  int compress;
  MAPSIZE Map;
  MAPDUMP DMap;
  int hour, month,jday;
  int y, x;
  int n;
  unsigned char ***ShadowMap;

  if(argc < 11) {
    printf("usage is: average_shadow:  \n");
    printf("inputfile, outputfile, # in, # out, nrows, ncols dx XOrig YOrig month\n");
    printf("the last 4 variable should all be entered as integers \n");
    exit(-1);
  }

  strcpy(infilename, argv[1]);    /* name of the nc_uchar input file - no header */
  strcpy(outfilename, argv[2]);   /* name of the nc_uchar output file        */
  nIn = GetNumber(argv[3]);       /* number of input images (typically 24 for hourly) */
  nOut = GetNumber(argv[4]);      /* number of output images (typically 8 for 3 hourly) */
  nRows = GetNumber(argv[5]);     /* number of rows in a single image */  
  nCols = GetNumber(argv[6]);     /* number of columns in a single image */
  dx = GetNumber(argv[7]);        /* cell size in a single image */
  /* extreme west coordinatee */
  if (!(CopyDouble(&Map.Xorig, argv[8], 1)))
	  exit (-1);;
  /* exterme north coordinate */
  if (!(CopyDouble(&Map.Yorig, argv[9], 1)))
	  exit (-1);;
  month = GetNumber(argv[10]);

  if(fmod((float)nIn,(float)nOut)!=0.0) {
    printf("Number of input images not wholly divisible by number of output images \n");
    exit(-1);
  }

  /* Assign values to map variables */
  strcpy(VarName, "Shade.Factor"); 
  Map.X = 0;
  Map.Y = 0;
  Map.OffsetX = 0;
  Map.OffsetY = 0;
  Map.NX = nCols;
  Map.NY = nRows;
  Map.DX = dx;
  Map.DY = dx;
  Map.DXY = (float) sqrt(Map.DX * Map.DX + Map.DY * Map.DY);

  strcpy(DMap.FileName, outfilename);
  DMap.ID = 304;			
  DMap.Layer = 1;
  DMap.Resolution = MAP_OUTPUT;	/* Full resolution maps */
  strcpy(DMap.Name, "Shade.Factor");
  strcpy(DMap.LongName, "Shade Factor");
  strcpy(DMap.Format, "%d");
  strcpy(DMap.FileLabel, "Shade Factor");
  strcpy(DMap.Units, "");
  DMap.NumberType = NC_BYTE;         
  DMap.MaxVal = 0;
  DMap.MinVal = 0;

  compress=nIn/nOut;
  
  if (!(input = (unsigned char *) calloc(nRows*nCols, sizeof(unsigned char))))
    exit(-1);

  if (!(input2 = (unsigned char *) calloc(nRows*nCols, sizeof(unsigned char))))
    exit(-1);

  if (!(output = (unsigned char *) calloc(nRows*nCols, sizeof(unsigned char))))
    exit(-1);

  if (!(temp = (float *) calloc(nRows*nCols, sizeof(float))))
    exit(-1);

  CreateMapFileNetCDF(DMap.FileName, DMap.FileLabel, &Map);
  /* netcdf map properties */
  DMap.N = nOut;					
  if (!(DMap.DumpDate = (DATE *) calloc(DMap.N, sizeof(DATE))))
      exit(-1);

  for(i = 0; i < nOut; i++){  
	hour = 24 / nOut;
	DMap.DumpDate[i].Year = 2000;
    DMap.DumpDate[i].Month = month;
    DMap.DumpDate[i].Day = 15;
	jday = DayOfYear(DMap.DumpDate[i].Year, month, DMap.DumpDate[i].Day);
    DMap.DumpDate[i].JDay = jday;
    DMap.DumpDate[i].Hour = hour;

    for(k = 0; k < nRows * nCols; k++) 
      temp[k] = 0.0;
	
	for(j = 0; j < compress; j++){
		q = i*compress + j;
		// fread(input, sizeof(unsigned char), nCols*nRows, infile);
		flag = Read2DMatrixNetCDF(infilename, input, NC_BYTE, Map.NY, Map.NX, q, VarName, q); 
		if (flag == 0) {
			for(k = 0; k < nCols * nRows; k++) {
				temp[k] += (float)input[k]/(float)compress;
			}
		}
		else if (flag == 1){
			ReverseUCHARMatrix(input, input2, nRows, nCols);
			for(k = 0; k < nCols * nRows; k++) {
				temp[k] += (float)input2[k]/(float)compress;
			}
		}
		else 
			exit (-1);
	}
	for(k = 0;k < nCols*nRows; k++){
		output[k] = (unsigned char)temp[k];   
		if(temp[k] > 255.0) 
			output[k] = 255;
	}
	Write2DMatrixNetCDF(DMap.FileName, (void *)output, DMap.NumberType, Map.NY, Map.NX, &DMap, i);
  }
  return EXIT_SUCCESS;
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
    exit(-1);

  return number;
}
/*****************************************************************************
  CopyDouble()
*****************************************************************************/
int CopyDouble(double *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = strtod(Str, &EndPtr);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}
/*****************************************************************************
  Function: ReverseFloatMatrix()
  Usage: Revert the matrix Array2 by 180 degrees
  Author: Ning, 1/31/2013
******************************************************************************/
void ReverseUCHARMatrix(unsigned char *Array1, unsigned char *Array2, 
						int nrow, int ncol)
{
	int i;
	int j;

	for (i = 1; i <= nrow; i++) {
		for (j = 1; j <= ncol; j ++) {
			Array2[ncol*(nrow-i) + j - 1] = 
				Array1[(i-1)*ncol + j - 1] ;
		}
	}
	free(Array1);
}
/*****************************************************************************
  DayOfYear()
*****************************************************************************/
int DayOfYear(int Year, int Month, int Day) 
{
  int i;
  int Jday;
  int DaysPerMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if (IsLeapYear(Year))
    DaysPerMonth[1] = 29;
  else
    DaysPerMonth[1] = 28;

  for (i = 0, Jday = 0; i < (Month - 1); i++)
    Jday += DaysPerMonth[i];

  Jday += Day;

  return Jday;
}
/*****************************************************************************
  IsLeapYear()
*****************************************************************************/
unsigned char IsLeapYear(int Year) 
{
  if ((Year % 4 == 0 && Year % 100 != 0) || Year % 400 == 0) 
    return TRUE;
  return FALSE;
}

