/*
 * SUMMARY:      calculate skyview factor from a dem
 * USAGE:        
 *
 * AUTHOR:       Ning Sun
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       ning@hydro.washington.edu
 * ORIG-DATE:    02/04/2013
 *
 * Last Change:  
 * Modify:       
 * $Id:          skyviewNetCDF.c, v 3.1.1  2013/2/5   Ning Exp $  
 * COMMENTS:
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fifoNetCDF.h"
#include <math.h>
#include "sizeofNetCDF.h"
#include "data.h"
#include "settings.h"

int GetNumber(char *numberStr);
int CopyDouble(double *Value, char *Str, const int NValues);

int main(int argc, char **argv)
{
  FILE  *demfile,*outfile;
  char   demfilename[255],outfilename[255];
  char   VarName[255];	
  int    nRows;                    /* Number of rows */
  int    nCols;                    /* Number of columns */
  float *temp;
  float **elev;
  float **skyview;
  int    i;
  int    ny,nx;
  float  lx,ly;
  double max_angle,angle;
  int    nLook;
  double theta;
  float  dx;
  float  x,y,sx,sy,dz,dist;
  float  mx,my;
  float  start_elev;
  float  max_elev;
  int    stop_flag;
  float *Array;
  char FileLabel[BUFSIZE + 1];
  int eflag = 0;
  int flag;
  MAPSIZE Map;
  MAPDUMP DMap;

  /****************************************************************************/
  /*                             INITIALIZATION                               */
  /****************************************************************************/
  /* Fill the Map structure */
  strcpy(Map.System, "Coordinate system");

  if(argc < 9) {
    printf("usage is: skyview:  \n");
    printf("demfilename, outfilename, # of look direction, nrows, ncols, cellsize, XOrigin, YOrigina\n");
    exit(-1);
  }

  strcpy(demfilename, argv[1]);			/* name of the binary float dem input file - no header */
  strcpy(outfilename, argv[2]);			/* name of the binary float skyview output file        */
  nLook = GetNumber(argv[3]);			/* number of directions to look in 8 or 16 are adequate */
  nRows = GetNumber(argv[4]);			/* number of rows in the input dem */
  nCols = GetNumber(argv[5]);			/* number of columns in the input dem */
  dx    = (float)GetNumber(argv[6]);	/* the cellsize of the dem (program assumes that */
										/* x and y are the same and that the units of dx */
										/* are the same units as in the dem)*/
  /* extreme west coordinatee */
  if (!(CopyDouble(&Map.Xorig, argv[7], 1)))
	  exit (-1);;
  /* exterme north coordinate */
  if (!(CopyDouble(&Map.Yorig, argv[8], 1)))
	  exit (-1);;

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
  DMap.ID = 305;				/* Snow water equivalent */
  DMap.Layer = 1;
  DMap.Resolution = MAP_OUTPUT;	/* Full resolution maps */
  DMap.N = 1;					/* Dump for a single timesteps */
  strcpy(DMap.Name, "SkyView.Factor");
  strcpy(DMap.LongName, "SkyView.Factor");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "-");
  strcpy(DMap.FileLabel, "SkyView Factor");
  DMap.NumberType = NC_FLOAT;          /* NC_FLOAT */
  DMap.DumpDate = (DATE *) calloc(DMap.N, sizeof(DATE));
  if (DMap.DumpDate == NULL)
    exit(-1);
  DMap.DumpDate[0].Year = 1999;
  DMap.DumpDate[0].Month = 12;
  DMap.DumpDate[0].Day = 31;
  DMap.DumpDate[0].JDay = 365;
  DMap.DumpDate[0].Hour = 23;
  DMap.MaxVal = 0;
  DMap.MinVal = 0;

  /* allocate memory */
  if (!((elev) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((elev)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }

  if (!((skyview) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((skyview)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }
  
  if (!(Array = (float*) calloc(nCols*nRows, sizeof(float))))
     exit(-1);

  /* read input file dem.nc */
  if (!(temp = (float *) calloc(nRows*nCols, sizeof(float))))
    exit(-1);
  strcpy(VarName, "Basin.DEM");
  flag = Read2DMatrixNetCDF(demfilename, temp, NC_FLOAT, Map.NY, Map.NX, 0,
	       VarName, 0);
  
  if (flag == 0){
	  for (ny = 0, i = 0; ny < Map.NY; ny++) {
		  for (nx = 0; nx < Map.NX; nx++, i++) {
			  elev[ny][nx] = temp[i]; }
	  }
  }
  else if (flag == 1){
	  for (ny = Map.NY - 1, i = 0; ny >= 0; ny--) {
		  for (nx = 0; nx < Map.NX; nx++, i++) {
			  elev[ny][nx] = temp[i]; }
	  }
  }
  else exit (-1);
  free(temp);

  max_elev = 0.0;
  for (ny = 0; ny < nRows; ny++) {
    for (nx = 0; nx < nCols; nx++) {
      if(elev[ny][nx]>max_elev) 
		  max_elev = elev[ny][nx];
    }
  }

  ly=(float)(nRows*dx-dx);
  lx=(float)(nCols*dx-dx);

  printf("beginning skyview calculations \n");
   
  for (ny = 0; ny < nRows; ny++) {
    for (nx = 0; nx < nCols; nx++) {
      skyview[ny][nx]=0.0;
      start_elev=elev[ny][nx];

      if(start_elev>0) {
		  for (i = 0; i < nLook; i++) {
		  theta=6.283185/((double)nLook)*(double)i;
		  sx=(float)nx*dx+0.5*dx;
		  sy=(float)ny*dx+0.5*dx;
		  x=sx;
	      y=sy;
	      max_angle=0.0;
		  
		  while(x>dx && x < lx && y>dx && y<ly) {
			  x=x+((float)cos(theta))*dx;
			  y=y+((float)sin(theta))*dx;
			  dz=elev[(int)(y/dx)][(int)(x/dx)]-start_elev;
	          dist=sqrt((x-sx)*(x-sx)+(y-sy)*(y-sy));
	          if(dz>0) {
				  angle=atan((double)(dz/dist));
				  if(angle > max_angle) 
					  max_angle = angle;
			  }
		  }
		  skyview[ny][nx] +=(cos(max_angle)*cos(max_angle));
		  }
		  skyview[ny][nx]=skyview[ny][nx]/(float)nLook;
	  }
	  ((float *) Array)[ny * nCols + nx] = skyview[ny][nx];
	}
  }

  CreateMapFileNetCDF(DMap.FileName, DMap.FileLabel, &Map);
  Write2DMatrixNetCDF(DMap.FileName, (void *)Array, DMap.NumberType, Map.NY, Map.NX, &DMap, 0);

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