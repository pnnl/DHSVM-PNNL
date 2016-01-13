/*
 * SUMMARY:      MakeModelStateNetCDF.c - Create initial model state for DHSVM
 * USAGE:        MakeModelStateBin <infofile>
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD:     02/06/2013
 * $Id:          MakeModelStateNetCDF.c, v 3.1.1  2013/2/6  Ning Exp $
 * DESCRIPTION:  Create initial model state for DHSVM (netcdf I/O format)
 * DESCRIP-END.
 * FUNCTIONS:    ()
 * COMMENTS:
 */

static const char vcid[] = "$Id: MakeModelStateNetCDF.c,v 3.1";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include "sizeofNetCDF.h"
#include "fileio.h"
#include "settings.h"
#include "fifoNetCDF.h"
#include "Calendar.h"
#include "DHSVMerror.h"
#include "data.h"
#include "init.h"
#define MAXLAYERS 10
#define MONTHSPYR 12

void StoreModelState(char *Path, DATE Current, int NY, int NX, 
                     int NVegLayers, float *RainInt, float *SnowInt,
                     float TempIntStorage, unsigned char SnowMask, 
                     unsigned short LastSnow, float Swq, float LwBottom, 
                     float TBottom, float LwTop, float TTop, float Cold, 
                     int NSoilLayers, float *Moist, float SoilTSurf,
                     float *Temp, float GroundHeat, float Runoff,
					 MAPSIZE Map, MAPDUMP DMap);
float GetFloat(char *numberStr);

int CopyDouble(double *Value, char *Str, const int NValues);

/*****************************************************************************
  MakeModelStateNetCDF.c()

  Make a model state file to use as an initialization for DHSVM

  The state variables for DHSVM include the following variables:

    - Rain interception for each vegetation layer
    - Snow interception for top vegetation layer
    - Temporary interception (used for mass release algorithm)

    - Snow pack conditions:
      - presence/absence
      - number of days since last snowfall (used in albedo calculation)
      - snow water equivalent
      - for each layer of the snow pack:
        - liquid water content
        - temperature
      - cold content

    - Soil conditions:
      - for each soil layer:
        - soil moisture (also for the layer below the deepest root zone)
        - temperature
      - surface temperature
      - ground heat storage
*****************************************************************************/
int main(int argc, char *argv[])
{
  char InfoFileName[MAXSTRING+1];
  char Path[MAXSTRING+1];
  DATE Day;
  int NY;
  int NX;
  float dx; /* cell size */
  int NVegLayers;
  int i;
  float RainInt[MAXLAYERS];
  float SnowInt[MAXLAYERS];
  float TempIntStorage;
  unsigned char SnowMask;
  unsigned short LastSnow;
  float Swq;
  float LwBottom;
  float TBottom;
  float LwTop;
  float TTop;
  float Cold;
  int NSoilLayers;
  float Moist[MAXLAYERS];
  float SoilTSurf;
  float Temp[MAXLAYERS];
  float GroundHeat;
  float Runoff;
  int Junk;
  double Xorig;					
  double Yorig;					
  FILE *InfoFile;
  MAPSIZE *Map;
  MAPDUMP *DMap;
  
  if (!(DMap = (MAPDUMP *) calloc(1, sizeof(MAPDUMP))))
    exit(-1);
  if (!(Map = (MAPSIZE *) calloc(1, sizeof(MAPSIZE))))
    exit(-1);

  if (argc != 5) {
    fprintf(stderr, "Usage: MakeModelState <infofile> <cellsize> <Xorig> <Yorig>\n");
	fprintf(stderr, "The cellsize (in the same units as the DEM elevation\n");
	fprintf(stderr, "Xorig is the extreme west and Yorig is the extreme north");
    fprintf(stderr, "The info file MUST contain the following information:\n");
    fprintf(stderr, " - path for output file\n");
    fprintf(stderr, " - date for the model state, in mm/dd/yyyy-hh\n");
    fprintf(stderr, " - number of rows (ny) and number of columns (nx)\n");
    fprintf(stderr, " - maximum number of vegetation layers\n");
    fprintf(stderr, " - rain interception in m for each vegetation layer\n");
    fprintf(stderr, " - snow interception in m for top vegetation layer\n");
    fprintf(stderr, " - snow cover mask\n");
    fprintf(stderr, " - number of days since last snow fall\n");
    fprintf(stderr, " - snow water equivalent in m\n");
    fprintf(stderr, " - liquid water content in m of bottom layer of snowpack\n");
    fprintf(stderr, " - temperature in C of bottom layer of snow pack\n");
    fprintf(stderr, " - liquid water content in m of top layer of snowpack\n");
    fprintf(stderr, " - temperature in C of top layer of snow pack\n");
    fprintf(stderr, " - cold content of snow pack\n");
    fprintf(stderr, " - maximum number of root zone layers\n");
    fprintf(stderr, " - volumetric soil moisture content for each layer\n");
    fprintf(stderr, "   (including the layer below the lowest root zone layer)\n");
    fprintf(stderr, " - temperature in C at soil surface\n");
    fprintf(stderr, " - soil temperature in C for each root zone layer\n");
    fprintf(stderr, " - ground heat storage\n");
    fprintf(stderr, " - runoff\n");
    exit(1);
  }
  strcpy(InfoFileName, argv[1]); 
  if (!(InfoFile = fopen(InfoFileName, "r"))) {
    fprintf(stderr, "Canot open info file %s\n", InfoFileName);
    exit(1);
  }
  
  dx = GetFloat(argv[2]); /* the cellsize of the dem */   
  
  /* extreme west coordinatee */
  if (!(CopyDouble(&(Map->Xorig), argv[3], 1)))
	  exit (-1);;
  /* exterme north coordinate */
  if (!(CopyDouble(&(Map->Yorig), argv[4], 1)))
	  exit (-1);;

  fscanf(InfoFile, "%s", Path);
  ScanDate(InfoFile, &Day);
  fscanf(InfoFile, "%d %d", &NY, &NX);
  fscanf(InfoFile, "%d", &NVegLayers);
  for (i = 0; i < NVegLayers; i++) 
    fscanf(InfoFile, "%f", &RainInt[i]);
  fscanf(InfoFile, "%f", &SnowInt[0]);
  for (i = 1; i < NVegLayers; i++) 
    SnowInt[i] = 0;
  TempIntStorage = 0.0;
  fscanf(InfoFile, "%d", &Junk);
  SnowMask = (unsigned char) Junk;
  fscanf(InfoFile, "%d", &Junk);
  LastSnow = (unsigned short) Junk;
  fscanf(InfoFile, "%f", &Swq);
  fscanf(InfoFile, "%f", &LwBottom);
  fscanf(InfoFile, "%f", &TBottom);
  fscanf(InfoFile, "%f", &LwTop);
  fscanf(InfoFile, "%f", &TTop);
  fscanf(InfoFile, "%f", &Cold);
  fscanf(InfoFile, "%d", &NSoilLayers);
  for (i = 0; i <= NSoilLayers; i++)
    fscanf(InfoFile, "%f", &Moist[i]);
  fscanf(InfoFile, "%f", &SoilTSurf);
  for (i = 0; i < NSoilLayers; i++)
    fscanf(InfoFile, "%f", &Temp[i]);
  fscanf(InfoFile, "%f", &GroundHeat); 
  fscanf(InfoFile, "%f", &Runoff);

  Map->X = 0;
  Map->Y = 0;
  Map->OffsetX = 0;
  Map->OffsetY = 0;
  Map->NX = NX;
  Map->NY = NY;
  Map->DX = dx;
  Map->DY = dx;
  Map->DXY = (float) sqrt(Map->DX * Map->DX + Map->DY * Map->DY);

  StoreModelState(Path, Day, NY, NX, NVegLayers, 
                  RainInt, SnowInt,TempIntStorage, SnowMask,
                  LastSnow, Swq, LwBottom, TBottom, LwTop, TTop, Cold,
                  NSoilLayers, Moist, SoilTSurf, Temp, GroundHeat, Runoff,
				  *Map, *DMap);

  return EXIT_SUCCESS;
}

/*****************************************************************************
  StoreModelState()

  Store the current state of the model.

  The state variables for DHSVM include the following variables:

    - Canopy interception for each vegetation layer

    - Snow pack conditions:
      - presence/absence
      - number of days since last snowfall (used in albedo calculation)
      - snow water equivalent
      - for each layer of the snow pack:
        - liquid water content
        - temperature
      - cold content

    - Soil conditions:
      - for each soil layer:
        - soil moisture (also for the layer below the deepest root zone)
        - temperature
      - surface temperature
      - ground heat storage
*****************************************************************************/
void StoreModelState(char *Path, DATE Current, int NY, int NX, 
                     int NVegLayers, float *RainInt, float *SnowInt,
                     float TempIntStorage, unsigned char SnowMask, 
                     unsigned short LastSnow, float Swq, float LwBottom, 
                     float TBottom, float LwTop, float TTop, float Cold, 
                     int NSoilLayers, float *Moist, float SoilTSurf,
                     float *Temp, float GroundHeat, float Runoff, 
					 MAPSIZE Map, MAPDUMP DMap)
{
  const char *Routine = "StoreModelState";
  char Str[NAMESIZE+1];
  char FileLabel[MAXSTRING+1];
  char FileName[NAMESIZE+1];
  char VarName[BUFSIZE + 1], LongName[BUFSIZE + 1];
  char Units[MAXSTRING+1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  void *Array;

  /* print a message to stdout that state is being stored */
  printf("Storing model state\n");
 
  /* Store the canopy interception */
  sprintf(Str, "%02d.%02d.%04d.%02d.00.00", Current.Month, Current.Day, 
	  Current.Year, Current.Hour);
  MakeFileNameNetCDF(Path, "Interception.State.", Str, FileName); 
  strcpy(FileLabel, "Interception storage for each vegetation layer");
  CreateMapFileNetCDF(FileName, FileLabel, &Map);
  if (!(Array = (float *) calloc(NY * NX, sizeof(float)))) {
    perror(Routine);
    exit(1);
  }
  for (i = 0; i < NVegLayers; i++) {
    for (y = 0; y < NY; y++) {
      for (x = 0; x < NX; x++) { 
		  ((float *) Array)[y * NX + x] = RainInt[i];}
	}
	DMap.ID = 202;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	strcpy(VarName, "Precip.IntRain");
	sprintf(DMap.Name, "%d.%s", i, VarName);
	strcpy(LongName, "Interception Storage (liquid)");
	sprintf(DMap.LongName, "%s (Layer %d)", LongName, i);
	strcpy(DMap.Format, "%.4g");
	strcpy(DMap.Units, "m");
	sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
	strncpy(DMap.FileName, Str, BUFSIZE);
	strcpy(DMap.FileLabel, "Interception Storage (liquid)");
	DMap.NumberType = NC_FLOAT;
	Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);
  }
  
  for (i = 0; i < NVegLayers; i++) {
    for (y = 0; y < NY; y++) {
      for (x = 0; x < NX; x++) {
		  ((float *)Array)[y*NX + x] = SnowInt[i];}
	}
	DMap.ID = 203;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	strcpy(VarName, "Precip.IntSnow");
	sprintf(DMap.Name, "%d.%s", i, VarName);
	strcpy(LongName, "Interception Storage (frozen)");
	sprintf(DMap.LongName, "%s (Layer %d)", LongName, i);
	strcpy(DMap.Format, "%.4g");
	strcpy(DMap.Units, "m");
	sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
	strncpy(DMap.FileName, Str, BUFSIZE);
	strcpy(DMap.FileLabel, "Interception storage (frozen)");
	DMap.NumberType = NC_FLOAT;
	Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);
  }

  for (y = 0; y < NY; y++) 
    for (x = 0; x < NX; x++) 
      ((float *)Array)[y * NX + x] = TempIntStorage;
  DMap.ID = 204;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Temp.Instor");
  strcpy(DMap.LongName, "Temporary interception storage for top vegetation layer");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "m");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Temporary interception storage for top vegetation layer");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);
  free(Array);

  /* Store the snow pack conditions */ 
  sprintf(Str, "%02d.%02d.%04d.%02d.00.00", Current.Month, Current.Day, 
	  Current.Year, Current.Hour);
  MakeFileNameNetCDF(Path, "Snow.State.", Str, FileName); 
  strcpy(FileLabel, "Snow pack moisture and temperature state");
  CreateMapFileNetCDF(FileName, FileLabel, &Map);

  if (!(Array = (float *) calloc(NY * NX, sizeof(float)))) {
    perror(Routine);
    exit(1);
  }
  
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = (float)SnowMask;
  DMap.ID = 401;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.HasSnow");
  strcpy(DMap.LongName, "Snow Presence/Absence");
  strcpy(DMap.Format, "%1d");
  strcpy(DMap.Units, "");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Snow cover flag");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = (float)LastSnow;
  DMap.ID = 403;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.LastSnow");
  strcpy(DMap.LongName, "Last Snowfall");
  strcpy(DMap.Format, "%4d");
  strcpy(DMap.Units, "days");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Days since last snowfall");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = Swq;
  DMap.ID = 404;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.Swq");
  strcpy(DMap.LongName, "Snow Water Equivalent");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "m");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Snow water equivalent");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = LwBottom;
  DMap.ID = 406;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.PackWater");
  strcpy(DMap.LongName, "Liquid Water Content (Deep Layer)");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "m");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Liquid water content of snow pack");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);
    
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = TBottom;
  DMap.ID = 407;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.TPack");
  strcpy(DMap.LongName, "Snow Temperature (Deep Layer)");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "C");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Temperature of snow pack");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = LwTop;
  DMap.ID = 408;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.SurfWater");
  strcpy(DMap.LongName, "Liquid Water Content (Surface Layer)");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "m");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Liquid water content of surface layer");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = TTop;
  DMap.ID = 409;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.TSurf");
  strcpy(DMap.LongName, "Snow Temperature (Surface Layer)");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "C");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Temperature of snow pack surface layer");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = Cold;;
  DMap.ID = 410;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Snow.ColdContent");
  strcpy(DMap.LongName, "Snow Cold Content");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "J");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Cold content of snow pack");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  free(Array);

  /* Store the soil conditions */   
  sprintf(Str, "%02d.%02d.%04d.%02d.00.00", Current.Month, Current.Day, 
	  Current.Year, Current.Hour);
  MakeFileNameNetCDF(Path, "Soil.State.", Str, FileName); 
  strcpy(FileLabel, "Soil moisture and temperature state");
  CreateMapFileNetCDF(FileName, FileLabel, &Map);

  if (!(Array = (float *) calloc(NY * NX, sizeof(float)))) {
    perror(Routine);
    exit(1);
  }
  
  for (i = 0; i < NSoilLayers+1; i++) {
    for (y = 0; y < NY; y++) 
      for (x = 0; x < NX; x++) 
        ((float *)Array)[y*NX + x] = Moist[i];
	DMap.ID = 501;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	strcpy(VarName, "Soil.Moist");
	sprintf(DMap.Name, "%d.%s", i, VarName);
	strcpy(LongName, "Soil Moisture Content");
	sprintf(DMap.LongName, "%s (Layer %d)", LongName, i);
	strcpy(DMap.Format, "%.4g");
	strcpy(DMap.Units, "");
	sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
	strncpy(DMap.FileName, Str, BUFSIZE);
	strcpy(DMap.FileLabel, "Soil moisture");
	DMap.NumberType = NC_FLOAT;
	Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);
  }
  
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = SoilTSurf;
  DMap.ID = 505;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Soil.TSurf");
  strcpy(DMap.LongName, "Surface Temperature");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "C");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Soil surface temperature");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (i = 0; i < NSoilLayers; i++) {
    for (y = 0; y < NY; y++) 
      for (x = 0; x < NX; x++) 
        ((float *)Array)[y*NX + x] = Temp[i];
	DMap.ID = 511;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	strcpy(VarName, "Soil.Temp");
	sprintf(DMap.Name, "%d.%s", i, VarName);
	strcpy(LongName, "Soil Temperature");
	sprintf(DMap.LongName, "%s (Layer %d)", LongName, i);
	strcpy(DMap.Format, "%.4g");
	strcpy(DMap.Units, "C");
	sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
	strncpy(DMap.FileName, Str, BUFSIZE);
	strcpy(DMap.FileLabel, "Soil Temperature");
	DMap.NumberType = NC_FLOAT;
	Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);
  }
  
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = GroundHeat;
  DMap.ID = 510;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Soil.Qst");
  strcpy(DMap.LongName, "Ground Heat Storage");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "W/m2");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Ground heat storage");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = Runoff;
  DMap.ID = 512;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  strcpy(DMap.Name, "Soil.Runoff");
  strcpy(DMap.LongName, "Surface Ponding");
  strcpy(DMap.Format, "%.4g");
  strcpy(DMap.Units, "m");
  sprintf(Str, "%sMap.%s.nc", DMap.FileName, DMap.Name);
  strncpy(DMap.FileName, Str, BUFSIZE);
  strcpy(DMap.FileLabel, "Surface Ponding");
  DMap.NumberType = NC_FLOAT;
  Write2DMatrixNetCDF(FileName, Array, DMap.NumberType, NY, NX, &DMap, 0);

  free(Array);
}

/*****************************************************************************
  GetFloat()
*****************************************************************************/
float GetFloat(char *numberStr) 
{
  char *endPtr;
  float number = 0;

  number = (float) strtod(numberStr, &endPtr);
  if (*endPtr != '\0'){
    printf("problem extracting float from %s \n",numberStr);
    exit(-1);
  }

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











