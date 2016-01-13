/*
 * SUMMARY:      InitNewMonth.c - Initialize new time periods
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialization functions that have to be executed at the
 *               beginning of certain timestep
 * DESCRIP-END.
 * FUNCTIONS:    InitNewMonth()
 *               InitNewDay()
 *               InitNewStep()
 * COMMENTS:
 * $Id: InitNewMonth.c,v 3.1 2013/02/06 ning Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "fifobin.h"
#include "fileio.h"
#include "rad.h"
#include "slopeaspect.h"
#include "sizeofnt.h"
#include "varid.h"

/*****************************************************************************
  InitNewMonth()
  At the start of a new month, read the new radiation files 
  (diffuse and direct beam), and potentially a new LAI value.
*****************************************************************************/
void InitNewMonth(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
		  TOPOPIX **TopoMap, float **PrismMap,
		  unsigned char ***ShadowMap, RADCLASSPIX **RadMap, 
		  INPUTFILES *InFiles, int NVegs, VEGTABLE *VType, int NStats,
		  METLOCATION *Stat, char *Path)
{
  const char *Routine = "InitNewMonth";
  char FileName[MAXSTRING + 1];
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;
  int j;
  int y, x;
  float a, b, l;
  int NumberType;
  float *Array = NULL;
  unsigned char *Array1 = NULL;
  int flag;

  if (DEBUG)
    printf("Initializing new month\n");

  /* If PRISM precipitation fields are being used to interpolate the 
     observed precipitation fields, then read in the new months field */

  if (Options->Prism == TRUE) {
    printf("reading in new PRISM field for month %d \n", Time->Current.Month);
    sprintf(FileName, "%s.%02d.%s", Options->PrismDataPath, 
	    Time->Current.Month, Options->PrismDataExt);
    GetVarName(205, 0, VarName);
    GetVarNumberType(205, &NumberType);
    if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
      ReportError((char *) Routine, 1);
    flag = Read2DMatrix(FileName, Array, NumberType, Map->NY, Map->NX, 0, VarName, 0);
	
	if ((Options->FileFormat == NETCDF && flag == 0) 
	  || (Options->FileFormat == BIN)) {
      for (y = 0, i = 0; y < Map->NY; y++) 
        for (x = 0; x < Map->NX; x++, i++) 
          PrismMap[y][x] = Array[i];
	}
    else if (Options->FileFormat == NETCDF && flag == 1){
	  for (y = Map->NY - 1, i = 0; y >= 0; y--) 
		for (x = 0; x < Map->NX; x++, i++) 
		  PrismMap[y][x] = Array[i];
  }
  else ReportError((char *) Routine, 57);
  
  free(Array);
  }

  if (Options->Shading == TRUE) {
    printf("reading in new shadow map for month %d \n", Time->Current.Month);
    sprintf(FileName, "%s.%02d.%s", Options->ShadingDataPath, 
	    Time->Current.Month, Options->ShadingDataExt);
    GetVarName(304, 0, VarName);
    GetVarNumberType(304, &NumberType);
    if (!(Array1 = (unsigned char *) calloc(Map->NY * Map->NX, sizeof(unsigned char))))
      ReportError((char *) Routine, 1);
    for (i = 0; i < Time->NDaySteps; i++) {
      Read2DMatrix(FileName, Array1, NumberType, Map->NY, Map->NX, i, VarName, i);
	  for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
		  ShadowMap[i][y][x] = Array1[y * Map->NX + x];
		}
      }
    }
    free(Array1);
  }

  printf("changing LAI, albedo and diffuse transmission parameters\n");
  for (i = 0; i < NVegs; i++) {
    for (j = 0; j < VType[i].NVegLayers; j++) {
      VType[i].LAI[j] = VType[i].LAIMonthly[j][Time->Current.Month - 1];
      VType[i].MaxInt[j] = VType[i].LAI[j] * VType[i].Fract[j] *
	LAI_WATER_MULTIPLIER;
      VType[i].Albedo[j] = VType[i].AlbedoMonthly[j][Time->Current.Month - 1];
    }
    if (VType[i].OverStory) {
      a = VType[i].LeafAngleA;
      b = VType[i].LeafAngleB;
      l = VType[i].LAI[0] / VType[i].ClumpingFactor;
      if (l == 0)
	VType[i].Taud = 1.0;
      else
	VType[i].Taud =
	  exp(-b * l) * ((1 - a * l) * exp(-a * l) +
			 (a * l) * (a * l) * evalexpint(1, a * l));
    }
    else {
      VType[i].Taud = 0.0;
    }
  }
  

}

/*****************************************************************************
  Function name: InitNewDay()

  Purpose      : Initialize the Earth-Sun geometry variables at the beginning
                 of each day

  Required     : 
    int DayOfYear           - day of year (January 1 = 1)
    SOLARGEOMETRY *SolarGeo - structure with information about Earth-Sun 
                              geometry

  Returns      : void

  Modifies     : 
    SOLARGEOMETRY *SolarGeo

  Comments     : To be excuted at the beginning of each new day
*****************************************************************************/
void InitNewDay(int DayOfYear, SOLARGEOMETRY * SolarGeo)
{
  SolarDay(DayOfYear, SolarGeo->Longitude, SolarGeo->Latitude,
	   SolarGeo->StandardMeridian, &(SolarGeo->NoonHour),
	   &(SolarGeo->Declination), &(SolarGeo->HalfDayLength),
	   &(SolarGeo->Sunrise), &(SolarGeo->Sunset),
	   &(SolarGeo->TimeAdjustment), &(SolarGeo->SunEarthDistance));
}

/*****************************************************************************
  Function name: InitNewStep()

  Purpose      : Initialize Earth-Sun geometry and meteorological data at the
                 beginning of each timestep

  Required     :
    MAPSIZE Map              - Structure with information about location
    TIMESTRUCT Time          - Structure with time information
    int PrecipType           - Type of precipitation input, RADAR, STATION or
                               OROGRAPHIC
    int FlowGradient         - Type of FlowGradient calculation
    int NStats               - Number of meteorological stations
    METLOCATION *Stat        - Structure with information about the
                               meteorological stations in or near the study
                               area 
    char *RadarFileName      - Name of file with radar images
    MAPSIZE Radar            - Structure with information about the
                               precipitation radar coverage
    RADCLASSPIX **RadMap     - Structure with radiation data for each pixel
    RADARPIX **RadarMap      - Structure with precipitation information for
                               each radar pixel
    SOLARGEOMETRY *SolarGeo  - structure with information about Earth-Sun 
                               geometry
    SOILPIX **SoilMap        - structure with soil information
    float ***MM5Input        - MM5 input maps
    float ***WindModel       - Wind model maps
                           
  Returns      : void

  Modifies     :

  Comments     : To be executed at the beginning of each time step
*****************************************************************************/
void InitNewStep(INPUTFILES *InFiles, MAPSIZE *Map, TIMESTRUCT *Time,
		 int NSoilLayers, OPTIONSTRUCT *Options, int NStats,
		 METLOCATION *Stat, char *RadarFileName, MAPSIZE *Radar,
		 RADARPIX **RadarMap, SOLARGEOMETRY *SolarGeo, 
		 TOPOPIX **TopoMap, RADCLASSPIX **RadMap, SOILPIX **SoilMap,
		 float ***MM5Input, float ***WindModel, MAPSIZE *MM5Map)
{
  const char *Routine = "InitNewStep";
  int i;			/* counter */
  int j;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type in MM5 input */
  int Step;			/* Step in the MM5 Input */
  float *Array = NULL;
  int MM5Y, MM5X;

  /*printf("current time is %4d-%2d-%2d-%2d\n", Time->Current.Year,Time->Current.Month, Time->Current.Day, Time->Current.Hour);*/

  /* Calculate variables related to the position of the sun above the
     horizon, this is only necessary if shading is TRUE */

  SolarHour(SolarGeo->Latitude,
	    (Time->DayStep + 1) * ((float) Time->Dt) / SECPHOUR,
	    ((float) Time->Dt) / SECPHOUR, SolarGeo->NoonHour,
	    SolarGeo->Declination, SolarGeo->Sunrise, SolarGeo->Sunset,
	    SolarGeo->TimeAdjustment, SolarGeo->SunEarthDistance,
	    &(SolarGeo->SineSolarAltitude), &(SolarGeo->DayLight),
	    &(SolarGeo->SolarTimeStep), &(SolarGeo->SunMax),
	    &(SolarGeo->SolarAzimuth));

/*printf("SunMax is %f\n",SolarGeo->SunMax);*/
  if (Options->MM5 == TRUE) {

    /* Read the data from the MM5 files */

    if (!(Array = (float *) calloc(MM5Map->NY * MM5Map->NX, sizeof(float))))
      ReportError((char *) Routine, 1);
    NumberType = NC_FLOAT;

    Step = NumberOfSteps(&(Time->StartMM5), &(Time->Current), Time->Dt);

    Read2DMatrix(InFiles->MM5Temp, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
		  MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
		  MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
		  MM5Input[MM5_temperature - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
      }

    Read2DMatrix(InFiles->MM5Humidity, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);

    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
		  MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
		  MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
		  MM5Input[MM5_humidity - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
      }

    Read2DMatrix(InFiles->MM5Wind, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
		  MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
		  MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
		  MM5Input[MM5_wind - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
      }

    Read2DMatrix(InFiles->MM5ShortWave, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
	MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
	MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
	MM5Input[MM5_shortwave - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
      }

    Read2DMatrix(InFiles->MM5LongWave, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
	MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
	MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
	MM5Input[MM5_longwave - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
      }

    Read2DMatrix(InFiles->MM5Precipitation, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
	MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
	MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
	MM5Input[MM5_precip - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
	if (MM5Input[MM5_precip - 1][y][x] < 0.0) {
	  printf("Warning: MM5 precip is less than zero %f\n",
		 MM5Input[MM5_precip - 1][y][x]);
	  MM5Input[MM5_precip - 1][y][x] = 0.0;
	}
      }
    Read2DMatrix(InFiles->MM5Terrain, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
	MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
	MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
	MM5Input[MM5_terrain - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
      }
    Read2DMatrix(InFiles->MM5Lapse, Array, NumberType, MM5Map->NY,
		 MM5Map->NX, Step);
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++) {
	MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
	MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
	MM5Input[MM5_lapse - 1][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
      }

    if (Options->HeatFlux == TRUE) {

      for (i = 0, j = MM5_lapse; i < NSoilLayers; i++, j++) {
	Read2DMatrix(InFiles->MM5SoilTemp[i], Array, NumberType, MM5Map->NY,
		     MM5Map->NX, Step);
	for (y = 0; y < Map->NY; y++)
	  for (x = 0; x < Map->NX; x++) {
	    MM5Y = (int) ((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
	    MM5X = (int) ((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
	    MM5Input[j][y][x] = Array[MM5Y * MM5Map->NX + MM5X];
	  }
      }
    }
    free(Array);
  }
/*end if MM5*/

  /* if the flow gradient is based on the water table, recalculate the water
     table gradients.  Flow directions are now calculated in RouteSubSurface*/
  if (Options->FlowGradient == WATERTABLE) {
    /* Calculate the WaterLevel, i.e. the height of the water table above 
       some datum */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  SoilMap[y][x].WaterLevel =
	    TopoMap[y][x].Dem - SoilMap[y][x].TableDepth;
	}
      }
    }
/*     HeadSlopeAspect(Map, TopoMap, SoilMap); */
  }

  if ((Options->MM5 == TRUE && Options->QPF == TRUE) || Options->MM5 == FALSE)
    GetMetData(Options, Time, NSoilLayers, NStats, SolarGeo->SunMax, Stat,
	       Radar, RadarMap, RadarFileName);
}
