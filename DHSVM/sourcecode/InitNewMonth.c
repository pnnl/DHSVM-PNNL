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
 *               InitNewWaterYear()
 * COMMENTS:
 * $Id: InitNewMonth.c,v 3.1 2013/02/06 ning Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  TOPOPIX **TopoMap, float **PrismMap, unsigned char ***ShadowMap, 
  INPUTFILES *InFiles, int NVegs, VEGTABLE *VType, int NStats,
  METLOCATION *Stat, char *Path, VEGPIX ***VegMap)
{
  const char *Routine = "InitNewMonth";
  char FileName[MAXSTRING + 1];
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;
  int j, jj;
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
    if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(FileName, Array, NumberType, Map, 0, VarName, 0);

    if ((Options->FileFormat == NETCDF && flag == 0)
      || (Options->FileFormat == BIN)) {
      for (y = 0, i = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++, i++)
          PrismMap[y][x] = Array[i];
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--)
        for (x = 0; x < Map->NX; x++, i++)
          PrismMap[y][x] = Array[i];
    }
    else ReportError((char *)Routine, 57);

    free(Array);
  }

  if (Options->Shading == TRUE) {
    printf("reading in new shadow map for month %d \n", Time->Current.Month);
    sprintf(FileName, "%s.%02d.%s", Options->ShadingDataPath,
      Time->Current.Month, Options->ShadingDataExt);
    GetVarName(304, 0, VarName);
    GetVarNumberType(304, &NumberType);
    if (!(Array1 = (unsigned char *)calloc(Map->NY * Map->NX, sizeof(unsigned char))))
      ReportError((char *)Routine, 1);
    for (i = 0; i < Time->NDaySteps; i++) {
	  /* if computational time step is finer than hourly, make the shade factor equal within
	  the hourly interval */
	  if (Time->NDaySteps > 24) {
		jj = round(i / (Time->NDaySteps / 24));
		Read2DMatrix(FileName, Array1, NumberType, Map, jj, VarName, jj);
	  }
	  else   
      Read2DMatrix(FileName, Array1, NumberType, Map, i, VarName, i);
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          ShadowMap[i][y][x] = Array1[y * Map->NX + x];
        }
      }
    }
    free(Array1);
  }

  printf("changing LAI, albedo and diffuse transmission parameters\n");

  for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        for (j = 0; j < VType[(*VegMap)[y][x].Veg - 1].NVegLayers; j++) {
          (*VegMap)[y][x].LAI[j] = (*VegMap)[y][x].LAIMonthly[j][Time->Current.Month - 1];
          /*Due to LAI and FC change, have to change MaxInt to spatial as well*/
          (*VegMap)[y][x].MaxInt[j] = (*VegMap)[y][x].LAI[j] * (*VegMap)[y][x].Fract[j] * LAI_WATER_MULTIPLIER;
        }
      }
		}
	}
  for (i = 0; i < NVegs; i++) {
    if (Options->ImprovRadiation) {
      if (VType[i].OverStory == TRUE) {
        VType[i].ExtnCoeff = VType[i].MonthlyExtnCoeff[Time->Current.Month - 1];
      }
      else
        VType[i].ExtnCoeff = 0.;
    }
    for (j = 0; j < VType[i].NVegLayers; j++) {
      VType[i].Albedo[j] = VType[i].AlbedoMonthly[j][Time->Current.Month - 1];
    }
	if (Options->CanopyRadAtt == VARIABLE) {
      if (VType[i].OverStory) {
        a = VType[i].LeafAngleA;
        b = VType[i].LeafAngleB;
        l = VType[i].LAI[0] / VType[i].ClumpingFactor;
        if (l == 0)
          VType[i].Taud = 1.0;
        else
          VType[i].Taud = exp(-b * l) * ((1 - a * l) * exp(-a * l) +
            (a * l) * (a * l) * evalexpint(1, a * l));
      }
      else {
        VType[i].Taud = 0.0;
      }
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


/******************************************************************************/
/* UpdateMM5Field                                                              */
/******************************************************************************/
static void
UpdateMM5Field(char *input, int Step, MAPSIZE *Map, MAPSIZE *MM5Map,
               float *Array, float **MM5InputField)
{
  int x;
  int y;
  const int NumberType = NC_FLOAT;
  int MM5Y, MM5X;

  Read2DMatrix(input, Array, NumberType, MM5Map, Step, "", 0);
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      MM5Y = (int)((y + MM5Map->OffsetY) * Map->DY / MM5Map->DY);
      MM5X = (int)((x - MM5Map->OffsetX) * Map->DX / MM5Map->DY);
      MM5InputField[y][x] = Array[MM5Y * MM5Map->NX + MM5X];
    }
  }
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
                 TOPOPIX **TopoMap, SOILPIX **SoilMap,
                 float ***MM5Input, float **PrecipLapseMap, 
                 float ***WindModel, MAPSIZE *MM5Map)
{
  const char *Routine = "InitNewStep";
  int i;			/* counter */
  int j;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int Step;			/* Step in the MM5 Input */
  float *Array = NULL;
  int MM5Y, MM5X;
  int rdprecip, rdstep;
  uchar first;
  const int NumberType = NC_FLOAT;

  first = IsEqualTime(&(Time->Current), &(Time->Start));

  /*printf("current time is %4d-%2d-%2d-%2d\n", Time->Current.Year,Time->Current.Month, Time->Current.Day, Time->Current.Hour);*/

  /* Calculate variables related to the position of the sun above the
     horizon, this is only necessary if shading is TRUE */

  SolarHour(SolarGeo->Latitude,
            (Time->DayStep + 1) * ((float)Time->Dt) / SECPHOUR,
            ((float)Time->Dt) / SECPHOUR, SolarGeo->NoonHour,
            SolarGeo->Declination, SolarGeo->Sunrise, SolarGeo->Sunset,
            SolarGeo->TimeAdjustment, SolarGeo->SunEarthDistance,
            &(SolarGeo->SineSolarAltitude), &(SolarGeo->DayLight),
            &(SolarGeo->SolarTimeStep), &(SolarGeo->SunMax),
            &(SolarGeo->SolarAzimuth));

  if (Options->MM5 == TRUE) {
    /* Read the data from the MM5 files */
    if (!(Array = (float *)calloc(MM5Map->NY * MM5Map->NX, sizeof(float))))
      ReportError((char *)Routine, 1);

    Step = NumberOfSteps(&(Time->StartMM5), &(Time->Current), Time->Dt);

    UpdateMM5Field(InFiles->MM5Temp, Step, Map, MM5Map, Array,
                   MM5Input[MM5_temperature - 1]);
    UpdateMM5Field(InFiles->MM5Humidity, Step, Map, MM5Map, Array,
                   MM5Input[MM5_humidity - 1]);
    UpdateMM5Field(InFiles->MM5Wind, Step, Map, MM5Map, Array,
                   MM5Input[MM5_wind - 1]);
    UpdateMM5Field(InFiles->MM5ShortWave, Step, Map, MM5Map, Array,
                   MM5Input[MM5_shortwave - 1]);
    UpdateMM5Field(InFiles->MM5LongWave, Step, Map, MM5Map, Array,
                   MM5Input[MM5_longwave - 1]);
    UpdateMM5Field(InFiles->MM5Precipitation, Step, Map, MM5Map, Array,
                   MM5Input[MM5_precip - 1]);

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (MM5Input[MM5_precip - 1][y][x] < 0.0) {
          printf("Warning: MM5 precip is less than zero %f\n",
                 MM5Input[MM5_precip - 1][y][x]);
          MM5Input[MM5_precip - 1][y][x] = 0.0;
        }
      }
    }

    /* Terrain does not change during the simulation, so only read it
       at step 0 */
    if (first) {
      rdstep = 0;
      UpdateMM5Field(InFiles->MM5Terrain, rdstep, Map, MM5Map, Array,
                     MM5Input[MM5_terrain - 1]);
    }

    if (strlen(InFiles->MM5Lapse) > 0) {
      rdprecip = 0;
      rdstep = 0;

      switch (InFiles->MM5LapseFreq) {
      case (FreqSingle):
        if (first) {
          rdprecip = 1;
          rdstep = 0;
        }
        break;
      case (FreqMonth):
        rdstep = Time->Current.Month - 1;
        rdprecip = 1;
        break;
      case (FreqContinous):
        /* Step unchanged */
        rdprecip = 1;
        rdstep = Step;
        break;
      default:
        ReportError("InitNewStep", 15);
      }
      if (rdprecip) {
        UpdateMM5Field(InFiles->MM5Lapse, rdstep, Map, MM5Map, Array,
                       MM5Input[MM5_lapse - 1]);
      }
      
    } else if (first) {
      
      /* If a MM5 temperature lapse map is not specified, fill the map
         with the domain-wide temperature lapse rate (which must be
         specified). Only need to do this once. */
      
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          MM5Input[MM5_lapse - 1][y][x] = TEMPLAPSE;
        }
      }
    }

    if (Options->HeatFlux == TRUE) {


      for (i = 0, j = MM5_lapse; i < NSoilLayers; i++, j++) {
        UpdateMM5Field(InFiles->MM5SoilTemp[i], Step, Map, MM5Map, Array,
                       MM5Input[j]);
      }
    }
    free(Array);

    /* MM5 precip lapse rate is at the DEM resolution, so needs to be
       read differently */

    if (strlen(InFiles->PrecipLapseFile) > 0) {

      rdprecip = 0;
      

      switch (InFiles->MM5PrecipDistFreq) {
      case (FreqSingle):
        if (first) {
          rdprecip = 1;
          rdstep = 0;
        }
        break;
      case (FreqMonth):
        rdstep = Time->Current.Month - 1;
        rdprecip = 1;
        break;
      case (FreqContinous):
        /* Step unchanged */
        rdprecip = 1;
        rdstep = Step;
        break;
      default:
        ReportError("InitNewStep", 15);
      }

      if (rdprecip) {

        if (!(Array = (float *)calloc(Map->NY * Map->NX, sizeof(float))))
          ReportError((char *)Routine, 1);
        
        
        Read2DMatrix(InFiles->PrecipLapseFile, Array, NumberType, Map, rdstep, "", 0);
        for (y = 0; y < Map->NY; y++) {
          for (x = 0; x < Map->NX; x++) {
            PrecipLapseMap[y][x] = Array[y * Map->NX + x];
          }
        }
        free(Array);
      }
    }

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

/*****************************************************************************
   InitNewWaterYear()
   At the start of a new water year, re-initiate the SWE stats maps 
 *****************************************************************************/
void InitNewWaterYear(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
                TOPOPIX **TopoMap, SNOWPIX **SnowMap)
{
  const char *Routine = "InitNewYear";
  int y, x;
  if (DEBUG)
    printf("Initializing new water year \n");

  /* If PRISM precipitation fields are being used to interpolate the
     observed precipitation fields, then read in the new months field */

  if (Options->SnowStats == TRUE) {
    printf("resetting SWE stats map %d \n", Time->Current.Year);
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          SnowMap[y][x].MaxSwe = 0.0;
          SnowMap[y][x].MaxSweDate = 0;
          SnowMap[y][x].MeltOutDate = 0;
        }
      }
    }
  }
}
