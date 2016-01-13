/*
 * SUMMARY:      StoreModelState.c - Store the state of the model
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Store the state of the model.  This allows restarts of the
 *               model with the correct initial conditions
 * DESCRIP-END.
 * FUNCTIONS:    StoreModelState()
 * COMMENTS:
 * $Id: StoreModelState.c,v 1.8 2004/08/16 18:26:38 colleen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"
#include "varid.h"

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
void StoreModelState(char *Path, DATE * Current, MAPSIZE * Map,
		     OPTIONSTRUCT * Options, TOPOPIX ** TopoMap,
		     PRECIPPIX ** PrecipMap, SNOWPIX ** SnowMap,
		     MET_MAP_PIX ** MetMap, RADCLASSPIX ** RadMap,
		     VEGPIX ** VegMap, LAYER * Veg, SOILPIX ** SoilMap,
		     LAYER * Soil, ROADSTRUCT ** Network, 
		     UNITHYDRINFO * HydrographInfo, float *Hydrograph,
		     CHANNEL * ChannelData)
{
  const char *Routine = "StoreModelState";
  char Str[NAMESIZE + 1];
  char FileLabel[MAXSTRING + 1];
  char FileName[NAMESIZE + 1];
  FILE *HydroStateFile;
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */
  MAPDUMP DMap;			/* Dump Info */
  void *Array;
  float RoadIExcess = 0.0;
 
  /* print a message to stdout that state is being stored */

  printf("Storing model state\n");
  PrintDate(Current, stdout);
  printf("\n");

  if (MetMap != NULL) {

    sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Current->Month, Current->Day,
	    Current->Year, Current->Hour, Current->Min, Current->Sec);
    sprintf(FileName, "%sMet.State.%s%s", Path, Str, fileext);
    strcpy(FileLabel, "Basic Meteorology at time step");

    CreateMapFile(FileName, FileLabel, Map);

    if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
      ReportError((char *) Routine, 1);

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {

	  ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].Precip;

	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 201;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {

	  ((float *) Array)[y * Map->NX + x] = MetMap[y][x].accum_precip;

	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 701;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {

	  ((float *) Array)[y * Map->NX + x] = MetMap[y][x].air_temp;

	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 702;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {

	  ((float *) Array)[y * Map->NX + x] = MetMap[y][x].wind_speed;

	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 703;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {

	  ((float *) Array)[y * Map->NX + x] = MetMap[y][x].humidity;

	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 704;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {

	  ((float *) Array)[y * Map->NX + x] = RadMap[y][x].Beam +
	    RadMap[y][x].Diffuse;

	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 303;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

    free(Array);
  }
  /* Store the canopy interception */

  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Current->Month, Current->Day,
	  Current->Year, Current->Hour, Current->Min, Current->Sec);
  sprintf(FileName, "%sInterception.State.%s%s", Path, Str, fileext);
  strcpy(FileLabel, "Interception storage for each vegetation layer");

  CreateMapFile(FileName, FileLabel, Map);

  if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *) Routine, 1);

  for (i = 0; i < Veg->MaxLayers; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	  if (i < NVeg)
	    ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].IntRain[i];
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 202;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);
  }

  for (i = 0; i < Veg->MaxLayers; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	  if (i < NVeg)
	    ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].IntSnow[i];
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 203;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].TempIntStorage;
      }
      else {
	((float *) Array)[y * Map->NX + x] = NA;
      }
    }
  }
  DMap.ID = 204;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  free(Array);

  /* Store the snow pack conditions */

  sprintf(FileName, "%sSnow.State.%s%s", Path, Str, fileext);
  strcpy(FileLabel, "Snow pack moisture and temperature state");
  CreateMapFile(FileName, FileLabel, Map);

  if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *) Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = (float) SnowMap[y][x].HasSnow;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 401;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = (float) SnowMap[y][x].LastSnow;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 403;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SnowMap[y][x].Swq;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 404;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SnowMap[y][x].PackWater;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 406;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TPack;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 407;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SnowMap[y][x].SurfWater;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 408;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TSurf;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 409;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SnowMap[y][x].ColdContent;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 410;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  free(Array);

  /* Store the soil conditions */

  sprintf(FileName, "%sSoil.State.%s%s", Path, Str, fileext);
  strcpy(FileLabel, "Soil moisture and temperature state");
  CreateMapFile(FileName, FileLabel, Map);

  if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError((char *) Routine, 1);

  for (i = 0; i < Soil->MaxLayers + 1; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	  if (i <= NSoil)
	    ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Moist[i];
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 501;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TSurf;
      else
	((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TSurf;
    }
  }
  DMap.ID = 505;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (i = 0; i < Soil->MaxLayers; i++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  NSoil = Soil->NLayers[SoilMap[y][x].Soil - 1];
	  if (i < NSoil)
	    ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Temp[i];
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
	else
	  ((float *) Array)[y * Map->NX + x] = NA;
      }
    }
    DMap.ID = 511;
    DMap.Layer = i;
    DMap.Resolution = MAP_OUTPUT;
    strcpy(DMap.FileName, "");
    GetVarAttr(&DMap);
    Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
	((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qst;
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 510;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)){
	RoadIExcess = 0.0; 
	if(Options->RoadRouting){
	  if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
	    for (i = 0; i < CELLFACTOR; i++)
	      RoadIExcess += (Network[y][x].h[i] * Network[y][x].RoadArea)/
		((float)CELLFACTOR * (Map->DX*Map->DY));
	  } 
	}
	((float *) Array)[y * Map->NX + x] = SoilMap[y][x].IExcess + RoadIExcess;
      }
      else
	((float *) Array)[y * Map->NX + x] = NA;
    }
  }
  DMap.ID = 512;
  DMap.Resolution = MAP_OUTPUT;
  strcpy(DMap.FileName, "");
  GetVarAttr(&DMap);
  Write2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, &DMap, 0);

  free(Array);

  /* If the unit hydrograph is used for flow routing, store the unit 
     hydrograph array */

  if (Options->Extent == BASIN && Options->HasNetwork == FALSE) {
    sprintf(FileName, "%sHydrograph.State.%s", Path, Str);
    OpenFile(&HydroStateFile, FileName, "w", FALSE);
    for (i = 0; i < HydrographInfo->TotalWaveLength; i++)
      fprintf(HydroStateFile, "%f\n", Hydrograph[i]);
    fclose(HydroStateFile);
  }
}
