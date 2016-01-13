/*
 * SUMMARY:      InitTerrainMaps() - Initialize terrain coverages
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize terrain coverages
 * DESCRIP-END.
 * FUNCTIONS:    InitTerrainMaps()
 *               InitTopoMap()
 *               InitSoilMap()
 *               InitVegMap()
 * COMMENTS:
 * $Id: InitTerrainMaps.c,v 3.1 2013/2/3 00:08:33 Ning Exp $     
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
#include "getinit.h"
#include "sizeofnt.h"
#include "slopeaspect.h"
#include "varid.h"

/*****************************************************************************
  InitTerrainMaps()
*****************************************************************************/
void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
		     LAYER * Soil, TOPOPIX *** TopoMap, SOILPIX *** SoilMap, VEGPIX *** VegMap)
{
  printf("\nInitializing terrain maps\n");

  InitTopoMap(Input, Options, Map, TopoMap);
  InitSoilMap(Input, Options, Map, Soil, *TopoMap, SoilMap);
  InitVegMap(Options,Input, Map, VegMap);
}

/*****************************************************************************
  InitTopoMap()
*****************************************************************************/
void InitTopoMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
		 TOPOPIX *** TopoMap)
{
  const char *Routine = "InitTopoMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  int flag;         /* either or not reverse the matrix */
  int NumberType;		/* Number type of data set */
  unsigned char *Mask = NULL;	/* Basin mask */
  float *Elev;			/* Surface elevation */
  STRINIENTRY StrEnv[] = {
    {"TERRAIN", "DEM FILE", "", ""},
    {"TERRAIN", "BASIN MASK FILE", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Process the [TERRAIN] section in the input file */
  if (!(*TopoMap = (TOPOPIX **) calloc(Map->NY, sizeof(TOPOPIX *))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*TopoMap)[y] = (TOPOPIX *) calloc(Map->NX, sizeof(TOPOPIX))))
      ReportError((char *) Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the elevation data from the DEM dataset */
  GetVarName(001, 0, VarName);
  GetVarNumberType(001, &NumberType);
  if (!(Elev = (float *) calloc(Map->NX * Map->NY,
				SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);

  flag = Read2DMatrix(StrEnv[demfile].VarStr, Elev, NumberType, Map->NY, Map->NX, 0,
	       VarName, 0);

  /* Assign the attributes to the map pixel */
  /* Reverse the matrix is flag = 1 & netcdf option is selected */
  if ((Options->FileFormat == NETCDF && flag == 0) || (Options->FileFormat == BIN)){
	for (y = 0, i = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++, i++) {
	    (*TopoMap)[y][x].Dem = Elev[i]; }
	}
  }
  else if (Options->FileFormat == NETCDF && flag == 1){
	for (y = Map->NY - 1, i = 0; y >= 0; y--) {
	  for (x = 0; x < Map->NX; x++, i++) {
	    (*TopoMap)[y][x].Dem = Elev[i]; }
	}
  }
  else ReportError((char *) Routine, 57);
  free(Elev);
  
  /* find out the minimum grid elevation of the basin */
  MINELEV = 9999;
  for (y = 0, i = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++, i++) {
	  if ((*TopoMap)[y][x].Dem < MINELEV) {
		MINELEV = (*TopoMap)[y][x].Dem;
	  }
	}
  }

  /* Read the mask */
  GetVarName(002, 0, VarName);
  GetVarNumberType(002, &NumberType);
  if (!(Mask = (unsigned char *) calloc(Map->NX * Map->NY,
					SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  flag = Read2DMatrix(StrEnv[maskfile].VarStr, Mask, NumberType, Map->NY, Map->NX, 0,
	       VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0) 
	  || (Options->FileFormat == BIN))
  {
	  for (y = 0, i = 0; y < Map->NY; y++) {


		  for (x = 0; x < Map->NX; x++, i++) {
			  (*TopoMap)[y][x].Mask = Mask[i]; }
	  }
  }
  else if (Options->FileFormat == NETCDF && flag == 1){
	  for (y = Map->NY - 1, i = 0; y >= 0; y--) {
		  for (x = 0; x < Map->NX; x++, i++) {
			  (*TopoMap)[y][x].Mask = Mask[i]; }
	  }
  }
  else ReportError((char *) Routine, 57);
  free(Mask);
  
  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and 
     fraction of flow flowing in each direction based on the land surface 
     slope. */
  ElevationSlopeAspect(Map, *TopoMap);

  
  /* After calculating the slopes and aspects for all the points, reset the 
     mask if the model is to be run in point mode */
  if (Options->Extent == POINT) {
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
	(*TopoMap)[y][x].Mask = OUTSIDEBASIN;
    (*TopoMap)[Options->PointY][Options->PointX].Mask = (1 != OUTSIDEBASIN);
  }
}

/*****************************************************************************
  InitSoilMap()
*****************************************************************************/
void InitSoilMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
		 LAYER * Soil, TOPOPIX ** TopoMap, SOILPIX *** SoilMap)
{
  const char *Routine = "InitSoilMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  unsigned char *Type;		/* Soil type */
  float *Depth;			/* Soil depth */
  int flag;
  STRINIENTRY StrEnv[] = {
    {"SOILS", "SOIL MAP FILE", "", ""},
    {"SOILS", "SOIL DEPTH FILE", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Process the filenames in the [SOILS] section in the input file */
  /* Assign the attributes to the correct map pixel */
  if (!(*SoilMap = (SOILPIX **) calloc(Map->NY, sizeof(SOILPIX *))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SoilMap)[y] = (SOILPIX *) calloc(Map->NX, sizeof(SOILPIX))))
      ReportError((char *) Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the soil type */
  GetVarName(003, 0, VarName);
  GetVarNumberType(003, &NumberType);
  if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY,
					SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  flag = Read2DMatrix(StrEnv[soiltype_file].VarStr, Type, NumberType, Map->NY,
	       Map->NX, 0, VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0) 
	  || (Options->FileFormat == BIN))
  {
	  for (y = 0, i = 0; y < Map->NY; y++) {
		  for (x = 0; x < Map->NX; x++, i++) {
			  if (((int) Type[i]) > Soil->NTypes)
				  ReportError(StrEnv[soiltype_file].VarStr, 32);
			  (*SoilMap)[y][x].Soil = Type[i]; }
	  }
  }
  else if (Options->FileFormat == NETCDF && flag == 1){
	  for (y = Map->NY - 1, i = 0; y >= 0; y--) {
		  for (x = 0; x < Map->NX; x++, i++) {
			  if (((int) Type[i]) > Soil->NTypes)
				  ReportError(StrEnv[soiltype_file].VarStr, 32);
			  (*SoilMap)[y][x].Soil = Type[i]; }
	  }
  }
  else ReportError((char *) Routine, 57);

  /* Read the total soil depth  */
  GetVarName(004, 0, VarName);
  GetVarNumberType(004, &NumberType);
  if (!(Depth = (float *) calloc(Map->NX * Map->NY,
				 SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  flag = Read2DMatrix(StrEnv[soildepth_file].VarStr, Depth, NumberType, Map->NY,
	       Map->NX, 0, VarName, 0);

  /* Assign the attributes to the correct map pixel */
  if ((Options->FileFormat == NETCDF && flag == 0) 
	  || (Options->FileFormat == BIN))
  {
	  for (y = 0, i = 0; y < Map->NY; y++) {
		  for (x = 0; x < Map->NX; x++, i++) {
			  (*SoilMap)[y][x].Depth = Depth[i]; }
	  }
  }
  else if (Options->FileFormat == NETCDF && flag == 1){
	  for (y = Map->NY - 1, i = 0; y >= 0; y--) {
		  for (x = 0; x < Map->NX; x++, i++) {
			  (*SoilMap)[y][x].Depth = Depth[i]; }
	  }
  }
  else ReportError((char *) Routine, 57);

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (Options->Infiltration == DYNAMIC)
		  (*SoilMap)[y][x].InfiltAcc = 0.; 
	  (*SoilMap)[y][x].MoistInit = 0.; 

     /* allocate memory for the number of root layers, plus an additional 
	  layer below the deepest root layer */
	  if (INBASIN(TopoMap[y][x].Mask)) {
		  if (!((*SoilMap)[y][x].Moist = 
			  (float *) calloc((Soil->NLayers[Type[i] - 1] + 1), sizeof(float))))
			  ReportError((char *) Routine, 1);
		  if (!((*SoilMap)[y][x].Perc =
			  (float *) calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
			  ReportError((char *) Routine, 1);
		  if (!((*SoilMap)[y][x].Temp =
			  (float *) calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
			  ReportError((char *) Routine, 1);
	  }
      else {
		  (*SoilMap)[y][x].Moist = NULL;
		  (*SoilMap)[y][x].Perc = NULL;
		  (*SoilMap)[y][x].Temp = NULL;
      }
    }
  }
  free(Type);
  free(Depth);
}

/*****************************************************************************
  InitVegMap()
*****************************************************************************/
void InitVegMap(OPTIONSTRUCT * Options, LISTPTR Input, MAPSIZE * Map, VEGPIX *** VegMap)
{
  const char *Routine = "InitVegMap";
  char VarName[BUFSIZE + 1];
  char VegMapFileName[BUFSIZE + 1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int flag;
  int NumberType;		/* number type */
  unsigned char *Type;		/* Vegetation type */

  /* Get the map filename from the [VEGETATION] section */
  GetInitString("VEGETATION", "VEGETATION MAP FILE", "", VegMapFileName,
		(unsigned long) BUFSIZE, Input);
  if (!VegMapFileName)
    ReportError("VEGETATION MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(005, 0, VarName);
  GetVarNumberType(005, &NumberType);
  if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY,
					SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  flag = Read2DMatrix(VegMapFileName, Type, NumberType, Map->NY, Map->NX, 0, VarName, 0);

  /* Assign the attributes to the correct map pixel */
  if (!(*VegMap = (VEGPIX **) calloc(Map->NY, sizeof(VEGPIX *))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*VegMap)[y] = (VEGPIX *) calloc(Map->NX, sizeof(VEGPIX))))
      ReportError((char *) Routine, 1);
  }

  if ((Options->FileFormat == NETCDF && flag == 0) 
	  || (Options->FileFormat == BIN))
  {
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      (*VegMap)[y][x].Veg = Type[i];
      (*VegMap)[y][x].Tcanopy = 0.0;
		  }
	  }
  }
  else if (Options->FileFormat == NETCDF && flag == 1){
	  for (y = Map->NY - 1, i = 0; y >= 0; y--) {
		  for (x = 0; x < Map->NX; x++, i++) {
			  (*VegMap)[y][x].Veg = Type[i]; 
			  (*VegMap)[y][x].Tcanopy = 0.0;
		  }
	  }
  }
  else ReportError((char *) Routine, 57);

  free(Type);
}










