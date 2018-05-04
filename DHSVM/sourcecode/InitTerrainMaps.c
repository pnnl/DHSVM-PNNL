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

#include <ga.h>
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
#include "ParallelDHSVM.h"

 /*****************************************************************************
   InitTerrainMaps()
 *****************************************************************************/
void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE * GMap, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, TOPOPIX ***TopoMap, SOILTABLE *SType, SOILPIX ***SoilMap, 
  VEGTABLE *VType, VEGPIX ***VegMap)

{
  if (ParallelRank() == 0) printf("\nInitializing terrain maps\n");

  InitTopoMap(Input, Options, GMap, Map, TopoMap);
  InitSoilMap(Input, Options, Map, Soil, *TopoMap, SoilMap);
  InitVegMap(Options, Input, Map, VegMap);
  if (Options->CanopyGapping)
    InitCanopyGapMap(Options, Input, Map, Soil, Veg, VType, VegMap, SType, SoilMap);
}

/*****************************************************************************
  InitTopoMap()
*****************************************************************************/
void InitTopoMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * GMap, MAPSIZE * Map,
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
  int dodump;                   /* Flag to dump topography */
  int masked_decomposition;     /* Flag to do domain decomposition using the mask */
  int striped;                  /* Flag to do just stripe the masked domain */
  MAPSIZE TMap;                 /* temporary local domain */
  STRINIENTRY StrEnv[] = {
    {"TERRAIN", "DEM FILE", "", ""},
    {"TERRAIN", "BASIN MASK FILE", "", ""},
    {"TERRAIN", "DUMP TOPO", "", "FALSE"},
    {"TERRAIN", "DECOMPOSITION", "", "STRIPED"},
    {NULL, NULL, "", NULL}
  };

  /* Process the [TERRAIN] section in the input file */

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* determine how to do domain decomposition, then do it */
  if (strncmp(StrEnv[decompose].VarStr, "SIMPLE", 6) == 0) {
    masked_decomposition = FALSE;
  } else if (strncmp(StrEnv[decompose].VarStr, "MASKED", 6) == 0) {
    masked_decomposition = TRUE;
    striped = 0;
  } else if (strncmp(StrEnv[decompose].VarStr, "STRIPED", 7) == 0) {
    masked_decomposition = TRUE;
    striped = 1;
  } else if (strncmp(StrEnv[decompose].VarStr, "STRIPEX", 7) == 0) {
    masked_decomposition = TRUE;
    striped = 2;
  } else if (strncmp(StrEnv[decompose].VarStr, "STRIPEY", 7) == 0) {
    masked_decomposition = TRUE;
    striped = 3;
  } else {
    ReportError(StrEnv[decompose].KeyName, 51);
  }

  /* let GA decide */
  SimpleDomainDecomposition(GMap, &TMap);

  /* if called for, use the mask to adjust the simple decomposition to
     hopefully produce a better load balance */ 

  if (masked_decomposition && ParallelSize() > 1) {
  
    /* read the mask into an array using the default, simple decomposition */

    GetVarName(002, 0, VarName);
    GetVarNumberType(002, &NumberType);
    if (!(Mask = (unsigned char *)calloc(TMap.NX * TMap.NY,
                                         SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(StrEnv[maskfile].VarStr, Mask, NumberType, &TMap, 0,
                        VarName, 0);

    MaskedDomainDecomposition(GMap, &TMap, Map, striped, Mask);

    free(Mask);
  } else {
    memcpy(Map, &TMap, sizeof(MAPSIZE));
  }


  /* now allocate the topography data structures with appropriate
     decomposition */

  if (!(*TopoMap = (TOPOPIX **)calloc(Map->NY, sizeof(TOPOPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*TopoMap)[y] = (TOPOPIX *)calloc(Map->NX, sizeof(TOPOPIX))))
      ReportError((char *)Routine, 1);
  }

  /* Read the elevation data from the DEM dataset */
  GetVarName(001, 0, VarName);
  GetVarNumberType(001, &NumberType);
  if (!(Elev = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);

  flag = Read2DMatrix(StrEnv[demfile].VarStr, Elev, NumberType, Map, 0,
    VarName, 0);

  /* Assign the attributes to the map pixel */
  /* Reverse the matrix is flag = 1 & netcdf option is selected */
  if ((Options->FileFormat == NETCDF && flag == 0) || (Options->FileFormat == BIN)) {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Dem = Elev[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Dem = Elev[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  free(Elev);

  /* Read the mask */
  GetVarName(002, 0, VarName);
  GetVarNumberType(002, &NumberType);
  if (!(Mask = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[maskfile].VarStr, Mask, NumberType, Map, 0,
    VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Mask = Mask[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Mask = Mask[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  free(Mask);

  /* get flag to dump topography */
  if (strncmp(StrEnv[dumptopo].VarStr, "TRUE", 4) == 0)
    dodump = TRUE;
  else 
    dodump = FALSE;



  /* find out the minimum grid elevation of the basin (using the basin mask) */
  MINELEV = DHSVM_HUGE;
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INBASIN((*TopoMap)[y][x].Mask)) {
        if ((*TopoMap)[y][x].Dem < MINELEV) {
          MINELEV = (*TopoMap)[y][x].Dem;
        }
      }
    }
  }
  /* printf("%d: local MINELEV = %.3f\n", ParallelRank(), MINELEV); */
  GA_Fgop(&MINELEV, 1, "min");
  if (ParallelRank() == 0) 
    printf("global MINELEV = %.3f\n", MINELEV);

  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and
     fraction of flow flowing in each direction based on the land surface
     slope. */
  ElevationSlopeAspect(Map, *TopoMap);
  GMap->NumCells = Map->AllCells;
  GMap->AllCells = Map->AllCells;

  /* After calculating the slopes and aspects for all the points, reset the
     mask if the model is to be run in point mode */
  if (Options->Extent == POINT) {
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
        (*TopoMap)[y][x].Mask = OUTSIDEBASIN;
    (*TopoMap)[Options->PointY][Options->PointX].Mask = (1 != OUTSIDEBASIN);
  }

  /* for debugging
  if (dodump) DumpTopo(Map, *TopoMap);
  */
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
  if (!(*SoilMap = (SOILPIX **)calloc(Map->NY, sizeof(SOILPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SoilMap)[y] = (SOILPIX *)calloc(Map->NX, sizeof(SOILPIX))))
      ReportError((char *)Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the soil type */
  GetVarName(003, 0, VarName);
  GetVarNumberType(003, &NumberType);
  if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[soiltype_file].VarStr, Type, NumberType, 
	Map, 0, VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (((int)Type[i]) > Soil->NTypes)
          ReportError(StrEnv[soiltype_file].VarStr, 32);
        (*SoilMap)[y][x].Soil = Type[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (((int)Type[i]) > Soil->NTypes)
          ReportError(StrEnv[soiltype_file].VarStr, 32);
        (*SoilMap)[y][x].Soil = Type[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);

  /* Read the total soil depth  */
  GetVarName(004, 0, VarName);
  GetVarNumberType(004, &NumberType);
  if (!(Depth = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[soildepth_file].VarStr, Depth, NumberType, 
	Map, 0, VarName, 0);

  /* Assign the attributes to the correct map pixel */
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*SoilMap)[y][x].Depth = Depth[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*SoilMap)[y][x].Depth = Depth[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (Options->Infiltration == DYNAMIC)
        (*SoilMap)[y][x].InfiltAcc = 0.;
      (*SoilMap)[y][x].MoistInit = 0.;

      /* allocate memory for the number of root layers, plus an additional
       layer below the deepest root layer */
      if (INBASIN(TopoMap[y][x].Mask)) {
        if (!((*SoilMap)[y][x].Moist =
          (float *)calloc((Soil->NLayers[Type[i] - 1] + 1), sizeof(float))))
          ReportError((char *)Routine, 1);
        if (!((*SoilMap)[y][x].Perc =
          (float *)calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
          ReportError((char *)Routine, 1);
        if (!((*SoilMap)[y][x].Temp =
          (float *)calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
          ReportError((char *)Routine, 1);
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
    (unsigned long)BUFSIZE, Input);
  if (!VegMapFileName)
    ReportError("VEGETATION MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(005, 0, VarName);
  GetVarNumberType(005, &NumberType);
  if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(VegMapFileName, Type, NumberType, Map, 0, VarName, 0);

  /* Assign the attributes to the correct map pixel */
  if (!(*VegMap = (VEGPIX **)calloc(Map->NY, sizeof(VEGPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*VegMap)[y] = (VEGPIX *)calloc(Map->NX, sizeof(VEGPIX))))
      ReportError((char *)Routine, 1);
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
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Veg = Type[i];
        (*VegMap)[y][x].Tcanopy = 0.0;
      }
    }
  }
  else ReportError((char *)Routine, 57);

  free(Type);
}


/*****************************************************************************
InitCanopyGapMap()
*****************************************************************************/
void InitCanopyGapMap(OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, VEGTABLE *VType, VEGPIX ***VegMap, 
  SOILTABLE *SType, SOILPIX ***SoilMap)
{
  const char *Routine = "InitCanopyGapMap";
  char VarName[BUFSIZE + 1];
  char CanopyMapFileName[BUFSIZE + 1];
  int i, j;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int flag;
  int NVeg;
  int NSoil;
  int NumberType;		/* number type */
  unsigned char *Gap;		/* presence of gap */

  /* Get the canopy gap map filename from the [VEGETATION] section */
  GetInitString("VEGETATION", "CANOPY GAP MAP FILE", "", CanopyMapFileName,
    (unsigned long)BUFSIZE, Input);
  if (!CanopyMapFileName)
    ReportError("CANOPY GAP MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(007, 0, VarName);
  GetVarNumberType(007, &NumberType);
  if (!(Gap = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(CanopyMapFileName, Gap, NumberType, Map, 0, VarName, 0);

  /* if NetCDF, may need to reverse the matrix */
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Gapping = Gap[i];
        /* set gapping to false for cells with no overstory */
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
          (*VegMap)[y][x].Gapping = 0;
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Gapping = Gap[i];
        /* set gapping to false for cells with no overstory */
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
          (*VegMap)[y][x].Gapping = 0;
        /* set gapping to false given glacier cell */
        if (VType[(*VegMap)[y][x].Veg - 1].Index == GLACIER)
          (*VegMap)[y][x].Gapping = 0;
      }
    }
  }
  else ReportError((char *)Routine, 57);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      NVeg = Veg->MaxLayers;
      NSoil = Soil->MaxLayers;
      if (Options->CanopyGapping) {
        if (!((*VegMap)[y][x].Type = (CanopyGapStruct *)calloc(2, sizeof(CanopyGapStruct))))
          ReportError((char *)Routine, 1);
        for (i = 0; i < CELL_PARTITION; i++) {
          if (!((*VegMap)[y][x].Type[i].IntRain = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].IntSnow = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].Moist = (float *)calloc(NSoil+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EPot = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EAct = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EInt = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].ESoil = (float **)calloc(NVeg, sizeof(float *))))
            ReportError((char *)Routine, 1);

          for (j = 0; j < NVeg; j++) {
            if (!((*VegMap)[y][x].Type[i].ESoil[j] = (float *)calloc(NSoil, sizeof(float))))
              ReportError((char *)Routine, 1);
          }
        }
      }
    }
  }
  free(Gap);
}







