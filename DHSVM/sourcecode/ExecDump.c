/*
 * SUMMARY:      ExecDump.c - Write selected output
 * USAGE:        Part of DHSVM
 *
 * DESCRIPTION:  Write selected output files
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIP-END.
 * FUNCTIONS:    ExecDump()
 *               DumpMap()
 *               DumpPix()
 *               DumpPixSed()
 * COMMENTS:
 * $Id: ExecDump.c, v 4.0  2013/1/5   Ning Exp $       
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "fileio.h"
#include "sizeofnt.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  ExecDump()
*****************************************************************************/
void ExecDump(MAPSIZE *Map, DATE *Current, DATE *Start, OPTIONSTRUCT *Options,
	      DUMPSTRUCT *Dump, TOPOPIX **TopoMap, EVAPPIX **EvapMap, 
		  PIXRAD **RadiationMap, PRECIPPIX ** PrecipMap, RADCLASSPIX **RadMap, 
		  SNOWPIX **SnowMap, MET_MAP_PIX **MetMap, VEGPIX **VegMap, LAYER *Veg, 
		  SOILPIX **SoilMap, SEDPIX **SedMap, ROADSTRUCT **Network, 
		  CHANNEL *ChannelData, FINEPIX ***FineMap, LAYER *Soil, AGGREGATED *Total, 
	      UNITHYDRINFO *HydrographInfo, float *Hydrograph)
{
  int i;			/* counter */
  int j;			/* counter */
  int x;
  int y;
  int ii;			/* FineMap counter */
  int jj;			/* FineMap counter */
  int xx;			/* FineMap x-coordinate */
  int yy;			/* FineMap y-coordinate */
  float overlandinflow;          /* Hillslope erosion that enters the channel network */
  float overroadinflow;          /* Road surface erosion that enters the channel network */
  FINEPIX PixAggFineMap;	/* FineMap quanitities aggregated over a pixel */

  /* dump the aggregated basin values for this timestep */
  DumpPix(Current, IsEqualTime(Current, Start), &(Dump->Aggregate),
	  &(Total->Evap),&(Total->Precip), &(Total->RadClass), &(Total->Snow),
	  &(Total->Soil), Soil->MaxLayers, Veg->MaxLayers, Options);
  fprintf(Dump->Aggregate.FilePtr, " %lu", Total->Saturated);
  fprintf(Dump->Aggregate.FilePtr, "\n");

  if (Options->Sediment)
  DumpPixSed(Current, IsEqualTime(Current, Start), &(Dump->AggregateSediment),
          &(Total->Sediment), &(Total->Road), Total->SedimentOverlandInflow, 
          Total->SedimentOverroadInflow, &(Total->Fine) );

  if (Options->Extent != POINT) {
    /* check whether the model state needs to be dumped at this timestep, and
       dump state if needed */
    if (Dump->NStates < 0) {
      StoreModelState(Dump->Path, Current, Map, Options, TopoMap, PrecipMap,
		      SnowMap, MetMap, RadMap, VegMap, Veg, SoilMap, Soil,
		      Network, HydrographInfo, Hydrograph, ChannelData);
      if (Options->HasNetwork)
	StoreChannelState(Dump->Path, Current, ChannelData->streams);
    }
    else {
      for (i = 0; i < Dump->NStates; i++) {
	if (IsEqualTime(Current, &(Dump->DState[i]))) {
	  StoreModelState(Dump->Path, Current, Map, Options, TopoMap,
			  PrecipMap, SnowMap, MetMap, RadMap, VegMap, Veg,
			  SoilMap, Soil, Network, HydrographInfo, Hydrograph,
			  ChannelData);
	  if (Options->HasNetwork)
	    StoreChannelState(Dump->Path, Current, ChannelData->streams);
	}
	  }
	}

    /* check which pixels need to be dumped, and dump if needed */
    for (i = 0; i < Dump->NPix; i++) {
      y = Dump->Pix[i].Loc.N;
      x = Dump->Pix[i].Loc.E;
      // If we include any FineMap quantities, we need to aggregate them over the
      // pixel, which is on the coarse grid
      if (Options->MassWaste) {
        // Initialize PixAggFineMap
        PixAggFineMap.SatThickness = 0.0;
        PixAggFineMap.DeltaDepth = 0.0;
        PixAggFineMap.Probability = 0.0;
        PixAggFineMap.MassWasting = 0.0;
        PixAggFineMap.MassDeposition = 0.0;
        PixAggFineMap.SedimentToChannel = 0.0;
        // Sum up PixAggFineMap quantities
        for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
          for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
            yy = (int) y*Map->DY/Map->DMASS + ii;
            xx = (int) x*Map->DX/Map->DMASS + jj;
            PixAggFineMap.SatThickness += (*FineMap[yy][xx]).SatThickness;
            PixAggFineMap.DeltaDepth += (*FineMap[yy][xx]).DeltaDepth;
            PixAggFineMap.Probability += (*FineMap[yy][xx]).Probability;
            PixAggFineMap.MassWasting += (*FineMap[yy][xx]).MassWasting;
            PixAggFineMap.MassDeposition += (*FineMap[yy][xx]).MassDeposition;
            PixAggFineMap.SedimentToChannel += (*FineMap[yy][xx]).SedimentToChannel;
          }
        }
        // Normalize PixAggFineMap quantities by # FineMap cells in a pixel
        PixAggFineMap.SatThickness /= Map->DMASS*Map->DMASS;
        PixAggFineMap.DeltaDepth /= Map->DMASS*Map->DMASS;
        PixAggFineMap.Probability /= Map->DMASS*Map->DMASS;
      }

      else {
	PixAggFineMap.SatThickness = -999.;
        PixAggFineMap.DeltaDepth = -999.;
        PixAggFineMap.Probability = -999.;
        PixAggFineMap.MassWasting = -999.;
        PixAggFineMap.MassDeposition = -999.;
        PixAggFineMap.SedimentToChannel = -999.;
      }

      if (Options->InitSedFlag && channel_grid_has_channel(ChannelData->stream_map, x, y)){
	overlandinflow = 0.;
	for(i=0;i<NSEDSIZES;i++) {
	  overlandinflow += ChannelData->stream_map[x][y]->channel->sediment.overlandinflow[i];
	}
      }
      else overlandinflow = -999.;

      if (Options->RoadRouting && channel_grid_has_channel(ChannelData->road_map, x, y) ){
	overroadinflow = 0.;
	for(i=0;i<NSEDSIZES;i++) {
	    overroadinflow += ChannelData->road_map[x][y]->channel->sediment.overroadinflow[i];
	}
      }
      else overroadinflow = -999.;
      
      /* output sediment-related variable at the pixel */
      if ( Options->Sediment) 
        DumpPixSed(Current, IsEqualTime(Current, Start), &(Dump->Pix[i].OutFileSediment),
              &(SedMap[y][x]), &(Network[y][x]), overlandinflow, overroadinflow, &PixAggFineMap);

      /* output variable at the pixel */
      DumpPix(Current, IsEqualTime(Current, Start), &(Dump->Pix[i].OutFile),
              &(EvapMap[y][x]), &(PrecipMap[y][x]),&(RadMap[y][x]), &(SnowMap[y][x]),
              &(SoilMap[y][x]), Soil->NLayers[(SoilMap[y][x].Soil - 1)],
              Veg->NLayers[(VegMap[y][x].Veg - 1)], Options);
      fprintf(Dump->Pix[i].OutFile.FilePtr, "\n");
    }

    /* check which maps need to be dumped at this timestep, and dump maps if needed */
    for (i = 0; i < Dump->NMaps; i++) {
      for (j = 0; j < Dump->DMap[i].N; j++) {
	    if (IsEqualTime(Current, &(Dump->DMap[i].DumpDate[j]))) {
	      fprintf(stdout, "Dumping Maps at ");
	      PrintDate(Current, stdout);
	      fprintf(stdout, "\n");
	      DumpMap(Map, Current, &(Dump->DMap[i]), TopoMap, EvapMap,
		      PrecipMap, RadiationMap, SnowMap, SoilMap, SedMap, FineMap,
		      Soil, VegMap, Veg, Network, Options);
		}
      }
    }
  }
}

/*****************************************************************************
  DumpMap()
*****************************************************************************/
void DumpMap(MAPSIZE *Map, DATE *Current, MAPDUMP *DMap, TOPOPIX **TopoMap,
	     EVAPPIX **EvapMap, PRECIPPIX **PrecipMap, PIXRAD **RadMap,
	     SNOWPIX **SnowMap, SOILPIX **SoilMap, SEDPIX **SedMap, 
		 FINEPIX ***FineMap, LAYER *Soil, VEGPIX **VegMap, LAYER *Veg, 
		 ROADSTRUCT **Network, OPTIONSTRUCT *Options)
{
  const char *Routine = "DumpMap";
  char DataLabel[MAXSTRING + 1];
  float Offset;
  float Range;
  int Index;
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int ii, jj, yy, xx; 		/* counters for FineMap variables */
  void *Array;
  int numPoints;
  char VarIDStr[4];		/* stores VarID for sending to ReportError */

  sprintf(DataLabel, "%02d.%02d.%04d.%02d.%02d.%02d", Current->Month,
	  Current->Day, Current->Year, Current->Hour, Current->Min,
	  Current->Sec);

  /* find out what date we are dumping */
  for (Index = 0; Index < DMap->N; Index++) {
    if (IsEqualTime(Current, &(DMap->DumpDate[Index])))
      break;
  }

  sprintf(VarIDStr, "%d", DMap->ID);

  if(DMap->ID >= 800 && DMap->ID < 900) {
    numPoints = Map->NX * Map->NY * Map->DX * Map->DY / (Map->DMASS*Map->DMASS);
  }
  else {
    numPoints = Map->NX * Map->NY;
  }

  switch (DMap->NumberType) {
  case NC_BYTE:
    if (!(Array = calloc(numPoints, SizeOfNumberType(NC_BYTE))))
      ReportError((char *) Routine, 1);
    break;
  case NC_CHAR:
    if (!(Array = calloc(numPoints, SizeOfNumberType(NC_CHAR))))
      ReportError((char *) Routine, 1);
    break;
  case NC_SHORT:
    if (!(Array = calloc(numPoints, SizeOfNumberType(NC_SHORT))))
      ReportError((char *) Routine, 1);
    break;
  case NC_INT:
    if (!(Array = calloc(numPoints, SizeOfNumberType(NC_INT))))
      ReportError((char *) Routine, 1);
    break;

  case NC_FLOAT:
    if (!(Array = calloc(numPoints, SizeOfNumberType(NC_FLOAT))))
      ReportError((char *) Routine, 1);
    break;
  case NC_DOUBLE:
    if (!(Array = calloc(numPoints, SizeOfNumberType(NC_DOUBLE))))
      ReportError((char *) Routine, 1);
    break;
  default:
    Array = NULL;
    ReportError((char *) Routine, 40);
    break;
  }

  Offset = DMap->MinVal;
  Range = DMap->MaxVal - DMap->MinVal;

  switch (DMap->ID) {

  case 101:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].ETot;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((EvapMap[y][x].ETot - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 102:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      /* soil */
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EPot[NVeg];
	    else if (DMap->Layer <= NVeg)
	      /* vegetation layer */
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EPot[DMap->Layer - 1];
	    else
	      /* vegetation layer not present at this pixel */
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EPot[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EPot[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 103:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EInt[NVeg];
	    else if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EInt[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EInt[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EInt[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 104:
    /* NETCDFWORK: This does not work for NETCDF.  Fix */
    if (DMap->Resolution == MAP_OUTPUT) {
      for (i = 0; i < Soil->MaxLayers; i++) {
	    for (y = 0; y < Map->NY; y++) {
	      for (x = 0; x < Map->NX; x++) {
	        if (INBASIN(TopoMap[y][x].Mask)) {
	          NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	          if (DMap->Layer <= NVeg)
		        ((float *) Array)[y * Map->NX + x] = 
				EvapMap[y][x].ESoil[DMap->Layer - 1][i];
			  else
		        ((float *) Array)[y * Map->NX + x] = NA;
			}
	        else
	          ((float *) Array)[y * Map->NX + x] = NA;
		  }
		}
	    Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY,
		      Map->NX, DMap, Index);
	  }
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (i = 0; i < Soil->MaxLayers; i++) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
	    if (INBASIN(TopoMap[y][x].Mask)) {
	      NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	      if (DMap->Layer <= NVeg)
		((unsigned char *) Array)[y * Map->NX + x] =
		  (unsigned char) ((EvapMap[y][x].ESoil[DMap->Layer - 1][i] -
				    Offset) / Range * MAXUCHAR);
	      else
		((unsigned char *) Array)[y * Map->NX + x] = 0;
	    }
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	}
	Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		      Index);
      }
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 105:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EAct[NVeg];
	    else if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EAct[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer > NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EAct[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EAct[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 201:
    if (DMap->Resolution == MAP_OUTPUT) {
     for (y = 0; y < Map->NY; y++) {
	   for (x = 0; x < Map->NX; x++) {
	     ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].Precip;
	   }
	 }
	 Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		           DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	    for (x = 0; x < Map->NX; x++)
	      ((unsigned char *) Array)[y * Map->NX + x] =
	       (unsigned char) ((PrecipMap[y][x].Precip - Offset) /
			     Range * MAXUCHAR);
	  Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 202:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	    for (x = 0; x < Map->NX; x++) {
	      if (INBASIN(TopoMap[y][x].Mask)) {
	        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	        if (DMap->Layer <= NVeg)
			  ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].IntRain[DMap->Layer - 1];
			else
	          ((float *) Array)[y * Map->NX + x] = NA;
		  }
	      else
	        ((float *) Array)[y * Map->NX + x] = NA;
		}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	    for (x = 0; x < Map->NX; x++) {
	      if (INBASIN(TopoMap[y][x].Mask)) {
	        NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
			if (DMap->Layer <= NVeg)
			  ((unsigned char *) Array)[y * Map->NX + x] = 
			      (unsigned char) ((PrecipMap[y][x].IntRain[DMap->Layer - 1] -
				  Offset) / Range * MAXUCHAR);
			else
	          ((unsigned char *) Array)[y * Map->NX + x] = 0;
		  }
		  else
	        ((unsigned char *) Array)[y * Map->NX + x] = 0;
		}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 203:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		PrecipMap[y][x].IntSnow[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((PrecipMap[y][x].IntSnow[DMap->Layer - 1] -
				  Offset) / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 206:
    if (DMap->Resolution == MAP_OUTPUT) {
     for (y = 0; y < Map->NY; y++) {
	   for (x = 0; x < Map->NX; x++) {
		 ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].SumPrecip;
	   }
	 }
	 Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		           DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	    for (x = 0; x < Map->NX; x++)
	      ((unsigned char *) Array)[y * Map->NX + x] =
	       (unsigned char) ((PrecipMap[y][x].Precip - Offset) /
			     Range * MAXUCHAR);
	  Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 301:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	    for (x = 0; x < Map->NX; x++)
	      ((float *) Array)[y * Map->NX + x] = RadMap[y][x].ObsShortIn;
	  Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	    for (x = 0; x < Map->NX; x++)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		     (unsigned char) ((RadMap[y][x].ObsShortIn - Offset) / Range * MAXUCHAR);
	  Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 302:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = RadMap[y][x].RBMNetShort;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((RadMap[y][x].RBMNetShort - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

    case 303:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = RadMap[y][x].PixelBeam;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((RadMap[y][x].PixelBeam - Offset) /Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 401:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] = SnowMap[y][x].HasSnow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] = SnowMap[y][x].HasSnow;
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 402:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    SnowMap[y][x].SnowCoverOver;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    SnowMap[y][x].SnowCoverOver;
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 403:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned short *) Array)[y * Map->NX + x] = SnowMap[y][x].LastSnow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) (((float) SnowMap[y][x].LastSnow - Offset) / Range
			     * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 404:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].Swq;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].Swq - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 405:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].Melt;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].Melt - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 406:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].PackWater;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].PackWater - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 407:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TPack;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].TPack - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 408:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].SurfWater;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].SurfWater - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 409:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TSurf;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].TSurf - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 410:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].ColdContent;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].ColdContent - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 501:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((float *) Array)[y * Map->NX + x] =
		SoilMap[y][x].Moist[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((SoilMap[y][x].Moist[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 502:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((float *) Array)[y * Map->NX + x] =
		SoilMap[y][x].Perc[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((SoilMap[y][x].Perc[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 503:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TableDepth;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].TableDepth - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 504:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].SatFlow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].SatFlow - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 505:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TSurf;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].TSurf - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 506:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qnet;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qnet - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 507:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qs;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qs - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 508:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qe;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qe - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 509:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qg;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qg - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 510:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qst;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qst - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 513:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].IExcess;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].IExcess - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 514:
    if (Options->Infiltration != DYNAMIC) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].InfiltAcc;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].InfiltAcc - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 801:
    if (!Options->MassWaste) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = (*FineMap[yy][xx]).Dem;
	      }
	    }
	  }
 	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((unsigned char *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] =
		  (unsigned char) (((*FineMap[yy][xx]).Dem - Offset) / Range * MAXUCHAR);
	      }
	    }
	  }
	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 803:
    if (!Options->MassWaste) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = (*FineMap[yy][xx]).SatThickness;
	      }
	    }
	  }
 	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((unsigned char *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] =
		  (unsigned char) (((*FineMap[yy][xx]).SatThickness - Offset) / Range * MAXUCHAR);
	      }
	    }
	  }
	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 804:
    if (!Options->MassWaste) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = (*FineMap[yy][xx]).DeltaDepth;
	      }
	    }
	  }
 	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((unsigned char *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] =
		  (unsigned char) (((*FineMap[yy][xx]).DeltaDepth - Offset) / Range * MAXUCHAR);
	      }
	    }
	  }
	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 805:
    if (!Options->MassWaste) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = (*FineMap[yy][xx]).Probability;
	      }
	    }
	  }
 	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((unsigned char *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] =
		  (unsigned char) (((*FineMap[yy][xx]).Probability - Offset) / Range * MAXUCHAR);
	      }
	    }
	  }
	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 806:
    if (!Options->MassWaste) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = (*FineMap[yy][xx]).SedimentToChannel;
	      }
	    }
	  }
 	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((unsigned char *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] =
		  (unsigned char) (((*FineMap[yy][xx]).SedimentToChannel - Offset) / Range * MAXUCHAR);
	      }
	    }
	  }
	  else {
	    for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	      for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		yy = (int) y*Map->DY/Map->DMASS + ii;
		xx = (int) x*Map->DX/Map->DMASS + jj;
		((float *) Array)[yy * (int)(Map->NX*Map->DX/Map->DMASS) + xx] = NA;
	      }
	    }
          }
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, (int)(Map->NY*Map->DY/Map->DMASS),
		    (int)(Map->NX*Map->DX/Map->DMASS), DMap, Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 901:
    if (!Options->InitSedFlag) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SedMap[y][x].SedFluxOut;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SedMap[y][x].SedFluxOut - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 902:
    if (!Options->InitSedFlag) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SedMap[y][x].Erosion;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SedMap[y][x].Erosion - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;

  case 903:
    if (!Options->RoadRouting) {
      ReportError(VarIDStr, 67);
    }
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = Network[y][x].Erosion;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((Network[y][x].Erosion - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError(VarIDStr, 66);
    break;


  default:
    ReportError(VarIDStr, 66);
    break;
  }

  free(Array);
}

/*****************************************************************************
  DumpPix()
*****************************************************************************/
  /*DumpPix(Current, IsEqualTime(Current, Start), &(Dump->Aggregate),
	  &(Total->Evap),&(Total->Precip), &(Total->RadClass), &(Total->Snow),
	  &(Total->Soil), Soil->MaxLayers, Veg->MaxLayers, Options);*/

void DumpPix(DATE * Current, int first, FILES * OutFile, EVAPPIX * Evap,
             PRECIPPIX * Precip, RADCLASSPIX * Rad, SNOWPIX * Snow,
	         SOILPIX * Soil, int NSoil, int NVeg, OPTIONSTRUCT *Options)
{
  int i, j;			/* counter */

  if (first == 1) {
	  
	  // Main Aggregate Values File
	  fprintf(OutFile->FilePtr, "         Date        ");
	  fprintf(OutFile->FilePtr, "HasSnow SnowCover LastSnow    Swq       Melt   ");
	  fprintf(OutFile->FilePtr, "PackWater TPack ");
	  fprintf(OutFile->FilePtr, " TotEvap  "); /*total evapotranspiration*/
	  for (i = 0; i < NVeg + 1; i++)
		  fprintf(OutFile->FilePtr, "EPot%d ", i);
	  for (i = 0; i < NVeg + 1; i++)
		  fprintf(OutFile->FilePtr, "EAct%d ", i);
	  for (i = 0; i < NVeg; i++)
		  fprintf(OutFile->FilePtr, "EInt%d ", i);
	  for (i = 0; i < NVeg; i++)
		  for (j = 0; j < NSoil; j++)
			  fprintf(OutFile->FilePtr, "ESoil%d%d ", i, j);
	  fprintf(OutFile->FilePtr, "   ESoil   ");

	  fprintf(OutFile->FilePtr, "  Precip(m) ");
	  fprintf(OutFile->FilePtr, " Snow(m) ");
	  
	  for (i = 0; i < NVeg; i++)
		  fprintf(OutFile->FilePtr, "IntRain%d ", i);
	  for (i = 0; i < NVeg; i++)
		  fprintf(OutFile->FilePtr, "IntSnow%d ", i);

	  for (i = 0; i < NSoil; i++)
		  fprintf(OutFile->FilePtr, "  SoilMoist%d", (i+1));
	  for (i = 0; i < NSoil; i++)
		  fprintf(OutFile->FilePtr, "     Perc%d   ", (i+1));
	  fprintf(OutFile->FilePtr, "  TableDepth   SatFlow   Runoff    IMP-DS    IExcess  ");
	  fprintf(OutFile->FilePtr, "SoilTemp Qnet Qs Qe Qg Qst Ra"); 

	  if (Options->Infiltration == DYNAMIC)
		  fprintf(OutFile->FilePtr, " InfiltAcc"); 
	  
	  fprintf(OutFile->FilePtr, " RadBeam    RadDiff  ");
	  fprintf(OutFile->FilePtr, "\n"); 

  }

  /* All variables are dumped in the case of a pixel dump */
  // Main Aggregate Values File

  // Date
  PrintDate(Current, OutFile->FilePtr);

  /* Snow */
  fprintf(OutFile->FilePtr, " %1d %1d %4d %g %g %g %g",
	  Snow->HasSnow, Snow->SnowCoverOver, Snow->LastSnow, Snow->Swq,
	  Snow->Melt, Snow->PackWater, Snow->TPack);
    
  /* fprintf(OutFile->FilePtr, " %7d %5d   %9.3E %9.3E",
		Snow->HasSnow, Snow->LastSnow, Snow->Swq, Snow->Melt); */
  fprintf(OutFile->FilePtr, " %9.3E", Evap->ETot);
  for (i = 0; i < NVeg + 1; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EPot[i]);           /* Potential transpiration */
  for (i = 0; i < NVeg + 1; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EAct[i]);           /* Actual transpiration */
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EInt[i]);
  for (i = 0; i < NVeg; i++)
    for (j = 0; j < NSoil; j++)
      fprintf(OutFile->FilePtr, " %g", Evap->ESoil[i][j]);    /*transpiration from each veg layer*/
  fprintf(OutFile->FilePtr, " %9.3E", Evap->EvapSoil);

  fprintf(OutFile->FilePtr, " %9.3E", Precip->Precip);
  fprintf(OutFile->FilePtr, " %9.3E", Precip->SnowFall);

  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Precip->IntRain[i]);
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Precip->IntSnow[i]);


  for (i = 0; i < NSoil; i++)
    fprintf(OutFile->FilePtr, " %9.3E ", Soil->Moist[i]);
  for (i = 0; i < NSoil; i++)
    fprintf(OutFile->FilePtr, " %9.3E ", Soil->Perc[i]);
  fprintf(OutFile->FilePtr, " %9.2E %9.2E %9.2E %9.2E %9.2E", Soil->TableDepth,
	  Soil->SatFlow, Soil->Runoff, Soil->DetentionStorage, Soil->IExcess);
  
  fprintf(OutFile->FilePtr, " %g %g %g %g %g %g %g", 
	  Soil->TSurf, Soil->Qnet, Soil->Qs, Soil->Qe, Soil->Qg, Soil->Qst,
	  Soil->Ra);

  if (Options->Infiltration == DYNAMIC)
    fprintf(OutFile->FilePtr, " %g", Soil->InfiltAcc);

  fprintf(OutFile->FilePtr, " %9.3E %9.3E", Rad->Beam, Rad->Diffuse);
}
/*****************************************************************************
  DumpPixSed()
*****************************************************************************/
void DumpPixSed(DATE * Current, int first, FILES * OutFileSediment,
             SEDPIX * SedMap, ROADSTRUCT * Network, 
             float SedimentOverlandInflow, float SedimentOverroadInflow,
             FINEPIX * FineMap)
{
  if (first == 1) {
    // Sediment Values File
    fprintf(OutFileSediment->FilePtr, "Date ");
    fprintf(OutFileSediment->FilePtr, "SatThick DeltaDepth Probability TotMassWasting TotMassDeposition TotSedToChannel ");
    fprintf(OutFileSediment->FilePtr, "Erosion SedFluxOut OverlandInflow ");
    fprintf(OutFileSediment->FilePtr, "RdErosion RdSedToHill OverroadInflow \n");
  }

  // Date
  PrintDate(Current, OutFileSediment->FilePtr);

  // Sediment
  fprintf(OutFileSediment->FilePtr, " %g %g %g %g %g %g %g %g %g %g %g %g\n",
            FineMap->SatThickness, FineMap->DeltaDepth, FineMap->Probability,
            FineMap->MassWasting, FineMap->MassDeposition, FineMap->SedimentToChannel,
            SedMap->Erosion, SedMap->SedFluxOut, SedimentOverlandInflow,
            Network->Erosion, SedMap->RoadSed, SedimentOverroadInflow);

}
