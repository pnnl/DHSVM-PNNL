/*
 * SUMMARY:      Aggregate.c - calculate basin-wide values
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the average values for the different fluxes and
 *               state variables over the basin.
 * DESCRIP-END.
 * FUNCTIONS:    Aggregate()
 * COMMENTS:
 * $Id: Aggregate.c,v 1.17 2004/08/18 01:01:25 colleen Exp $
 */

#include <ga.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"


/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  
  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().
  
  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void Aggregate(MAPSIZE *Map, OPTIONSTRUCT *Options, TOPOPIX **TopoMap,
	       LAYER *Soil, LAYER * Veg, VEGPIX **VegMap, EVAPPIX **Evap,
	       PRECIPPIX **Precip, PIXRAD **RadMap, SNOWPIX **Snow,
	       SOILPIX **SoilMap, AGGREGATED *Total, VEGTABLE *VType,
	       ROADSTRUCT **Network, CHANNEL *ChannelData, float *roadarea)
{
  int NPixels, gNPixels;	/* Number of pixels in the basin */
  int NSoilL;			/* Number of soil layers for current pixel */
  int NVegL;			/* Number of vegetation layers for current pixel */
  int i;				/* counter */
  int j;				/* counter */
  int x;
  int y;
  float DeepDepth;		/* depth to bottom of lowest rooting zone */

  NPixels = 0;
  *roadarea = 0.;

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        NPixels++;
        NSoilL = Soil->NLayers[SoilMap[y][x].Soil - 1];
        NVegL = Veg->NLayers[VegMap[y][x].Veg - 1];
		  
        /* aggregate the evaporation data */
        Total->Evap.ETot += Evap[y][x].ETot;
        for (i = 0; i < NVegL; i++) {
          Total->Evap.EPot[i] += Evap[y][x].EPot[i];
          Total->Evap.EAct[i] += Evap[y][x].EAct[i];
          Total->Evap.EInt[i] += Evap[y][x].EInt[i];
        }
        Total->Evap.EPot[Veg->MaxLayers] += Evap[y][x].EPot[NVegL];
        Total->Evap.EAct[Veg->MaxLayers] += Evap[y][x].EAct[NVegL];
		  
        for (i = 0; i < NVegL; i++) {
          for (j = 0; j < NSoilL; j++) {
            Total->Evap.ESoil[i][j] += Evap[y][x].ESoil[i][j];
          }
        }
        Total->Evap.EvapSoil += Evap[y][x].EvapSoil;
		  
        /* aggregate precipitation data */
        Total->Precip.Precip += Precip[y][x].Precip;
        Total->Precip.SnowFall += Precip[y][x].SnowFall;
        for (i = 0; i < NVegL; i++) {
          Total->Precip.IntRain[i] += Precip[y][x].IntRain[i];
          Total->Precip.IntSnow[i] += Precip[y][x].IntSnow[i];
          Total->CanopyWater += Precip[y][x].IntRain[i] +
            Precip[y][x].IntSnow[i];
        }

	/* aggregate radiation data */
	if (Options->MM5 == TRUE) {
	  Total->Rad.BeamIn = NOT_APPLICABLE;
	  Total->Rad.DiffuseIn = NOT_APPLICABLE;
	}
	else {
          Total->Rad.Tair += RadMap[y][x].Tair;
          Total->Rad.ObsShortIn += RadMap[y][x].ObsShortIn;
	  Total->Rad.BeamIn += RadMap[y][x].BeamIn;
	  Total->Rad.DiffuseIn += RadMap[y][x].DiffuseIn;
          Total->Rad.PixelNetShort += RadMap[y][x].PixelNetShort;
          Total->NetRad += RadMap[y][x].NetRadiation[0] + RadMap[y][x].NetRadiation[1];
	}

	/* aggregate snow data */
	if (Snow[y][x].HasSnow)
          Total->Snow.HasSnow = TRUE;
        Total->Snow.Swq += Snow[y][x].Swq;
        Total->Snow.Glacier += Snow[y][x].Glacier;
        /* Total->Snow.Melt += Snow[y][x].Melt; */
        Total->Snow.Melt += Snow[y][x].Outflow;
        Total->Snow.PackWater += Snow[y][x].PackWater;
        Total->Snow.TPack += Snow[y][x].TPack;
        Total->Snow.SurfWater += Snow[y][x].SurfWater;
        Total->Snow.TSurf += Snow[y][x].TSurf;
        Total->Snow.ColdContent += Snow[y][x].ColdContent;
        Total->Snow.Albedo += Snow[y][x].Albedo;
        Total->Snow.Depth += Snow[y][x].Depth;
        Total->Snow.VaporMassFlux += Snow[y][x].VaporMassFlux;
        Total->Snow.CanopyVaporMassFlux += Snow[y][x].CanopyVaporMassFlux;

        /* aggregate soil moisture data */
        Total->Soil.Depth += SoilMap[y][x].Depth;
        DeepDepth = 0.0;

        for (i = 0; i < NSoilL; i++) {
          Total->Soil.Moist[i] += SoilMap[y][x].Moist[i];
          assert(SoilMap[y][x].Moist[i] >= 0.0);
          Total->Soil.Perc[i] += SoilMap[y][x].Perc[i];
          Total->Soil.Temp[i] += SoilMap[y][x].Temp[i];
          Total->SoilWater += SoilMap[y][x].Moist[i] * VType[VegMap[y][x].Veg - 1].RootDepth[i] * Network[y][x].Adjust[i]; 
          DeepDepth += VType[VegMap[y][x].Veg - 1].RootDepth[i];
        }

        Total->Soil.Moist[Soil->MaxLayers] += SoilMap[y][x].Moist[NSoilL];
        Total->SoilWater += SoilMap[y][x].Moist[NSoilL] * (SoilMap[y][x].Depth - DeepDepth) * Network[y][x].Adjust[NSoilL];
        Total->Soil.TableDepth += SoilMap[y][x].TableDepth;

        if (SoilMap[y][x].TableDepth <= 0)
          (Total->Saturated)++;
		
        Total->Soil.WaterLevel += SoilMap[y][x].WaterLevel;
        Total->Soil.SatFlow += SoilMap[y][x].SatFlow;
        Total->Soil.TSurf += SoilMap[y][x].TSurf;
        Total->Soil.Qnet += SoilMap[y][x].Qnet;
        Total->Soil.Qs += SoilMap[y][x].Qs;
        Total->Soil.Qe += SoilMap[y][x].Qe;
        Total->Soil.Qg += SoilMap[y][x].Qg;
        Total->Soil.Qst += SoilMap[y][x].Qst;
        Total->Soil.IExcess += SoilMap[y][x].IExcess;
        Total->Soil.DetentionStorage += SoilMap[y][x].DetentionStorage;
		
        if (Options->Infiltration == DYNAMIC)
          Total->Soil.InfiltAcc += SoilMap[y][x].InfiltAcc;
		
        Total->Soil.Runoff += SoilMap[y][x].Runoff;
        Total->ChannelInt += SoilMap[y][x].ChannelInt;
        SoilMap[y][x].ChannelInt = 0.0;
        Total->RoadInt += SoilMap[y][x].RoadInt;
        SoilMap[y][x].RoadInt = 0.0;
      }
    }
  }

  /* Aggreation over processors is a series of all-reduce
     operations. This could get slow. The message sizes are really
     small */

  GA_Igop(&NPixels, 1, "+");

  /* divide road area by pixel area so it can be used to calculate depths
     over the road surface in FinalMassBalancs */
  *roadarea /= Map->DX * Map->DY * NPixels;
  GA_Fgop(roadarea, 1, "+");

  /* calculate average values for all quantities except the surface flow */

#define ALL_REDUCE_AVG_FLOAT(x, npix)           \
  x /= (float)npix;  GA_Fgop(&(x), 1, "+");

  /* average evaporation data */
  ALL_REDUCE_AVG_FLOAT(Total->Evap.ETot, NPixels);
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    ALL_REDUCE_AVG_FLOAT(Total->Evap.EPot[i], NPixels);
    ALL_REDUCE_AVG_FLOAT(Total->Evap.EAct[i], NPixels);
  }
  for (i = 0; i < Veg->MaxLayers; i++) {
    ALL_REDUCE_AVG_FLOAT(Total->Evap.EInt[i], NPixels);
  }
  for (i = 0; i < Veg->MaxLayers; i++) {
    for (j = 0; j < Soil->MaxLayers; j++) {
      ALL_REDUCE_AVG_FLOAT(Total->Evap.ESoil[i][j], NPixels);
    }
  }
  ALL_REDUCE_AVG_FLOAT(Total->Evap.EvapSoil, NPixels);

  /* average precipitation data */
  ALL_REDUCE_AVG_FLOAT(Total->Precip.Precip, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Precip.SnowFall, NPixels);
  for (i = 0; i < Veg->MaxLayers; i++) {
    ALL_REDUCE_AVG_FLOAT(Total->Precip.IntRain[i], NPixels);
    ALL_REDUCE_AVG_FLOAT(Total->Precip.IntSnow[i], NPixels);
  }
  ALL_REDUCE_AVG_FLOAT(Total->CanopyWater, NPixels);

  /* average radiation data */
  ALL_REDUCE_AVG_FLOAT(Total->Rad.Tair, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Rad.ObsShortIn, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Rad.PixelNetShort, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->NetRad, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Rad.BeamIn, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Rad.DiffuseIn, NPixels);
  for (i = 0; i <= 2; i++) {
    ALL_REDUCE_AVG_FLOAT(Total->Rad.NetShort[i], NPixels);
  }

  /* average snow data */
  ALL_REDUCE_AVG_FLOAT(Total->Snow.Swq, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.Melt, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.PackWater, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.TPack, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.SurfWater, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.TSurf, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.ColdContent, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.Albedo, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.Depth, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.VaporMassFlux, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Snow.CanopyVaporMassFlux, NPixels);

  /* average soil moisture data */
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Depth, NPixels);
  for (i = 0; i < Soil->MaxLayers; i++) {
    ALL_REDUCE_AVG_FLOAT(Total->Soil.Moist[i], NPixels);
    ALL_REDUCE_AVG_FLOAT(Total->Soil.Perc[i], NPixels);
    ALL_REDUCE_AVG_FLOAT(Total->Soil.Temp[i], NPixels);
  }
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Moist[Soil->MaxLayers], NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.TableDepth, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.WaterLevel, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.SatFlow, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.TSurf, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Qnet, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Qs, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Qe, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Qg, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Qst, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.IExcess, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.DetentionStorage, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Road.IExcess, NPixels);
  
  if (Options->Infiltration == DYNAMIC) {
    ALL_REDUCE_AVG_FLOAT(Total->Soil.InfiltAcc, NPixels);
  }

  ALL_REDUCE_AVG_FLOAT(Total->SoilWater, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->Soil.Runoff, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->ChannelInt, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->RoadInt, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->CulvertReturnFlow, NPixels);
  ALL_REDUCE_AVG_FLOAT(Total->CulvertToChannel, NPixels);

#undef ALL_REDUCE_AVG_FLOAT
}
