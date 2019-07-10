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
 * $Id: Aggregate.c,v 1.17 2018/02/18 ning Exp $
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
#include "timing.h"

/* This macro is defined so several different operations can be
   applied to each field, including summation over processors, without
   repeating essentially the same code for each operation. Another
   macro called "MACRO" must be #define'd prior to using this and
   should be #undef'd afterward. If a new field is added to
   AGGREGATED, it would have to appear in cell-by-cell summation, in
   Aggregate(), and here in this macro. */

#define APPLY_MACRO(aggregate) \
  { \
    int i; \
  MACRO(aggregate->Evap.ETot); \
  for (i = 0; i < Veg->MaxLayers + 1; i++) { \
    MACRO(aggregate->Evap.EPot[i]); \
    MACRO(aggregate->Evap.EAct[i]); \
  } \
  for (i = 0; i < Veg->MaxLayers; i++) { \
    MACRO(aggregate->Evap.EInt[i]); \
  } \
  for (i = 0; i < Veg->MaxLayers; i++) { \
    for (j = 0; j < Soil->MaxLayers; j++) { \
      MACRO(aggregate->Evap.ESoil[i][j]); \
    } \
  } \
  MACRO(aggregate->Evap.EvapSoil); \
  MACRO(aggregate->Precip.Precip); \
  MACRO(aggregate->Precip.SnowFall); \
  for (i = 0; i < Veg->MaxLayers; i++) { \
    MACRO(aggregate->Precip.IntRain[i]); \
    MACRO(aggregate->Precip.IntSnow[i]); \
  } \
  MACRO(aggregate->CanopyWater); \
  MACRO(aggregate->Rad.Tair); \
  MACRO(aggregate->Rad.ObsShortIn); \
  MACRO(aggregate->Rad.PixelNetShort); \
  MACRO(aggregate->NetRad); \
  MACRO(aggregate->Rad.BeamIn); \
  MACRO(aggregate->Rad.DiffuseIn); \
  for (i = 0; i < 2; i++) { \
    MACRO(aggregate->Rad.NetShort[i]); \
    MACRO(aggregate->Rad.LongIn[i]); \
    MACRO(aggregate->Rad.LongOut[i]); \
  } \
  MACRO(aggregate->Snow.Swq); \
  MACRO(aggregate->Snow.Melt); \
  MACRO(aggregate->Snow.PackWater); \
  MACRO(aggregate->Snow.TPack); \
  MACRO(aggregate->Snow.SurfWater); \
  MACRO(aggregate->Snow.TSurf); \
  MACRO(aggregate->Snow.ColdContent); \
  MACRO(aggregate->Snow.Albedo); \
  MACRO(aggregate->Snow.Depth); \
  MACRO(aggregate->Snow.Qe); \
  MACRO(aggregate->Snow.Qs); \
  MACRO(aggregate->Snow.Qsw); \
  MACRO(aggregate->Snow.Qlw); \
  MACRO(aggregate->Snow.Qp); \
  MACRO(aggregate->Snow.MeltEnergy); \
  MACRO(aggregate->Snow.VaporMassFlux); \
  MACRO(aggregate->Snow.CanopyVaporMassFlux); \
  MACRO(aggregate->Soil.Moist[Soil->MaxLayers]); \
  MACRO(aggregate->Soil.TableDepth); \
  MACRO(aggregate->Soil.WaterLevel); \
  MACRO(aggregate->Soil.SatFlow); \
  MACRO(aggregate->Soil.TSurf); \
  MACRO(aggregate->Soil.Qnet); \
  MACRO(aggregate->Soil.Qs); \
  MACRO(aggregate->Soil.Qe); \
  MACRO(aggregate->Soil.Qg); \
  MACRO(aggregate->Soil.Qst); \
  MACRO(aggregate->Soil.IExcess); \
  MACRO(aggregate->Soil.DetentionStorage); \
  MACRO(aggregate->Road.IExcess); \
  if (Options->Infiltration == DYNAMIC) { \
    MACRO(aggregate->Soil.InfiltAcc); \
  } \
  MACRO(aggregate->SoilWater); \
  MACRO(aggregate->Soil.Runoff); \
  MACRO(aggregate->ChannelInt); \
  MACRO(aggregate->RoadInt); \
  MACRO(aggregate->CulvertReturnFlow); \
  MACRO(aggregate->CulvertToChannel); \
  if (TotNumGap > 0) { \
    MACRO(aggregate->Veg.Type[Opening].Qsw); \
    MACRO(aggregate->Veg.Type[Opening].Qlin); \
    MACRO(aggregate->Veg.Type[Opening].Qlw); \
    MACRO(aggregate->Veg.Type[Opening].Qe); \
    MACRO(aggregate->Veg.Type[Opening].Qs); \
    MACRO(aggregate->Veg.Type[Opening].Qp); \
    MACRO(aggregate->Veg.Type[Opening].Swq); \
    MACRO(aggregate->Veg.Type[Opening].MeltEnergy); \
  } \
  MACRO(aggregate->Soil.Depth); \
  for (i = 0; i < Soil->MaxLayers; i++) { \
    MACRO(aggregate->Soil.Moist[i]); \
    MACRO(aggregate->Soil.Perc[i]); \
    MACRO(aggregate->Soil.Temp[i]); \
  } \
  } 


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
	       LAYER *Soil, LAYER *Veg, VEGPIX **VegMap, EVAPPIX **Evap,
	       PRECIPPIX **Precip, PIXRAD **RadMap, SNOWPIX **Snow,
	       SOILPIX **SoilMap, AGGREGATED *Total, VEGTABLE *VType,
	       ROADSTRUCT **Network, CHANNEL *ChannelData, float *roadarea,
               int Dt)
{
  int NPixels;                  /* Number of pixels in the basin */
  int NSoilL;                   /* Number of soil layers for current pixel */
  int NVegL;                    /* Number of vegetation layers for current pixel */
  int i;                                /* counter */
  int j;                                /* counter */
  int x;
  int y;
  int n;
  float DeepDepth;              /* depth to bottom of lowest rooting zone */
  float *sums;

  TIMING_TASK_START("Aggregate", 2);

  NPixels = Map->AllCells;
  *roadarea = 0.;

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
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
        Total->Snow.Qe += Snow[y][x].Qe;
        Total->Snow.Qs += Snow[y][x].Qs;
        Total->Snow.Qsw += Snow[y][x].Qsw;
        Total->Snow.Qlw += Snow[y][x].Qlw;
        Total->Snow.Qp += Snow[y][x].Qp;
        Total->Snow.MeltEnergy += Snow[y][x].MeltEnergy;
        Total->Snow.VaporMassFlux += Snow[y][x].VaporMassFlux;
        Total->Snow.CanopyVaporMassFlux += Snow[y][x].CanopyVaporMassFlux;

        if (VegMap[y][x].Gapping) {
          Total->Veg.Type[Opening].Qsw += VegMap[y][x].Type[Opening].Qsw;
          Total->Veg.Type[Opening].Qlin += VegMap[y][x].Type[Opening].Qlin;
          Total->Veg.Type[Opening].Qlw += VegMap[y][x].Type[Opening].Qlw;
          Total->Veg.Type[Opening].Qe += VegMap[y][x].Type[Opening].Qe;
          Total->Veg.Type[Opening].Qs += VegMap[y][x].Type[Opening].Qs;
          Total->Veg.Type[Opening].Qp += VegMap[y][x].Type[Opening].Qp;
          Total->Veg.Type[Opening].Swq += VegMap[y][x].Type[Opening].Swq;
          Total->Veg.Type[Opening].MeltEnergy += VegMap[y][x].Type[Opening].MeltEnergy;
        }

        if (VegMap[y][x].Gapping > 0.0 ) {
          Total->Veg.Type[Opening].Qsw += VegMap[y][x].Type[Opening].Qsw;
          Total->Veg.Type[Opening].Qlin += VegMap[y][x].Type[Opening].Qlin;
          Total->Veg.Type[Opening].Qlw += VegMap[y][x].Type[Opening].Qlw;
          Total->Veg.Type[Opening].Qe += VegMap[y][x].Type[Opening].Qe;
          Total->Veg.Type[Opening].Qs += VegMap[y][x].Type[Opening].Qs;
          Total->Veg.Type[Opening].Qp += VegMap[y][x].Type[Opening].Qp;
          Total->Veg.Type[Opening].Swq += VegMap[y][x].Type[Opening].Swq;
          Total->Veg.Type[Opening].MeltEnergy += VegMap[y][x].Type[Opening].MeltEnergy;
        }
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

  /* At this point each Total field needs to be summed over processors
     (all-reduce).  This is a very expensive communication
     operation. In order to reduce the number of messages, it's
     necessary to do the all-reduce for all fields at once. */

                                /* count the number of slots we need (roadarea is 0) */
  n = 1;

#define MACRO(x) ++n;
  APPLY_MACRO(Total);
#undef MACRO


                                /* make an array to store the
                                   aggregate fields from this
                                   processor */

  if (!(sums = (float *) calloc(n, sizeof(float)))) {
    ReportError("Aggregate", 1);
  }

  n = 0;
  sums[n++] = *roadarea;

                                /* sum aggregated fields across processors */

#define MACRO(x) sums[n++] = x;
  APPLY_MACRO(Total);
#undef MACRO

  GA_Fgop(sums, n, "+");

                                /* put aggregated fields back in Total */
  
  n = 0;
  *roadarea = sums[n++] / (Map->DX * Map->DY * NPixels);

#define MACRO(x) x = sums[n++]; x /= (float)NPixels;
  APPLY_MACRO(Total);
#undef MACRO
  

  /* These are averaged differently than other aggregate fields */
  if (TotNumGap > 0) {
    Total->Veg.Type[Opening].Qsw /= TotNumGap/NPixels;
    Total->Veg.Type[Opening].Qlin /= TotNumGap/NPixels;
    Total->Veg.Type[Opening].Qlw /= TotNumGap/NPixels;
    Total->Veg.Type[Opening].Qe /= TotNumGap/NPixels;
    Total->Veg.Type[Opening].Qs /= TotNumGap/NPixels;
    Total->Veg.Type[Opening].Qp /= TotNumGap/NPixels;
    Total->Veg.Type[Opening].Swq /= TotNumGap/NPixels;
    Total->Veg.Type[Opening].MeltEnergy /= TotNumGap/NPixels;
  }
  
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    /* convert EPot from m/s to m */
    Total->Evap.EPot[i] *= Dt;
  }


  free(sums);
  TIMING_TASK_END("Aggregate", 2);
}
