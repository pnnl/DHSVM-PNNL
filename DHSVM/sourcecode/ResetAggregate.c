/*
 * SUMMARY:      ResetAggregate.c - Reset basin-wide values to zero
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Reset basin-wide values to zero
 * DESCRIP-END.
 * FUNCTIONS:    ResetAggregate.()
 * COMMENTS:
 * $Id: ResetAggregate.c,v 1.12 2004/05/03 03:28:46 colleen Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "constants.h"

 /*****************************************************************************
   ResetAggregate()

   Set all the area averages to zero
 *****************************************************************************/
void ResetAggregate(LAYER * Soil, LAYER * Veg, AGGREGATED * Total,
  OPTIONSTRUCT *Options)
{
  int i;			/* counter */
  int j;			/* counter */

  if (DEBUG)
    printf("Resetting the aggregate values\n");

  /* initialize evaporation data */
  Total->Evap.ETot = 0.0;
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    Total->Evap.EPot[i] = 0.0;
    Total->Evap.EAct[i] = 0.0;
  }
  for (i = 0; i < Veg->MaxLayers; i++) {
    Total->Evap.EInt[i] = 0.0;
    for (j = 0; j < Soil->MaxLayers; j++)
      Total->Evap.ESoil[i][j] = 0.0;
  }
  Total->Evap.EvapSoil = 0.0;

  /* initialize precipitation data */
  Total->Precip.Precip = 0.0;
  Total->Precip.SnowFall = 0.0;
  for (i = 0; i < Veg->MaxLayers; i++) {
    Total->Precip.IntRain[i] = 0.0;
    Total->Precip.IntSnow[i] = 0.0;
  }

  /* initialize radiation data */
  Total->Rad.Tair = 0.;
  Total->Rad.PixelNetShort = 0.0;
  Total->Rad.BeamIn = 0.;
  Total->Rad.DiffuseIn = 0.;
  Total->NetRad = 0.;
  Total->Rad.ObsShortIn = 0.;
  for (i = 0; i < 2; i++) {
    Total->Rad.NetShort[i] = 0.;
    Total->Rad.LongIn[i] = 0.;
  }

  /* initialize snow data */
  Total->Snow.HasSnow = FALSE;
  Total->Snow.SnowCoverOver = FALSE;
  Total->Snow.LastSnow = 0;
  Total->Snow.Swq = 0.0;
  Total->Snow.Iwq = 0.0;
  Total->Snow.GlMelt = 0.0;
  Total->Snow.IceRemoved = 0.0;
  Total->Snow.Melt = 0.0;
  Total->Snow.PackWater = 0.0;
  Total->Snow.TPack = 0.0;
  Total->Snow.SurfWater = 0.0;
  Total->Snow.TSurf = 0.0;
  Total->Snow.ColdContent = 0.0;
  Total->Snow.Albedo = 0.0;
  Total->Snow.Depth = 0.0;
  Total->Snow.Qe = 0.0;
  Total->Snow.Qs = 0.0;
  Total->Snow.Qsw = 0.0;
  Total->Snow.Qlw = 0.0;
  Total->Snow.Qp = 0.0;
  Total->Snow.MeltEnergy = 0.0;
  Total->Snow.VaporMassFlux = 0.0;
  Total->Snow.CanopyVaporMassFlux = 0.0;

  /* initialize soil moisture data.  The total amount of runoff is calculated
     in the RouteSurface() routine */
  Total->Soil.Soil = 0;
  Total->Soil.Depth = 0.0;
  for (i = 0; i < Soil->MaxLayers + 1; i++)
    Total->Soil.Moist[i] = 0.0;
  for (i = 0; i < Soil->MaxLayers; i++) {
    Total->Soil.Perc[i] = 0.0;
    Total->Soil.Temp[i] = 0.0;
  }
  Total->Soil.TableDepth = 0.0;
  Total->Soil.WaterLevel = 0.0;
  Total->Soil.SatFlow = 0.0;
  Total->Soil.TSurf = 0.0;
  Total->Soil.Qnet = 0.0;
  Total->Soil.Qs = 0.0;
  Total->Soil.Qe = 0.0;
  Total->Soil.Qg = 0.0;
  Total->Soil.Qst = 0.0;
  Total->Soil.IExcess = 0.0;
  Total->Road.IExcess = 0.0;
  Total->Soil.DetentionStorage = 0.0;

  if (Options->Infiltration == DYNAMIC)
    Total->Soil.InfiltAcc = 0.0;
  Total->SoilWater = 0.0;
  Total->CanopyWater = 0.0;
  Total->Runoff = 0.0;
  Total->ChannelInt = 0.0;
  Total->RoadInt = 0.0;
  Total->Saturated = 0;
  Total->CulvertReturnFlow = 0;
  Total->CulvertToChannel = 0;

  if (Options->CanopyGapping && TotNumGap > 0) {
	for (i = 0; i < CELL_PARTITION; i++) {
	  Total->Veg.Type[i].Qsw = 0;
	  Total->Veg.Type[i].Qlin = 0;
	  Total->Veg.Type[i].Qlw = 0;
	  Total->Veg.Type[i].Qe = 0;
	  Total->Veg.Type[i].Qs = 0;
	  Total->Veg.Type[i].Qp = 0;
	  Total->Veg.Type[i].Swq = 0;
	  Total->Veg.Type[i].MeltEnergy = 0;
	}
  }
}
