/*
 * SUMMARY:      InitAggregated.c - Initialize basin-wide values
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize basin-wide values
 * DESCRIP-END.
 * FUNCTIONS:    InitAggregated()
 * COMMENTS:
 * $Id: InitAggregated.c,v 1.4 2003/07/01 21:26:15 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  InitAggregated()

  Allocates memory for the structure that will hold basin total and/or basin 
  average values
*****************************************************************************/
void InitAggregated(int MaxVegLayers, int MaxSoilLayers, AGGREGATED * Total)
{
  int i;			/* counter */

  if (!(Total->Evap.EPot = (float *) calloc(MaxVegLayers + 1, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Evap.EAct = (float *) calloc(MaxVegLayers + 1, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Evap.EInt = (float *) calloc(MaxVegLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Evap.ESoil = (float **) calloc(MaxVegLayers, sizeof(float *))))
    ReportError("InitAggregated()", 1);

  for (i = 0; i < MaxVegLayers; i++) {
    if (!(Total->Evap.ESoil[i] =
	  (float *) calloc(MaxSoilLayers, sizeof(float))))
      ReportError("InitAggregated()", 1);
  }

  if (!(Total->Precip.IntRain = (float *) calloc(MaxVegLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Precip.IntSnow = (float *) calloc(MaxVegLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Soil.Moist = (float *) calloc(MaxSoilLayers + 1, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Soil.Perc = (float *) calloc(MaxSoilLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Soil.Temp = (float *) calloc(MaxSoilLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);
}
