/*
 * SUMMARY:      CalcAvailableWater.c - Calculate moisture in soil column
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Mark Wigmosta (*)
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculates the amount of soil moisture available for
 *               saturated flow below the groundwater table.
 * DESCRIP-END.
 * FUNCTIONS:    CalcAvailableWater()
 * COMMENTS:     Mark Wigmosta, Batelle Pacific Northwest Laboratories,
 *               ms_wigmosta@pnl.gov
 * $Id: CalcAvailableWater.c,v 1.4 2003/07/01 21:26:09 olivier Exp $      
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "soilmoisture.h"

/****************************************************************************
  Function name: CalcAvailableWater()

  Purpose      : This routine calculates the amount of soil moisture
                 available for saturated flow below the groundwater table.
                 No flow is allowed to occur if the moisture content falls
                 below the field capacity.

  Required     :
    int NRootLayers  - Number of soil layers
    float TotalDepth - Total depth of the soil profile (m)
    float *RootDepth - Depth of each of the soil layers (m)
    float *Porosity  - Porosity of each soil layer
    float *FCap      - Field capacity of each soil layer
    float TableDepth - Depth of the water table below the surface (m)
    float *Adjust    - Correction for each layer for loss of soil storage
                       due to channel/road-cut.  Multiplied with RootDepth
                       to give the layer thickness for use in calculating
                       soil moisture

  Returns      :
    float AvailableWater - Total water available for saturated flow

  Modifies     : void

  Comments     :
*****************************************************************************/
float CalcAvailableWater(int NRootLayers, float TotalDepth, float *RootDepth,
			 float *Porosity, float *FCap, float TableDepth,
			 float *Adjust)
{
  float AvailableWater;		/* amount of water available for movement
				   (m) */
  float DeepFCap;		/* field capacity of the layer below the
				   deepest root layer */
  float DeepLayerDepth;		/* depth of layer below deepest root zone
				   layer */
  float DeepPorosity;		/* porosity of the layer below
				   the deepest root layer */
  float Depth;			/* depth below the ground surface (m) */
  int i;			/* counter */

  AvailableWater = 0;

  Depth = 0.0;
  for (i = 0; i < NRootLayers && Depth < TotalDepth; i++) {
    if (RootDepth[i] < (TotalDepth - Depth))
      Depth += RootDepth[i];
    else
      Depth = TotalDepth;
    if (Depth > TableDepth) {
      if ((Depth - TableDepth) > RootDepth[i])
	AvailableWater += (Porosity[i] - FCap[i]) * RootDepth[i] * Adjust[i];
      else
	AvailableWater += (Porosity[i] - FCap[i]) * (Depth - TableDepth) *
	  Adjust[i];
    }
  }

  if (Depth < TotalDepth) {

    DeepPorosity = Porosity[NRootLayers - 1];
    DeepFCap = FCap[NRootLayers - 1];

    DeepLayerDepth = TotalDepth - Depth;
    Depth = TotalDepth;

    if ((Depth - TableDepth) > DeepLayerDepth)
      AvailableWater += (DeepPorosity - DeepFCap) * DeepLayerDepth *
	Adjust[NRootLayers];
    else
      AvailableWater += (DeepPorosity - DeepFCap) * (Depth - TableDepth) *
	Adjust[NRootLayers];
  }

  assert(AvailableWater >= 0.0);
  return AvailableWater;
}
