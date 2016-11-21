/*
 * SUMMARY:      CalcTotalWater.c - Calculate total water in soil column
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Batelle Pacific Northwest Laboratories
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Jul-96
 * DESCRIPTION:  Calculates the total amount of moisture in the soil column 
 * DESCRIP-END.
 * FUNCTIONS:    CalcTotalWater()
 * COMMENTS:
 * $Id: CalcTotalWater.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "soilmoisture.h"

/*****************************************************************************
  Function name: CalcTotalWater()

  Purpose      : Calculate the total amount of water in the soil column,
                 corrected for the effects of road cut-banks and channels in
                 grid cells
                 
  Required     :
    int NSoilLayers  - Number of soil layers
    float TotalDepth - Total thickness of the soil column (m)
    float *RootDepth - Array with thicknesses of root layers (m)
    float *Moist     - Array with soil moisture in each layer
    float *Adjust    - Array with coefficients to correct for loss of soil
                       storage due to channel/road-cut for each soil layer.
                       Multiplied with RootDepth to give the zone thickness
                       for use in calculating soil moisture 

  Returns      : 
    float TotalWater - Total amount of water in the soil profile (m)

  Modifies     : void

  Comments     :
*****************************************************************************/
float CalcTotalWater(int NSoilLayers, float TotalDepth, float *RootDepth,
		     float *Moist, float *Adjust)
{
  float TotalWater;		/* Total water in the soil column (m) */
  float DeepLayerDepth;		/* depth of layer below deepest root zone
				   layer */
  float Depth;			/* depth below the ground surface (m) */
  int i;			/* counter */

  TotalWater = 0.0;
  Depth = 0.0;

  for (i = 0; i < NSoilLayers; i++) {
    Depth += RootDepth[i];
    TotalWater += Moist[i] * RootDepth[i] * Adjust[i];
  }

  DeepLayerDepth = TotalDepth - Depth;

  TotalWater += Moist[NSoilLayers] * DeepLayerDepth * Adjust[NSoilLayers];

  return TotalWater;
}
