/*
 * SUMMARY:      WaterTableDepth.c - Calculate the depth of the water table
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculate the depth of the water table below the ground
 *               surface corrected for road and channel effects
 * DESCRIP-END.
 * FUNCTIONS:    WaterTableDepth()
 * COMMENTS:
 * $Id: WaterTableDepth.c,v 1.4 2003/07/01 21:26:26 olivier Exp $     
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "soilmoisture.h"

/*****************************************************************************
  Function name: WaterTableDepth()

  Purpose      :   This function calculates the depth of the water table
                   below the ground surface, based on the amount of soil
                   moisture in the root zone layers. 

  Required     :
    int NRootLayers  - Number of soil layers
    float TotalDepth - Total depth of the soil profile (m)
    float *RootDepth - Depth of each of the soil layers (m)
    float *Porosity  - Porosity of each soil layer
    float *FCap      - Field capacity of each soil layer
    float *Adjust    - Correction for each layer for loss of soil storage
                       due to channel/road-cut.  Multiplied with RootDepth
                       to give the layer thickness for use in calculating
                       soil moisture 

  Returns      :
    float TableDepth - Depth of the water table below the soil surface

  Modifies     :
    float *Moist - Amount of soil moisture in each layer

  Comments     :   
    Since no unsaturated flow is allowed to occur when the moisture content
    is below the field capacity, and because lateral saturated flow is the
    only mechanism by which water can disappear from the soil below the
    deepest root layer, the soil moisture content in the soil below the
    deepest root layer can never fall below field capacity.  The water
    immediately above the water table is assumed to be at field capacity. 

    Changes have been made to account for the potential loss of soil storage
    in a grid cell due to a road-cut or channel.
*****************************************************************************/
float WaterTableDepth(int NRootLayers, float TotalDepth, float *RootDepth,
		      float *Porosity, float *FCap, float *Adjust, float *Moist)
{
  float DeepFCap;		/* field capacity of the layer below the
				   deepest root layer */
  float DeepLayerDepth;		/* depth of layer below deepest root zone
				   layer */
  float DeepPorosity;		/* porosity of the layer below the deepest
				   root layer */
  float TableDepth;		/* depth of the water table (m) */
  float MoistureTransfer;	/* amount of soil moisture transferred from
				   the current layer to the layer above (m)
				 */
  int i;			/* counter */
  float TotalStorage = 0.0;
  float ExcessFCap;
  float TotalExcessFCap = 0.0;

  MoistureTransfer = 0.0;
  DeepLayerDepth = TotalDepth;
  DeepPorosity = Porosity[NRootLayers - 1];
  DeepFCap = FCap[NRootLayers - 1];
  for (i = 0; i < NRootLayers; i++)
    DeepLayerDepth -= RootDepth[i];

  /* Redistribute soil moisture.  I.e. water from supersaturated layers is
     transferred to the layer immediately above */

  if (Moist[NRootLayers] >= DeepPorosity) {
    MoistureTransfer = (Moist[NRootLayers] - DeepPorosity) *
      DeepLayerDepth * Adjust[NRootLayers];
    Moist[NRootLayers] = DeepPorosity;

    for (i = NRootLayers - 1; i >= 0; i--) {
      Moist[i] += MoistureTransfer / (RootDepth[i] * Adjust[i]);
      if (Moist[i] >= Porosity[i]) {
	MoistureTransfer = (Moist[i] - Porosity[i]) * RootDepth[i] * Adjust[i];
	Moist[i] = Porosity[i];
      }
      else {
	MoistureTransfer = 0.0;
	break;
      }
    }
  }
  else {
    MoistureTransfer = 0.0;
  }

  if (MoistureTransfer > 0) {
    /* Surface ponding occurs */
    TableDepth = -MoistureTransfer;
  }
  else {
    /* Warning added by Pascal Storck, 08/15/2000 */
    /* Based on a single bad parameter in a DHSVM input file (a third layer
       vertical hydraulic conductivity that was 10 times smaller than the layer
       above it), it was noted that DHSVM can develop what are basically perched
       water tables.  These perched water tables greatly complicate the calculation
       of the pixel water table depth because the soil below the perched table
       is not completely saturated above field capacity.
       For example, if we have three soil layers and a deep layer, all 1 meter thick,
       and we saturate the second layer from the surface, what is, or should be, the water
       table depth. Should we allow subsurface flow to occur, should we include the saturation
       of disconnected overlying layers in the calculation of the hydraulic gradient.
       At this point, just be cautious.  Using any combination of soil parameters or 
       intial water states which can cause the lower layers of the soil profile
       to drain more quickly than water can flow down through the matrix will
       result in mass balance problems.  I.e. water will be forced out of the cell
       to the downslope, this water will be taken from the deepest soil layer, which
       can cause the deep layer soil moisture to go negative. */

    TotalStorage += DeepLayerDepth * Adjust[NRootLayers] * (DeepPorosity - DeepFCap);
   
	ExcessFCap = DeepLayerDepth * Adjust[NRootLayers] * (Moist[NRootLayers] - DeepFCap);
    
	if (ExcessFCap < 0.0)
		ExcessFCap = 0.0;
	TotalExcessFCap += ExcessFCap;

    for (i = 0; i < NRootLayers; i++) {
		TotalStorage += RootDepth[i] * Adjust[i] * (Porosity[i] - FCap[i]);
		ExcessFCap = RootDepth[i] * Adjust[i] * (Moist[i] - FCap[i]);
		if (ExcessFCap < 0.0)
			ExcessFCap = 0.0;
		TotalExcessFCap += ExcessFCap;
    }
    TableDepth = TotalDepth * (1 - TotalExcessFCap / TotalStorage);
  
    if (TableDepth < 0)
		TableDepth = -(TotalExcessFCap - TotalStorage);
  }

  if (TableDepth > TotalDepth)
	  printf("TableDepth = %.2f, TotalDepth = %.2f\n", TableDepth, TotalDepth);
  if (TableDepth != TableDepth)
	  printf("TableDepth = %.2f", TableDepth);
  
  assert(TableDepth <= TotalDepth);
  return TableDepth;
}
