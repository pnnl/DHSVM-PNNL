/*
 * SUMMARY:      UnsaturatedFlow.c - Calculate the unsaturated flow
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Mark Wigmosta (*)
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculate the unsaturated flow in the soil column (vertical
 *               flow) 
 * DESCRIP-END.
 * FUNCTIONS:    UnsaturatedFlow()
 * COMMENTS: (*) Mark Wigmosta, Batelle Pacific Northwest Laboratories, 
 *               ms_wigmosta@pnl.gov
 * $Id: UnsaturatedFlow.c,v 1.12 2004/05/04 19:38:59 colleen Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "functions.h"
#include "soilmoisture.h"

/*****************************************************************************
  Function name: UnsaturatedFlow()

  Purpose      : Calculate the unsaturated flow in the soil column, and
                 adjust the moisture in each soil layer

  Required     :
    int Dt             - Time step (seconds)
    float DX           - Grid cell width (m)
    float DY           - Grid cell width (m)
    float Infiltration - Amount of infiltration entering the top of the soil
                         column (m)
    float RoadbedInfiltration - Amount of infiltration through the roadbed(m)
    float SatFlow      - Amount of saturated flow entering the soil column
                         from neighbouring pixels (m) 
    int NSoilLayers    - Number of soil layers
    float TotalDepth   - Total depth of the soil profile (m)
    float Area         - Area of channel or road surface (m)
    float *RootDepth   - Depth of each of the soil layers (m)
    float *Ks          - Vertical saturated hydraulic conductivity in each
                         soil layer (m/s)
    float *PoreDist    - Pore size distribution index for each soil layer
    float *Porosity    - Porosity of each soil layer
    float *FCap        - Field capacity of each soil layer
    float *Perc        - Amount of water percolating from each soil layer to
                         the layer below (m)
    float *PercArea    - Area of the bottom of each soil layer as a fraction
                         of the gridcell area DX*DY
    float *Adjust      - Correction for each layer for loss of soil storage
                         due to channel/road-cut.  Multiplied with RootDepth
                         to give the layer thickness for use in calculating
                         soil moisture 
    int CutBankZone    - Number of the soil layer containing the bottom of
                         the cut-bank 
    float BankHeight   - Distance from ground surface to channel bed or
                         bottom of road-cut (m) 

  Returns      : void

  Modifies     :
    float *TableDepth - Depth of the water table below the ground surface (m)
    float *RunOff     - Amount of runoff produced at the pixel (m)
    float *Moist      - Moisture content in each soil layer

  Comments     :
    Sources: 
    Bras, R. A., Hydrology, an introduction to hydrologic science, Addisson 
        Wesley, Inc., Reading, etc., 1990.
    Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
        hydrology-vegetation model for complex terrain, Water Resour. Res.,
        30(6), 1665-1679, 1994.

    This function is based on Wigmosta et al [1994], and assumes a unit 
    hydraulic gradient in the unsaturated zone.  This implies a 
    steady-state situation and uniform moisture distribution.  This is not
    unreasonable for situations where the groundwater is fairly deep 
    [Bras, 1990], but may be unrealistic for the boreal forest system, or 
    other similar situations where the groundwater table is relatively close
    to the surface. The current formulation does not allow for upward
    movement of water in the unsaturated zone, by capillary and matrix
    forces.  No unsaturated flux is assumed to occur if the water content
    drops below the field capacity.

    The unsaturated hydraulic conductivity is calculated using the
    Brooks-Corey equation see for example Wigmosta et al. [1994], eq.41

    The calculated amount of drainage is averaged with the amount calculated
    for the previous timestep, see eq. 42 [Wigmosta et al., 1994].
  
    The residual moisture content is assumed to be zero (easy to adapt if 
    needed).

    If the amount of soil moisture in a layer after drainage still 
    exceeds the porosity for that layer, this additional amount of water is
    added to the drainage to the next layer.
  
    CHANGES:
  
    Changes have been made to account for the potental loss of soil storage
    in a grid cell due to a road-cut or channel.  Correction coefficents are
    calculated in AdjustStorage() and CutBankGeometry()
*****************************************************************************/
void UnsaturatedFlow(int Dt, float DX, float DY, float Infiltration, 
		     float RoadbedInfiltration, float SatFlow, int NSoilLayers, 
		     float TotalDepth, float Area, float *RootDepth, float *Ks, 
		     float *PoreDist, float *Porosity, float *FCap, 
		     float *Perc, float *PercArea, float *Adjust, 
		     int CutBankZone, float BankHeight, float *TableDepth, 
		     float *Runoff, float *Moist, int RoadRouteOption,
		     int InfiltOption, float *RoadIExcess)
{
  float DeepDrainage;		/* amount of drainage from the lowest root 
				   zone to the layer below it (m) */
  float DeepLayerDepth;		/* depth of the layer below the deepest
				   root layer */
  float Drainage;		/* amount of water drained from each soil 
				   layer during the current timestep */
  float Exponent;		/* Brooks-Corey exponent */
  float FieldCapacity;		/* amount of water in soil at field capacity
				   (m) */
  float MaxSoilWater;		/* maximum allowable amount of soil moiture
				   in each layer (m) */
  float SoilWater;		/* amount of water in each soil layer (m) */
  int i;			/* counter */

  DeepLayerDepth = TotalDepth;
  for (i = 0; i < NSoilLayers; i++)
    DeepLayerDepth -= RootDepth[i];

  /* first take care of infiltration through the roadbed/channel, then through the
     remaining surface */
  if (*TableDepth <= BankHeight) { /* watertable above road/channel surface */
    if(RoadRouteOption)
      *RoadIExcess += RoadbedInfiltration;
    else
      *Runoff += RoadbedInfiltration;
  }
  
  else {
    if (CutBankZone == NSoilLayers) {
      Moist[NSoilLayers] += RoadbedInfiltration / 
	(DeepLayerDepth * Adjust[NSoilLayers]);
    }
    else if (CutBankZone >= 0){
      Moist[CutBankZone] += RoadbedInfiltration / 
	(RootDepth[CutBankZone] * Adjust[CutBankZone]);
    }
  }
  if (*TableDepth <= 0) { /* watertable above surface */
    *Runoff += Infiltration;
    if(InfiltOption == DYNAMIC) Infiltration = 0.; 
  }
  else {
    Moist[0] += Infiltration / (RootDepth[0] * Adjust[0]);
  }
  
  
  for (i = 0; i < NSoilLayers; i++) {
    
    /* No movement if soil moisture is below field capacity */
    if (Moist[i] > FCap[i]) {
      Exponent = 2.0 / PoreDist[i] + 3.0;
      
      if (Moist[i] > Porosity[i])
	/* this can happen because the moisture content can exceed the 
	   porosity the way the algorithm is implemented */
	Drainage = Ks[i];
      else
	Drainage =
	  Ks[i] * pow((double) (Moist[i] / Porosity[i]), (double) Exponent);
      
      Drainage *= Dt;
      Perc[i] = 0.5 * (Perc[i] + Drainage) * PercArea[i];
      
      MaxSoilWater = RootDepth[i] * Porosity[i] * Adjust[i];
      SoilWater = RootDepth[i] * Moist[i] * Adjust[i];
      FieldCapacity = RootDepth[i] * FCap[i] * Adjust[i];
      
      /* No unsaturated flow if the moisture content drops below field 
	 capacity */
      
      if ((SoilWater - Perc[i]) < FieldCapacity)
	Perc[i] = SoilWater - FieldCapacity;
      
      /* WORK IN PROGRESS */
      /* If the moisture content is greater than the porosity add the 
	 additional soil moisture to the percolation */
      
      SoilWater -= Perc[i];
      if (SoilWater > MaxSoilWater)
	Perc[i] += SoilWater - MaxSoilWater;
      
      /* Adjust the moisture content in the current layer, and the layer 
	 immediately below it */
      
      Moist[i] -= Perc[i] / (RootDepth[i] * Adjust[i]);
      if (i < (NSoilLayers - 1))
	Moist[i + 1] += Perc[i] / (RootDepth[i + 1] * Adjust[i + 1]);
    }
    else
      Perc[i] = 0.0;
    
    /* convert back to straight 1-d flux */
    Perc[i] /= PercArea[i];
  }
  
  DeepDrainage = (Perc[NSoilLayers - 1] * PercArea[NSoilLayers - 1]) + SatFlow;
  
  Moist[NSoilLayers] += DeepDrainage / (DeepLayerDepth * Adjust[NSoilLayers]);
  
  /* added 8/16/2000 by Pascal Storck */
  /* this following statement will force a trap of out of bounds 
     soil moisture in the lowest layer in the mass balance calculation */
  if (Moist[NSoilLayers] < FCap[NSoilLayers - 1]) {
    /*    Moist[NSoilLayers] = FCap[NSoilLayers - 1]; */
   /*  printf("Warning: Deep layer soil moisture is less than field capacity.\n"); */

  /*   if (Moist[NSoilLayers] < 0. ) */
/*       printf("Warning: Deep layer soil moisture is negative.\n"); */
  }
  
  /* Calculate the depth of the water table based on the soil moisture 
     profile and adjust the soil moisture profile, to assure that the soil 
     moisture is never more than the maximum allowed soil moisture amount,
     i.e. the porosity.  A negative water table depth means that the water is 
     ponding on the surface.  This amount of water becomes surface Runoff */
  
  *TableDepth = WaterTableDepth(NSoilLayers, TotalDepth, RootDepth, Porosity,
				FCap, Adjust, Moist);
  
  if (*TableDepth < 0.0) {
    *Runoff += -(*TableDepth);
    
    if(InfiltOption == DYNAMIC){
      if (Infiltration > -(*TableDepth))
	Infiltration += *TableDepth;
      else
	Infiltration = 0.;
    }
    
    *TableDepth = 0.0;
  }
}
