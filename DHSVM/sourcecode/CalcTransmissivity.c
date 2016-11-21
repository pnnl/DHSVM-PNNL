/*
 * SUMMARY:      CalcTransmissivity.c - Calculate saturated conductivity
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculates the transmissivity through the saturated part of
 *               the soil profile
 * DESCRIP-END.
 * FUNCTIONS:    CalcTransmissivity()
 * COMMENTS: Modified by Ted Bohn on 10/1/2013
             Implemented 2-part transmissivity v depth function. бн
             Introduced a new parameter, DEPTH_THRESHOLD.  When water table depth 
			 is shallower than this threshold, transmissivity decays exponentially 
			 with depth as before.  When water table depth is deeper than this threshold, 
			 transmissivity decays linearly until it reaches 0.  DEPTH_THRESHOLD must be 
			 specified in the config file. 

 * $Id: CalcTransmissivity.c,v 1.4 2003/07/01 21:26:11 olivier Exp $   
 * $Id: CalcTransmissivity.c,v 3.1.1 2013/10/01 13:12 Ted Bohn Exp $   
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "functions.h"

/*****************************************************************************
  Function name: CalcTransmissivity()

  Purpose      : Calculates the transmissivity through the saturated part of
                 the soil profile
                 
  Required     : 
    float SoilDepth  - Total soil depth in m
    float WaterTable - Depth of the water table below the soil surface in m
    float LateralKs  - Lateral hydraulic conductivity in m/s
    float KsExponent - Exponent that describes exponential decay of LateralKs
                       with depth below the soil surface

  Returns      : Transmissivity in m2/s

  Modifies     : NA

  Comments     :
    Source:
    Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
      hydrology-vegetation model for complex terrain, Water Resour. Res.,
      30(6), 1665-1679, 1994.

    Based on:
    Beven, K. J., On subsurface stormflow:  An analysis of response times,
      Hydrolog. Sci. J., 4, 505-521, 1982.

    The hydraulic conductivity is assumed exponentially with depth, based on
    material in Beven [1982].
*****************************************************************************/
float CalcTransmissivity(float SoilDepth, float WaterTable, float LateralKs,
			 float KsExponent, float DepthThresh)
{
  float Transmissivity;		/* Transmissivity (m^2/s) */
  float TransThresh;

  if (fequal(KsExponent, 0.0))
    Transmissivity = LateralKs * (SoilDepth - WaterTable);
  else {
	/* a smaller value of WaterTable variables indicates a higher actual water table depth */
	if (WaterTable < DepthThresh) {
	  Transmissivity = (LateralKs / KsExponent) * (exp(-KsExponent * WaterTable) - exp(-KsExponent * SoilDepth));
	}
    else  {
	  TransThresh = (LateralKs / KsExponent) * (exp(-KsExponent * DepthThresh) - exp(-KsExponent * SoilDepth));
	  if(SoilDepth < DepthThresh) {
		printf("Warning: Soil DepthThreshold (%.2f) > the soil depth (%.2f)!\n", DepthThresh, SoilDepth);
		printf("Transmissivity is set to zero!");
	  }
	  Transmissivity = (SoilDepth-WaterTable)/(SoilDepth-DepthThresh)*TransThresh;
	}
  }

  return Transmissivity;
}
