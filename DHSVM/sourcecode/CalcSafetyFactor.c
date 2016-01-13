/*
 * SUMMARY:      CalcSafetyFactor.c - Calculate the factor of safety
 * USAGE:        Part of MWM
 *
 * AUTHOR:       Laura Bowling and Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * DESCRIPTION:  Calculate the factor of safety
 * DESCRIP-END.
 * FUNCTIONS:    CalcSafetyFactor()
 * COMMENTS:
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "DHSVMerror.h"
#include "settings.h"
#include "constants.h"
#include "data.h"

float FindValue(STATSTABLE Stats, int iter);

/*****************************************************************************
  Function name: CalcSafetyFactor()

  Purpose      : Calculate the factor of safety for mass wasting failure.
 
                 
  Required     : 
    float Swq   - Snow water equivalent (m) 
    float Depth - Snow depth (m)


  Returns      : float, values between 0 and 1 for failure:
                        > 1 for stable
                        -0.1 for unconditionally unstable
                        and -999 for not in basin/invalid slope.

  Modifies     : none
   
  Comments     :
*****************************************************************************/

float CalcSafetyFactor(float Slope, int Soil, float SoilDepth, int Veg, 
		       SEDTABLE *SedType, VEGTABLE *VType, float M, 
		       SOILTABLE *SType, float Swq, float Depth,
		       int iter)
{
  double angle_int_frict_rad, soil_cohes_kg, slope_angle_rad, fc_soil_density;
  double root_cohes_kg;
  float FrictionAngle;                 /* soil parameter for infinite slope model (deg)*/
  float SoilCohesion;                  /* soil parameter for infinite slope model (kPa)*/
  float RootCohesion;                  /* veg parameter for infinite slope model (kPa) */ 
  float Surcharge;                     /* surcharge from snow and vegetation (kg/m2) */
  double term1;
  float loading;
  float safetyfactor;
  float SnowDensity;                   /* Density of snow (kg/m2) */ 

  if (Slope >= 0.) { 

    if(SoilDepth<=0.0) SoilDepth=0.001;

    M /= SoilDepth;
    if(M>=1.0) M=0.99;

    /* Get stochastic parameter values. */
    /* Need to check for valid soil and vegetation types ! */
    RootCohesion = FindValue(VType[Veg - 1].RootCoh, iter);
    FrictionAngle = FindValue(SedType[Soil - 1].Friction, iter);
    SoilCohesion = FindValue(SedType[Soil - 1].Cohesion, iter);
    Surcharge = FindValue(VType[Veg - 1].VegSurcharge, iter);

    /* Depth is not calculated anywhwere, so SnowDensity
       is not included in the Factor of Safety calculation */
    SnowDensity = (Swq - Depth) * WATER_DENSITY;
 /*    Surcharge += SnowDensity; */

/*     if(SnowDensity>0) printf("Surcharge w/  snow %4.3f SnowDensity %4.3f\n",  */
/* 			     Surcharge, SnowDensity);  */

    /* converting cohesion from kPa to kg/m2 and angles from degrees to radians */
    root_cohes_kg = (RootCohesion * 1000.) / G;
    soil_cohes_kg = (SoilCohesion * 1000.) / G;
    angle_int_frict_rad = RADPDEG * FrictionAngle;
    slope_angle_rad = RADPDEG * Slope;
    
    fc_soil_density = SType[Soil - 1].Dens[0] + 
      (SType[Soil - 1].FCap[0] * WATER_DENSITY);
    
    loading = (Surcharge / (WATER_DENSITY * SoilDepth)) + 
      ((M * SedType[Soil - 1].SatDensity) / WATER_DENSITY) +
      ((1 - M) * (fc_soil_density / WATER_DENSITY));
    
    /* Check to see if slope is unconditionally unstable.*/
    term1 = ((soil_cohes_kg + root_cohes_kg)/
	     ((Surcharge + SoilDepth * fc_soil_density)*
	      cos(slope_angle_rad)*cos(slope_angle_rad))) + 
      tan(angle_int_frict_rad);
		
    if(term1 <= tan(slope_angle_rad)) {
      /* Slope is unconditionally unstable for these values. */
      safetyfactor = -.1;
    }
    else {
      safetyfactor = (((2. * (soil_cohes_kg + root_cohes_kg)) / (WATER_DENSITY * SoilDepth * (sin(2. * slope_angle_rad)))) + ((loading - M) * ((tan(angle_int_frict_rad)) / (tan(slope_angle_rad))))) / loading; 
    }

  }    /* End of factor of safety calculation loop */
  else if (Slope <  0.)
    safetyfactor = -999.;
  else  
    safetyfactor = -999.;
  
  return safetyfactor;
}
