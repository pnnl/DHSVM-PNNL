/*
 * SUMMARY:      SoilEvaporation.c - Calculate evaporation from the soil 
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculate evaporation from the soil 
 * DESCRIP-END.
 * FUNCTIONS:    SoilEvaporation()
 * COMMENTS:
 * $Id: SoilEvaporation.c,v 1.4 2003/07/01 21:26:25 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "settings.h"

/*****************************************************************************
  SoilEvaporation()
*****************************************************************************/
float SoilEvaporation(int Dt, float Temp, float Slope, float Gamma, float Lv,
		      float AirDens, float Vpd, float NetRad, float RaSoil,
		      float Transpiration, float Porosity, float Ks,
		      float Press, float m, float RootDepth,
		      float *MoistContent, float Adjust)
{
  float DesorptionVolume;	/* Amount of water the soil can deliver to the
				   atmosphere during a timestep (mm) */
  float EPot;			/* Potential evaporation from soil during
				   timestep (mm) */
  float SoilEvap;		/* Amount of evaporation directly from the soil
				   (mm) */
  float SoilMoisture;		/* Amount of water in surface soil layer 
				   (mm) */

  DesorptionVolume = Desorption(Dt, *MoistContent, Porosity, Ks, Press, m);

  /* Eq.4 Wigmosta et al [1994] */

  /* Calculate the density of pure water as a function of temperature.
     Thiesen, Scheel-Diesselhorst Equation (in Handbook of hydrology, fig
     11.1.1) */

  EPot = (Slope * NetRad + AirDens * CP * Vpd / RaSoil) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;

  /* The potential evaporation rate accounts for the amount of moisture that
     the atmosphere can absorb.  If we do not account for the amount of
     evaporation from overlying evaporation, we can end up with the situation
     that all vegetation layers and the soil layer transpire/evaporate at the
     potential rate, resulting in an overprediction of the actual evaporation
     rate.  Thus we subtract the amount of evaporation that has already
     been calculated for overlying layers from the potential evaporation.
     Another mechanism that could be used to account for this would be to 
     decrease the vapor pressure deficit while going down through the canopy
     (not implemented here) */

  EPot -= Transpiration;
  if (EPot < 0.0)
    EPot = 0.0;

  /* Eq.8 Wigmosta et al [1994] */

  SoilEvap = MIN(EPot, DesorptionVolume);
  SoilEvap *= Adjust;

  SoilMoisture = *MoistContent * RootDepth * Adjust;

  if (SoilEvap > SoilMoisture) {
    SoilEvap = SoilMoisture;
    *MoistContent = 0.0;
  }
  else {
    SoilMoisture -= SoilEvap;
    *MoistContent = SoilMoisture / (RootDepth * Adjust);
  }

  return SoilEvap;
}
