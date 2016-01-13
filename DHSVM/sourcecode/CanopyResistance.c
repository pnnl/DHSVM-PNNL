/*
 * SUMMARY:      CanopyResistance.c - Calculate canopy resistance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate canopy resistance
 * DESCRIP-END.
 * FUNCTIONS:    CanopyResistance()
 * COMMENTS:
 * $Id: CanopyResistance.c,v 1.4 2003/07/01 21:26:11 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  CanopyResistance()
*****************************************************************************/
float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
		       float VpdThres, float MoistThres, float WP,
		       float TSoil, float SoilMoisture, float Vpd, float Rp)
{
  float MoistFactor;		/* multiplier for resistance due to soil
				   moisture feed-back */
  float Resistance;		/* Canopy resistance (s/m) */
  float RpFactor;		/* multiplier for resistance due to light
				   level feed-back */
  float TFactor;		/* multiplier for resistance due to soil 
				   temperaure feed-back */
  float VpdFactor;		/* multiplier for resistance due to vapor
				   pressure deficit feed-back */

  if (TSoil <= 0) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }

  /* for OBS */
  TFactor = 1.0 / (0.176 + 0.0770 * TSoil - 0.0018 * TSoil * TSoil);

  /* for OJP */
/*   TFactor = 1.0/(.0705 * TSoil - 0.0013 * pow(TSoil, (double) 2.0)); */

  if (TFactor <= 0) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }

  /* equation 14, Wigmosta et al [1994] */

  if (Vpd >= VpdThres) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else
    VpdFactor = 1.0 / (1 - Vpd / VpdThres);

  /* equation 15, Wigmosta et al [1994 */

  RpFactor = 1.0 / ((RsMin / RsMax + Rp / Rpc) / (1 + Rp / Rpc));

  /* equation 16, Wigmosta et al [1994] */

  if (SoilMoisture <= WP) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else if (SoilMoisture < MoistThres)
    MoistFactor = (MoistThres - WP) / (SoilMoisture - WP);
  else
    MoistFactor = 1.0;

  Resistance = TFactor * VpdFactor * RpFactor * MoistFactor * RsMin / LAI;

  return Resistance;

}
