/*
 * SUMMARY:      SurfaceEnergyBalance.c - Calculate surface energy balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculate surface energy balance.  This group of functions
 *               is used by the iterative Brent method to determine the
 *               surface temperature 
 * DESCRIP-END.
 * FUNCTIONS:    SurfaceEnergyBalance()
 * COMMENTS:
 * $Id: SurfaceEnergyBalance.c,v 1.4 2003/07/01 21:26:26 olivier Exp $     
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  Function name: SurfaceEnergyBalance()

  Purpose      : Calculate the surface energy balance in the absence of snow

  Required     :
    float TSurf           - new estimate of effective surface temperature
    va_list ap            - Argument list initialized by va_start().  For
                            elements of list and order, see beginning of
                            routine

  Returns      :
    float RestTerm        - Rest term in the energy balance

  Modifies     : none

  Comments     :
*****************************************************************************/
float SurfaceEnergyBalance(float TSurf, va_list ap)
{
  /* start of list of arguments in variable argument list */

  int Dt;			/* Model time step (seconds) */
  float Ra;			/* Aerodynamic resistance (s/m) */
  float Z;			/* Reference height (m) */
  float Displacement;		/* Displacement height (m) */
  float Z0;			/* Surface roughness (m) */
  float Wind;			/* Wind speed (m/s) */
  float ShortRad;		/* Net incident shortwave radiation (W/m2) */
  float LongRadIn;		/* Incoming longwave radiation (W/m2) */
  float AirDens;		/* Density of air (kg/m3) */
  float Lv;			/* Latent heat of vaporization (J/kg3) */
  float ETot;			/* Total evapotranspiration (m) */
  float Kt;			/* Effective soil thermal conductivity 
				   (W/(m*K)) */
  float ChSoil;			/* Soil thermal capacity (J/(kg*K)) */
  float Porosity;		/* Porosity of upper soil layer */
  float MoistureContent;	/* Moisture content of upper soil layer */
  float Depth;			/* Depth of soil heat profile (m) */
  float Tair;			/* Air temperature (C) */
  float TSoilUpper;		/* Soil temperature in upper layer (C) */
  float TSoilLower;		/* Soil temperature at Depth (C) */
  float OldTSurf;		/* Surface temperature during previous time
				   step */
  float MeltEnergy;		/* Energy used to melt/refreeze snow pack 
				   (W/m2) */

  /* end of list of arguments in variable argument list */

  float GroundHeat;		/* ground heat exchange at surface (W/m2) */
  float HeatCapacity;		/* soil heat capacity (J/(m3*C) */
  float HeatStorageChange;	/* change in ground heat storage (W/m2) */
  float LatentHeat;		/* latent heat exchange at surface (W/m2) */
  float LongRadOut;		/* long wave radiation emitted by surface
				   (W/m2) */
  float NetRad;			/* net radiation exchange at surface (W/m2) */
  float RestTerm;		/* rest term in surface energy balance
				   (W/m2) */
  float SensibleHeat;		/* sensible heat exchange at surface (W/m2) */
  float TMean;			/* Mean temperature during interval (C) */
  double Tmp;			/* temporary variable */

  /* Assign the elements of the array to the appropriate variables.  The list
     is traversed as if the elements are doubles, because:

     In the variable-length part of variable-length argument lists, the old
     ``default argument promotions'' apply: arguments of type float are
     always promoted (widened) to type double, and types char and short int
     are promoted to int. Therefore, it is never correct to invoke
     va_arg(argp, float); instead you should always use va_arg(argp,
     double). 

     (quoted from the comp.lang.c FAQ list)
   */

  Dt = va_arg(ap, int);
  Ra = (float) va_arg(ap, double);
  Z = (float) va_arg(ap, double);
  Displacement = (float) va_arg(ap, double);
  Z0 = (float) va_arg(ap, double);
  Wind = (float) va_arg(ap, double);
  ShortRad = (float) va_arg(ap, double);
  LongRadIn = (float) va_arg(ap, double);
  AirDens = (float) va_arg(ap, double);
  Lv = (float) va_arg(ap, double);
  ETot = (float) va_arg(ap, double);
  Kt = (float) va_arg(ap, double);
  ChSoil = (float) va_arg(ap, double);
  Porosity = (float) va_arg(ap, double);
  MoistureContent = (float) va_arg(ap, double);
  Depth = (float) va_arg(ap, double);
  Tair = (float) va_arg(ap, double);
  TSoilUpper = (float) va_arg(ap, double);
  TSoilLower = (float) va_arg(ap, double);
  OldTSurf = (float) va_arg(ap, double);
  MeltEnergy = (float) va_arg(ap, double);

  /* In this routine transport of energy to the surface is considered 
     positive */

  TMean = 0.5 * (OldTSurf + TSurf);

  /* Apply the stability correction to the aerodynamic resistance */

  if (Wind > 0.0)
    Ra /= StabilityCorrection(Z, Displacement, TMean, Tair, Wind, Z0);
  else
    Ra = DHSVM_HUGE;

  /* Calculate the longwave exchange and net radiation, assuming black body */

  Tmp = TMean + 273.15;
  LongRadOut = STEFAN * (Tmp * Tmp * Tmp * Tmp);
  NetRad = ShortRad + LongRadIn - LongRadOut;

  /* Calculate the sensible heat flux */

  SensibleHeat = AirDens * CP * (Tair - TMean) / Ra;

  /* Calculate the latent heat flux */

  LatentHeat = -(Lv * ETot) / Dt * WATER_DENSITY;

  /* Calculate the ground heat flux */

  GroundHeat = Kt * (TSoilLower - TMean) / Depth;

  /* Calculate the change in the ground heat storage in the upper 
     0.1 m of the soil */

  HeatCapacity = (1 - Porosity) * ChSoil;
  if (TSoilUpper >= 0.0)
    HeatCapacity += MoistureContent * CH_WATER;
  else
    HeatCapacity += MoistureContent * CH_ICE;

  HeatStorageChange = (HeatCapacity * (OldTSurf - TMean) * DZ_TOP) / Dt;

  /* Calculate the net energy exchange at the surface.  The left hand side of 
     the equation should go to zero for the balance to close, so we want to 
     minimize the absolute value of the left hand side */

  RestTerm =
    MeltEnergy + NetRad + SensibleHeat + LatentHeat +
    GroundHeat + HeatStorageChange;

  return RestTerm;
}
