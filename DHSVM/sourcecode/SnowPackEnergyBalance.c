/*
 * SUMMARY:      SnowPackEnergyBalance.c - Calculate snow pack energy balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 09:09:29
 * DESCRIPTION:  Calculate snow pack energy balance
 * DESCRIP-END.
 * FUNCTIONS:    SnowPackEnergyBalance()
 * COMMENTS:
 * $Id: SnowPackEnergyBalance.c,v 1.4 2003/07/01 21:26:25 olivier Exp $     
 */

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "settings.h"
#include "constants.h"
#include "massenergy.h"
#include "snow.h"
#include "functions.h"

/*****************************************************************************
  Function name: SnowPackEnergyBalance()

  Purpose      : Calculate the surface energy balance for the snow pack

  Required     :
    float TSurf           - new estimate of effective surface temperature
    va_list ap            - Argument list initialized by va_start().  For
                            elements of list and order, see beginning of
                            routine

  Returns      :
    float RestTerm        - Rest term in the energy balance

  Modifies     : 
    float *RefreezeEnergy - Refreeze energy (W/m2) 
    float *VaporMassFlux  - Mass flux of water vapor to or from the
                            intercepted snow 

  Comments     :
    Reference:  Bras, R. A., Hydrology, an introduction to hydrologic
                science, Addisson Wesley, Inc., Reading, etc., 1990.
*****************************************************************************/
float SnowPackEnergyBalance(float TSurf, va_list ap)
{
  /* start of list of arguments in variable argument list */

  int Dt;			/* Model time step (hours) */
  float Ra;			/* Aerodynamic resistance (s/m) */
  float Z;			/* Reference height (m) */
  float Displacement;		/* Displacement height (m) */
  float Z0;			/* Roughness length (m) */
  float Wind;			/* Wind speed (m/s) */
  float ShortRad;		/* Net incident shortwave radiation (W/m2) */
  float LongRadIn;		/* Incoming longwave radiation (W/m2) */
  float AirDens;		/* Density of air (kg/m3) */
  float Lv;			/* Latent heat of vaporization (J/kg3) */
  float Tair;			/* Air temperature (C) */
  float Press;			/* Air pressure (Pa) */
  float Vpd;			/* Vapor pressure deficit (Pa) */
  float EactAir;		/* Actual vapor pressure of air (Pa) */
  float Rain;			/* Rain fall (m/timestep) */
  float SweSurfaceLayer;	/* Snow water equivalent in surface layer (m)
				 */
  float SurfaceLiquidWater;	/* Liquid water in the surface layer (m) */
  float OldTSurf;		/* Surface temperature during previous time
				   step */
  float *RefreezeEnergy;	/* Refreeze energy (W/m2) */
  float *VaporMassFlux;		/* Mass flux of water vapor to or from the
				   intercepted snow */

  /* end of list of arguments in variable argument list */

  float AdvectedEnergy;		/* Energy advected by precipitation (W/m2) */
  float DeltaColdContent;	/* Change in cold content (W/m2) */
  float EsSnow;			/* saturated vapor pressure in the snow pack
				   (Pa)  */
  float LatentHeat;		/* Latent heat exchange at surface (W/m2) */
  float LongRadOut;		/* long wave radiation emitted by surface
				   (W/m2) */
  float Ls;			/* Latent heat of sublimation (J/kg) */
  float NetRad;			/* Net radiation exchange at surface (W/m2) */
  float RestTerm;		/* Rest term in surface energy balance
				   (W/m2) */
  float SensibleHeat;		/* Sensible heat exchange at surface (W/m2) */
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
  Tair = (float) va_arg(ap, double);
  Press = (float) va_arg(ap, double);
  Vpd = (float) va_arg(ap, double);
  EactAir = (float) va_arg(ap, double);
  Rain = (float) va_arg(ap, double);
  SweSurfaceLayer = (float) va_arg(ap, double);
  SurfaceLiquidWater = (float) va_arg(ap, double);
  OldTSurf = (float) va_arg(ap, double);
  RefreezeEnergy = (float *) va_arg(ap, double *);
  VaporMassFlux = (float *) va_arg(ap, double *);

  /* Calculate active temp for energy balance as average of old and new  */

  TMean = 0.5 * (OldTSurf + TSurf);

  /* Correct aerodynamic conductance for stable conditions
     Note: If air temp >> snow temp then aero_cond -> 0 (i.e. very stable)
     velocity (vel_2m) is expected to be in m/sec */

  /* Apply the stability correction to the aerodynamic resistance 
     NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
     think that it is more correct to calculate ALL fluxes at the same
     reference level */

  if (Wind > 0.0)
    Ra /= StabilityCorrection(2.0f, 0.f, TMean, Tair, Wind, Z0);
  else
    Ra = DHSVM_HUGE;

  /* Calculate longwave exchange and net radiation */

  Tmp = TMean + 273.15;
  LongRadOut = STEFAN * (Tmp * Tmp * Tmp * Tmp);
  NetRad = ShortRad + LongRadIn - LongRadOut;

  /* Calculate the sensible heat flux */

  SensibleHeat = AirDens * CP * (Tair - TMean) / Ra;

  /* Calculate the mass flux of ice to or from the surface layer */

  /* Calculate the saturated vapor pressure in the snow pack, 
     (Equation 3.32, Bras 1990) */

  EsSnow = SatVaporPressure(TMean);

  *VaporMassFlux = AirDens * (EPS / Press) * (EactAir - EsSnow) / Ra;
  *VaporMassFlux /= WATER_DENSITY;
  if (fequal(Vpd, 0.0) && *VaporMassFlux < 0.0)
    *VaporMassFlux = 0.0;

  /* Calculate latent heat flux */

  if (TMean >= 0.0) {
    /* Melt conditions: use latent heat of vaporization */
    LatentHeat = Lv * *VaporMassFlux * WATER_DENSITY;
  }
  else {
    /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
    Ls = (677. - 0.07 * TMean) * JOULESPCAL * GRAMSPKG;
    LatentHeat = Ls * *VaporMassFlux * WATER_DENSITY;
  }

  /* Calculate advected heat flux from rain 
     WORK IN PROGRESS:  Should the following read (Tair - Tsurf) ?? */

  AdvectedEnergy = (CH_WATER * Tair * Rain) / Dt;

  /* Calculate change in cold content */

  DeltaColdContent = CH_ICE * SweSurfaceLayer * (TSurf - OldTSurf) / Dt;

  /* Calculate net energy exchange at the snow surface */

  RestTerm = NetRad + SensibleHeat + LatentHeat + AdvectedEnergy -
    DeltaColdContent;

  *RefreezeEnergy = (SurfaceLiquidWater * LF * WATER_DENSITY) / Dt;

  if (fequal(TSurf, 0.0) && RestTerm > -(*RefreezeEnergy)) {
    *RefreezeEnergy = -RestTerm;	/* available energy input over cold content
					   used to melt, i.e. Qrf is negative value
					   (energy out of pack) */
    RestTerm = 0.0;
  }
  else {
    RestTerm += *RefreezeEnergy;	/* add this positive value to the pack */
  }

  return RestTerm;
}
