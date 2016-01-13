/*
 * SUMMARY:      SnowMelt.c - Calculate snow accumulation and melt
 * USAGE:        
 *
 * AUTHOR:       Mark Wigmosta and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 08:50:06
 * DESCRIPTION:  Calculate snow accumulation and melt using an energy balance
 *               approach for a two layer snow model
 * DESCRIP-END.
 * FUNCTIONS:    SnowMelt()
 * COMMENTS:
 * $Id: SnowMelt.c,v 1.4 2003/07/01 21:26:25 olivier Exp $     
 */

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "brent.h"
#include "constants.h"
#include "settings.h"
#include "massenergy.h"
#include "functions.h"
#include "snow.h"

static float CalcSnowPackEnergyBalance(float Tsurf, ...);

/*****************************************************************************
  Function name: SnowMelt()

  Purpose      : Calculate snow accumulation and melt using an energy balance
                 approach for a two layer snow model

  Required     :
    int y                  - Row counter
    int x                  - Column counter
    int Dt                 - Model timestep (seconds)
    float Z                - Reference height (m) 
    float Displacement     - Displacement height (m)
    float Z0               - Surface roughness (m)
    float BaseRa           - Aerodynamic resistance (uncorrected for
                             stability) (s/m)
    float AirDens          - Density of air (kg/m3)
    float EactAir          - Actual vapor pressure of air (Pa) 
    float Lv               - Latent heat of vaporization (J/kg3)
    float ShortRad         - Net exchange of shortwave radiation (W/m2)
    float LongRadIn        - Incoming long wave radiation (W/m2)
    float Press            - Air pressure (Pa)
    float RainFall         - Amount of rain (m)
    float Snowfall         - Amount of snow (m)
    float Tair             - Air temperature (C)
    float Vpd              - Vapor pressure deficit (Pa)
    float Wind             - Wind speed (m/s)
    float *PackWater       - Liquid water content of snow pack 
    float *SurfWater	   - Liquid water content of surface layer 
    float *Swq             - Snow water equivalent at current pixel (m)
    float *VaporMassFlux;  - Mass flux of water vapor to or from the
                             intercepted snow (m)
    float *TPack           - Temperature of snow pack (C)
    float *TSurf           - Temperature of snow pack surface layer (C)
    float *MeltEnergy      - Energy used for melting and heating of snow pack
                             (W/m2)

  Returns      :
    float Outflow          - Amount of snowpack outflow (m)

  Modifies     :
    float *PackWater       - Liquid water content of snow pack 
    float *SurfWater	   - Liquid water content of surface layer 
    float *Swq             - Snow water equivalent at current pixel (m)
    float *VaporMassFlux;  - Mass flux of water vapor to or from the
                             intercepted snow (m)
    float *TPack           - Temperature of snow pack (C)
    float *TSurf           - Temperature of snow pack surface layer (C)
    float *MeltEnergy      - Energy used for melting and heating of snow pack
                             (W/m2)

  Comments     :
*****************************************************************************/
float SnowMelt(int y, int x, int Dt, float Z, float Displacement, float Z0,
	       float BaseRa, float AirDens, float EactAir, float Lv,
	       float ShortRad, float LongRadIn, float Press, float RainFall,
	       float SnowFall, float Tair, float Vpd, float Wind,
	       float *PackWater, float *SurfWater, float *Swq,
	       float *VaporMassFlux, float *TPack, float *TSurf,
	       float *MeltEnergy)
{
  float DeltaPackCC;		/* Change in cold content of the pack */
  float DeltaPackSwq;		/* Change in snow water equivalent of the
				   pack (m) */
  float Ice;			/* Ice content of snow pack (m) */
  float InitialSwq;		/* Initial snow water equivalent (m) */
  float MassBalanceError;	/* Mass balance error (m) */
  float MaxLiquidWater;		/* Maximum liquid water content of pack (m) */
  float OldTSurf;		/* Old snow surface temperature (C) */
  float Outflow;		/* Amount water flowing out of the snow pack
				   during the time interval (m) */
  float PackCC;			/* Cold content of snow pack (J) */
  float PackSwq;		/* Snow pack snow water equivalent (m) */
  float Qnet;			/* Net energy exchange at the surface (W/m2) */
  float RefreezeEnergy;		/* refreeze energy (W/m2) */
  float RefrozenWater;		/* Amount of refrozen water (m) */
  float SnowFallCC;		/* Cold content of new snowfall (J) */
  float SnowMelt;		/* Amount of snow melt during time interval
				   (m water equivalent) */
  float SurfaceCC;		/* Cold content of snow pack (J) */
  float SurfaceSwq;		/* Surface layer snow water equivalent (m) */

  InitialSwq = *Swq;
  OldTSurf = *TSurf;

  /* Initialize snowpack variables */

  Ice = *Swq - *PackWater - *SurfWater;

  /* Reconstruct snow pack */
  if (Ice > MAX_SURFACE_SWE)
    SurfaceSwq = MAX_SURFACE_SWE;
  else
    SurfaceSwq = Ice;
  PackSwq = Ice - SurfaceSwq;

  /* Calculate cold contents */
  SurfaceCC = CH_ICE * SurfaceSwq * *TSurf;
  PackCC = CH_ICE * PackSwq * *TPack;
  if (Tair > 0.0)
    SnowFallCC = 0.0;
  else
    SnowFallCC = CH_ICE * SnowFall * Tair;

  /* Distribute fresh snowfall */
  if (SnowFall > (MAX_SURFACE_SWE - SurfaceSwq)) {
    DeltaPackSwq = SurfaceSwq + SnowFall - MAX_SURFACE_SWE;
    if (DeltaPackSwq > SurfaceSwq)
      DeltaPackCC = SurfaceCC + (SnowFall - MAX_SURFACE_SWE) / SnowFall *
	SnowFallCC;
    else
      DeltaPackCC = DeltaPackSwq / SurfaceSwq * SurfaceCC;
    SurfaceSwq = MAX_SURFACE_SWE;
    SurfaceCC += SnowFallCC - DeltaPackCC;
    PackSwq += DeltaPackSwq;
    PackCC += DeltaPackCC;
  }
  else {
    SurfaceSwq += SnowFall;
    SurfaceCC += SnowFallCC;
  }
  if (SurfaceSwq > 0.0)
    *TSurf = SurfaceCC / (CH_ICE * SurfaceSwq);
  else
    *TSurf = 0.0;
  if (PackSwq > 0.0)
    *TPack = PackCC / (CH_ICE * PackSwq);
  else
    *TPack = 0.0;

  /* Adjust ice and *SurfWater */
  Ice += SnowFall;
  *SurfWater += RainFall;

  /* Calculate the surface energy balance for snow_temp = 0.0 */

  Qnet = CalcSnowPackEnergyBalance((float) 0.0, Dt, BaseRa, Z, Displacement,
				   Z0, Wind, ShortRad, LongRadIn, AirDens,
				   Lv, Tair, Press, Vpd, EactAir, RainFall,
				   SurfaceSwq, *SurfWater, OldTSurf,
				   &RefreezeEnergy, VaporMassFlux);

  /* If Qnet == 0.0, then set the surface temperature to 0.0 */
  if (fequal(Qnet, 0.0)) {
    *TSurf = 0.0;
    if (RefreezeEnergy >= 0.0) {
      RefrozenWater = RefreezeEnergy / (LF * WATER_DENSITY) * Dt;
      if (RefrozenWater > *SurfWater) {
	RefrozenWater = *SurfWater;
	RefreezeEnergy = (RefrozenWater * LF * WATER_DENSITY) / Dt;
      }
      *MeltEnergy += RefreezeEnergy;
      SurfaceSwq += RefrozenWater;
      Ice += RefrozenWater;
      *SurfWater -= RefrozenWater;
      assert(*SurfWater >= 0.0);
      SnowMelt = 0.0;

    }
    else {

      /* Calculate snow melt */
      SnowMelt = fabs(RefreezeEnergy) / (LF * WATER_DENSITY) * Dt;
      *MeltEnergy += RefreezeEnergy;
    }

    /* Convert vapor mass flux to a depth per timestep and adjust *SurfWater */
    *VaporMassFlux *= Dt;

    if (*SurfWater < -(*VaporMassFlux)) {
      *VaporMassFlux = -(*SurfWater);
      *SurfWater = 0.0;
    }
    else
      *SurfWater += *VaporMassFlux;

    /* If SnowMelt < Ice, there was incomplete melting of the pack */

    if (SnowMelt < Ice) {
      if (SnowMelt <= PackSwq) {
	*SurfWater += SnowMelt;
	PackSwq -= SnowMelt;
	Ice -= SnowMelt;
      }
      else {
	*SurfWater += SnowMelt + *PackWater;
	*PackWater = 0.0;
	PackSwq = 0.0;
	Ice -= SnowMelt;
	SurfaceSwq = Ice;
      }
    }

    /* Else, SnowMelt > Ice and there was complete melting of the pack */
    else {
      SnowMelt = Ice;
      *SurfWater += Ice;
      SurfaceSwq = 0.0;
      *TSurf = 0.0;
      PackSwq = 0.0;
      *TPack = 0.0;
      Ice = 0.0;
    }
  }

  /* Else, SnowPackEnergyBalance(T=0.0) <= 0.0 */
  else {
    /* Calculate surface layer temperature using "Brent method" */

    *TSurf = RootBrent(y, x, (float) (*TSurf - DELTAT), (float) 0.0,
		       SnowPackEnergyBalance, Dt, BaseRa, Z, Displacement,
		       Z0, Wind, ShortRad, LongRadIn, AirDens, Lv, Tair,
		       Press, Vpd, EactAir, RainFall, SurfaceSwq, *SurfWater,
		       OldTSurf, &RefreezeEnergy, VaporMassFlux);

    /* since we iterated, the surface layer is below freezing and no snowmelt
     */

    SnowMelt = 0.0;

    /* Since updated snow_temp < 0.0, all of the liquid water in the surface
       layer has been frozen */

    SurfaceSwq += *SurfWater;
    Ice += *SurfWater;
    *SurfWater = 0.0;
    *MeltEnergy += (*SurfWater * LF * WATER_DENSITY) / Dt;

    /* Convert mass flux to a depth per timestep and adjust SurfaceSwq */

    *VaporMassFlux *= Dt;

    if (SurfaceSwq < -(*VaporMassFlux)) {
      *VaporMassFlux = -SurfaceSwq;
      SurfaceSwq = 0.0;
      Ice = PackSwq;
    }
    else {
      SurfaceSwq += *VaporMassFlux;
      Ice += *VaporMassFlux;
    }
  }

  /* Done with iteration etc, now Update the liquid water content of the
     surface layer */

  MaxLiquidWater = LIQUID_WATER_CAPACITY * SurfaceSwq;
  if (*SurfWater > MaxLiquidWater) {
    Outflow = *SurfWater - MaxLiquidWater;
    *SurfWater = MaxLiquidWater;
  }
  else
    Outflow = 0.0;

  /* Refreeze liquid water in the pack.                                   
     variable 'RefreezeEnergy' is the heat released to the snow pack            
     if all liquid water were refrozen.                                   
     if RefreezeEnergy < PackCC then all water IS refrozen           
     PackCC always <=0.0 

     WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this does 
     not involve energy transported to the pixel.  Instead heat from the snow 
     pack is used to refreeze water */

  *PackWater += Outflow;	/* add surface layer outflow to pack liquid water */
  RefreezeEnergy = *PackWater * LF * WATER_DENSITY;

  /* calculate energy released to freeze */

  if (PackCC < -RefreezeEnergy) {	/* cold content not fully depleted */
    PackSwq += *PackWater;	/* refreeze all water and update */
    Ice += *PackWater;
    *PackWater = 0.0;
    if (PackSwq > 0.0)
      *TPack = (PackCC + RefreezeEnergy) / (CH_ICE * PackSwq);
    else
      *TPack = 0.0;
  }
  else {
    /* cold content has been either exactly satisfied or exceeded. If
       PackCC = refreeze then pack is ripe and all pack water is
       refrozen, else if energy released in refreezing exceeds PackCC 
       then exactly the right amount of water is refrozen to satify PackCC.
       The refrozen water is added to PackSwq and Ice */

    *TPack = 0.0;
    DeltaPackSwq = -PackCC / (LF * WATER_DENSITY);
    *PackWater -= DeltaPackSwq;
    PackSwq += DeltaPackSwq;
    Ice += DeltaPackSwq;
  }

  /* Update the liquid water content of the pack */

  MaxLiquidWater = LIQUID_WATER_CAPACITY * PackSwq;
  if (*PackWater > MaxLiquidWater) {
    Outflow = *PackWater - MaxLiquidWater;
    *PackWater = MaxLiquidWater;
  }
  else
    Outflow = 0.0;

  /* Update snow properties */

  Ice = PackSwq + SurfaceSwq;

  if (Ice > MAX_SURFACE_SWE) {
    SurfaceCC = CH_ICE * *TSurf * SurfaceSwq;
    PackCC = CH_ICE * *TPack * PackSwq;
    if (SurfaceSwq > MAX_SURFACE_SWE) {
      PackCC += SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq;
      SurfaceCC -= SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq;
      PackSwq += SurfaceSwq - MAX_SURFACE_SWE;
      SurfaceSwq -= SurfaceSwq - MAX_SURFACE_SWE;
    }
    else if (SurfaceSwq < MAX_SURFACE_SWE) {
      PackCC -= PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / PackSwq;
      SurfaceCC += PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / PackSwq;
      PackSwq -= MAX_SURFACE_SWE - SurfaceSwq;
      SurfaceSwq += MAX_SURFACE_SWE - SurfaceSwq;
    }
    *TPack = PackCC / (CH_ICE * PackSwq);
    *TSurf = SurfaceCC / (CH_ICE * SurfaceSwq);
  }
  else {
    PackSwq = 0.0;
    PackCC = 0.0;
    *TPack = 0.0;
  }

  *Swq = Ice + *PackWater + *SurfWater;

  if (fequal(*Swq, 0.0)) {
    *TSurf = 0.0;
    *TPack = 0.0;
  }

  /* Mass balance test */

  MassBalanceError = (InitialSwq - *Swq) + (RainFall + SnowFall) - Outflow +
    *VaporMassFlux;

  return (Outflow);
}

/*****************************************************************************
  Function name: CalcSnowPackEnergyBalance()

  Purpose      : Dummy function to make a direct call to
                 SnowEnergyBalance() possible.

  Required     : 
    float TSurf - SnowPack surface temperature (C)
    other arguments required by SnowPackEnergyBalance()

  Returns      :
    float Qnet - Net energy exchange at the SnowPack snow surface (W/m^2)

  Modifies     : none

  Comments     : function is local to this module
*****************************************************************************/
static float CalcSnowPackEnergyBalance(float Tsurf, ...)
{
  va_list ap;			/* Used in traversing variable argument list
				 */
  float Qnet;			/* Net energy exchange at the SnowPack snow
				   surface (W/m^2) */

  va_start(ap, Tsurf);
  Qnet = SnowPackEnergyBalance(Tsurf, ap);
  va_end(ap);

  return Qnet;
}
