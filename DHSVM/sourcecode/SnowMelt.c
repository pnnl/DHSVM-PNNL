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
 * $Id: SnowMelt.c,v 3.2 2018/5/21 ning Exp $
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
  float *MeltEnergy, float *Iwq, float *GlMelt, float *depth,
  float *density, float *glwater, float *Qold,
  OPTIONSTRUCT * Options, float *IceRemoved)
{
  float DeltaPackCC;		  /* Change in cold content of the pack */
  float DeltaPackSwq;		  /* Change in snow water equivalent of the pack (m) */
  float Ice;			        /* Ice content of glacier pack (m) */
  float InitialSwq;		    /* Initial snow water equivalent (m) */
  float MassBalanceError;	/* Mass balance error (m) */
  float MaxLiquidWater;		/* Maximum liquid water content of pack (m) */
  float OldTSurf;		      /* Old snow surface temperature (C) */
  float Outflow;		      /* Amount water flowing out of the snow pack during the time interval (m) */
  float PackCC;			      /* Cold content of snow pack (J) */
  float PackSwq;		      /* Snow pack snow water equivalent (m) */
  float Qnet;			        /* Net energy exchange at the surface (W/m2) */
  float RefreezeEnergy;		/* refreeze energy (W/m2) */
  float RefrozenWater;		/* Amount of refrozen water (m) */
  float SnowFallCC;		    /* Cold content of new snowfall (J) */
  float SnowMelt;		      /* Amount of snow melt during time interval (m water equivalent) */
  float SurfaceCC;		    /* Cold content of snow pack (J) */
  float SurfaceSwq;		    /* Surface layer snow water equivalent (m) */

  /* needed for third ice layer */
  float SnowIce;          /* Ice content of snow pack (m)*/
  float GlacierIce;
  float InitialIwq;
  float GLIceMelt;
  //Snow Density Variables
  float new_snow;
  float Tavg;
  float ddz2;
  float Ps;
  float delta_depth;
  float density_new;
  float depth_new;


  GLIceMelt = 0.0;
  InitialSwq = *Swq;
  InitialIwq = *Iwq;
  OldTSurf = *TSurf;

  /* Initialize snowpack variables */
  if (isnan(*TSurf))
    printf("OldTSurf= %f\n", *TSurf);

  SnowIce = *Swq - *PackWater - *SurfWater;

  GlacierIce = InitialIwq;
  Ice = SnowIce + GlacierIce;

  /* Reconstruct snow pack */
  if (SnowIce > MAX_SURFACE_SWE) {
    SurfaceSwq = MAX_SURFACE_SWE;
    PackSwq = SnowIce - SurfaceSwq;

  }
  else
  {
    SurfaceSwq = SnowIce;
    PackSwq = 0.0;

  }

  if (SurfaceSwq <= SnowIce)
    PackSwq = SnowIce - SurfaceSwq;
  else
    PackSwq = 0.;

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
  SnowIce += SnowFall;
  Ice += SnowFall;
  *SurfWater += RainFall;

  /* Calculate the surface energy balance for snow_temp = 0.0 */
  if (isnan(*TSurf))
    printf("before snowpack energy= %f\n", *TSurf);

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
      SnowIce += RefrozenWater;
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
    if (SnowIce == 0.0) {
      /*No Snowpack present, handle vapor fluxes on glacier surface */
      if (GlacierIce > 0.0) {
        if (GlacierIce < -(*VaporMassFlux)) {
          *VaporMassFlux = -GlacierIce;
          GlacierIce = 0.0;
          Ice = 0.0;
        }
        else {
          GlacierIce += *VaporMassFlux;
          Ice = GlacierIce;
        }
      }
    }
    else {

      if (*SurfWater < -(*VaporMassFlux)) {
        *VaporMassFlux = -(*SurfWater);
        *SurfWater = 0.0;
      }
      else
        *SurfWater += *VaporMassFlux;
    }
    /* If SnowMelt < Ice, there was incomplete melting of the pack */
    if (SnowMelt <= SnowIce) {
      if (SnowMelt <= PackSwq) {
        *SurfWater += SnowMelt;
        PackSwq -= SnowMelt;
        Ice -= SnowMelt;
        SnowIce -= SnowMelt;
      }
      else { /* Melt all of pack layer and part of surface layer. */
     //else {
        *SurfWater += SnowMelt + *PackWater;
        *PackWater = 0.0;
        SurfaceSwq -= (SnowMelt - PackSwq);
        PackSwq = 0.0;
        SnowIce -= SnowMelt;
        Ice -= SnowMelt;
      }
    }

    else { /* Snowmelt > SnowIce: Melt snow pack completely and also part of the ice */
      // printf("Melt snow pack completely and also part of the ice\n");
      if (SnowMelt < Ice) {
        *SurfWater += SnowIce + *PackWater;
        *PackWater = 0.0;
        PackSwq = 0.0;
        Ice -= SnowMelt;
        GLIceMelt = SnowMelt - SnowIce;
        GlacierIce -= GLIceMelt;
        SnowMelt = SnowIce;
        SurfaceSwq = 0.0;
        SnowIce = 0.0;
        *TSurf = 0.0;
        *TPack = 0.0;
      }
      else { /* Else, SnowMelt > Ice and there was complete melting of the glacier */
        SnowMelt = Ice;
        GLIceMelt = GlacierIce;
        GlacierIce = 0.0;
        Ice = 0.0;
        *SurfWater += SnowIce + *PackWater;
        SurfaceSwq = 0.0;
        *TSurf = 0.0;
        PackSwq = 0.0;
        *TPack = 0.0;
        SnowIce = 0.0;
        /*NOTE: GLIceMelt is added outflow at the end of the code */
        /* readjust melt energy to account for melt only of available snow (W.I.P. NOT SURE ABOUT THIS)*/
        //*MeltEnergy -= RefreezeEnergy;
        //RefreezeEnergy = RefreezeEnergy / fabs(RefreezeEnergy) * SnowMelt  * LF * WATER_DENSITY / Dt;
        //*MeltEnergy  += RefreezeEnergy;
      }
    }
  }

  /* Else, SnowPackEnergyBalance(T=0.0) <= 0.0 */
  else {
    /* Calculate surface layer temperature using "Brent method" */
    //if(isnan(*TSurf))
    //if(*TSurf<-50.0)
    //printf("before root brent call TSurf = %f Swq = %f Iwq = %f\n", *TSurf, *Swq, *Iwq);
    *TSurf = RootBrent(y, x, (float)(*TSurf - DELTAT), (float)0.0,
      SnowPackEnergyBalance, Dt, BaseRa, Z, Displacement,
      Z0, Wind, ShortRad, LongRadIn, AirDens, Lv, Tair,
      Press, Vpd, EactAir, RainFall, SurfaceSwq, *SurfWater,
      OldTSurf, &RefreezeEnergy, VaporMassFlux);

    /* since we iterated, the surface layer is below freezing and no snowmelt */
    if (fabs(*TSurf) <= 1e-6)
      *TSurf = 0.0;

    SnowMelt = 0.0;
    GLIceMelt = 0.0;
    /* Since updated snow_temp < 0.0, all of the liquid water in the surface
       layer has been frozen */
    SurfaceSwq += *SurfWater;
    Ice += *SurfWater;
    SnowIce += *SurfWater;
    *SurfWater = 0.0;
    *MeltEnergy += (*SurfWater * LF * WATER_DENSITY) / Dt;

    /* Convert mass flux to a depth per timestep and adjust SurfaceSwq */
    *VaporMassFlux *= Dt;
    if (SnowIce == 0.0) {
      /*No Snowpack present, handle vapor fluxes on glacier surface */
      if (GlacierIce > 0.0) {
        if (GlacierIce < -(*VaporMassFlux)) {
          *VaporMassFlux = -GlacierIce;
          GlacierIce = 0.0;
          Ice = 0.0;
        }
        else {
          GlacierIce += *VaporMassFlux;
          Ice = GlacierIce;
        }
      }
    }
    else {
      if (SurfaceSwq < -(*VaporMassFlux)) {
        *VaporMassFlux = -SurfaceSwq;
        SurfaceSwq = 0.0;
        SnowIce = PackSwq;
        Ice = PackSwq + GlacierIce;
      }
      else {
        SurfaceSwq += *VaporMassFlux;
        SnowIce += *VaporMassFlux;
      }
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
    SnowIce += *PackWater;
    *PackWater = 0.0;
    if (PackSwq > 0.0) {
      *TPack = (PackCC + RefreezeEnergy) / (CH_ICE * (PackSwq));
      if (*TPack > 0.)
        *TPack = 0.0;
      else
        *TPack = 0.0;
    }
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
    SnowIce += DeltaPackSwq;
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

  Ice = GlacierIce + PackSwq + SurfaceSwq;

  if (SnowIce > MAX_SURFACE_SWE) {
    SurfaceCC = CH_ICE * *TSurf * SurfaceSwq;
    PackCC = CH_ICE * *TPack * PackSwq;
    if (SurfaceSwq > MAX_SURFACE_SWE) {
      PackCC += SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq;
      SurfaceCC -= SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq;
      PackSwq += SurfaceSwq - MAX_SURFACE_SWE;
      SurfaceSwq -= SurfaceSwq - MAX_SURFACE_SWE;
    }
    else if (SurfaceSwq < MAX_SURFACE_SWE) {
      PackCC -= PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / (PackSwq);
      SurfaceCC += PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / (PackSwq);
      PackSwq -= MAX_SURFACE_SWE - SurfaceSwq;
      SurfaceSwq += MAX_SURFACE_SWE - SurfaceSwq;
    }
    *TPack = PackCC / (CH_ICE * (PackSwq));
    *TSurf = SurfaceCC / (CH_ICE * SurfaceSwq);
  }
  else {
    PackSwq = 0.0;
    PackCC = 0.0;
    *TPack = 0.0;
  }

  *Swq = SnowIce + *PackWater + *SurfWater;
  *Iwq = GlacierIce;
  Outflow += (GLIceMelt);
  *GlMelt = GLIceMelt;

  if (fequal(*Swq, 0.0)) {
    *TSurf = 0.0;
    *TPack = 0.0;
  }
  if (*Iwq <= 0.0) {
    *Iwq = 0.0;
  }

  /* Calculate Snow Density following algorithm using in the VIC model taken from BRAS/SNTHERM89 */
  new_snow = (*Swq - InitialSwq) * 1000;
  if (new_snow <= 0.0) {
    new_snow = 0.0;
  }
  /* Estimate density of new snow based on air temperature */
  if (new_snow > 0.)
    density_new = 67.9 + 51.3 * exp(Tair / 2.6);
  else
    density_new = 0.0;

  /* Estimate average snowpack temperature */
  Tavg = *TSurf + KELVIN;
  if (new_snow > 0.0) {
    if (*depth > 0.) {
      /* Compact current snowpack by weight of new snowfall */
      delta_depth = (((new_snow / 25.4) * (*depth / 0.0254)) / (InitialSwq / 0.0254)
        * pow((*depth / 0.0254) / 10., 0.35)) * 0.0254;
      if (delta_depth > 0.9 * *depth) delta_depth = 0.9 * *depth;
      depth_new = new_snow / density_new;
      *depth = *depth - delta_depth + depth_new;
      *density = 1000. * *Swq / *depth;
    }
    else {
      /* no snowpack present, so snow density equals that of new snow */
      *density = density_new;
      *depth = 1000. * *Swq / *density;
    }
  }
  else *density = 1000. * *Swq / *depth;

  /** Densification of the snow pack due to aging **/
  /** based on SNTHRM89 R. Jordan 1991**/
  if (*depth > 0.) {
    Ps = 0.5 * G * RHO_W * *Swq;
    ddz2 = SNDENS_ETA0 * exp(-SNDENS_C5*(Tavg - KELVIN) + SNDENS_C6* *density);
    delta_depth = Ps / ddz2 * *depth * Dt * SECPHOUR;
    if (delta_depth > 0.9 * *depth) delta_depth = 0.9 * *depth;
    *depth -= delta_depth;
    *density = 1000. * *Swq / *depth;
  }
  /* Adjust PackSwq and Iwq for ice formation from pack to glacier */
  /* W.I.P. Due to certain instabilities in the density formulation */
  /* Conversion to ice is restricted to when snowpack is greater than*/
  /* 3m w.e., and 1m w.e. is left after conversion, this ensures that */
  /* Snow does not convert erroneous to ice too soon at high elevations*/
  /* Resulting in falsely exposing low albedo ice at high elevations */
  /* These values can be adjused regionally, depending on accumulation */
  if (PackSwq > 5.0) {
    if (*density > 850.0) {
      //printf("SNOW TO ICE Pack SWQ = %f Denisty = %f\n", PackSwq, *density);
      *Iwq += PackSwq - 4.0;
      PackCC *= (4.0 / PackSwq);
      PackSwq = 4.0;
      //*PackWater = 0.0;
      *Swq = SurfaceSwq + *SurfWater + PackSwq + *PackWater;
      *density = 537.098;
      SnowIce = SurfaceSwq;
      *depth = *Swq * 1000 / *density;
      GlacierIce = *Iwq;
      Ice = SnowIce + GlacierIce;
      *TPack = PackCC / (CH_ICE * (PackSwq));
    }
  }
  /*Delete Ice From Simulation if No Glaciers are too be simulated */
  if (Options->Glacier == NO_GLACIER)
  {
    *IceRemoved += *Iwq;
    *Iwq = 0.0;
  }

  /* Mass balance test */
  MassBalanceError = (InitialSwq - *Swq) + (InitialIwq - *Iwq) + (RainFall + SnowFall) - Outflow + *VaporMassFlux;
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
  va_list ap;			/* Used in traversing variable argument list */
  float Qnet;			/* Net energy exchange at the SnowPack snow surface (W/m^2) */

  va_start(ap, (float)Tsurf);
  Qnet = SnowPackEnergyBalance((float)Tsurf, ap);
  va_end(ap);

  return Qnet;
}
