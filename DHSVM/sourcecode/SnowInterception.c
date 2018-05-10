/*
* SUMMARY:      SnowInterception.c - simulates snow interception and release
* USAGE:
*
* AUTHOR:       Pascal Storck
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       pstorck@u.washington.edu
* ORIG-DATE:    29-Aug-1996 at 13:42:17
* DESCRIPTION:  Calculates the interception and subsequent release of
*               by the forest canopy using an energy balance approach
* DESCRIP-END.
* FUNCTIONS:    SnowInterception()
* COMMENTS:
* $Id: SnowInterception.c,v 1.5 2003/11/12 20:01:52 colleen Exp $
*/

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "brent.h"
#include "constants.h"
#include "settings.h"
#include "massenergy.h"
#include "snow.h"
#include "functions.h"
#include "timing.h"

/*****************************************************************************
Function name: SnowInterception()

Purpose      : Calculate snow interception and release by the canopy

Required     :
int y                  - Row counter
int x                  - Column counter
int Dt                 - Model timestep (seconds)
float F                - Fractional coverage
float LAI              - Leaf Area Index
float MaxInt           - Maximum rainfall interception storage (m)
float BaseRa           - Aerodynamic resistance (uncorrected for
stability) (s/m)
float AirDens          - Density of air (kg/m3)
float EactAir          - Actual vapor pressure of air (Pa)
float Lv               - Latent heat of vaporization (J/kg3)
PIXRAD *LocalRad       - Components of radiation balance for current pixel
(W/m2)
float Press            - Air pressure (Pa)
float Tair             - Air temperature (C)
float Vpd	           - Vapor pressure deficit (Pa)
float Wind             - Wind speed (m/s)
float *RainFall        - Amount of rain (m)
float *Snowfall        - Amount of snow (m)
float *IntRain         - Intercepted rain (m)
float *IntSnow         - Snow water equivalent of intercepted snow (m)
float *TempIntStorage  - Temporary storage for snowmelt and rainfall
involved in mass release calculations (m)
float *VaporMassFlux   - Vapor mass flux to/from intercepted snow
(m/timestep)
float *Tcanopy         - Canopy temperature (C)
float *MeltEnergy      - Energy used in heating and melting of the snow
(W/m2)
float *Height          - Height of vegetation (m)

Returns      : none

Modifies     :
float *RainFall        - Amount of rain (m)
float *Snowfall        - Amount of snow (m)
float *IntRain         - Intercepted rain (m)
float *IntSnow         - Snow water equivalent of intercepted snow (m)
float *TempIntStorage  - Temporary storage for snowmelt and rainfall
involved in mass release calculations (m)
float *VaporMassFlux   - Vapor mass flux to/from intercepted snow
(m/timestep)
float *Tcanopy         - Canopy temperature (C)

Comments     : Only the top canopy layer is taken into account for snow
interception.  Snow interception by lower canopy is
disregarded.  Rain water CAN be intercepted by lower canopy
layers (similar to InterceptionStorage()).
Of course:  NO vegetation -> NO interception
*****************************************************************************/
void SnowInterception(OPTIONSTRUCT *Options, int y, int x, int Dt, float F,
  float Vf, float LAI, float MaxInt, float MaxSnowIntCap, float MDRatio,
  float SnowIntEff, float Ra, float AirDens, float EactAir,
  float Lv, PIXRAD *LocalRad, float Press, float Tair,
  float Vpd, float Wind, float *RainFall, float *SnowFall,
  float *IntRain, float *IntSnow, float *TempIntStorage,
  float *VaporMassFlux, float *Tcanopy, float *MeltEnergy,
  float *Height, unsigned char Understory)
{
  float AdvectedEnergy;		/* Energy advected by the rain (W/m2) */
  float DeltaSnowInt;		/* Change in the physical swe of snow
                            interceped on the branches. (m) */
  float Drip;			    /* Amount of drip from intercepted snow as a
                            result of snowmelt (m) */
  float ExcessSnowMelt;		/* Snowmelt in excess of the water holding
                            capacity of the tree (m) */
  float EsSnow;			    /* saturated vapor pressure in the snow pack (Pa)  */
  float InitialSnowInt;		/* Initial intercepted snow (m) */
  float InitialWaterInt;	/* Initial intercepted water (snow and rain) (m) */
  float LatentHeat;		    /* Latent heat flux (W/m2) */
  float LongOut;		    /* Longwave radiation emitted by canopy (W/m2) */
  float Ls;			        /* Latent heat of sublimation (J/(kg K) */
  float MassBalanceError;	/* Mass blalnce to make sure no water is
                            being destroyed/created (m) */
  float MaxWaterInt;		/* Water interception capacity (m) */
  float MaxSnowInt;		    /* Snow interception capacity (m) - multiplier w/ temp */
  float NetRadiation;
  float PotSnowMelt;		/* Potential snow melt (m) */
  float RainThroughFall;	/* Amount of rain reaching to the ground (m) */
  float RefreezeEnergy;		/* Energy available for refreezing or melt */
  float ReleasedMass;		/* Amount of mass release of intercepted snow (m) */
  float SensibleHeat;		/* Sensible heat flux (W/m2) */
  float SnowThroughFall;	/* Amount of snow reaching to the ground (m) */
  float Tmp;			    /* Temporary variable */
  float MaxIntercept;		/* max snow interception - regardless of temp */
  float overload;		    /* overload of intercepted snow due to rainfall
                            or condensation */
  float intrainfrac;		/* fraction of intercepted water which is liquid */
  float intsnowfrac;		/*fraction of intercepted water which is solid */
  float OriginalRainfall;

  TIMING_TASK_START("Snow interception", 3);

  /* Initialize Drip, H2O balance, and mass release variables. */
  OriginalRainfall = *RainFall;
  InitialWaterInt = *IntSnow + *IntRain;

  /* Convert from pixel depth to physical depth*/
  *IntSnow /= F;
  *IntRain /= F;

  InitialSnowInt = *IntSnow;

  Drip = 0.0;
  ReleasedMass = 0.0;

  /* Determine the maximum snow interception water equivalent.
  Kobayashi, D., 1986, Snow Accumulation on a Narrow Board,
  Cold Regions Science and Technology, (13), pp. 239-245. Figure 4. */
  if (Tair > -5.0)
    MaxSnowInt = 1.0;
  else
    MaxSnowInt = 0.25;

  /* Adjust snow interception capacity (unit: m) to the user-provided input
  parameter, maximum snow interception capacity (%), which varies with
  vegetation type */
  MaxSnowInt *= MaxSnowIntCap;
  /* max snow interception - regardless of temp */
  MaxIntercept = MaxSnowIntCap;

  /* Calculate snow interception from snowfall of current timestep limited
  by max snow interception capacity */
  DeltaSnowInt = SnowIntEff * *SnowFall;
  if (DeltaSnowInt + *IntSnow > MaxSnowInt)
    DeltaSnowInt = MaxSnowInt - *IntSnow;
  if (DeltaSnowInt < 0.0)
    DeltaSnowInt = 0.0;

  /* pixel snow fall depth -- compposed of throughfall through canopy
  and snowfall on open space */
  SnowThroughFall = (*SnowFall - DeltaSnowInt) * F + (*SnowFall) * (1 - F);

  /* update total intercepted snow depth */
  *IntSnow += DeltaSnowInt;

  /* Calculate amount of water/rain intercepted on branches and stored in
  intercepted snow.
  MaxInt = LAI * F * Rain LAI Multiplier (input constant -- around 0.0001)
  LAI Multiplier for rain interception */
  MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;

  if ((*IntRain + *RainFall) <= MaxWaterInt) {
    *IntRain += *RainFall;
    RainThroughFall = *RainFall * (1 - F);
  }
  else {
    RainThroughFall = (*IntRain + *RainFall - MaxWaterInt) * F +
      (*RainFall * (1 - F));
    *IntRain = MaxWaterInt;
  }

  /* Now that total intercepted water has been calculated, allow for structural
  unloading of branches i.e. if absolute maximum capacity is reached then
  allow sliding due to branch bending. Of course, if chunks of snow are
  falling, they can contain both ice and liquid water - Let both of these
  come off in the correct proportions */
  if (*IntRain + *IntSnow > MaxIntercept) {
    overload = (*IntRain + *IntSnow) - MaxIntercept;
    intsnowfrac = *IntSnow / (*IntSnow + *IntRain);
    intrainfrac = *IntRain / (*IntSnow + *IntRain);
    *IntRain = *IntRain - overload * intrainfrac;
    *IntSnow = *IntSnow - overload * intsnowfrac;
    SnowThroughFall = SnowThroughFall + overload * intsnowfrac * F;
    RainThroughFall = RainThroughFall + overload * intrainfrac * F;
  }

  /* The canopy temperature is assumed to be equal to the air temperature if
  the air temperature is below 0 deg C, otherwise the canopy temperature is
  equal to 0 deg C */
  if (Tair > 0.)
    *Tcanopy = 0.;
  else
    *Tcanopy = Tair;

  /* Calculate the net radiation at the canopy surface, using the canopy
  temperature. The outgoing longwave is subtracted twice, because the
  canopy radiates in two directions */
  Tmp = *Tcanopy + 273.15;
  LongOut = STEFAN * (Tmp * Tmp * Tmp * Tmp);

  /* If the improved radiation scheme, the forested cell is considered continuous with
  small canopy gaps. This assumption is valid for the typical DHSVM implmentation on the
  spatial scale of 30-150 m. As such,  F is not needed for weighting radiation  */
  if (Options->ImprovRadiation) {
    NetRadiation = LocalRad->NetShort[0] + LocalRad->LongIn[0] - 2 * Vf*LongOut;
    NetRadiation /= F;
  }
  else {
    NetRadiation = LocalRad->NetShort[0] + LocalRad->LongIn[0] - 2 * F*LongOut;
    NetRadiation /= F;
  }

  /* Calculate the vapor mass flux between the canopy and the surrounding air mass
  - snow covered aerodynamic resistance is assumed to increase by an order
  of magnitude based on Lunderg et al 1998, Journal of Hydrological Processes */
  EsSnow = SatVaporPressure(*Tcanopy);
  *VaporMassFlux = AirDens * (EPS / Press) * (EactAir - EsSnow) / (Ra * 10.0);
  *VaporMassFlux /= WATER_DENSITY;
  if (fequal(Vpd, 0.0) && *VaporMassFlux < 0.0)
    *VaporMassFlux = 0.0;

  /* Calculate the latent heat flux */
  Ls = (677. - 0.07 * *Tcanopy) * JOULESPCAL * GRAMSPKG;
  LatentHeat = Ls * *VaporMassFlux * WATER_DENSITY;

  /* Calculate the sensible heat flux */
  SensibleHeat = AirDens * CP * (Tair - *Tcanopy) / (Ra * 10.0);

  /* Calculate the advected energy */
  AdvectedEnergy = (CH_WATER * Tair * *RainFall) / Dt;

  /* Calculate the amount of energy available for refreezing */
  RefreezeEnergy = SensibleHeat + LatentHeat + NetRadiation + AdvectedEnergy;
  RefreezeEnergy *= Dt;

  /* if RefreezeEnergy is positive it means energy is available to melt the
  intercepted snow in the canopy.  If it is negative, it means that
  intercepted water will be refrozen */

  /* Update maximum water interception storage */
  MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;

  /* Convert the vapor mass flux from a flux to a depth per interval */
  *VaporMassFlux *= Dt;

  if (RefreezeEnergy > 0.0) {	/*we've got melt */
    if (-(*VaporMassFlux) > *IntRain) {
      *VaporMassFlux = -(*IntRain);
      *IntRain = 0.;
    }
    else
      *IntRain += *VaporMassFlux;

    PotSnowMelt = MIN((RefreezeEnergy / (LF * WATER_DENSITY)), *IntSnow);

    *MeltEnergy -= (LF * PotSnowMelt * WATER_DENSITY) / Dt;

    /* if the intercepted rain and potential snowmelt is less than the
    liquid water holding capacity of the intercepted snowpack, then simply
    add the total potential snowmelt to the liquid water content of the
    intercepted snowpack. */
    if ((*IntRain + PotSnowMelt) <= MaxWaterInt) {
      *IntSnow -= PotSnowMelt;
      *IntRain += PotSnowMelt;
      PotSnowMelt = 0.0;
    }
    else {
      ExcessSnowMelt = PotSnowMelt + *IntRain - MaxWaterInt;

      *IntSnow -= MaxWaterInt - (*IntRain);
      *IntRain = MaxWaterInt;
      if (*IntSnow < 0.0)
        *IntSnow = 0.0;

      if (SnowThroughFall > 0.0 && InitialSnowInt <= MIN_INTERCEPTION_STORAGE) {
        /* Water in excess of MaxWaterInt has been generated.  If it is
        snowing and there was little intercepted snow at the beginning of the
        time step ( <= MIN_INTERCEPTION_STORAGE), then allow the snow to melt
        as it is intercepted.  Also, enforce that the following holds true:
        if Intercpeted snow is below minimum thresold then it can only be
        removed via melting */
        Drip += ExcessSnowMelt;
        *IntSnow -= ExcessSnowMelt;
        if (*IntSnow < 0.0)
          *IntSnow = 0.0;
        ExcessSnowMelt = 0.0;
      }
      else
        /* Else, SnowThroughFall = 0.0 or SnowThroughFall > 0.0 and there is a
        substantial amount of intercepted snow at the beginning of the time
        step ( > MIN_INTERCEPTION_STORAGE).  Snow melt may generate mass release. */
        *TempIntStorage += ExcessSnowMelt;
      MassRelease(IntSnow, TempIntStorage, &ReleasedMass, &Drip, MDRatio);
    }
    /* If intercepted snow has melted, add the water it held to drip */
    MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;
    if (*IntRain > MaxWaterInt) {
      Drip += *IntRain - MaxWaterInt;
      *IntRain = MaxWaterInt;
    }
  }
  else {			/* else (RefreezeEnergy <= 0.0) */

                    /* Reset *TempIntStorage to 0.0 when energy balance is negative */
    *TempIntStorage = 0.0;

    /* Refreeze as much surface water as you can */
    if (RefreezeEnergy > -(*IntRain) * LF) {
      *IntSnow += fabs(RefreezeEnergy) / LF;
      *IntRain -= fabs(RefreezeEnergy) / LF;

      *MeltEnergy += (fabs(RefreezeEnergy) * WATER_DENSITY) / Dt;
      RefreezeEnergy = 0.0;
    }
    else {
      /* All of the water in the surface layer has been frozen. */
      *IntSnow += *IntRain;

      /* Energy released by freezing of intercepted water is added
      to the MeltEnergy */
      *MeltEnergy += (LF * *IntRain * WATER_DENSITY) / Dt;
      *IntRain = 0.0;
    }
    if (-(*VaporMassFlux) > *IntSnow) {
      *VaporMassFlux = -(*IntSnow);
      *IntSnow = 0.0;
    }
    else
      *IntSnow += *VaporMassFlux;
  }

  /* Convert drip, mass, *IntSnow, *IntRain, *MeltEnergy and
  *int_vapor_flux from physical depths to pixel depths.  Update p0 and
  snow_fract. */
  *IntSnow *= F;
  *IntRain *= F;
  *MeltEnergy *= F;
  *VaporMassFlux *= F;
  Drip *= F;
  ReleasedMass *= F;

  /* Calculate intercepted H2O balance. */
  MassBalanceError = (InitialWaterInt - (*IntSnow + *IntRain)) + (*SnowFall + *RainFall) -
    (SnowThroughFall + RainThroughFall + Drip + ReleasedMass) + *VaporMassFlux;

  *RainFall = RainThroughFall + Drip;
  *SnowFall = SnowThroughFall + ReleasedMass;

  TIMING_TASK_END("Snow interception", 3);
}
