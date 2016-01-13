/*
 * SUMMARY:      SensibleHeatFlux.c - Calculate sensible heat flux
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate sensible heat flux
 * DESCRIP-END.
 * FUNCTIONS:    SensibleHeatFlux()
 *               NoSensibleHeatFlux()
 * COMMENTS:
 * $Id: SensibleHeatFlux.c,v 1.4 2003/07/01 21:26:24 olivier Exp $     
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "brent.h"
#include "functions.h"

/*****************************************************************************
  SensibleHeatFlux()
*****************************************************************************/
void SensibleHeatFlux(int y, int x, int Dt, float Ra, float ZRef,
		      float Displacement, float Z0, PIXMET * LocalMet,
		      float NetShort, float LongIn, float ETot,
		      int NSoilLayers, float *SoilDepth, SOILTABLE * SoilType,
		      float MeltEnergy, SOILPIX * LocalSoil)
{
  float FluxDepth;		/* Lower boundary for soil heat flux (m) */
  float HeatCapacity;		/* Soil heat capacity */
  float MaxTSurf;		/* Upper bracket for effective surface 
				   temperature (C) */
  float MinTSurf;		/* Lower bracket for effective surface
				   temperature (C) */
  float OldTSurf;		/* Effective surface temperature at the end of
				   the last timestep (C) */
  float KhEff;			/* Thermal conductivity of the soil (W/(m*C)
				 */
  float TMean;			/* Average surface temperature (C) */
  float TSoilLower;		/* Temperature os the soil at FluxDepth (C)
				 */
  float TSoilUpper;		/* Temperature os the soil in top layer (C)
				 */
  double Tmp;			/* Temporary value */

  OldTSurf = LocalSoil->TSurf;
  MaxTSurf = 0.5 * (LocalSoil->TSurf + LocalMet->Tair) + DELTAT;
  MinTSurf = 0.5 * (LocalSoil->TSurf + LocalMet->Tair) - DELTAT;

  /* WORK IN PROGRESS */
  /* the lower boundary for the soil heat flux is currently fixed at a depth 
     of 1.0 m */

  FluxDepth = 1.0;
  TSoilLower = LocalSoil->Temp[NSoilLayers - 1];
  TSoilUpper = LocalSoil->Temp[0];

  /* Calculate the effective thermal conductivity of the soil above between
     FluxDepth and DZ_TOP */

  KhEff = CalcEffectiveKh(NSoilLayers, DZ_TOP, FluxDepth, SoilDepth,
			  SoilType->KhDry, SoilType->KhSol, LocalSoil->Moist,
			  SoilType->Porosity, LocalSoil->Temp);

  /*   KhEff = 1; */

  /* Calculate the effective surface temperature that makes sure that the 
     sum of the terms of the energy balance equals 0 */

  LocalSoil->TSurf =
    RootBrent(y, x, MinTSurf, MaxTSurf, SurfaceEnergyBalance, Dt, Ra, ZRef,
	      Displacement, Z0, LocalMet->Wind, NetShort, LongIn,
	      LocalMet->AirDens, LocalMet->Lv, ETot, KhEff,
	      SoilType->Ch[0], SoilType->Porosity[0], LocalSoil->Moist[0],
	      FluxDepth, LocalMet->Tair, TSoilUpper,
	      TSoilLower, OldTSurf, MeltEnergy);

  /* Calculate the terms of the energy balance.  This is similar to the
     code in SurfaceEnergyBalance.c */

  TMean = 0.5 * (OldTSurf + LocalSoil->TSurf);

  if (LocalMet->Wind > 0.0)
    Ra /= StabilityCorrection(ZRef, Displacement, TMean, LocalMet->Tair,
			      LocalMet->Wind, Z0);
  else
    Ra = DHSVM_HUGE;

  LocalSoil->Ra = Ra;

  Tmp = TMean + 273.15;
  LocalSoil->Qnet = NetShort + LongIn - STEFAN * (Tmp * Tmp * Tmp * Tmp);

  LocalSoil->Qs = LocalMet->AirDens * CP * (LocalMet->Tair - TMean) / Ra;

  LocalSoil->Qe = -(LocalMet->Lv * ETot) / Dt * WATER_DENSITY;

  LocalSoil->Qg = KhEff * (TSoilLower - TMean) / FluxDepth;

  HeatCapacity = (1 - SoilType->Porosity[0]) * SoilType->Ch[0];

  if (TSoilUpper >= 0.0)
    HeatCapacity += LocalSoil->Moist[0] * CH_WATER;
  else
    HeatCapacity += LocalSoil->Moist[0] * CH_ICE;

  LocalSoil->Qst = (HeatCapacity * (OldTSurf - TMean) * DZ_TOP) / Dt;

  LocalSoil->Qrest = LocalSoil->Qnet + LocalSoil->Qs + LocalSoil->Qe +
    LocalSoil->Qg + LocalSoil->Qst + MeltEnergy;
}

/*****************************************************************************
  Function name: NoSensibleHeatFlux()

  Purpose      : Calculate latent heat flux in W/m2

  Required     : 
    int Dt             - Model timestep (seconds)
    PIXMET LocalMet    - Met data for the current pixel 
    float ETot         - Total vapor flux (mm/timestep) 
    SOILPIX *LocalSoil - Structure with soil moisture data and energy data
                         for the current pixel

  Returns      :
    void

  Modifies     :
    members of LocalSoil

  Comments     : This function sets all the energy fluxes at the pixel level 
                 to 0.0, but calculates the evapotranspiration in W/m2
*****************************************************************************/
void NoSensibleHeatFlux(int Dt, PIXMET * LocalMet, float ETot,
			SOILPIX * LocalSoil)
{
  LocalSoil->TSurf = 0.0;

  LocalSoil->Ra = 0.0;

  LocalSoil->Qnet = 0.0;

  LocalSoil->Qs = 0.0;

  LocalSoil->Qe = -(LocalMet->Lv * ETot) / Dt * WATER_DENSITY;

  LocalSoil->Qg = 0.0;

  LocalSoil->Qst = 0.0;

  LocalSoil->Qrest = 0.0;
}
