/*
 * SUMMARY:      MassEnergyBalance.c - Calculate mass and energy balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate mass and energy balance at each pixel
 * DESCRIP-END.
 * FUNCTIONS:    MassEnergyBalance()
 * COMMENTS:
 * $Id: MassEnergyBalance.c,v3.1.2 2013/08/18 ning Exp $     
 */
#ifdef SNOW_ONLY
#define NO_ET
#define NO_SOIL
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "snow.h"
#include "constants.h"
#include "soilmoisture.h"
#include "Calendar.h"

/*****************************************************************************
  Function name: MassEnergyBalance()

  Purpose      : Calculate mass and energy balance

  Required     :

  Returns      : void

  Modifies     :

  Comments     :

  Reference    :
    Epema, G.F. and H.T. Riezbos, 1983, Fall Velocity of waterdrops at different
heights as a factor influencing erosivity of simulated rain. Rainfall simulation,
Runoff and Soil Erosion. Catena suppl. 4, Braunschweig. Jan de Ploey (Ed), 1-17. 

  Laws, J.o., and D.A. Parsons, 1943, the relation of raindrop size to intensity. 
Trans. Am. Geophys. Union, 24: 452-460.

  Wicks, J.M. and J.C. Bathurst, 1996, SHESED: a physically based, distributed
erosion and sediment yield component for the SHE hydrological modeling system,
Journal of Hydrology, 175, 213-238.

*****************************************************************************/
void MassEnergyBalance(OPTIONSTRUCT *Options, int y, int x, float SineSolarAltitude, 
		       float DX, float DY, int Dt, int HeatFluxOption, 
		       int CanopyRadAttOption, int RoadRouteOption,
		       int InfiltOption, int MaxVegLayers, PIXMET *LocalMet,
		       ROADSTRUCT *LocalNetwork, PRECIPPIX *LocalPrecip,
		       VEGTABLE *VType, VEGPIX *LocalVeg, SOILTABLE *SType,
		       SOILPIX *LocalSoil, SNOWPIX *LocalSnow,
		       EVAPPIX *LocalEvap, PIXRAD *TotalRad,
		       CHANNEL *ChannelData, float **skyview)
{
  PIXRAD LocalRad;			/* Radiation balance components (W/m^2) */
  float SurfaceWater;		/* Pixel average depth of water before infiltration is calculated (m) */

  float RoadWater;          /* Average depth of water on the road surface 
							   (normalized by grid cell area) 
							   before infiltration is calculated (m) */
  float ChannelWater;       /* Precip that hits the channel */
  float Infiltration;		/* Infiltration into the top soil layer (m) */
  float Infiltrability;     /* Dynamic infiltration capacity (m/s)*/
  float B;                  /* Capillary drive and soil saturation deficit used 
							   in dynamic infiltration calculation*/
  float LowerRa;		    /* Aerodynamic resistance for lower layer (s/m) */

  float LowerWind;			/* Wind for lower layer (m/s) */
  float MaxInfiltration;	/* Maximum infiltration into the top soil layer (m) */

  float MaxRoadbedInfiltration;	/* Maximum infiltration through the road bed soil layer (m) */

  float MeltEnergy;			/* Energy used to melt snow and  change of cold content of snow pack */

  float MoistureFlux;		/* Amount of water transported from the pixel 
							   to the atmosphere (m/timestep) */
  float NetRadiation;		/* Total Net long- and shortwave radiation for each veg layer (W/m2) */
  float PercArea;           /* Surface area of percolation corrected for 
							   channel and road area, divided by the grid cell area (0-1)  */

  float Reference;			/* Reference height for sensible heat calculation (m) */

  float RoadbedInfiltration;/* Infiltration through the road bed (m) */
  float Roughness;			/* Roughness length (m) */
  float Rp;					/* radiation flux in visible part of the spectrum (W/m^2) */

  float UpperRa;		    /* Aerodynamic resistance for upper layer (s/m) */

  float UpperWind;			/* Wind for upper layer (m/s) */
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */



  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */
  float Tsurf;				/* Surface temperature used in LongwaveBalance() (C) */
  float RainfallIntensity;  /* Rainfall intensity (mm/h) */
  float MS_Rainfall;        /* Momentum squared for rain throughfall((kg* m/s)^2 /(m^2 * s)) */
  int MS_Index;             /* Index for determining alpha and beta cooresponding to RainfallIntensity*/
  int   NVegLActual;		/* Number of vegetation layers above snow */
  float alpha[4]={2.69e-8,3.75e-8,6.12e-8,11.75e-8}; /* empirical coefficient
				for rainfall momentum after Wicks and Bathurst (1996) */ 
  float beta[4]={1.6896,1.5545,1.4242,1.2821};       /* empirical coefficient
				for rainfall momentum after Wicks and Bathurst (1996) */ 
  int i;
  float CanopyHeight[18] = {0.5,1,1.5,2,3,4,5,6,7,8,9,10,11,12,
			    13,14,15,16};        /* Canopy height at which drip fall 
			          velocity is prescribed after Epema and Riezebos (1983) (m)*/           
  float FallVelocity[18] = {2.96,4.12,5.12,5.82,6.84,7.54,8.05,8.36,
			    8.54,8.66,8.75,8.82,8.87,8.91,8.96,9.02,9.07,9.13};
                                     /* Drip fall velocity corresponding to
                                     CanopyHeight after Epema and Riezebos (1983) (m/s)*/
  float LD_FallVelocity;             /* Leaf drip fall velocity corresponding to the
                                     canopy height in vegetation map (m/s) */

  /* Calculate the number of vegetation layers above the snow */

  NVegLActual = VType->NVegLayers;
  if (LocalSnow->HasSnow == TRUE && VType->UnderStory == TRUE)
    --NVegLActual;

  /* initialize the total amount of evapotranspiration, and MeltEnergy */

  LocalEvap->ETot = 0.0;
  MeltEnergy = 0.0;
  MoistureFlux = 0.0;

  /* calculate the radiation balance for the ground/snow surface and the
     vegetation layers above that surface */
  RadiationBalance(Options, HeatFluxOption, CanopyRadAttOption, SineSolarAltitude, 
	       LocalMet->VICSin, LocalMet->Sin, LocalMet->SinBeam, LocalMet->SinDiffuse, 
		   LocalMet->Lin, LocalMet->Tair, LocalVeg->Tcanopy, 
		   LocalSoil->TSurf, SType->Albedo, VType, LocalSnow, &LocalRad);

  /* calculate the actual aerodynamic resistances and wind speeds */
  UpperWind = VType->U[0] * LocalMet->Wind;
  UpperRa = VType->Ra[0] / LocalMet->Wind;
  if (VType->OverStory == TRUE) {
    LowerWind = VType->U[1] * LocalMet->Wind;
    LowerRa = VType->Ra[1] / LocalMet->Wind;
  }
  else {
    LowerWind = UpperWind;
    LowerRa = UpperRa;
  }

  /* Leaf drip impact*/
  /* Find corresponding fall velocity for overstory and understory heights
     by weighting scheme */

  if (VType->OverStory){
    /* staring at 1 assumes the overstory height > 0.5 m */
    for (i = 1; i <= 17; i++ ) {           
      if (VType->Height[0] < CanopyHeight[i]) {
        LD_FallVelocity = ((VType->Height[0] - CanopyHeight[i-1])
			   *FallVelocity[i] +
			   (CanopyHeight[i] - VType->Height[0])*FallVelocity[i-1]) /
	  (CanopyHeight[i] - CanopyHeight[i-1]);
      }
    }
    if (VType->UnderStory) {                 
      /* ending at 16 assumes the understory height < 16 m */
      for (i = 0; i <= 16; i++) {
	if (VType->Height[1] < CanopyHeight[i]) {
	  LD_FallVelocity = ((VType->Height[1] - CanopyHeight[i])*FallVelocity[i] +
			     (CanopyHeight[i+1] - VType->Height[1])*FallVelocity[i-1]) /
	    (CanopyHeight[i+1] - CanopyHeight[i]);
	}
      }
    }
  }
  else if (VType->UnderStory) {                 
    /* ending at 16 assumes the understory height < 16 m */
    for (i = 0; i <= 16; i++) {
      if (VType->Height[0] < CanopyHeight[i]) {
	LD_FallVelocity = ((VType->Height[0] - CanopyHeight[i])*FallVelocity[i] +
			   (CanopyHeight[i+1] - VType->Height[0])*FallVelocity[i-1]) /
	  (CanopyHeight[i+1] - CanopyHeight[i]);
      }
    }
  }
  else LD_FallVelocity = 0;

  
  /* RainFall impact */
  /* 3600 is conversion factor (number of seconds per hour) */
  if (LocalPrecip->RainFall > 0.) {
    RainfallIntensity = LocalPrecip->RainFall * (1./MMTOM) * (3600./Dt);
    
    /* Momentum is later weighted with the overstory/understory fraction */
    if (RainfallIntensity < 10.)
      MS_Index = 0;
    else if (RainfallIntensity >= 10. && RainfallIntensity < 100.) 
      MS_Index = floor((RainfallIntensity + 49)/50);
    else 
      MS_Index = 3;
    
    /* Eq. 1, Wicks and Bathurst (1996) */
    MS_Rainfall = alpha[MS_Index] * pow(RainfallIntensity, beta[MS_Index]);
    
    /* Calculating mediam raindrop diameter after Laws and Parsons (1943) */
    LocalPrecip->Dm =  0.00124 * pow((double)RainfallIntensity, 0.182); 
  }
  else {
    MS_Rainfall = 0;
    LocalPrecip->Dm = LEAF_DRIP_DIA;
  }

  /* calculate the amount of interception storage, and the amount of 
     throughfall.  Of course this only needs to be done if there is
     vegetation present. */

#ifndef NO_SNOW

  if (VType->OverStory == TRUE &&
      (LocalPrecip->IntSnow[0] || LocalPrecip->SnowFall > 0.0)) {
    SnowInterception(y, x, Dt, VType->Fract[0], VType->LAI[0],
		     VType->MaxInt[0], VType->MaxSnowInt, VType->MDRatio,
		     VType->SnowIntEff, UpperRa, LocalMet->AirDens,
		     LocalMet->Eact, LocalMet->Lv, &LocalRad, LocalMet->Press,
		     LocalMet->Tair, LocalMet->Vpd, UpperWind,
		     &(LocalPrecip->RainFall), &(LocalPrecip->SnowFall),
		     &(LocalPrecip->IntRain[0]), &(LocalPrecip->IntSnow[0]),
		     &(LocalPrecip->TempIntStorage),
		     &(LocalSnow->CanopyVaporMassFlux), &(LocalVeg->Tcanopy),
		     &MeltEnergy, &(LocalPrecip->MomentSq), VType->Height, 
		     VType->UnderStory, MS_Rainfall, LD_FallVelocity);

    MoistureFlux -= LocalSnow->CanopyVaporMassFlux;

    /* Because we now have a new estimate of the canopy temperature we can
       recalculate the longwave balance */
    if (LocalSnow->HasSnow == TRUE)
      Tsurf = LocalSnow->TSurf;
    else if (HeatFluxOption == TRUE)
      Tsurf = LocalSoil->TSurf;
    else
      Tsurf = LocalMet->Tair;
    LongwaveBalance(Options, VType->OverStory, VType->Fract[0], LocalMet->Lin, 
	               LocalVeg->Tcanopy, Tsurf, &LocalRad);
  }
  else if (VType->NVegLayers > 0) {
    LocalVeg->Tcanopy = LocalMet->Tair;
    LocalSnow->CanopyVaporMassFlux = 0.0;
    LocalPrecip->TempIntStorage = 0.0;
    InterceptionStorage(VType->NVegLayers, NVegLActual, VType->MaxInt,
			VType->Fract, LocalPrecip->IntRain,
			&(LocalPrecip->RainFall), &(LocalPrecip->MomentSq), 
			VType->Height, VType->UnderStory, Dt, MS_Rainfall,
			LD_FallVelocity);
  }
  else {
    /* If no vegetation, kinetic energy is all due to direct precipitation. */
    if(LocalPrecip->RainFall > 0.0)
    LocalPrecip->MomentSq = MS_Rainfall;
  }

  /* If snow on the ground, assume no overland flow erosion. */
  if(LocalSnow->HasSnow)
    LocalPrecip->MomentSq = 0.0;
						    
  /* if snow is present, simulate the snow pack dynamics */

  if (LocalSnow->HasSnow || LocalPrecip->SnowFall > 0.0) {

    if (VType->OverStory == TRUE) {
      SnowLongIn = LocalRad.LongIn[1];
      SnowNetShort = LocalRad.NetShort[1];
    }
    else {
      SnowLongIn = LocalRad.LongIn[0];
      SnowNetShort = LocalRad.NetShort[0];
    }

    SnowWind = VType->USnow * LocalMet->Wind;
    SnowRa = VType->RaSnow / LocalMet->Wind;

    LocalSnow->Outflow =
      SnowMelt(y, x, Dt, 2. + Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
	       LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
	       LocalMet->Press, LocalPrecip->RainFall, LocalPrecip->SnowFall,
	       LocalMet->Tair, LocalMet->Vpd, SnowWind,
	       &(LocalSnow->PackWater), &(LocalSnow->SurfWater),
	       &(LocalSnow->Swq), &(LocalSnow->VaporMassFlux),
	       &(LocalSnow->TPack), &(LocalSnow->TSurf), &MeltEnergy);

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */

    LocalPrecip->RainFall = 0.0;
    MoistureFlux -= LocalSnow->VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
       can recalculate the longwave balance */

    Tsurf = LocalSnow->TSurf;
    LongwaveBalance(Options, VType->OverStory, VType->Fract[0], LocalMet->Lin,
		    LocalVeg->Tcanopy, Tsurf, &LocalRad);
  }
  else {
    LocalSnow->Outflow = 0.0;
    LocalSnow->VaporMassFlux = 0.0;
  }

  /* Determine whether a snow pack is still present, or whether everything
     has melted */

  if (LocalSnow->Swq > 0.0)
    LocalSnow->HasSnow = TRUE;
  else
    LocalSnow->HasSnow = FALSE;

  /*do the glacier add */
  if (LocalSnow->Swq < 1.0 && VType->Index == GLACIER) {
    printf("resetting glacier swe of %f to 5.0 meters\n", LocalSnow->Swq);
    LocalSnow->Glacier += (5.0 - LocalSnow->Swq);
    LocalSnow->Swq = 5.0;
    LocalSnow->TPack = 0.0;
    LocalSnow->TSurf = 0.0;
  }
#endif

#ifndef NO_ET
  /* calculate the amount of evapotranspiration from each vegetation layer 
     above the ground/soil surface.  Also calculate the total amount of 
     evapotranspiration from the vegetation */

  if (VType->OverStory == TRUE) {
    Rp = VISFRACT * LocalRad.NetShort[0];
    NetRadiation =
      LocalRad.NetShort[0] +
      LocalRad.LongIn[0] - 2 * VType->Fract[0] * LocalRad.LongOut[0];
    EvapoTranspiration(0, Dt, LocalMet, NetRadiation, Rp, VType, SType,
		       MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[0]),
		       LocalEvap, LocalNetwork->Adjust, UpperRa);
    MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];

    if (LocalSnow->HasSnow != TRUE && VType->UnderStory == TRUE) {
      Rp = VISFRACT * LocalRad.NetShort[1];
      NetRadiation =
	LocalRad.NetShort[1] +
	LocalRad.LongIn[1] - VType->Fract[1] * LocalRad.LongOut[1];
      EvapoTranspiration(1, Dt, LocalMet, NetRadiation, Rp, VType, SType,
			 MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[1]),
			 LocalEvap, LocalNetwork->Adjust, LowerRa);
      MoistureFlux += LocalEvap->EAct[1] + LocalEvap->EInt[1];
    }
    else if (VType->UnderStory == TRUE) {
      LocalEvap->EAct[1] = 0.;
      LocalEvap->EInt[1] = 0.;
    }
  }				/* end if(VType->OverStory == TRUE) */
  else if (LocalSnow->HasSnow != TRUE && VType->UnderStory == TRUE) {
    Rp = VISFRACT * LocalRad.NetShort[0];
    NetRadiation =
      LocalRad.NetShort[0] +
      LocalRad.LongIn[0] - VType->Fract[0] * LocalRad.LongOut[0];
    EvapoTranspiration(0, Dt, LocalMet, NetRadiation, Rp, VType, SType,
		       MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[0]),
		       LocalEvap, LocalNetwork->Adjust, LowerRa);
    MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];
  }
  else if (VType->UnderStory == TRUE) {
    LocalEvap->EAct[0] = 0.;
    LocalEvap->EInt[0] = 0.;
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is 
     present and there is no understory */

  if (LocalSnow->HasSnow != TRUE && VType->UnderStory != TRUE) {

    if (VType->OverStory == TRUE)
      NetRadiation =
	LocalRad.NetShort[1] + LocalRad.LongIn[1] - LocalRad.LongOut[1];
    else
      NetRadiation =
	LocalRad.NetShort[0] + LocalRad.LongIn[0] - LocalRad.LongOut[0];

    LocalEvap->EvapSoil =
      SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
		      LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
		      NetRadiation, LowerRa, MoistureFlux, SType->Porosity[0],
		      SType->Ks[0], SType->Press[0], SType->PoreDist[0],
		      VType->RootDepth[0], &(LocalSoil->Moist[0]),
		      LocalNetwork->Adjust[0]);
  }
  else
    LocalEvap->EvapSoil = 0.0;

  MoistureFlux += LocalEvap->EvapSoil;
  LocalEvap->ETot += LocalEvap->EvapSoil;

#endif

  /* add the water that was not intercepted to the upper soil layer */

#ifndef NO_SOIL

  /* This has been modified so that PercArea for infiltration is calculated
     locally to account for the fact that some cells have roads and streams.
     I am not sure if the old PercArea (which does not account for the fact 
     that some cells have roads and streams) needs to remain the same. 
     Currently, the old PercArea is passed to UnsaturatedFlow */
  
  MaxRoadbedInfiltration = 0.;
  MaxInfiltration = 0.;
  ChannelWater = 0.;
  RoadWater = 0.;
  SurfaceWater = 0.;
  PercArea = 1.;
  RoadbedInfiltration = 0.;
  
  /* ChannelWater is precipitation falling on the channel */
  /* (if there is no road, LocalNetwork->RoadArea = 0) */
  if (channel_grid_has_channel(ChannelData->stream_map, x, y)){
    PercArea = 1. - (LocalNetwork->Area + LocalNetwork->RoadArea)/(DX*DY);
    ChannelWater = LocalNetwork->Area/(DX*DY) * LocalPrecip->RainFall;
  }
  /* If there is a road and no channel, the PercArea is 
     based on the road only */
  else if (channel_grid_has_channel(ChannelData->road_map, x, y)){
    PercArea = 1. - (LocalNetwork->RoadArea)/(DX*DY);
    MaxRoadbedInfiltration = (1. - PercArea) * 
      LocalNetwork->MaxInfiltrationRate * Dt; 
  }
  
  /* SurfaceWater is rain falling on the hillslope + 
     snowmelt on the hillslope (there is no snowmelt on the channel) +
     existing IExcess */
  SurfaceWater = (PercArea * LocalPrecip->RainFall) +
    ((1. - (LocalNetwork->RoadArea)/(DX*DY)) * LocalSnow->Outflow) + 
    LocalSoil->IExcess;
  
  /* RoadWater is rain falling on the road surface +
     snowmelt on the road surface + existing Road IExcess 
     (Existing road IExcess = 0). WORK IN PROGRESS*/
  RoadWater = (LocalNetwork->RoadArea/(DX*DY) * 
	       (LocalPrecip->RainFall + LocalSnow->Outflow)) + 
    LocalNetwork->IExcess;
  
  
  if(InfiltOption == STATIC)
    MaxInfiltration = (1. - VType->ImpervFrac) * PercArea * SType->MaxInfiltrationRate * Dt; 
  
  else { /* InfiltOption == DYNAMIC 
	    Dynamic Infiltration Capacity after Parlange and Smith 1978, 
	    as used in KINEROS and THALES */
    Infiltration = 0.0;
    
    if (SurfaceWater > 0.) {
      /* Infiltration is a function of the amount of water infiltrated since 
	 the storm started */
      if (LocalPrecip->PrecipStart){
	LocalSoil->MoistInit = LocalSoil->Moist[0];
	LocalSoil->InfiltAcc = 0.0;
      }
      
      /* Check that the B parameter > 0 */
      if ((LocalSoil->InfiltAcc > 0.) && (SType->Porosity[0] > LocalSoil->MoistInit)) {
	
	B = (SType->Porosity[0] - LocalSoil->MoistInit) * 
	  (SType->G_Infilt + SurfaceWater);
	Infiltrability = SType->Ks[0] * exp((LocalSoil->InfiltAcc)/B) / 
	  (exp((LocalSoil->InfiltAcc)/B) - 1.);
      }
      
      else 
	Infiltrability = SurfaceWater/Dt ; 
      
      MaxInfiltration = Infiltrability * PercArea *
	(1. - VType->ImpervFrac) *  Dt;
      
      LocalPrecip->PrecipStart = FALSE; 
    }/* end  if (SurfaceWater > 0.) */
    else 
      LocalPrecip->PrecipStart = TRUE;
    
  } /* end Dynamic MaxInfiltration calculation */ 
  
  Infiltration = (1. - VType->ImpervFrac) * SurfaceWater;
  
  if (Infiltration > MaxInfiltration) 
    Infiltration = MaxInfiltration;
  
  RoadbedInfiltration = RoadWater;
  
  if (RoadbedInfiltration > MaxRoadbedInfiltration) 
    RoadbedInfiltration = MaxRoadbedInfiltration;
  
  if (RoadRouteOption == FALSE)
    LocalSoil->IExcess = SurfaceWater - Infiltration + 
      RoadWater - RoadbedInfiltration;
  else {
    LocalSoil->IExcess = SurfaceWater - Infiltration;
    LocalNetwork->IExcess = RoadWater - RoadbedInfiltration;
    if (LocalNetwork->IExcess < 0.){
      LocalNetwork->IExcess = 0.;
      printf("MEB: NetIExcess(%f), reset to 0\n", LocalNetwork->IExcess);
    }
  }

  if (LocalSoil->IExcess < 0.){
    printf("MEB: SoilIExcess(%f), reset to 0\n", LocalSoil->IExcess);
    LocalSoil->IExcess = 0.;
  }

  /*Add water that hits the channel network to the channel network */
  if (ChannelWater > 0.){
    channel_grid_inc_inflow(ChannelData->stream_map, x, y, ChannelWater * DX * DY);
    LocalSoil->ChannelInt += ChannelWater;
   }
  
  /* Calculate unsaturated soil water movement, and adjust soil water 
     table depth */

  UnsaturatedFlow(Dt, DX, DY, Infiltration, RoadbedInfiltration,
		  LocalSoil->SatFlow, SType->NLayers,
		  LocalSoil->Depth, LocalNetwork->Area, VType->RootDepth,
		  SType->Ks, SType->PoreDist, SType->Porosity, SType->FCap,
		  LocalSoil->Perc, LocalNetwork->PercArea,
		  LocalNetwork->Adjust, LocalNetwork->CutBankZone,
		  LocalNetwork->BankHeight, &(LocalSoil->TableDepth),
		  &(LocalSoil->IExcess), LocalSoil->Moist, RoadRouteOption,
		  InfiltOption, &(LocalNetwork->IExcess));
  
  /* Infiltration is updated in UnsaturatedFlow and accumulated 
     below */
  if ((InfiltOption == DYNAMIC) && (SurfaceWater > 0.)) 
    LocalSoil->InfiltAcc += Infiltration; 
  
  if (HeatFluxOption == TRUE) {
    if (LocalSnow->HasSnow == TRUE) {
      Reference = 2. + Z0_SNOW;
      Roughness = Z0_SNOW;
    }
    else {
      Reference = 2. + Z0_GROUND;
      Roughness = Z0_GROUND;
    }

    SensibleHeatFlux(y, x, Dt, LowerRa, Reference, 0.0f, Roughness,
		     LocalMet, LocalRad.PixelNetShort, LocalRad.PixelLongIn,
		     MoistureFlux, SType->NLayers, VType->RootDepth,
		     SType, MeltEnergy, LocalSoil);
    Tsurf = LocalSoil->TSurf;
    LongwaveBalance(Options, VType->OverStory, VType->Fract[0], LocalMet->Lin,
		    LocalVeg->Tcanopy, Tsurf, &LocalRad);
  }
  else
    NoSensibleHeatFlux(Dt, LocalMet, MoistureFlux, LocalSoil);

#endif

  /* add the components of the radiation balance for the current pixel to 
     the total */
  AggregateRadiation(MaxVegLayers, VType->NVegLayers, &LocalRad, TotalRad);
  /* For RBM model, save the energy fluxes for outputs */
  if (Options->StreamTemp) {
    if (channel_grid_has_channel(ChannelData->stream_map, x, y))
      channel_grid_inc_other(ChannelData->stream_map, x, y, &LocalRad, LocalMet, skyview[y][x]);
  }
}
