/*
* SUMMARY:      CanopyGapEnergyBalance.c
* USAGE:        Part of DHSVM
*
* AUTHOR:       Ning Sun
* E-MAIL:       ning.sun@pnnl.gov
* ORIG-DATE:    Jul-18
* DESCRIPTION:  Calculate snow balance under an idealized cylindrical canopy
*               gap/openning
* DESCRIP-END.
* FUNCTIONS:    CanopyGapSnow()
* Reference:
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "constants.h"
#include "snow.h"

/********************************************************************************
Function Name: CanopyGapInterception()

Purpose      : Calculate snow/rain interception
Returns      :
Comments     :
********************************************************************************/
void CanopyGapInterception(OPTIONSTRUCT *Options, CanopyGapStruct **Gap,
  int HeatFluxOption, int y, int x, int Dt, int NVegLActual,
  float DX, float DY, float UpperRa, float UpperWind, VEGTABLE *VType,
  SOILPIX *LocalSoil, VEGPIX *LocalVeg, SNOWPIX *LocalSnow,
  PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, PIXMET *LocalMet) {

  if ((*Gap)[Opening].UnderStory == TRUE) {
    (*Gap)[Opening].Tcanopy = LocalMet->Tair;
    (*Gap)[Opening].CanopyVaporMassFlux = 0.0;
    (*Gap)[Opening].TempIntStorage = 0.0;

    /* calculate rain interception by the gap*/
    CanopyGapInterceptionStorage((*Gap)[Opening].NVegLActual, VType->MaxInt,
      LocalVeg->Fract, (*Gap)[Opening].IntRain, &((*Gap)[Opening].RainFall));
  }
}

/********************************************************************************
Function Name: CanopyGapSnowMelt()

Purpose      : Calculate snow melt
Returns      :
Comments     :
********************************************************************************/
void CanopyGapSnowMelt(OPTIONSTRUCT *Options, int y, int x, int Dt,
  CanopyGapStruct **Gap, float DX, float DY, VEGTABLE *VType, VEGPIX *LocalVeg,
  SNOWPIX *LocalSnow, PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, PIXMET *LocalMet)
{
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */
  float Tsurf;              /* Surface temperature */
  float OldTSurf;           /* Effective surface temperature at the end of the last timestep (C)*/
  float Tmean;              /* Average snow surface temperature*/
  float tmp;                /* temporary variable */
  float Tmp;                /* temporary variable */
  float Ls;			        /* Latent heat of sublimation (J/kg) */

  /********************* calculate gap snow melt *********************/
  if ((*Gap)[Opening].HasSnow || (*Gap)[Opening].SnowFall > 0.0) {
    SnowLongIn = (*Gap)[Opening].LongIn[1];
    SnowNetShort = (*Gap)[Opening].NetShort[1];

    SnowWind = (*Gap)[Opening].USnow * LocalMet->Wind;  
    SnowRa = (*Gap)[Opening].RaSnow / LocalMet->Wind; 
    
    /* adjust the wind and Ra values for gap so they fall betwene open 
    and forested values */
    tmp = VType->USnow*LocalMet->Wind; //forested snow wind
    SnowWind = tmp + (SnowWind- tmp)*GAPWIND_FACTOR;

    tmp = VType->RaSnow/LocalMet->Wind;
    SnowRa = tmp - (tmp-SnowRa)*GAPWIND_FACTOR;
    
    OldTSurf = (*Gap)[Opening].TSurf;
    (*Gap)[Opening].SnowPackOutflow =
      SnowMelt(y, x, Dt, 2.+Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, (*Gap)[Opening].RainFall, (*Gap)[Opening].SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &((*Gap)[Opening].PackWater), &((*Gap)[Opening].SurfWater),
        &((*Gap)[Opening].Swq), &((*Gap)[Opening].VaporMassFlux),
        &((*Gap)[Opening].TPack), &((*Gap)[Opening].TSurf), &((*Gap)[Opening].MeltEnergy));

    /* Calculate the terms of the snow energy balance.  This is similar to the
    code in SnowPackEnergyBalance.c */
    Tmean = 0.5 * (OldTSurf + (*Gap)[Opening].TSurf);

    /* Apply the stability correction to the aerodynamic resistance */
    if (SnowWind > 0.0)
      SnowRa /= StabilityCorrection(2.0f, 0.f, Tmean, LocalMet->Tair, SnowWind, Z0_SNOW);
    else
      SnowRa = DHSVM_HUGE;

    /* convert snow surface temperature from C to K */
    Tmp = Tmean + 273.15;
    /* net shortwave radiation */
    (*Gap)[Opening].Qsw = SnowNetShort;
    /* net longwave radiation */
	(*Gap)[Opening].Qlin = SnowLongIn;
    (*Gap)[Opening].Qlw = SnowLongIn - STEFAN * (Tmp * Tmp * Tmp * Tmp);
    /* sensible heat */
    (*Gap)[Opening].Qs = LocalMet->AirDens * CP * (LocalMet->Tair - Tmean) / SnowRa;

    /* Calculate latent heat flux */
    if (Tmean >= 0.0) {
      /* Melt conditions: use latent heat of vaporization */
      (*Gap)[Opening].Qe = LocalMet->Lv * (*Gap)[Opening].VaporMassFlux * WATER_DENSITY;
    }
    else {
      /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
      Ls = (677. - 0.07 * Tmean) * JOULESPCAL * GRAMSPKG;
      (*Gap)[Opening].Qe = Ls * (*Gap)[Opening].VaporMassFlux * WATER_DENSITY;
    }
    (*Gap)[Opening].Qe /= Dt;

    /* Calculate advected heat flux from rain */
    (*Gap)[Opening].Qp = (CH_WATER * LocalMet->Tair * (*Gap)[Opening].RainFall) / Dt;

    /* end of snow energy tems */

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    (*Gap)[Opening].RainFall = 0.0;
    (*Gap)[Opening].MoistureFlux -= (*Gap)[Opening].VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
    can recalculate the longwave balance */
    Tsurf = Gap[Opening]->TSurf;
    CanopyGapLongRadiation(&((*Gap)[Opening]), VType->Height[0],
      LocalVeg->Gapping, LocalMet->Lin, LocalVeg->Tcanopy, LocalVeg->Fract[0]);
    (*Gap)[Opening].LongIn[0] = 0.;
  }
  else {
    (*Gap)[Opening].SnowPackOutflow = 0.0;
    (*Gap)[Opening].VaporMassFlux = 0.0;

	/* snow is all melt */
	(*Gap)[Opening].Qs = 0.;
	(*Gap)[Opening].Qe = 0.;
	(*Gap)[Opening].Qp = 0.;
	(*Gap)[Opening].Qsw = 0;
	(*Gap)[Opening].Qlin = 0;
	(*Gap)[Opening].Qlw = 0; /* not used in calculation */
	(*Gap)[Opening].MeltEnergy = 0;
  }

  if ((*Gap)[Opening].Swq > 0.0)
    (*Gap)[Opening].HasSnow = TRUE;
  else
    (*Gap)[Opening].HasSnow = FALSE;
}

/*****************************************************************************
Function name: CalcCanopyGapAerodynamic()

Purpose      : Calculate the aerodynamic resistance for each vegetation
layer, and the wind 2m above the layer boundary.
*****************************************************************************/
void CalcCanopyGapAerodynamic(CanopyGapStruct **Gap, int NVegLayers,
  float *Height)
{
  float K2;
  float Z0_Lower;    /* roughness length for understory (m) */
  float d_Lower;     /* displacement height (m) */

  K2 = VON_KARMAN * VON_KARMAN;

  /* bare soil */
  if ((*Gap)[Opening].UnderStory == FALSE) {
    Z0_Lower = Z0_GROUND;
    d_Lower = 0;
  }
  /* understory present*/
  else {
    Z0_Lower = Z0_MULTIPLIER * Height[1];
    d_Lower = D0_MULTIPLIER * Height[1];
  }

  /* No snow: get wind speed & aerodynamic resistence value */
  (*Gap)[Opening].U[1] =
    log((2.+Z0_Lower)/Z0_Lower) / log((Zref-d_Lower)/Z0_Lower);
  (*Gap)[Opening].Ra[1] =
    log((2.+Z0_Lower)/Z0_Lower) * log((Zref-d_Lower)/Z0_Lower)/K2;

  /* get the wind speed and aerodynamic resistence value for snow surface*/
  (*Gap)[Opening].USnow = log((2.+Z0_SNOW)/Z0_SNOW) / log(Zref/Z0_SNOW);
  (*Gap)[Opening].RaSnow = log((2.+Z0_SNOW)/Z0_SNOW) * log(Zref/Z0_SNOW)/K2;
}

/*****************************************************************************
Function name: CalcCanopyGapET()

Purpose      : Calculate the ET
*****************************************************************************/
void CalcCanopyGapET(CanopyGapStruct **Gap, int NSoil, VEGTABLE *VType,
  VEGPIX *LocalVeg, SOILTABLE *SType, SOILPIX *LocalSoil, PIXMET *LocalMet,
  EVAPPIX *LocalEvap, ROADSTRUCT *LocalNetwork, int Dt, float UpperRa,
  float LowerRa)
{
  float NetRadiation;		/* Total Net long- and shortwave radiation (W/m2) */
  float Rp;					/* radiation flux in visible part of the spectrum (W/m^2) */

  /********************** for gap **********************/
  if ((*Gap)[Opening].HasSnow != TRUE && VType->UnderStory == TRUE) {
    Rp = VISFRACT * (*Gap)[Opening].NetShort[1];
    NetRadiation =
      (*Gap)[Opening].NetShort[1] +
      (*Gap)[Opening].LongIn[1] - (*Gap)[Opening].LongOut[1];
      
    EvapoTranspiration(1, 1, Dt, 1, LocalMet, NetRadiation,
      Rp, VType, SType, (*Gap)[Opening].MoistureFlux, (*Gap)[Opening].Moist,
      LocalSoil->Temp, &((*Gap)[Opening].IntRain[0]),
      (*Gap)[Opening].EPot, (*Gap)[Opening].EInt, (*Gap)[Opening].ESoil,
      (*Gap)[Opening].EAct, &((*Gap)[Opening].ETot), LocalNetwork->Adjust, LowerRa, LocalVeg);

    (*Gap)[Opening].MoistureFlux += (*Gap)[Opening].EAct[1] + (*Gap)[Opening].EInt[1];

    (*Gap)[Opening].NetRadiation[1] = NetRadiation;
    (*Gap)[Opening].NetRadiation[0] = 0.;
  }
  else if (VType->UnderStory == TRUE) {
    (*Gap)[Opening].EAct[1] = 0.;
    (*Gap)[Opening].EInt[1] = 0.;
    (*Gap)[Opening].NetRadiation[0] = 0;
    (*Gap)[Opening].NetRadiation[1] = 0.;
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is
  present and there is no understory */
  if ((*Gap)[Opening].HasSnow != TRUE && VType->UnderStory != TRUE) {
    NetRadiation =
      (*Gap)[Opening].NetShort[1] + (*Gap)[Opening].LongIn[1] - (*Gap)[Opening].LongOut[1];
    (*Gap)[Opening].NetRadiation[1] = NetRadiation;
    (*Gap)[Opening].NetRadiation[0] = 0.;
    (*Gap)[Opening].EvapSoil =
    SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
        LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
        NetRadiation, UpperRa, (*Gap)[Opening].MoistureFlux, LocalSoil->Porosity[0],
        LocalSoil->FCap[0], SType->Ks[0], SType->Press[0], SType->PoreDist[0],
        VType->RootDepth[0], &((*Gap)[Opening].Moist[0]),
        LocalNetwork->Adjust[0]);
  }
  else
    (*Gap)[Opening].EvapSoil = 0.0;

  (*Gap)[Opening].MoistureFlux += (*Gap)[Opening].EvapSoil;
  (*Gap)[Opening].ETot += (*Gap)[Opening].EvapSoil;
}

/*****************************************************************************
Function name: CalcGapSurroudingIntercept()

Purpose      : .
*****************************************************************************/
void CalcGapSurroudingIntercept(OPTIONSTRUCT *Options, int HeatFluxOption,
  int y, int x, int Dt, int NVegLActual, CanopyGapStruct **Gap, VEGTABLE *VType,
  PIXRAD *LocalRad, PIXMET *LocalMet, float UpperRa, float UpperWind, VEGPIX *LocalVeg)
{
  float Tsurf;
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */

  if (((*Gap)[Forest].IntSnow[0] || (*Gap)[Forest].SnowFall > 0.0)) {
    SnowInterception(Options, y, x, Dt, LocalVeg->Fract[0], LocalVeg->Vf,
      LocalVeg->LAI[0], LocalVeg->MaxInt[0], VType->MaxSnowInt, VType->MDRatio,
      VType->SnowIntEff, UpperRa, LocalMet->AirDens,
      LocalMet->Eact, LocalMet->Lv, LocalRad, LocalMet->Press,
      LocalMet->Tair, LocalMet->Vpd, UpperWind,
      &((*Gap)[Forest].RainFall), &((*Gap)[Forest].SnowFall),
      &((*Gap)[Forest].IntRain[0]), &((*Gap)[Forest].IntSnow[0]),
      &((*Gap)[Forest].TempIntStorage), &((*Gap)[Forest].CanopyVaporMassFlux),
      &((*Gap)[Forest].Tcanopy), &((*Gap)[Forest].MeltEnergy), VType->Height,
      VType->UnderStory);
    (*Gap)[Forest].MoistureFlux -= (*Gap)[Forest].CanopyVaporMassFlux;

    /* Because we now have a new estimate of the canopy temperature we can
    recalculate the longwave balance */
    if ((*Gap)[Forest].HasSnow == TRUE)
      Tsurf = (*Gap)[Forest].TSurf;
    else if (HeatFluxOption == TRUE)
      Tsurf = (*Gap)[Forest].TSurf;
    else
      Tsurf = LocalMet->Tair;

    /* update longwave radiation */
    GapSurroundingLongRadiation(&((*Gap)[Forest]), LocalMet->Lin, LocalVeg->Vf,
      LocalVeg->Fract[0], (*Gap)[Forest].Tcanopy, Tsurf);
  }
  /* if no snow */
  else if (VType->NVegLayers > 0) {
    (*Gap)[Forest].Tcanopy = LocalMet->Tair;
    (*Gap)[Forest].CanopyVaporMassFlux = 0.0;
    (*Gap)[Forest].TempIntStorage = 0.0;
    InterceptionStorage(NVegLActual, LocalVeg->MaxInt, LocalVeg->Fract, (*Gap)[Forest].IntRain,
      &((*Gap)[Forest].RainFall));
  }

  /* if snow is present, simulate the snow pack dynamics */
  if ((*Gap)[Forest].HasSnow || (*Gap)[Forest].SnowFall > 0.0) {

    SnowLongIn = LocalRad->LongIn[1];
    SnowNetShort = LocalRad->NetShort[1];
    SnowWind = VType->USnow * LocalMet->Wind;
    SnowRa = VType->RaSnow / LocalMet->Wind;

    (*Gap)[Forest].SnowPackOutflow =
      SnowMelt(y, x, Dt, 2.+Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, (*Gap)[Forest].RainFall, (*Gap)[Forest].SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &((*Gap)[Forest].PackWater), &((*Gap)[Forest].SurfWater),
        &((*Gap)[Forest].Swq), &((*Gap)[Forest].VaporMassFlux),
        &((*Gap)[Forest].TPack), &((*Gap)[Forest].TSurf), &((*Gap)[Forest].MeltEnergy));

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    (*Gap)[Forest].RainFall = 0.0;
    (*Gap)[Forest].MoistureFlux -= (*Gap)[Forest].VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
    can recalculate the longwave balance */
    Tsurf = (*Gap)[Forest].TSurf;
    GapSurroundingLongRadiation(&((*Gap)[Forest]), LocalMet->Lin, LocalVeg->Vf,
      LocalVeg->Fract[0], (*Gap)[Forest].Tcanopy, Tsurf);
  }
  else {
    (*Gap)[Forest].SnowPackOutflow = 0.0;
    (*Gap)[Forest].VaporMassFlux = 0.0;
  }

  /* Determine whether a snow pack is still present, or whether everything
  has melted */
  if ((*Gap)[Forest].Swq > 0.0)
    (*Gap)[Forest].HasSnow = TRUE;
  else
    (*Gap)[Forest].HasSnow = FALSE;
}

/*****************************************************************************
Function name: CalcGapSurroudingET()

Purpose      :
*****************************************************************************/
void CalcGapSurroudingET(int Dt, CanopyGapStruct **Gap, 
  SOILTABLE *SType, VEGTABLE *VType, PIXRAD *LocalRad, PIXMET *LocalMet, 
  SOILPIX *LocalSoil, ROADSTRUCT *LocalNetwork, float UpperRa, float LowerRa,
  VEGPIX *LocalVeg)

{
  float Rp;
  float NetRadiation;

  if (VType->OverStory == TRUE) {
    Rp = VISFRACT * (*Gap)[Forest].NetShort[0];
    NetRadiation = (*Gap)[Forest].NetShort[0] +
      (*Gap)[Forest].LongIn[0] - 2 * LocalVeg->Vf * (*Gap)[Forest].LongOut[0];
    (*Gap)[Forest].NetRadiation[0] = NetRadiation;

    EvapoTranspiration(0, 1, Dt, LocalVeg->Fract[0], LocalMet, NetRadiation,
      Rp, VType, SType, (*Gap)[Forest].MoistureFlux, (*Gap)[Forest].Moist, LocalSoil->Temp,
      &((*Gap)[Forest].IntRain[0]), (*Gap)[Forest].EPot, (*Gap)[Forest].EInt, (*Gap)[Forest].ESoil,
      (*Gap)[Forest].EAct, &((*Gap)[Forest].ETot), LocalNetwork->Adjust, UpperRa, LocalVeg);
    (*Gap)[Forest].MoistureFlux += (*Gap)[Forest].EAct[0] + (*Gap)[Forest].EInt[0];

    if ((*Gap)[Forest].HasSnow != TRUE && VType->UnderStory == TRUE) {
      Rp = VISFRACT * LocalRad->NetShort[1];
      NetRadiation =
        LocalRad->NetShort[1] +
        LocalRad->LongIn[1] - LocalVeg->Fract[1] * LocalRad->LongOut[1];
      LocalRad->NetRadiation[1] = NetRadiation;
      EvapoTranspiration(1, 1, Dt, LocalVeg->Fract[1], LocalMet, NetRadiation,
        Rp, VType, SType, (*Gap)[Forest].MoistureFlux, (*Gap)[Forest].Moist, LocalSoil->Temp,
        &((*Gap)[Forest].IntRain[1]), (*Gap)[Forest].EPot, (*Gap)[Forest].EInt, (*Gap)[Forest].ESoil,
        (*Gap)[Forest].EAct, &((*Gap)[Forest].ETot), LocalNetwork->Adjust, LowerRa, LocalVeg);
      (*Gap)[Forest].MoistureFlux += (*Gap)[Forest].EAct[1] + (*Gap)[Forest].EInt[1];
    }
    else if (VType->UnderStory == TRUE) {
      (*Gap)[Forest].EAct[1] = 0.;
      (*Gap)[Forest].EInt[1] = 0.;
      (*Gap)[Forest].NetRadiation[1] = 0.;
    }
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is
  present and there is no understory */
  if ((*Gap)[Forest].HasSnow != TRUE && VType->UnderStory != TRUE) {
    if (VType->OverStory == TRUE) {
      NetRadiation =
        (*Gap)[Forest].NetShort[1] + (*Gap)[Forest].LongIn[1] - (*Gap)[Forest].LongOut[1];
      (*Gap)[Forest].NetRadiation[1] = NetRadiation;
    }
    (*Gap)[Forest].EvapSoil =
    SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
        LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
        NetRadiation, LowerRa, (*Gap)[Forest].MoistureFlux, LocalSoil->Porosity[0],
        LocalSoil->FCap[0], SType->Ks[0], SType->Press[0], SType->PoreDist[0],
        VType->RootDepth[0], &((*Gap)[Forest].Moist[0]),
        LocalNetwork->Adjust[0]);
  }
  else
    (*Gap)[Forest].EvapSoil = 0.0;

  (*Gap)[Forest].MoistureFlux += (*Gap)[Forest].EvapSoil;
  (*Gap)[Forest].ETot += (*Gap)[Forest].EvapSoil;
}

/*****************************************************************************
Function name: AggregateCanopyGap()

Purpose      : Aggregate the gap and non-gap mass balance variables based
on area weight.
*****************************************************************************/
void AggregateCanopyGap(CanopyGapStruct **Gap, VEGPIX *LocalVeg,
  SOILPIX *LocalSoil, SNOWPIX *LocalSnow, EVAPPIX *LocalEvap,
  PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, double weight, int NSoil,
  int NVeg, int NVegLayers)
{
  int i, j;

  LocalPrecip->RainFall =
    weight*(*Gap)[Opening].RainFall + (1-weight)*(*Gap)[Forest].RainFall;
  LocalPrecip->SnowFall = 
    weight*(*Gap)[Opening].SnowFall + (1-weight)*(*Gap)[Forest].SnowFall;
  LocalPrecip->Precip = 
    weight*(*Gap)[Opening].Precip + (1-weight)*(*Gap)[Forest].Precip;

  LocalSnow->Outflow = 
    weight*(*Gap)[Opening].SnowPackOutflow + (1-weight)*(*Gap)[Forest].SnowPackOutflow;
  LocalSnow->CanopyVaporMassFlux =
    weight*(*Gap)[Opening].CanopyVaporMassFlux + (1-weight)*(*Gap)[Forest].CanopyVaporMassFlux;
  LocalSnow->VaporMassFlux =
    weight*(*Gap)[Opening].VaporMassFlux + (1-weight)*(*Gap)[Forest].VaporMassFlux;

  for (i = 0; i < 2; i++) {
    LocalRad->NetShort[i] = 
      weight*(*Gap)[Opening].NetShort[i] + (1-weight)*(*Gap)[Forest].NetShort[i];
    LocalRad->LongIn[i] = 
      weight*(*Gap)[Opening].LongIn[i] + (1-weight)*(*Gap)[Forest].LongIn[i];
    LocalRad->LongOut[i] = 
      weight*(*Gap)[Opening].LongOut[i] + (1-weight)*(*Gap)[Forest].LongOut[i];
  }

  LocalSnow->Swq = 
    weight*(*Gap)[Opening].Swq + (1-weight)*(*Gap)[Forest].Swq;
  LocalSnow->TPack = 
    weight*(*Gap)[Opening].TPack + (1-weight)*(*Gap)[Forest].TPack;
  LocalSnow->PackWater = 
    weight*(*Gap)[Opening].PackWater + (1-weight)*(*Gap)[Forest].PackWater;
  LocalSnow->SurfWater = 
    weight*(*Gap)[Opening].SurfWater + (1-weight)*(*Gap)[Forest].SurfWater;

  for (j = 0; j <= NSoil; j++) {
    LocalSoil->Moist[j] = weight*(*Gap)[Opening].Moist[j] +
      (1-weight)*(*Gap)[Forest].Moist[j];
  }

  LocalVeg->MoistureFlux =
    weight*(*Gap)[Opening].MoistureFlux + (1-weight)*(*Gap)[Forest].MoistureFlux;
  LocalVeg->MeltEnergy =
    weight*(*Gap)[Opening].MeltEnergy + (1-weight)*(*Gap)[Forest].MeltEnergy;

  /* Intercepted rain/snow */
  for (i = 0; i < NVeg; i++) {
    LocalPrecip->IntRain[i] =
      weight*(*Gap)[Opening].IntRain[i] + (1-weight)*(*Gap)[Forest].IntRain[i];
    LocalPrecip->IntSnow[i] =
      weight*(*Gap)[Opening].IntSnow[i] + (1-weight)*(*Gap)[Forest].IntSnow[i];
  }
  /* ET */
  for (i = 0; i <= NVeg; i++) {
      LocalEvap->EPot[i] =
        weight*(*Gap)[Opening].EPot[i] + (1-weight)*(*Gap)[Forest].EPot[i];
      LocalEvap->EAct[i] =
        weight*(*Gap)[Opening].EAct[i] + (1-weight)*(*Gap)[Forest].EAct[i];
  }
  for (i = 0; i < NVeg; i++) {
    LocalEvap->EInt[i] =
      weight*(*Gap)[Opening].EInt[i] + (1-weight)*(*Gap)[Forest].EInt[i];
  }
  for (i = 0; i < NVeg; i++) {
    for (j = 0; j < NSoil; j++) {
        LocalEvap->ESoil[i][j] =
          weight*(*Gap)[Opening].ESoil[i][j] + (1-weight)*(*Gap)[Forest].ESoil[i][j];
    }
  }
  LocalEvap->ETot =
    weight*(*Gap)[Opening].ETot + (1-weight)*(*Gap)[Forest].ETot;
}
