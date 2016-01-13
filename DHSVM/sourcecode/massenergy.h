/*
 * SUMMARY:      massenergy.h - header file for MassEnergyBalance() part 
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for MassEnergyBalance() part of DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: massenergy.h,v 1.5 2003/11/12 20:01:53 colleen Exp $     
 */

#ifndef MASSENERGY_H
#define MASSENERGY_H

#include "data.h"
#include <stdarg.h>

void AggregateRadiation(int MaxVegLayers, int NVegL, PIXRAD * Rad,
			PIXRAD * TotalRad);

float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
		       float VpdThres, float MoistThres, float WP,
		       float TSoil, float SoilMoisture, float Vpd, float Rp);

float Desorption(int Dt, float Moisture, float Porosity, float Ks, 
			   float Press, float m);

void EvapoTranspiration(int Layer, int Dt, PIXMET * Met, float NetRad,
			   float Rp, VEGTABLE * VType, SOILTABLE * SType,
			   float MoistureFlux, SOILPIX * LocalSoil, float *Int,
			   EVAPPIX * LocalEvap, float *Adjust, float Ra);

void InitLocalRad(int HeatFluxOption, float Rs, float Ld,
		       float Tair, float Tcanopy, float Tsoil,
		       VEGTABLE * VType, SNOWPIX * LocalSnow, PIXRAD * LocalRad);

void InterceptionStorage(int NMax, int NAct, float *MaxInt, float *Fract,
			   float *Int, float *Precip, float *MomentSq, float *Height, 
			   unsigned char Understory, float Dt,
			   float MS_Rainfall, float LD_FallVelocity);

void LongwaveBalance(OPTIONSTRUCT *Options, unsigned char OverStory, 
			   float F, float Ld, float Tcanopy, float Tsurf, PIXRAD * LocalRad);

void NoEvap(int Layer, int NSoilLayers, EVAPPIX * LocalEvap);

void NoSensibleHeatFlux(int Dt, PIXMET * LocalMet, float ETot, SOILPIX * LocalSoil);

void RadiationBalance(OPTIONSTRUCT *Options, int HeatFluxOption, int CanopyRadAttOption, 
		      float SineSolarAltitude, float VICRs, float Rs,
		      float Rsd, float Rsb, float Ld, float Tair,
		      float Tcanopy, float Tsoil, float SoilAlbedo,
		      VEGTABLE *VType, SNOWPIX *LocalSnow, PIXRAD *LocalRad);

void SensibleHeatFlux(int y, int x, int Dt, float Ra, float ZRef,
		      float Displacement, float Z0, PIXMET * LocalMet,
		      float NetShort, float LongIn, float ETot,
		      int NSoilLayers, float *SoilDepth, SOILTABLE * SoilType,
		      float MeltEnergy, SOILPIX * LocalSoil);

void ShortwaveBalance(OPTIONSTRUCT *Options, unsigned char OverStory, 
					  float F, float Rs, float Rsb,
					  float Rsd, float Tau, float *Albedo, PIXRAD * LocalRad);

float SoilEvaporation(int Dt, float Temp, float Slope, float Gamma, float Lv,
		      float AirDens, float Vpd, float NetRad, float RaSoil,
		      float Transpiration, float Porosity, float Ks,
		      float Press, float m, float RootDepth,
		      float *MoistContent, float Adjust);

float StabilityCorrection(float Z, float d, float Tsurf, float Tair,
			  float Wind, float Z0);

float SurfaceEnergyBalance(float TSurf, va_list ap);

#endif
