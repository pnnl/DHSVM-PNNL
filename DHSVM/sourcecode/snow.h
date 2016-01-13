/*
 * SUMMARY:      snow.h - header file for DHSVM snow routines
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:    29-Aug-1996 at 16:03:11
 * DESCRIPTION:  header file for DHSVM snow routines
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: snow.h,v 1.5 2003/11/12 20:01:53 colleen Exp $     
 */

#ifndef SNOW_H
#define SNOW_H

#include <stdarg.h>

#define MAX_SURFACE_SWE          0.125	/* maximum depth of the surface layer
					   in water equivalent (m) */

void MassRelease(float *InterceptedSnow, float *TempInterceptionStorage,
		 float *ReleasedMass, float *Drip, float MDRatio);

void SnowInterception(int y, int x, int Dt, float F, float LAI,
		      float MaxInt, float MaxSnowIntCap, float MDRatio,
		      float SnowIntEff, float Ra, float AirDens, float EactAir,
		      float Lv, PIXRAD * LocalRad, float Press, float Tair,
		      float Vpd, float Wind, float *RainFall, float *SnowFall,
		      float *IntRain, float *IntSnow, float *TempIntStorage,
		      float *VaporMassFlux, float *Tcanopy, float *MeltEnergy, 
		      float *MomentSq, float *Height, unsigned char UnderStory,
		      float MS_Rainfall, float LD_FallVelocity);

float SnowMelt(int y, int x, int Dt, float Z, float Displacement, float Z0,
	       float BaseRa, float AirDens, float EactAir, float Lv,
	       float ShortRad, float LongRadIn, float Press, float RainFall,
	       float SnowFall, float Tair, float Vpd, float Wind,
	       float *PackWater, float *SurfWater, float *Swq,
	       float *VaporMassFlux, float *TPack, float *TSurf,
	       float *MeltEnergy);

float SnowPackEnergyBalance(float TSurf, va_list ap);

#endif
