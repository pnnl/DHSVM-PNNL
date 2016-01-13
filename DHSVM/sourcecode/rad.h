/*
 * SUMMARY:      rad.h header file for radiation functions for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for radiation functions for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: rad.h,v 1.4 2003/07/01 21:26:31 olivier Exp $     
 */

#ifndef RAD_H
#define RAD_H

#define ALBEDO  0.15		/* WORK IN PROGRESS, See InitNewStep() */

void SeparateRadiation(float TotalSolar, float ClearIndex,
		       float *Beam, float *Diffuse);

void SolarAngle(float Latitude, float Albedo, float Declination,
		float CellAspect, float CellSlope, float SunMax,
		float SineSolarAltitude, float SolarTimeStep, int DayLight,
		float SolarAzimuth, float Dt, float *Direct, float *Diffuse);

void SolarConst(float lat_deg, float lat_min, float lng_deg,
		float lng_min, float *StandardMeridian, float *Latitude,
		float *Longitude);

void SolarDay(int DayOfYear, float Longitude, float Latitude,
	      float StandardMeridian, float *NoonHour, float *Declination,
	      float *HalfDayLength, float *Sunrise, float *Sunset,
	      float *TimeAdjustment, float *SunEarthDist);

void SolarHour(float Latitude, float LocalHour, float Dt, float NoonHour,
	       float Declination, float Sunrise, float Sunset,
	       float TimeAdjustment, float SunEarthDist,
	       float *SineSolarAltitude, int *DayLight, float *SolarTimeStep,
	       float *SunMax, float *SolarAzimuth);

#endif
