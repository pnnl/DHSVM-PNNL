/*
 * SUMMARY:      CalcSolar.c - inline solar calculations
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Batelle Pacific Northwest Laboratories
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Jul-96
 * DESCRIPTION:  These functions make inline solar radiation calculations
 *               that take into slope and aspect of the pixel, but do not
 *               account for shadowing of neighbouring pixels. 
 * DESCRIP-END.
 * FUNCTIONS:    SolarDay()
 *               SolarHour()
 *               SolarAngle()
 *               SolarConst()
 * COMMENTS:
 * $Id: CalcSolar.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include <math.h>
#include "constants.h"
#include "settings.h"
#include "Calendar.h"
#include "functions.h"
#include "rad.h"

/*****************************************************************************
  Function name: SolarDay()

  Purpose:	 This subroutine calculates daily solar values 

  Required:
   int DayOfYear	  - day of year (January 1 is 1)
   float Longitude        - site longitude (rad)
   float Latitude         - site latitude (rad)
   float StandardMeridian - longitude of time zone of standard meridian (rad)

  Returns: void

  Modifies:
    float *NoonHour       - true solar noon (hr)
    float *Declination    - solar Declination (rad)
    float *HalfDayLength  - half-day length (hr)
    float *Sunrise        - time of Sunrise (hr)
    float *Sunset         - time of Sunset (hr)
    float *TimeAdjustment - required adjustment to local time (hr)
    float *SunEarthDist   - distance from sun to earth

  Comments     : EXECUTE AT START OF EACH DAY 
*****************************************************************************/
void SolarDay(int DayOfYear, float Longitude, float Latitude,
	      float StandardMeridian, float *NoonHour, float *Declination,
	      float *HalfDayLength, float *Sunrise, float *Sunset,
	      float *TimeAdjustment, float *SunEarthDist)
{
  float B;			/* coefficient for equation of time */
  float EqnOfTime;		/* adjustment for equation of time (min) */
  float LongitudeAdjust;	/* adjustment for longitude (min) */
  float CosineHalfDayLength;	/* cosine of the half-day length */

  /* note need to check if day light savings time calculate adjustment for 
     true solar time longitude adjustment add 4 min per degree away from 
     StandardMeridian (4 min/degree * 180 degree/pi radian) */

  LongitudeAdjust = (MINPDEG * DEGPRAD) * (StandardMeridian - Longitude);

  /* equation of time */
  B = (2.0 * PI * (DayOfYear - 81)) / 364.;
  EqnOfTime = 9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B);

  /* adjustment factor to convert local time to solar time
     solar time = local time + TimeAdjustment  */
  /* for example from GMT to the west coast of the us (PST) */
  /* is a -8 hour shift, i.e. PST = GMT-8 */
  *TimeAdjustment = -(LongitudeAdjust + EqnOfTime) / MINPHOUR;

  /* work in solar time  */
  *NoonHour = 12.0;

  /* solar Declinationation  */
  *Declination = .4098 * sin(2 * PI * (284 + DayOfYear) / DAYPYEAR);

  /* half-day length  */
  CosineHalfDayLength = -tan(Latitude) * tan(*Declination);
  if (CosineHalfDayLength >= 1.0)
    *HalfDayLength = PI;
  else
    *HalfDayLength = acos(CosineHalfDayLength);

  /* convert HalfDayLength from radians to Hours 
     1 radian = (180 deg / PI) * (1 hr / 15 degrees rotation) */
  *HalfDayLength = *HalfDayLength / RADPHOUR;

  /* solar time of Sunrise and Sunset  */
  *Sunrise = *NoonHour - *HalfDayLength;
  *Sunset = *NoonHour + *HalfDayLength;

  /* calculate the sun-earth distance */
  *SunEarthDist = 1.0 + 0.033 * cos(RADPDEG * (360. * DayOfYear / 365));
}

/*****************************************************************************
  Function name: SolarHour()  
  
  Purpose: This subroutine calculates position of the sun as a function of
           the time of day, the length of time the sun is above the horizon,
           and the maximum radiation.

  Required:
    float Latitude		- site laditude (rad)
    float LocalHour		- local time (hr)
    float Dt			- length of current timestep (hr)
    float NoonHour		- true solar noon (hr)
    float Declination		- solar Declination (rad)
    float Sunrise		- time of Sunrise (hr)
    float Sunset		- time of Sunset (hr)
    float TimeAdjustment	- required adjustment to convert local time
				  to solar time (hr)
    float SunEarthDist          - distance from Sun to Earth 

  Returns: void

  Modifies:
    float *SineSolarAltitude - sine of sun's SolarAltitude 
    int *DayLight	     - FALSE: measured solar radiation and the sun is
                                      below the horizon.  
			       TRUE: sun is above the horizon
    float *SolarTimeStep     - fraction of the timestep the sun is above the 
                               horizon 
    float *SunMax            - calculated solar radiation at the top of the 
                               atmosphere (W/m^2) 
                                                 
  Comments     : EXECUTE AT START OF EACH TIMESTEP 
*****************************************************************************/
void SolarHour(float Latitude, float LocalHour, float Dt, float NoonHour,
	       float Declination, float Sunrise, float Sunset,
	       float TimeAdjustment, float SunEarthDist,
	       float *SineSolarAltitude, int *DayLight, float *SolarTimeStep,
	       float *SunMax, float *SolarAzimuth)
{
  float SolarAltitude;		/* SolarAltitude of sun from horizon (rads) */
  float SolarZenith;		/* sun zenith angle (rads) */
  float StartHour = 0;		/* currect Hour in solar time (hr) */
  float EndHour = 0;		/* mid-point of current solar Hour (hr) */
  float Hour;			/* angle of current "halfhr" from solar noon
				   (rads) */

  /* NOTE THAT HERE Dt IS IN HOURS NOT IN SECONDS */

  *SunMax = 0.0;
  *SineSolarAltitude = 0.0;
  *SolarTimeStep = 1.0;

  /* all calculations based on hour, the solar corrected local time */

  Hour = LocalHour + TimeAdjustment;

  if (Hour < 0)
    Hour += 24;
  if (Hour > 24)
    Hour -= 24;

  *DayLight = FALSE;
  if ((Hour > Sunrise) && ((Hour - Dt) < Sunset))
    *DayLight = TRUE;

  if (*DayLight == TRUE) {
    /* sun is above the horizon */

    if (Dt > 0.0) {
      /* compute average solar SolarAltitude over the timestep */
      StartHour = (Hour - Dt > Sunrise) ? Hour - Dt : Sunrise;
      EndHour = (Hour < Sunset) ? Hour : Sunset;

      /*  convert to radians  */
      StartHour = RADPHOUR * (StartHour - NoonHour);
      EndHour = RADPHOUR * (EndHour - NoonHour);
      *SolarTimeStep = EndHour - StartHour;

      /*  determine the average geometry of the sun angle  */
      *SineSolarAltitude = sin(Latitude) * sin(Declination)
	+ (cos(Latitude) * cos(Declination) * (sin(EndHour) - sin(StartHour))
	   / *SolarTimeStep);
    }
    else {
      Hour = RADPHOUR * (Hour - NoonHour);
      *SineSolarAltitude = sin(Latitude) * sin(Declination)
	+ (cos(Latitude) * cos(Declination) * cos(Hour));
    }

    SolarAltitude = asin(*SineSolarAltitude);

    SolarZenith = PI / 2 - SolarAltitude;

    *SolarAzimuth = ((sin(Latitude) * (*SineSolarAltitude) - sin(Declination))
		     / (cos(Latitude) * sin(SolarZenith)));

    if (*SolarAzimuth > 1.)
      *SolarAzimuth = 1.;

    *SolarAzimuth = acos(-(*SolarAzimuth));

    if (Dt > 0.0) {
      if (fabs(EndHour) > fabs(StartHour)) {
	*SolarAzimuth = 2 * PI - (*SolarAzimuth);
      }
    }
    else {
      if (Hour > 0) {
	*SolarAzimuth = 2 * PI - (*SolarAzimuth);
      }
    }
/*     printf("Solar Azimuth = %g\n", *SolarAzimuth*180./PI); */

    *SunMax = SOLARCON * SunEarthDist * *SineSolarAltitude;
  }
}

/*****************************************************************************
  Function name: SolarAngle()  
  
  Purpose      : this subroutine uses solar radiation measured on a horizontal
                 surface to calculate direct, diffuse, and incoming reflected 
                 solar radiation for a sloping surface.

  Required     :
    float Latitude	        - site laditude (rad)
    float Albedo                - surface Albedo
    float Declination	        - solar Declination (rad)
    float CellAspect	        - cell aspect (rads eastward from north)
    float CellSlope	        - ground surface slope (rads)
    float SineSolarAltitude	- sine of sun's SolarAltitudede 
    float SunMax                - calculated solar radiation at the top of 
                                  the atmosphere (W/m^2) 
    float SolarTimeStep         - fraction of the timestep the sun is above 
                                  the horizon 
    int DayLight	        - FALSE: measured solar radiation and the sun
                                  is below the horizon.  
                                  TRUE: sun is above the horizon
     float SolarAzimuth    -  solar azimuth (rads eastward from north)
    float Dt			- length of current timestep (hr)

  Returns      : void

  Modifies     :
    float *Direct  - direct beam solar radiation (W/m^2)
    float *Diffuse - diffuse solar radiation (W/m^2)

  Comments     : EXECUTE EACH TIMESTEP FOR EACH GRID CELL 
                 Source:  Gates, D.M., "Biophysical ecology",
                 Springer-Verlag,  New York, etc.,  1980.
                 Especially Chapter 6, p. 136 ->
*****************************************************************************/
void SolarAngle(float Latitude, float Albedo, float Declination,
		float CellAspect, float CellSlope, float SunMax,
		float SineSolarAltitude, float SolarTimeStep, int DayLight,
		float SolarAzimuth, float Dt, float *Direct, float *Diffuse)
{
  float SolarAltitude;		/* SolarAltitude of sun from horizon (rads) */
  float CosineIncidenceAngle;	/* cosine of the incidence angle between
				   solar rays and the normal to the surface */
  float DiffuseSkyView;		/* sky view factor for diffuse radation 
				   (0.0 - 1.0)  */
  float ReflectedSkyView;	/* view factor for incoming refected solar 
				   radiation (0.0 - 1.0) */
  float Reflect;		/* incoming reflected solar radiation (W/m^2)
				 */

  /* NOTE THAT HERE Dt IS IN HOURS, NOT IN SECONDS */

  if (DayLight == TRUE) {

    /* calculate shortwave components horizontal surface */

    if (fequal(CellSlope, 0.0)) {
      *Direct = SunMax;
      *Diffuse = SunMax;
      Reflect = 0.0;
    }
    else {
      /*  sloping surface  */
      DiffuseSkyView = (PI - CellSlope) / PI;
      ReflectedSkyView = CellSlope / PI;
      SolarAltitude = asin(SineSolarAltitude);

      CosineIncidenceAngle = cos(SolarAltitude) * sin(CellSlope) *
	cos(SolarAzimuth - CellAspect) + cos(CellSlope) * sin(SolarAltitude);

      *Direct = SunMax * CosineIncidenceAngle / SineSolarAltitude;
      if (CosineIncidenceAngle <= 0.0)
	*Direct = 0.0;

      *Diffuse = SunMax * DiffuseSkyView;
      Reflect = Albedo * SunMax * ReflectedSkyView;
    }				/* end sloping surface */
  }				/* end daylight is true */
  else {
    /* sun is below the horizon */

    /* all measured solar radiation is diffuse, thus total rad on CellSlope 
       is (measured)*DiffuseSkyView */
    *Direct = 0.0;
    Reflect = 0.0;
    *Diffuse = (PI - CellSlope) / PI;
  }

  /* Average over timestep */
  if (Dt > 0.0) {

    /* The assumption is made that the incoming reflected radiation is diffuse 
       radiation for the "receiving" gridcell */

    *Diffuse += Reflect;
    *Direct *= SolarTimeStep / (Dt * RADPHOUR);
    *Diffuse *= SolarTimeStep / (Dt * RADPHOUR);
  }
}

/*****************************************************************************
  Function name: SolarConst()  
  
  Purpose      : this subroutine reads in site data, calculates constants,
                 and converts data from degrees to radians.

  Required     :
    float lat_deg           - latitude of site in degrees
    float lat_min           - latitude of site in minutes
    float lng_deg           - longitude of site in degrees
    float lng_min           - longitude of site in minutes
    float *StandardMeridian - standard meridian for local time zone in 
                              degrees
  
  Returns      : void

  Modifies     :
    float *Latitude         - latitude of site in radians
    float *Longitude        - longitude of site in radians
    float *StandardMeridian - standard meridian for local time zone (rads)

  Comments     : EXECUTE ONCE BEFORE TIME LOOP
*****************************************************************************/
void SolarConst(float lat_deg, float lat_min, float lng_deg,
		float lng_min, float *StandardMeridian, float *Latitude,
		float *Longitude)
{
  /*  convert to radians  */
  *Latitude = (lat_deg + (lat_min / 60.)) * PI / 180.;
  *Longitude = (lng_deg + (lng_min / 60.)) * PI / 180.;
  *StandardMeridian = *StandardMeridian * PI / 180.;
}
