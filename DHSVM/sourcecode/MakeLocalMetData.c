/*
 * SUMMARY:      MakeLocalMetData.c - Generates meteorological conditions
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Generates meteorological conditions for each individual cell
 * DESCRIP-END.
 * FUNCTIONS:    MakeLocalMetData()
 * COMMENTS:
 * $Id: MakeLocalMetData.c,v3.1.2 2014/01/1 ning Exp $     
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "snow.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "rad.h"

/*****************************************************************************
  Function name: MakeLocalMetData()

  Purpose      : Generates meteorological for each individual cell

  Required     :
    int y 
    int x
    MAPSIZE Map
    int DayStep
    unsigned char PrecipType
    int NStats
    METLOCATION *Stat
    uchar *MetWeights
    float LocalElev
    RADCLASSPIX *RadMap 
    PRECIPPIX *PrecipMap
    MAPSIZE Radar
    RADARPIX **RadarMap

  Returns      :
    PIXMET LocalMet

  Modifies     :

  Comments     :
    Reference: Shuttleworth, W.J., Evaporation,  In: Maidment, D. R. (ed.),
               Handbook of hydrology,  1993, McGraw-Hill, New York, etc..
*****************************************************************************/
PIXMET MakeLocalMetData(int y, int x, MAPSIZE * Map, int DayStep,
			OPTIONSTRUCT * Options, int NStats,
			METLOCATION * Stat, uchar * MetWeights,
			float LocalElev, RADCLASSPIX * RadMap,
			PRECIPPIX * PrecipMap, MAPSIZE * Radar,
			RADARPIX ** RadarMap, float **PrismMap,
			SNOWPIX * LocalSnow, SNOWTABLE * SnowAlbedo,
			float ***MM5Input, float ***WindModel,
			float **PrecipLapseMap, MET_MAP_PIX *** MetMap,
			int NGraphics, int Month, float skyview,
			unsigned char shadow, float SunMax,
			float SineSolarAltitude)
{
  float CurrentWeight;		/* weight for current station */
  float ScaleWind = 1;		/* Wind to be scaled by model factors if 
				   WindSource == MODEL */
  float Temp;			/* Temporary variable */
  float WeightSum;		/* sum of the weights */
  int i;			/* counter */
  int RadarX;			/* X coordinate of radar map coordinate */
  int RadarY;			/* Y coordinate of radar map coordinate */
  float TempLapseRate;
  int WindDirection = 0;	/* Direction of model wind */
  PIXMET LocalMet;		/* local met data */

  LocalMet.Tair = 0.0;
  LocalMet.Rh = 0.0;
  LocalMet.Wind = 0.0;
  LocalMet.Sin = 0.0;
  LocalMet.SinBeam = 0.0;
  LocalMet.SinDiffuse = 0.0;
  LocalMet.Lin = 0.0;
  TempLapseRate = 0.0;

  if (Options->MM5 == TRUE && Options->QPF == TRUE) {
    WeightSum = 0.0;
    for (i = 0; i < NStats; i++)
      WeightSum += (float) MetWeights[i];
  }

  if (Options->MM5 == TRUE) {
    LocalMet.Tair = MM5Input[MM5_temperature - 1][y][x] +
		(LocalElev - MM5Input[MM5_terrain - 1][y][x]) * 
		MM5Input[MM5_lapse - 1][y][x];
    LocalMet.Rh = MM5Input[MM5_humidity - 1][y][x];
    LocalMet.Wind = MM5Input[MM5_wind - 1][y][x];
    LocalMet.Sin = MM5Input[MM5_shortwave - 1][y][x];

    if (Options->Shading == TRUE) {
      if (SunMax > 0.0) {
		SeparateRadiation(LocalMet.Sin, LocalMet.Sin / SunMax,
		&(LocalMet.SinBeam), &(LocalMet.SinDiffuse)); 
	  }
      else {
		/* if sun is below horizon, the force all shortwave to zero */
		LocalMet.Sin = 0.0;
		LocalMet.SinBeam = 0.0;
		LocalMet.SinDiffuse = 0.0; 
	  }
    }
    LocalMet.Lin = MM5Input[MM5_longwave - 1][y][x];
    LocalMet.Press = 101300.0;
    PrecipMap->Precip = MM5Input[MM5_precip - 1][y][x];
  }
  else {			/* MM5 is false and we need to interpolate the basic met records */
    WeightSum = 0.0;
    for (i = 0; i < NStats; i++) {
      WeightSum += (float) MetWeights[i];
      if (Options->WindSource == MODEL && Stat[i].IsWindModelLocation) {
	    ScaleWind = Stat[i].Data.Wind;
	    WindDirection = Stat[i].Data.WindDirection;
      }
    }
    for (i = 0; i < NStats; i++) {
      CurrentWeight = ((float) MetWeights[i]) / WeightSum;
      LocalMet.Tair += CurrentWeight *
	  LapseT(Stat[i].Data.Tair, Stat[i].Elev, LocalElev,
	       Stat[i].Data.TempLapse);
      LocalMet.Rh += CurrentWeight * Stat[i].Data.Rh;
      if (Options->WindSource == STATION)
	  LocalMet.Wind += CurrentWeight * Stat[i].Data.Wind;
      LocalMet.Lin += CurrentWeight * Stat[i].Data.Lin;
      LocalMet.Sin += CurrentWeight * Stat[i].Data.Sin;
      if (Options->Shading == TRUE) {
	    LocalMet.SinBeam += CurrentWeight * Stat[i].Data.SinBeamObs;
	    LocalMet.SinDiffuse += CurrentWeight * Stat[i].Data.SinDiffuseObs;
      }
      TempLapseRate += CurrentWeight * Stat[i].Data.TempLapse;
    }
    if (Options->WindSource == MODEL)
      LocalMet.Wind = ScaleWind * WindModel[WindDirection - 1][y][x];

    if (Options->PrecipType == RADAR) {
      RadarY = (int) ((y + Radar->OffsetY) * Map->DY / Radar->DY);
      RadarX = (int) ((x - Radar->OffsetX) * Map->DX / Radar->DX);
      PrecipMap->Precip = RadarMap[RadarY][RadarX].Precip;
    }
    /* WORK IN PROGRESS, taken from old DHSVM version */
    /* Air pressure */
    /* In rare cases - i.e. when the lapse rate has a different sign for 
       different met stations - you can end up with a TemplapseRate of 0.0
       This will result in a crash, so a check was put in (Jul 28, 1997 - Bart
       Nijssen).  It is somewhat awkward to interpolate lapse rates anyway, so
       a better way of doing this would be welcome */
    if (TempLapseRate != 0.0) {
      Temp = 9.8067 / (TempLapseRate * 287.0);
      LocalMet.Press = 101300. * pow(((288.0 - TempLapseRate * LocalElev) / 288.0), Temp);
    }
    else
      LocalMet.Press = 101300.;

  }				/* end of else MM5==TRUE, i.e. all basic met, except for precip */
  /* has been interpolated */

  /* Here is how the following section works */
  /* Arc-Info (through use of the hillshade command) will give */
  /* an output file that ranges from 0 to 255 (the shade factor) */
  /* These correspond to the reflectance of the direct beam radiation */
  /* for a given sun position (altitude and azimuth) taking */
  /* into account the slope and aspect and topographic shading */
  /* of the local pixel.  If we wanted to use this value directly */
  /* then the correction to the observed beam w.r.t. a horizontal plane */
  /* would be   actual = horizontal*shadefactor/255/sin(solar_altitude) */
  /* the sin(solar_altitude) is necessary to convert horizontal into the maximum */
  /* possible flux */

  /* We can either have DHSVM make the solar_altitude calculation, which */
  /* is not all that hard, but is prone to user error (e.g. GMT time shifts, etc) */
  /* or we can simply include the solar_altitude info in the shade_factor */
  /* the question is how do we include the sin(solar_altitude) while */
  /* keeping new_shade_factor = shadefactor/255/sin(solar_altitude) defined */
  /* as a unsigned character */
  /* Answer:  At sal = 5 degrees max_new_shadefactor = 11.47 */
  /* i.e. the actual flux normal to sal is 11.47*observed_horizontal_flux */
  /* if we adopt this 5 degree value as a cutoff, we can then be assured that */
  /*    0<=newshadefactor<=11.47      and then scale it between 0 and 255 */
  /* the final calculation becomes,  */
  /*       actual = horizontal*(float)shadefactor/255.0*11.47 or simply
     acutal = horizontal*(float)shadefactor/22.23191               */
  /* thus radiation increases from 0 to 11.47 times the observed value in */
  /* increments of 4.5 percent */
  /* a finer resolution than this would require a higher min angle or more memory */

  if (Options->Shading == TRUE) {
    /* commented by Ning. the program script used to generate the shadow files
    are update to produce shadow factors ranging from 0 to 255 consistent with 
    arcinfo */
    LocalMet.SinBeam *= (float) shadow / 22.23191;
	
	/* if canopy shading is computed, then the skyview factor will be compared with
       riparian canopy openess */
	if (Options->CanopyShading && Options->StreamTemp)
	  LocalMet.SinDiffuse *= 1;
	else
	  LocalMet.SinDiffuse *= skyview;
    if (LocalMet.SinBeam + LocalMet.SinDiffuse > SOLARCON)
      LocalMet.SinBeam = SOLARCON - LocalMet.SinDiffuse;
  }
  else {
    LocalMet.SinBeam = LocalMet.Sin;
    LocalMet.SinDiffuse = 0;
  }
  RadMap->Beam = LocalMet.SinBeam;
  RadMap->Diffuse = LocalMet.SinDiffuse;

    /* Store the VIC incoming shortwave radiatio without topo or canopy shading */
  if (Options->StreamTemp)
	LocalMet.VICSin = LocalMet.Sin;
  LocalMet.Sin = RadMap->Beam + RadMap->Diffuse;

  if (Options->QPF == TRUE || Options->MM5 == FALSE) {
    if (Options->PrecipType == STATION && Options->Prism == FALSE) {
      PrecipMap->Precip = 0.0;
      for (i = 0; i < NStats; i++) {
	    CurrentWeight = ((float) MetWeights[i]) / WeightSum;
	    if (Options->PrecipLapse == MAP)
	      PrecipMap->Precip += CurrentWeight *
	      LapsePrecip(Stat[i].Data.Precip, 0, 1, PrecipLapseMap[y][x]);
		else
		  PrecipMap->Precip += CurrentWeight *
	      LapsePrecip(Stat[i].Data.Precip, Stat[i].Elev, LocalElev,
		  Stat[i].Data.PrecipLapse);
      }
    }
    else if (Options->PrecipType == STATION && Options->Prism == TRUE) {
      PrecipMap->Precip = 0.0;
      for (i = 0; i < NStats; i++) {
	    CurrentWeight = ((float) MetWeights[i]) / WeightSum;
		/* this is the real prism interpolation */
		/* note that X = position from left  boundary, ie # of columns */
		/* note that Y = position from upper boundary, ie # of rows   */
		if (Options->Outside == FALSE)
	      PrecipMap->Precip += CurrentWeight * Stat[i].Data.Precip /
		  PrismMap[Stat[i].Loc.N][Stat[i].Loc.E] * PrismMap[y][x];
		else
		  PrecipMap->Precip += CurrentWeight * Stat[i].Data.Precip /
		  Stat[i].PrismPrecip[Month - 1] * PrismMap[y][x];
		if (PrismMap[y][x] < 0){
		  printf("negative PrismMap value in MakeLocalMetData.c\n");
		  exit(0);
		}
      }
    }
  }

  /* due to the nature of the interpolation scheme in DHSVM and the */
  /* interpolation scheme to handle the mess of different formats of met stations */
  /* in the PRISM project */
  /* relative humidities can be quite low when precip is occuring */
  /* at times this will results in PET being greater than precip */
  /* allow an option in DHSVM to override RH if Precip is occuring */

  if (Options->Rhoverride == TRUE) {
    if (PrecipMap->Precip > 0.0)
      LocalMet.Rh = 100.0;
  }

  /* Separate precipitation into rainfall and snowfall */
  if (PrecipMap->Precip > 0.0 && LocalMet.Tair < MAX_SNOW_TEMP) {
    if (LocalMet.Tair > MIN_RAIN_TEMP)
      PrecipMap->SnowFall = PrecipMap->Precip *
	(MAX_SNOW_TEMP - LocalMet.Tair) / (MAX_SNOW_TEMP - MIN_RAIN_TEMP);
    else
      PrecipMap->SnowFall = PrecipMap->Precip;
  }
  else
    PrecipMap->SnowFall = 0.0;

  PrecipMap->RainFall = PrecipMap->Precip - PrecipMap->SnowFall;

  /* Local heat of vaporization, Eq. 4.2.1, Shuttleworth (1993) */
  LocalMet.Lv = 2501000 - 2361 * LocalMet.Tair;

  /* Psychrometric constant */
  LocalMet.Gamma = CP * LocalMet.Press / (EPS * LocalMet.Lv);

  /* Saturated vapor pressure, Eq. 4.2.2, Shuttleworth (1993) */
  LocalMet.Es = SatVaporPressure(LocalMet.Tair);

  /* Slope of vapor pressure curve, Eq. 4.2.3, Shuttleworth (1993) */
  LocalMet.Slope = 4098.0 * LocalMet.Es /
    ((237.3 + LocalMet.Tair) * (237.3 + LocalMet.Tair));

  /* Actual vapor pressure */
  LocalMet.Eact = LocalMet.Es * (LocalMet.Rh / 100.);

  /* Vapor pressure deficit */
  LocalMet.Vpd = LocalMet.Es - LocalMet.Eact;

  /* Air density, Eq. 4.2.4 Shuttleworth (1993) */
  LocalMet.AirDens = 0.003486 * LocalMet.Press / (275 + LocalMet.Tair);

  if (LocalSnow->HasSnow) {
    if (PrecipMap->SnowFall > 0.0)
      LocalSnow->LastSnow = 0;
    else
      LocalSnow->LastSnow++;
    LocalSnow->Albedo = CalcSnowAlbedo(LocalSnow->TSurf, LocalSnow->LastSnow,
				       SnowAlbedo);
  }
  else
    LocalSnow->LastSnow = 0;

  if (NGraphics > 0) {
    (*MetMap)[y][x].accum_precip =
    (*MetMap)[y][x].accum_precip + PrecipMap->Precip;
    (*MetMap)[y][x].air_temp = LocalMet.Tair;
    (*MetMap)[y][x].wind_speed = LocalMet.Wind;
    (*MetMap)[y][x].humidity = LocalMet.Rh;
  }

  return LocalMet;
}
