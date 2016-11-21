/*
 * SUMMARY:      constants.h - header file with constants for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD: Fri Jun  6 10:58:14 1997 by Bart Nijssen <nijssen@meter.ce.washington.edu>
 * DESCRIPTION:  header file with constants for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 */

/* 	$Id: constants.h,v 1.7 1997/04/23 22:55:20 nijssen Exp nijssen $	 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define CH_ICE     (2100.0e3)   /* Volumetric heat capacity (J/(m3*C) of ice 
                                   (0C) */

#define CH_WATER   (4186.8e3)   /* Volumetric heat capacity (J/(m3*C) of 
                                   water */

#define CP          1013.0      /* Specific heat of moist air at constant 
                                   pressure (J/(kg*C)) */

#define DAYSPYEAR    365.

#define DELTAT        25.       /* Used in SensibleHeatFlux to bracket the 
                                   effective surface temperature (C) */

#define DEGPRAD       57.29578  /* degree per radian */

#define D0_MULTIPLIER  0.63     /* Multiplier for vegetation height to get
                                   displacement height (m) */

#define DZ_TOP         0.1      /* Thickness of soil surface layer for which 
                                   heat stoarge change is calculated (m) */

#define EPS            0.622    /* ratio of molecular weight of water vapor to
                                   that for dry air */

#define G              9.81     /* gravitational accelleration (m/(s^2)) */

#define GRAMSPKG    1000.       /* grams per kilogram */

#define HOURSPDAY     24.

#define JOULESPCAL     4.1868   /* Joules per calorie */

#define LF           (333.7e3)  /* latent heat of fusion (J/kg) */

#define MINPDEG        4.       /* minutes per degree longitude */

#define MINPHOUR      60.       /* minutes per hour */

#define MMTOM          0.001    /* convert from mm to meter */

#define MONTHSPYR     12.

#undef PI
#define PI             3.14159265358979323846

#define RADPHOUR       0.2617994 /* radians per hour: Earth's Rotation 
                                    (2 PI rad/day) * (1 day/24 h) */
#define RADPDEG     (PI/180.0)   /* radians per degree */

#define SECPHOUR    3600.

#define SOLARCON    1360.       /* Solar constant (W/m^2) */

#define STEFAN    (5.6696e-8)   /* Stefan-Boltzmann constant (W/(M^2*C^4) */

#define VISFRACT       0.5      /* part of shortwave that is in the visible 
                                   range of the spectrum */

#define VON_KARMAN     0.4      /* Von Karman's constant */

#define WATER_DENSITY 1000.     /* Density of water in kg/m3 */

#define Z0_MULTIPLIER  0.13     /* Multiplier for vegetation height to get
                                   roughness length (m) */

/**************** extern constants - see globals.c ****************/

extern float LAI_SNOW_MULTIPLIER; /* multiplier to calculate the amount
				     of available snow interception as 
				     a function of LAI*/
extern float LAI_WATER_MULTIPLIER; /* multiplier to determine maximum
				      interception storage as a function of
				      LAI  */
extern float LIQUID_WATER_CAPACITY; /* water holding capacity of snow as a
				       fraction of snow-water-equivalent */ 
extern float MAX_SNOW_TEMP;	/* maximum temperature at which snow
				   can occur (C) */
extern float MIN_INTERCEPTION_STORAGE; /* the amount of snow on the canopy
                                          that can only be melted off. (m) */
extern float MIN_RAIN_TEMP;	/* minimum temperature at which rain
				   can occur (C) */
extern unsigned char OUTSIDEBASIN; /* Mask value indicating outside the basin */
extern float PRECIPLAPSE;       /* Precipitation lapse rate in m/timestep / m */
extern float TEMPLAPSE;         /* Temperature lapse rate in C/m */
extern int NWINDMAPS;		/* Number of wind maps in case the wind source
				   is MODEL */
extern float Z0_GROUND;		/* Roughness length for bare soil (m) */
extern float Z0_SNOW;		/* Roughness length for snow (m) */
extern float Zref;		/* Reference height (m) */

#endif

