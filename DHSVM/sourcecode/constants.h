/*
 * SUMMARY:      constants.h - header file with constants for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file with constants for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: constants.h,v 1.13 2004/08/18 01:01:34 colleen Exp $     
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define CELLFACTOR 3          /* For discretization of forest roads
				     for kinematic wave routing */

#define CH_ICE     (2100.0e3)	/* Volumetric heat capacity (J/(m3*C) of ice 
				   (0C) */

#define CH_WATER   (4186.8e3)	/* Volumetric heat capacity (J/(m3*C) of 
				   water */

#define CP          1013.0	/* Specific heat of moist air at constant 
				   pressure (J/(kg*C)) */

#define DELTAT        50.	/* Used in SensibleHeatFlux to bracket the 
				   effective surface temperature (C) */

#define DEGPRAD       57.29578	/* degree per radian */

#define D0_MULTIPLIER  0.63	/* Multiplier for vegetation height to get
				   displacement height (m) */

#define DZ_TOP         0.1	/* Thickness of soil surface layer for which 
				   heat stoarge change is calculated (m) */

#define EPS            0.622	/* ratio of molecular weight of water vapor to
				   that for dry air */

#define FS_CRITERIA    1.0       /* Criteria for determining is a pixel is failed. 
				    Pixels with factor of safety's less than this number,
				    are considered failed */

#define G              9.81	/* gravitational accelleration (m/(s^2)) */

#define GRAMSPKG    1000.	/* grams per kilogram */

#define JOULESPCAL     4.1868	/* Joules per calorie */

#define KhH2O 0.58		/* Thermal conductivity of water (W/(mk)) */

#define LEAF_DRIP_DIA  0.0055    /* Leaf drip diameter (m) */

#define LF           (333.7e3)	/* latent heat of fusion (J/kg) */

#define MINPDEG        4.	/* minutes per degree longitude */

#define MMTOM          0.001	/* convert from mm to meter */

#define MTHRESH 0.85            /* The 'critical' value of M for triggering the mass
				   wasting algorithm. */

#define PARTDENSITY 2685.         /* Particle density in kg/m3 */

#undef PI
#define PI             3.14159265358979323846

#define RADPHOUR       0.2617994	/* radians per hour: Earth's Rotation 
					   (2 PI rad/day) * (1 day/24 h) */

#define RADPDEG     (PI/180.0)	/* radians per degree */

#define ROADCROWN 0.02            /* This is the road crown slope, whether insloped,  
				     outsloped or crowned. This value was selected based 
				     on the Road Preconstruction Handbook */

#define SATPERCENT 0.2          /* Fraction of pixels which must exceed 
				   MTHRESH in order to call the mass wasting algorithm
				   when running in the old mode. */

#define SEDEXPONENT 2.0          /* Exponent for exponential decrease in rainfall detachment 
				    of soil particles with depth of overland flow, should 
				    be between 0.9 and 3.1 */

#define SETTLECRIT 0.004         /*JSL critical streampower threshold for sediment transport 
				   in m/s from KINEROS2 Documentation */

#define SOLARCON    1360.	        /* Solar constant (W/m^2) */

#define STEFAN    (5.6696e-8)	/* Stefan-Boltzmann constant (W/(M^2*C^4) */

#define TIMEWEIGHT 0.65          /* For kinematic wave surface routing */   

#define VISCOSITY 1.0            /* kinematic viscosity of water in mm2/s */

#define VISFRACT       0.5	/* part of shortwave that is in the visible 
				   range of the spectrum */

#define VON_KARMAN     0.4	/* Von Karman's constant */

#define WATER_DENSITY 1000.	/* Density of water in kg/m3 */

#define Z0_MULTIPLIER  0.13	/* Multiplier for vegetation height to get
				   roughness length (m) */

/**************** extern constants - see globals.c ****************/

extern float LAI_SNOW_MULTIPLIER;	/* multiplier to calculate the amount
					   of available snow interception as 
					   a function of LAI */
extern float LAI_WATER_MULTIPLIER;	/* multiplier to determine maximum
					   interception storage as a function of
					   LAI  */
extern float LIQUID_WATER_CAPACITY;	/* water holding capacity of snow as a
					   fraction of snow-water-equivalent */
extern float MAX_SNOW_TEMP;	/* maximum temperature at which snow
				   can occur (C) */
extern float MIN_INTERCEPTION_STORAGE;	/* the amount of snow on the canopy
					   that can only be melted off. (m) */
extern float MIN_RAIN_TEMP;	/* minimum temperature at which rain
				   can occur (C) */
extern unsigned char OUTSIDEBASIN;	/* Mask value indicating outside the basin */
extern float PRECIPLAPSE;	/* Precipitation lapse rate in m/timestep / m */
extern float TEMPLAPSE;		/* Temperature lapse rate in C/m */
extern int NWINDMAPS;		/* Number of wind maps in case the wind source
				   is MODEL */
extern float Z0_GROUND;		/* Roughness length for bare soil (m) */
extern float Z0_SNOW;		/* Roughness length for snow (m) */
extern float Zref;		/* Reference height (m) */
extern float MASSITER;           /* Maximum number of iterations. */
extern float DEBRISd50;          /* (mm) */
extern float DEBRISd90;          /* (mm) */ 
extern float CHANNELd50;         /* Currently not used */
extern float CHANNELd90;         /* Currently not used */
#endif
