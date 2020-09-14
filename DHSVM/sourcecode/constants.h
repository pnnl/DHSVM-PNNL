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
 * $Id: constants.h,v 3.1.1 2012/10/18 ning Exp $     
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define CELLFACTOR 3            /* For discretization of forest roads for kinematic wave routing */
#define CELL_PARTITION 2        /* Number of veg type in a grid cell */
#define CH_ICE     (2100.0e3)	/* Volumetric heat capacity (J/(m3*C) of ice (0C) */
#define CH_WATER   (4186.8e3)	/* Volumetric heat capacity (J/(m3*C) of water */
#define CP         1013.0		/* Specific heat of moist air at constant pressure (J/(kg*C)) */
#define DELTAT     50.		    /* Used in SensibleHeatFlux to bracket the effective surface temperature (C) */
#define DEGPRAD     57.29578	/* degree per radian */
#define D0_MULTIPLIER  0.63		/* Multiplier for vegetation height to get displacement height (m) */
#define DZ_TOP      0.1			/* Thickness of soil surface layer for which heat stoarge change is calculated (m) */
#define EPS         0.622		/* ratio of molecular weight of water vapor to that for dry air */
#define FS_CRITERIA 1.0			/* Criteria for determining is a pixel is failed. 
                                   Pixels with factor of safety's less than this number, are considered failed */
#define G           9.81		/* gravitational accelleration (m/(s^2)) */
#define GRAMSPKG    1000.	    /* grams per kilogram */
#define JOULESPCAL  4.1868	    /* Joules per calorie */
#define KhH2O       0.58		/* Thermal conductivity of water (W/(mk)) */
#define LEAF_DRIP_DIA  0.0055   /* Leaf drip diameter (m) */
#define LF            (333.7e3) /* latent heat of fusion (J/kg) */
#define MINPDEG        4.		/* minutes per degree longitude */
#define MMTOM          0.001	/* convert from mm to meter */
#define MTHRESH        0.85     /* The 'critical' value of M for triggering the mass wasting algorithm. */
#define PARTDENSITY    2685.    /* Particle density in kg/m3 */
#undef PI
#define PI   3.14159265358979323846
#define RADPHOUR    0.2617994   /* radians per hour: Earth's Rotation (2 PI rad/day) * (1 day/24 h) */
#define RADPDEG     (PI/180.0)	/* radians per degree */
#define ROADCROWN   0.02        /* This is the road crown slope, whether insloped,  
								outsloped or crowned. This value was selected based 
								on the Road Preconstruction Handbook */
#define SATPERCENT  0.2         /* Fraction of pixels which must exceed 
								MTHRESH in order to call the mass wasting algorithm
								when running in the old mode. */
#define SOLARCON    1360.	    /* Solar constant (W/m^2) */
#define STEFAN    (5.6696e-8)	/* Stefan-Boltzmann constant (W/(M^2*C^4) */
#define VISFRACT    0.5		    /* part of shortwave that is in the visible range of the spectrum */
#define VON_KARMAN  0.4		    /* Von Karman's constant */
#define WATER_DENSITY 1000.		/* Density of water in kg/m3 */
#define Z0_MULTIPLIER 0.13		/* Multiplier for vegetation height to get roughness length (m) */
#define MinDiff   (1.e-8)

/**************** extern constants - see globals.c ****************/

extern int NDIRS;               /* How many neighbors are used in surface/subsurface routing */
extern int xdirection4[];
extern int ydirection4[];
extern int xdirection8[];
extern int ydirection8[];
extern int *xdirection, *ydirection;

extern float LAI_SNOW_MULTIPLIER;		/* multiplier to calculate the amount of available
										snow interception as a function of LAI */
extern float LAI_WATER_MULTIPLIER;		/* multiplier to determine maximum
										interception storage as a function of LAI  */
extern float LIQUID_WATER_CAPACITY;		/* water holding capacity of snow as a
										fraction of snow-water-equivalent */
extern float MAX_SNOW_TEMP;				/* maximum temperature at which snow can occur (C) */
extern float MIN_INTERCEPTION_STORAGE;	/* the amount of snow on the canopy
										that can only be melted off. (m) */
extern float MIN_RAIN_TEMP;				/* minimum temperature at which rain can occur (C) */
extern unsigned char OUTSIDEBASIN;		/* Mask value indicating outside the basin */
extern float PRECIPLAPSE;				/* Precipitation lapse rate in m/timestep / m */
extern float MINELEV;
extern float TEMPLAPSE;					/* Temperature lapse rate in C/m */
extern int NWINDMAPS;					/* Number of wind maps in case the wind source is MODEL */
extern float Z0_GROUND;					/* Roughness length for bare soil (m) */
extern float Z0_SNOW;					/* Roughness length for snow (m) */
extern float Zref;						/* Reference height (m) */

/* snow albedo decay curve */
extern float ALB_MAX;                   /* fresh snow albedo */                                                               
extern float ALB_ACC_LAMBDA;            /* snow freeze albedo cruve control parameters */
extern float ALB_MELT_LAMBDA;           /* snow thaw albedo cruve control parameters */
extern float ALB_ACC_MIN;
extern float ALB_MELT_MIN;
extern float PRECIP_MULTIPLIER;        /* precipitatio multiplier */
extern float MAX_SURFACE_SWE; 	       /* maximum depth of the surface layer in snow water equivalent (m) */

extern float GAPWIND_FACTOR;
extern int TotNumGap;                  /* total number of grid cells with a gap structure */

extern float SNOWSLIDE1;               /* First Parameter in Snowslide equation */
extern float SNOWSLIDE2;               /* Second Parameter in Snowslide equation */
#endif
