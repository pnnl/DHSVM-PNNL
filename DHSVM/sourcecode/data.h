/*
 * SUMMARY:      data.h - header file with data structures
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file with data structures
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
* $Id: data.h, v 3.1.1  2012/10/31   Ning Exp $ 
 */

#ifndef DATA_H
#define DATA_H

#include "settings.h"
#include "Calendar.h"
#include "channel.h"

typedef struct {
  int N;			/* Northing */
  int E;			/* Easting */
} COORD;

typedef struct {
  char FileName[BUFSIZE + 1];
  FILE *FilePtr;
} FILES;

typedef struct {
  int ID;			/* Index for variable to dump */
  int Layer;			/* Layer for which to dump */
  char Name[BUFSIZE + 1];	/* Variable Name */
  char LongName[BUFSIZE + 1];	/* Long name */
  char Format[BUFSIZE + 1];	/* Output format (for netcdf files) */
  char Units[BUFSIZE + 1];	/* Units */
  uchar Resolution;		/* Resolution at which to dump */
  int N;			/* Number of timesteps for which to dump */
  float MinVal;			/* Lowest value for indexing low resolution */
  float MaxVal;			/* Highest value for indexing low resolution */
  char FileName[BUFSIZE + 1];	/* File to write dump to */
  char FileLabel[BUFSIZE + 1];	/* File label */
  int NumberType;		/* Number type of variable */
  DATE *DumpDate;		/* Date(s) at which to dump */
} MAPDUMP;

typedef struct {

  COORD Loc;			/* Location for which to dump */
  FILES OutFile;		/* Files in which to dump */
  FILES OutFileSediment;	/* Files in which to dump - Sediment values */
} PIXDUMP;

typedef struct {
  char Path[BUFSIZE + 1];			/* Path to dump to */
  char InitStatePath[BUFSIZE + 1];	/* Path for initial state */
  FILES Aggregate;					/* File with aggregated values for entire basin */
  FILES AggregateSediment;			/* File with aggregated sediment values for entire basin */
  FILES Balance;					/* File with summed mass balance values for entire basin */
  FILES FinalBalance;               /* File with summed mass balance values for the entire simulation period for entire basin */
  FILES SedBalance;					/* File with summed mass balance values for entire basin */
  FILES Stream;
  int NStates;						/* Number of model state dumps */
  DATE *DState;						/* Array with dates on which to dump state */
  int NPix;							/* Number of pixels for which to output timeseries */
  PIXDUMP *Pix;						/* Array with info on pixels for which to output timeseries */
  int NMaps;						/* Number of variables for which to output maps */
  MAPDUMP *DMap;					/* Array with info on each map to output */
} DUMPSTRUCT;

typedef struct {
  float ETot;			/* Total amount of evapotranspiration */
  float *EPot;			/* Potential transpiration from each vegetation/soil layer */
  float *EAct;			/* Actual transpiration from each vegetation soil layer */
  float *EInt;			/* Evaporation from interception for each vegetation layer */
  float **ESoil;		/* Transpiration for each vegetation layer from each soil zone */
  float EvapSoil;		/* Evaporation from the upper soil layer */
} EVAPPIX;

typedef struct {
  int TimeStep;
  float Fraction;
} UNITHYDR;

typedef struct {
  int MaxTravelTime;
  int TotalWaveLength;
  int *WaveLength;
} UNITHYDRINFO;

typedef struct {
  char Const[BUFSIZE + 1];	/* Filename for main input file  */
  char RadMapPath[BUFSIZE + 1];	/* Path and start of filename for rad files */
  char RadTablePath[BUFSIZE + 1];	/* Same for rad tables */
  char RadarFile[BUFSIZE + 1];	/* File with radar precipitation */
  char MM5Terrain[BUFSIZE + 1];	/* File with MM5 terrain (m) */
  char MM5Lapse[BUFSIZE + 1];	/* File with MM5 Lapse Rate (C/m) */
  char MM5Temp[BUFSIZE + 1];	/* File with MM5 temperature (C) */
  char MM5Humidity[BUFSIZE + 1];	/* File with MM5 humidity (%) */
  char MM5Wind[BUFSIZE + 1];	/* File with MM5 wind speed (m/s) */
  char MM5ShortWave[BUFSIZE + 1];	/* File with MM5 shortwave (W/m2) */
  char MM5LongWave[BUFSIZE + 1];	/* File with MM5 longwave (W/m2) */
  char MM5Precipitation[BUFSIZE + 1];	/* File with MM5 precipitation 
					   (m/timestep) */
  char **MM5SoilTemp;		/* Files with MM5 soil temperatures (C) */
  char PrecipLapseFile[BUFSIZE + 1];	/* File with precipitation 
					   lapse rate map */
  char WindMapPath[BUFSIZE + 1];	/* File with wind factors */
} INPUTFILES;

typedef struct {
  int NTypes;
  int *NLayers;
  int MaxLayers;
} LAYER;

typedef struct {
  float Tair;			/* Air temperature (C) */
  float Rh;			    /* Relative humidity (%) */
  float Wind;			/* Wind (m/s) */
  float VICSin;         /* Observed Incoming shortwave used for RBM only (W/m2) */
  float Sin;			/* Incoming shortwave (W/m^2) */
  float SinBeam;		/* Incoming beam radiation (W/m^2) */
  float SinDiffuse;		/* Incoming diffuse radiation (W/m^2) */
  float Lin;			/* Incoming longwave (W/m^2) */
  float AirDens;		/* Air density on kg/m^3 */
  float Lv;			    /* Latent heat of vaporization (J/kg) */
  float Press;			/* Atmospheric pressure (Pa) */
  float Gamma;			/* Psychrometric constant (Pa/C) */
  float Es;			    /* Saturated vapor pressure (Pa) */
  float Eact;			/* Actual vapor pressure (Pa) */
  float Slope;			/* Slope of vapor pressure curve (Pa/C) */
  float Vpd;			/* Vapor pressure deficit (Pa) */
} PIXMET;

typedef struct {
  float Rank;
  int   x;
  int   y;
} ITEM;

typedef struct {
  char System[BUFSIZE + 1];		 /* Coordinate system */
  double Xorig;					 /* X coordinate of Northwest corner */
  double Yorig;					 /* Y coordinate of Northwest corner */
  int X;						 /* Current x position */
  int Y;						 /* Current y position */
  int NX;						 /* Number of pixels in x direction */
  int NY;						 /* Number of pixels in y direction */
  float DX;						 /* Pixel spacing in x-direction */
  float DY;						 /* Pixel spacing in y-direction */
  float DXY;					 /* Pixel spacing in diagonal */
  int OffsetX;					 /* Offset in x-direction compared to basemap */
  int OffsetY;					 /* Offset in y-direction compared to basemap */
  int NumCells;                  /* Number of cells within the basin */
  int NYfine;                    /* Number of pixels for mass wasting algorithm in x direction */
  int NXfine;                    /* Number of pixels for mass wasting algorithm in y direction */
  float DMASS;					 /* Pixel spacing for mass wasting algorithm */
  int NumCellsfine;              /* Number of cells for mass wasting algorithm within the basin */
  int NumFineIn;                 /* Number of fine cells in one coarse cell */  
  ITEM *OrderedCells;            /* Structure array to hold the ranked elevations; NumCells in size */
} MAPSIZE;

typedef struct {
  float Tair;			/* Air temperature (C) */
  float TempLapse;		/* Temperature lapse rate (C/m) */
  float Rh;				/* Relative humidity (%) */
  float Wind;			/* Wind (m/s) */
  int WindDirection;    /* Wind direction, used when WindSource == MODEL  */

  float Sin;			/* Incoming shortwave (W/m^2) */
  float SinBeamObs;		/* Observed incoming beam radiation (W/m^2) */
  float SinDiffuseObs;	/* Observed incoming diffuse radiation (W/m^2) */

  float SinBeamMod;		/* Modeled incoming beam radiation (W/m^2) */
  float SinDiffuseMod;	/* Modeled incoming diffuse radiation (W/m^2) */
  float BeamRatio;		/* Ratio of observed beam to modeled beam */
  float DiffuseRatio;	/* Ratio of observed diffuse to modeled diffuse */

  float Lin;			/* Incoming longwave (W/m^2) */
  float ClearIndex;		/* Cloudiness index */
  /* The following is a hack, and needs to be replaced by a better method, WORK IN PROGRESS */


  float Precip;			/* Rainfall if available (m) */
  float Tsoil[3];		/* Soil temperature in upper three layers */
  float PrecipLapse;	/* Elevation Adjustment Factor for Precip */
} MET;

typedef struct {
  char Name[BUFSIZE + 1];		/* Station name */
  COORD Loc;					/* Station locations */
  float Elev;					/* Station elevations */
  float PrismPrecip[12];		/* MonthlyPrism Precip for each station if outside=TRUE */
  uchar IsWindModelLocation;	/* Only used in case the wind model option is
								specified.  In that case this field is TRUE
								for one (and only one) station, and FALSE for all others */

  FILES MetFile;				/* File with observations */
  MET Data;
} METLOCATION;

typedef struct {
  int FileFormat;				/* File format indicator, BIN or HDF */
  int HasNetwork;				/* Flag to indicate whether roads and/or channels are imposed on the model area,
								TRUE if NETWORK, FALSE if UNIT_HYDROGRAPH */
  int CanopyRadAtt;				/* Radiation attenuation through the canopy, either FIXED (old method) or VARIABLE (based
								on Nijssen and Lettenmaier) */
  int PrecipType;				/* Precipitation source indicator, either RADAR or STATION */
  int Prism;					/* If TRUE, user supplied PRISM maps will be  used to interpolate precipitation */
  int PrecipLapse;				/* Whether the precipitation lapse rate is CONSTANT or VARIABLE */
  int TempLapse;				/* Whether the temperature lapse rate is CONSTANT or VARIABLE */
  int CressRadius;
  int CressStations;
  int WindSource;				/* Wind source indicator, either MODEL or STATION */
  int HeatFlux;					/* Specifies whether a sensible heat flux 
								should be calculated, TRUE or FALSE */
  int Routing;                   /* Overland flow routing indicator, either CONVENTIONAL or KINEMATIC */
  int OldRouteFlag;              /* Initial Overland flow routing indicator, either 
								 CONVENTIONAL or KINEMATIC */
  int Sediment;                  /* Specifies whether sediment is run and variables 
									are output, TRUE or FALSE */
  int MassWaste;                 /* Specifies whether mass wasting model should be run
								 and variables should be output, TRUE or FALSE */
  int SurfaceErosion;            /* Specifies whether surface erosion model should be run
								variables should be output, TRUE or FALSE */
  int ErosionPeriod;             /* Specifies dates when erosion model should be run
								 and variables should be output, TRUE or FALSE */
  int OldSedFlag;                /* Surface erosion flag from previous timestep */
  int InitSedFlag;               /* Initial Surface erosion flag for dumping purposes */
  int RoadRouting;               /* Road Erosion indicator, either TRUE or FALSE */
  int ChannelRouting;            /* Specifies whether sediment should be routed through
								 the channel network, either TRUE or FALSE */
  int Infiltration;              /* Specifies static or dynamic maximum infiltration rate */
  int FlowGradient;				 /* Specifies whether the flow gradient is based
								 on the terrain elevation (TOPOGRAPHY) or the 
								 water table elevation (WATERTABLE).  The 
								 TOPOGRAPHY method is much faster, since the 
								 flow direction and gradient do not have to 
								 be recalculated every timestep */
  int Extent;					/* Specifies the extent of the model run, either POINT or BASIN */
  int Interpolation;
  int MM5;						/* TRUE if MM5 interface is to be used, FALSE otherwise */
  int QPF;						/* TRUE if QPF override, else FALSE */
  int PointX;					/* X-index of point to model in POINT mode */
  int PointY;					/* Y-index of point to model in POINT mode */
  int Snotel;					/* if TRUE then station veg = bare for output */
  int Outside;					/* if TRUE then all listed met stats are used */
  int Rhoverride;				/* if TRUE then RH=100% if Precip>0 */
  int Shading;					/* if TRUE then terrain shading for solar is on */
  int StreamTemp;
  int CanopyShading;
  char SedFile[BUFSIZE+1];		/* Filename for sediment input file  */
  char PrismDataPath[BUFSIZE + 1];
  char PrismDataExt[BUFSIZE + 1];
  char ShadingDataPath[BUFSIZE + 1];
  char ShadingDataExt[BUFSIZE + 1];
  char SkyViewDataPath[BUFSIZE + 1];
  char ImperviousFilePath[BUFSIZ + 1];
} OPTIONSTRUCT;

typedef struct {
  float Precip;					/* Total amount of precipitation at pixel (m) */
  float SumPrecip;              /* Accumulated precipitation at pixel (m) */
  float RainFall;		        /* Amount of rainfall (m) */
  float SnowFall;		        /* Amount of snowfall (m) */
  float MomentSq;               /* Momentum squared for rain, used in the sediment model (kg* m/s)^2 /m^2*s) */
  float *IntRain;		        /* Rain interception by each vegetation layer (m) */
  float *IntSnow;		        /* Snow interception by each vegetation layer (m) */
  float TempIntStorage;			/* Temporary snow and rain interception storage, used by MassRelease() */
  int PrecipStart;              /* TRUE if there was surface water in the last time step */ 
  float Dm;                     /* Median raindrop diameter (m) */
 } PRECIPPIX;

typedef struct {
  float Precip;			/* Radar precipitation for current bin */
} RADARPIX;

typedef struct {
  float Beam;			/* Beam value */
  float Diffuse;		/* Diffuse value */
} RADCLASSPIX;

typedef struct {
  float NetShort[2];    /* Shortwave radiation for vegetation surfaces and ground/snow surface W/m2 */

  float LongIn[2];		/* Incoming longwave radiation for vegetation surfaces and ground/snow surface W/m2 */

  float LongOut[2];		/* Outgoing longwave radiation for vegetation surfaces and ground/snow surface W/m2 */

  float PixelNetShort;	/* Net shortwave for the entire pixel W/m2 */
  float RBMNetLong;     /* Longwave radiation reaching the water surface W/m2 (for RBM) */
  float RBMNetShort;
  float PixelLongIn;	/* Incoming longwave for entire pixel W/m2 */
  float PixelLongOut;	/* Outgoing longwave for entire pixel W/m2 */
  float ObsShortIn;
  float PixelShortIn;   /* Incoming shortwave for entire pixel W/m2 */
  float PixelBeam;      /* Net direct beam radiation W/m2 */
  float PixelDiffuse;   /* Net diffuse beam radiation W/m2 */
} PIXRAD;

typedef struct {
  float Area;			/* Area of road or channel cut (m) */
  float BankHeight;		/* Height of road or channel cut (m) */
  int   CutBankZone;	/* Number of the soil layer that contains the bottom of the road/channel cut */
  float *PercArea;		/* Area of percolation zone for each soil layer, corrected for the road/channel cut,
						divided by the grid cell area (0-1)  */
  float *Adjust;		/* Array with coefficients to correct for loss of soil storage due to channel/road-cut for each soil layer.
						Multiplied with RootDepth to give the zone thickness for use in calculating soil moisture */
  float MaxInfiltrationRate;	 /* Area weighted infiltration rate through the road bed */
  uchar fraction;				 /* flow fraction intercepted by road channel */
  float RoadArea;                /* Road surface area (and area of percolation)*/
  float IExcess;                 /* Infiltration excess generated on road surface (m)*/
  float FlowLength;              /* Representative surface water flow length across the road surface (m) */
  float FlowSlope;               /* Representative road surface slope along the flow path (m/m) */
  ChannelClass *RoadClass;       /* Class of road with most area in the pixel */
  float *h;                      /* Infiltration excess on road grid cell (m)*/
  float *startRunoff;            /* Surface water flux from the previus (sub) time step. Used for road kinematic wave routing.*/
  float *startRunon;             /* Surface water flux from the previus (sub) time step. Used for road kinematic wave routing.*/
  float *OldSedIn;               /* Sediment inflow to road cell from previous time step (m3/m3). */
  float *OldSedOut;              /* Sediment outflow from road cell from previous time step (m3/m3). */
  float Erosion;                 /* Change in road elevation/call area due to erosion (m/timestep). */
} ROADSTRUCT;

typedef struct {
  float SolarAzimuth;		/* solar azimuth */
  float Latitude;		    /* Latitude of center of study area */
  float Longitude;			/* Longitude of center of study area */
  float StandardMeridian;	/* Standard meridian for current time zone */
  float NoonHour;		    /* Time at which solar noon occurs for current location */

  float Declination;		/* Solar declination */
  float HalfDayLength;		/* Length of half day in hours */
  float Sunrise;		    /* Hour of sunrise */
  float Sunset;				/* Hour of sunset */
  float TimeAdjustment;		/* Time adjustment to be made between center of study area and standard meridian */

  float SunEarthDistance;	/* Distance from Sun to Earth */
  float SineSolarAltitude;	/* Sine of sun's SolarAltitude  */
  int DayLight;				/* FALSE: measured solar radiation and the sun is below the horizon.  

							TRUE: sun is above the horizon */
  float SolarTimeStep;		/* Fraction of the timestep the sun is above the horizon  */

  float SunMax;				/* Calculated solar radiation at the top of the atmosphere (W/m^2) */
} SOLARGEOMETRY;

typedef struct {
  uchar HasSnow;			/* Snow cover flag */
  uchar SnowCoverOver;		/* Flag overstory can be covered */
  unshort LastSnow;			/* Days since last snowfall */
  float Swq;				/* Snow water equivalent */
  float Melt;				/* Snow Melt */
  float Outflow;		    /* Snow pack outflow (m) */
  float PackWater;			/* Liquid water content of snow pack */
  float TPack;				/* Temperature of snow pack */
  float SurfWater;			/* Liquid water content of surface layer */
  float TSurf;				/* Temperature of snow pack surface layer */
  float ColdContent;		/* Cold content of snow pack */
  float Albedo;				/* Albedo of snow pack */
  float Depth;				/* Snow depth; Does not appear to be calculated
							or used anywhere */
  float VaporMassFlux;		/* Vapor mass flux to/from snow pack
				   (m/timestep) */
  float CanopyVaporMassFlux;	/* Vapor mass flux to/from intercepted snow in
				   the canopy (m/timestep) */
  float Glacier;		/*amount of snow added to glacier during simulation */
} SNOWPIX;

typedef struct {
  int   Soil;			/* Soil type */
  float Depth;			/* Depth of total soil zone, including all root
						zone layers, and the saturated zone */
  float *Moist;			/* Soil moisture content in layers (0-1) */
  float *Perc;			/* Percolation from layers */
  float *Temp;			/* Temperature in each layer (C) */
  float TableDepth;		/* Depth of water table below ground surface (m) */
  float WaterLevel;		/* Absolute height of the watertable above datum (m), 
						i.e. corrected for terrain elevation */
  float SatFlow;		/* amount of saturated flow generated */
  float IExcess;		/* amount of surface runoff (m) generated from HOF and Return flow */
  float Runoff;         /* Surface water flux (m) from the grid cell. */
  float ChannelInt;		/* amount of subsurface flow intercepted by the channel */
  float RoadInt;		/* amount of water intercepted by the road */
  float TSurf;			/* Soil surface temperature */
  float Qnet;			/* Net radiation exchange at surface */
  float Qrest;			/* Rest term for energy balance (should be 0) */
  float Qs;				/* Sensible heat exchange */
  float Qe;				/* Latent heat exchange */
  float Qg;				/* Ground heat exchange */
  float Qst;			/* Ground heat storage */
  float Ra;				/* Soil surface aerodynamic resistance (s/m) */
  float InfiltAcc;               /* Accumulated water in the top layer (m) */
  float MoistInit;               /* Initial moisture content when ponding begins (0-1) */
  float startRunoff;             /* Surface water flux from the previus (sub) time step. Used for kinematic wave routing.*/
  float startRunon;              /* Surface water flux from the previus (sub) time step. Used for kinematic wave routing.*/
  float IExcessSed;              /* amount of surface runoff (m) generated from HOF and Return flow  - saved for sediment routing */
  float DetentionStorage;        /* amount of water kept in detention storage when impervious fraction > 0 */
  float DetentionIn;			 /* detention storage change in current time step */
  float DetentionOut;            /* water flow out of detention storage */
} SOILPIX;

typedef struct {
  char Desc[BUFSIZE + 1];	/* Soil type */
  int Index;
  int NLayers;				/* Number of soil layers */
  float Albedo;				/* Albedo of the soil surface */
  float Manning;		    /* Manning's roughness of the soil surface */ 
  float *Porosity;			/* Soil porosity for each layer */
  float *PoreDist;			/* Pore size distribution for each layer */
  float *Press;				/* Soil bubbling pressure for each layer */
  float *FCap;				/* Field capacity for each layer  */
  float *WP;				/* Wilting point for each layer */
  float *Dens;				/* Soil density (kg/m^3) */
  float *Ks;				/* Saturated hydraulic conductivity (vertical) for each layer */
  float KsLat;				/* Saturated hydraulic conductivity (lateral) */
  float KsLatExp;		    /* Exponent for vertical change of KsLat */
  float *KhDry;				/* Thermal conductivity for dry soil (W/(m*K)) */
  float *KhSol;				/* Effective solids thermal conductivity (W/(M*K)) */
  float *Ch;				/* Heat capacity for soil medium */
  float MaxInfiltrationRate;/* Maximum infiltration rate for upper layer (m/s) */
  float G_Infilt;                /* Mean capillary drive for dynamic maximum infiltration rate (m)   */
  float DepthThresh;    /* Threshold water table depth, beyond which transmissivity decays linearly with water table depth */
} SOILTABLE;

typedef struct {
  float Freeze;			/* albedo when surface temperature below 0 C */
  float Thaw;			/* albedo when surface temperature above 0 C */
} SNOWTABLE;

typedef struct {
  char Distribution[BUFSIZE+1];	/* Distribution type */
  float mean;
  float stdev;
  float min;
  float max;
  float mode;
} STATSTABLE;


typedef struct {
  float Dem;					/* Elevations */
  uchar Mask;					/* Mask for modeled area */
  unshort Travel;				/* Travel time */
  float Grad;					/* Sum of downslope slope-width products */
  float Slope;					/* Land surface slope */
  float Aspect;					/* Land surface slope direction */
  float FlowGrad;				/* Magnitude of subsurface flow gradient slope * width */
  unsigned char Dir[NDIRS];		/* Fraction of surface flux moving in each direction*/
  unsigned int TotalDir;	    /* Sum of Dir array */
  int drains_x;					/* x-loc of cell to which this impervious cell drains */
  int drains_y;					/* y-loc of cell to which this impervious cell drains */
  ITEM *OrderedTopoIndex;       /* Structure array to hold the ranked topoindex for fine pixels in a coarse pixel */
} TOPOPIX;

typedef struct {
  int Veg;			/* Vegetation type */
  float Tcanopy;		        /* Canopy temperature (C) */
} VEGPIX;

typedef struct {
  char Desc[BUFSIZE + 1];	/* Vegetation type */
  int Index;
  int NVegLayers;		/* Number of vegetation layers */
  int NSoilLayers;		/* Number of soil layers */
  unsigned char OverStory;	        /* TRUE if there is an overstory */
  unsigned char UnderStory;	/* TRUE if there is an understory */
  float *Height;		/* Height of vegetation (m) */
  float *Fract;			/* Fractional coverage */
  float *HemiFract;		/* used to calculated longwave radiation balance */
  float *LAI;			/* One Sided Leaf Area Index */
  float **LAIMonthly;	/* Monthly LAI (one-sided) */
  float *MaxInt;		/* Maximum interception storage (m) */
  float *RsMax;			/* Maximum stomatal resistance */
  float *RsMin;			/* Minimum stomatal resistance */
  float *MoistThres;	/* Soil moisture threshold above which soil 
						moisture does not restrict transpiration */
  float *VpdThres;		/* Vapor pressure deficit threshold above which
						stomatal closure occurs (Pa) */
  float **RootFract;	/* Fraction of roots in each soil layer */
  float *RootDepth;		/* Depth of root zones */
  float Atten;			/* Canopy attenuation for radiation, only used
						when the "canopy radiation attenuation" 
						option is set to fixed */
  float TotalDepth;		/* total depth of all root zones */
  float ClumpingFactor;	/* To convert LAI of overstory to Effective LAI
					    for canopy attenuation of shortwave radiation 
						taken after Chen and Black, 1991 */
  float Taud;			/* Transmission of Diffuse radiation through canopy */
						/* a function of the following two parameters and
						effective LAI (which can change monthly) */
  float LeafAngleA;		/* parameter describing the Leaf Angle Distribution */
  float LeafAngleB;		/* parameter describing the leaf Angle Distribution */
  float Scat;			/* scattering parameter (between 0.7 and 0.85) */
  float *Rpc;			/* fraction of radiaton that is photosynthetically active (PAR) */
  float *Albedo;		/* Albedo for each vegetation layer */
  float **AlbedoMonthly;
  float Cn;				/* Canopy attenuation coefficient for wind profile */
  float MaxSnowInt;		/* Maximum snow interception capacity for the overstory */
  float MDRatio;		/* Ratio of Mass Release to Meltwater drip from int snow */
  float SnowIntEff;		/* Efficiency of snow interception process */
  float ImpervFrac;		/* fraction of pixel that is impervious */
  float DetentionFrac;  /* fraction of flow on impervious area goes to flood detention storage */ 
  float DetentionDecay; /* Decay coefficient of linear reservoir storage for impervious surface detention facility.  */
  float Ra[2];			/* Aerodynamic resistance in the absence of snow  */
  float RaSnow;			/* Aerodynamic resistance for the lower boundary in the presence of snow */
  float Trunk;			/* Fraction of overstory height that identifies the top of the trunk space */
  float U[2];			/* Wind speed profile (m/s) */
  float USnow;			/* wind speed 2, above snow surface (m/s) */
  STATSTABLE RootCoh;   /* Used for the mass wasting model. */
  STATSTABLE VegSurcharge;      /* Used for the mass wasting model. */
} VEGTABLE;

typedef struct {
  float StartWaterStorage;
  float OldWaterStorage;
  float CumPrecipIn;
  float CumET;
  float CumIExcess;
  float CumChannelInt;
  float CumRoadInt;
  float CumSnowVaporFlux;
  float CumCulvertReturnFlow;
  float CumCulvertToChannel;
  float CumRunoffToChannel;
  float StartChannelSedimentStorage;
  float LastChannelSedimentStorage;
  float CumMassWasting;
  float CumSedimentToChannel;
  float CumMassDeposition;
  float CumSedimentErosion;
  float CumRoadErosion;
  float CumRoadSedHill;
  float CumDebrisInflow;
  float CumSedOverlandInflow;
  float CumCulvertSedToChannel;
  float CumSedOverroadInflow;
  float CumSedimentOutflow;
  float CumCulvertReturnSedFlow;
} WATERBALANCE;

typedef struct {
  float accum_precip;
  float air_temp;
  float wind_speed;
  float humidity;
} MET_MAP_PIX;

typedef struct {
  float SedFluxOut;              /* Time step total sediment flux from the grid cell (m3). */
  float OldSedIn;                /* Sediment inflow to grid cell from previous time step (m3/m3). */
  float OldSedOut;               /* Sediment outflow from grid cell from previous time step (m3/m3). */
  float Erosion;                 /* Change in grid cell elevation due to erosion (mm/timestep). */
  float RoadSed;                 /* Time step total sediment flux from the road surface to hillslope (m3). */
} SEDPIX;

typedef struct {
  char Desc[BUFSIZE + 1];	/* Soil type */
  float KIndex;                  /* Index of soil detachability via raindrop impact 
				    (1/J) */
  STATSTABLE Cohesion;		/* Soil cohesion (kPa)  */
  STATSTABLE Friction;		/* Angle of internal friction (degrees)*/	
  float SatDensity;	        /* Saturated density for each layer (kg/m3) */
  float d50;                     /* Median grainsize diameter for surface erosion (mm) */ 
} SEDTABLE;

typedef struct node node;
struct node {
  node *next;
  int x;
  int y;
};

typedef struct {
  float Dem;                     /* Elevations */
  uchar Mask;                   /* Mask for modeled area */
  float bedrock;                 /* Bedrock elevation (m) */
  float sediment;                /* Sediment thickness (m) */
  float SatThickness;            /* Water table thickness (m) */
  float DeltaDepth;              /* Change in sediment thickness (m) */
  float Probability;             /* Pixel failure probability. */
  float MassWasting;             /* Sediment (m3) lost due to mass wasting */
  float MassDeposition;          /* Sediment (m3) deposited in grid cell from mass wasting elsewhere */
  float SedimentToChannel;       /* Sediment (m3) deposited in channel from mass wasting */
  float TopoIndex;               /* Topographic Index used for soil moisture redistribution from coarse 
								grid to fine grid */
} FINEPIX; 
 
typedef struct {
  EVAPPIX Evap;
  PRECIPPIX Precip;
  PIXRAD Rad;
  RADCLASSPIX RadClass;
  ROADSTRUCT Road;
  SNOWPIX Snow;
  SOILPIX Soil;
  SEDPIX Sediment;
  FINEPIX Fine;
  float SoilWater;
  float CanopyWater;
  float Runoff;
  float ChannelInt;
  float RoadInt;
  unsigned long Saturated;
  float CulvertReturnFlow;
  float CulvertToChannel;
  float RunoffToChannel;
  float DebrisInflow;
  float SedimentOverlandInflow;
  float SedimentOverroadInflow;
  float ChannelSedimentStorage;
  float ChannelSuspendedSediment;
  float CulvertReturnSedFlow;
  float CulvertSedToChannel;
  float SedimentOutflow;
} AGGREGATED;

#endif
