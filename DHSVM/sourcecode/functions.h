/*
 * SUMMARY:      functions.h - header file for number of DHSVM functions 
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for large number of DHSVM functions 
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: functions.h, v 4.0  2013/1/10   Ning Exp $     
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "data.h"
#include "channel.h"
#include "DHSVMChannel.h"

void Aggregate(MAPSIZE *Map, OPTIONSTRUCT *Options, TOPOPIX **TopoMap,
	       LAYER *Soil, LAYER *Veg, VEGPIX **VegMap, EVAPPIX **Evap,
	       PRECIPPIX **Precip, RADCLASSPIX **RadMap, SNOWPIX **Snow,
	       SOILPIX **SoilMap, AGGREGATED *Total, VEGTABLE *VType,
	       ROADSTRUCT **Network, SEDPIX **SedMap, FINEPIX ***FineMap,
	       CHANNEL *ChannelData, float *roadarea);

void Alloc_Chan_Sed_Mem(float ** DummyVar);

void CalcAerodynamic(int NVegLayers, unsigned char OverStory,
		     float n, float *Height, float Trunk, float *U,
		     float *U2mSnow, float *Ra, float *RaSnow);

float CalcBagnold(float DS,TIMESTRUCT * Time, float, float, float, float);

double CalcDistance(COORD *LocA, COORD *LocB);

float CalcEffectiveKh(int NSoilLayers, float Top, float Bottom,
		      float *SoilDepth, float *KhDry, float *KhSol,
		      float *Moisture, float *Porosity, float *TSoil);

float CalcKhDry(float Density);

float CalcSafetyFactor(float Slope, int Soil, float SoilDepth, int Veg, 
		       SEDTABLE *SedType, VEGTABLE *VType, 
		       float M, SOILTABLE *SType, float Swq, float Depth,
		       int iter);

float CalcSnowAlbedo(float TSurf, unsigned short Last, SNOWTABLE *SnowAlbedo);

float CalcTransmissivity(float SoilDepth, float WaterTable, float LateralKs,
			 float KsExponent, float DepthThresh);

void CalcWeights(METLOCATION *Station, int NStats, int NX, int NY,
		 uchar **BasinMask, uchar ****WeightArray,
		 OPTIONSTRUCT *Options);

/* double cbrt (double x);*/ //compute the cubic root of a value

double ChannelCulvertSedFlow(int y, int x, CHANNEL * ChannelData, int i);

void CheckOut(int CanopyRadAttOption, LAYER Veg, LAYER Soil, 
	      VEGTABLE *VType, SOILTABLE *SType, MAPSIZE *Map, 
	      TOPOPIX **TopoMap, VEGPIX **VegMap, SOILPIX **SoilMap);

unsigned char dequal(double a, double b);

void draw(DATE *Day, int first, int DayStep, MAPSIZE *Map, int NGraphics,
	  int *which_graphics,
	  VEGTABLE *VType, SOILTABLE *SType,
	  SNOWPIX **SnowMap, SOILPIX **SoilMap, SEDPIX **SedMap, FINEPIX ***FineMap,
	  VEGPIX **VegMap, TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, float **PrismMap,
	  float **SkyViewMap, unsigned char ***ShadowMap, EVAPPIX **EvapMap,
	  RADCLASSPIX **RadMap, MET_MAP_PIX **MetMap, ROADSTRUCT **Network,
	  OPTIONSTRUCT *Options);

void DistributeSedimentDiams(float SedDiams[NSEDSIZES]);

void DumpMap(MAPSIZE *Map, DATE *Current, MAPDUMP *DMap, TOPOPIX **TopoMap,
	     EVAPPIX **EvapMap, PRECIPPIX **PrecipMap, RADCLASSPIX **RadMap,
	     SNOWPIX **Snowap, SOILPIX **SoilMap, SEDPIX **SedMap, FINEPIX ***FineMap,
	     LAYER *Soil, VEGPIX **VegMap, LAYER *Veg, ROADSTRUCT **Network,
	     OPTIONSTRUCT *Options);

void DumpPix(DATE *Current, int first, FILES *OutFile, EVAPPIX *Evap,
	     PRECIPPIX *Precip, RADCLASSPIX *Rad, SNOWPIX *Snow,
	     SOILPIX *Soil, int NSoil, int NVeg, OPTIONSTRUCT *Options);

void DumpPixSed(DATE *Current, int first, FILES *OutFileSediment, SEDPIX *SedMap,
             ROADSTRUCT *Network, float SedimentOverlandInflow, 
             float SedimentOverroadInflow, FINEPIX *FineMap);

void ExecDump(MAPSIZE * Map, DATE * Current, DATE * Start, OPTIONSTRUCT * Options,
	      DUMPSTRUCT * Dump, TOPOPIX ** TopoMap, EVAPPIX ** EvapMap,
	      PRECIPPIX ** PrecipMap, RADCLASSPIX ** RadMap, SNOWPIX ** SnowMap,
	      MET_MAP_PIX ** MetMap, VEGPIX ** VegMap, LAYER * Veg, SOILPIX ** SoilMap,
	      SEDPIX ** SedMap, ROADSTRUCT ** Network, CHANNEL * ChannelData, 
	      FINEPIX *** FineMap, LAYER * Soil, AGGREGATED * Total, 
	      UNITHYDRINFO * HydrographInfo, float *Hydrograph);

unsigned char fequal(float a, float b);

void FinalMassBalance(FILES *Out, AGGREGATED *Total, WATERBALANCE *Mass,
		      OPTIONSTRUCT * Options, float roadarea);

float FindDT(SOILPIX **SoilMap, MAPSIZE *Map, TIMESTRUCT *Time, 
	     TOPOPIX **TopoMap, SOILTABLE *SType); 

float FindDTRoad(ROADSTRUCT **Network, TIMESTRUCT *Time, int y, int x, 
		 float dx, float beta, float alpha);

void GenerateScales(MAPSIZE *Map, int NumberType, void **XScale,
		    void **YScale);

void GetMetData(OPTIONSTRUCT *Options, TIMESTRUCT *Time, int NSoilLayers,
		int NStats, float SunMax, METLOCATION *Stat, MAPSIZE *Radar,
		RADARPIX **RadarMap, char *RadarFileName);

uchar InArea(MAPSIZE *Map, COORD *Loc);

void InitAggregated(int MaxVegLayers, int MaxSoilLayers, AGGREGATED *Total);

int InitChannelSediment(Channel * Head, AGGREGATED *Total);

int InitChannelSedInflow(Channel * Head);

void InitCharArray(char *Array, int Size);

void InitConstants(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		   SOLARGEOMETRY *SolarGeo, TIMESTRUCT *Time);

void InitDump(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
	      int MaxSoilLayers, int MaxVegLayers, int Dt,
	      TOPOPIX **TopoMap, DUMPSTRUCT *Dump, int *NGraphics,
	      int **which_graphics);

void InitEvapMap(MAPSIZE *Map, EVAPPIX ***EvapMap, SOILPIX **SoilMap,
		 LAYER *Soil, VEGPIX **VegMap, LAYER *Veg,
		 TOPOPIX **TopoMap);

void InitFineMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, 
		     LAYER *Soil, TOPOPIX ***TopoMap, SOILPIX ***SoilMap, 
		  FINEPIX ****FineMap);

void InitImageDump(LISTPTR Input, int Dt, MAPSIZE *Map, int MaxSoilLayers,
		   int MaxVegLayers, char *Path, int NMaps, int NImages,
		   MAPDUMP **DMap);

void InitInFiles(INPUTFILES *InFiles);

void InitInterpolationWeights(MAPSIZE *Map, OPTIONSTRUCT *Options,
			      TOPOPIX **TopoMap, uchar ****MetWeights,
			      METLOCATION *Stats, int NStats);

void InitMapDump(LISTPTR Input, MAPSIZE *Map, int MaxSoilLayers,
		 int MaxVegLayers, char *Path, int TotalMapImages, int NMaps,
		 MAPDUMP **DMap);

void InitMassWaste(LISTPTR Input, TIMESTRUCT *Time);

void InitMetMaps(int NDaySteps, MAPSIZE *Map, MAPSIZE *Radar,
		 OPTIONSTRUCT *Options, char *WindPath, char *PrecipLapsePath,
		 float ***PrecipLapseMap, float ***PrismMap,
		 unsigned char ****ShadowMap, float ***SkyViewMap,
		 EVAPPIX ***EvapMap, PRECIPPIX ***PrecipMap,
		 RADARPIX ***RadarMap, RADCLASSPIX ***RadMap,
		 SOILPIX **SoilMap, LAYER *Soil, VEGPIX **VegMap,
		 LAYER *Veg, TOPOPIX **TopoMap, float ****MM5Input,
		 float ****WindModel);

void InitMetSources(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		    int NSoilLayers, TIMESTRUCT *Time, INPUTFILES *InFiles,
		    int *NStats, METLOCATION **Stat, MAPSIZE *Radar, 
		    MAPSIZE *MM5Map);

void InitMM5(LISTPTR Input, int NSoilLayers, TIMESTRUCT *Time,
	     INPUTFILES *InFiles, OPTIONSTRUCT *Options, MAPSIZE *MM5Map,
	     MAPSIZE *Map);

void InitMM5Maps(int NSoilLayers, int NY, int NX, float ****MM5Input,
		 RADCLASSPIX ***RadMap, OPTIONSTRUCT *Options);

void InitModelState(DATE *Start, MAPSIZE *Map, OPTIONSTRUCT *Options,
		    PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
		    SOILPIX **SoilMap, LAYER Soil, SOILTABLE *SType,
		    VEGPIX **VegMap, LAYER Veg, VEGTABLE *VType, char *Path,
		    SNOWTABLE *SnowAlbedo, TOPOPIX **TopoMap,
		    ROADSTRUCT **Network, UNITHYDRINFO *HydrographInfo,
		    float *Hydrograph);

void InitNetwork(int NY, int NX, float DX, float DY, TOPOPIX **TopoMap, 
		 SOILPIX **SoilMap, VEGPIX **VegMap, VEGTABLE *VType, 
		 ROADSTRUCT ***Network, CHANNEL *ChannelData, 
		 LAYER Veg, OPTIONSTRUCT *Options);

void InitNewDay(int DayOfYear, SOLARGEOMETRY *SolarGeo);

void InitNewMonth(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
		  TOPOPIX **TopoMap, float **PrismMap,
		  unsigned char ***ShadowMap, RADCLASSPIX **RadMap, 
		  INPUTFILES *InFiles, int NVegs, VEGTABLE *VType, int NStats,
		  METLOCATION *Stat, char *Path);

void InitNewStep(INPUTFILES *InFiles, MAPSIZE *Map, TIMESTRUCT *Time,
		 int NSoilLayers, OPTIONSTRUCT *Options, int NStats,
		 METLOCATION *Stat, char *RadarFileName, MAPSIZE *Radar,
		 RADARPIX **RadarMap, SOLARGEOMETRY *SolarGeo, 
		 TOPOPIX **TopoMap, RADCLASSPIX **RadMap, SOILPIX **SoilMap,
		 float ***MM5Input, float ***WindModel, MAPSIZE *MM5Map);

void InitParameters(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
		    ROADSTRUCT ***Network, CHANNEL *ChannelData, TOPOPIX **TopoMap,
		    TIMESTRUCT * Time, float *SedDiams);

int InitPixDump(LISTPTR Input, MAPSIZE *Map, uchar **BasinMask, char *Path,
		int NPix, PIXDUMP **Pix, OPTIONSTRUCT *Options);

void InitPrecipLapse(LISTPTR Input, INPUTFILES *InFiles);

void InitPrecipLapseMap(char *PrecipLapseFile, int NY, int NX,
			float ***PrecipLapseMap);

void InitPrismMap(int NY, int NX, float ***PrismMap);

void InitSurfaceSed(LISTPTR Input, TIMESTRUCT *Time);

void InitSedimentTables(int StepsPerDay, LISTPTR Input, SEDTABLE **SedType, SOILTABLE **SType,
			VEGTABLE **VType, LAYER *Soil, LAYER *Veg);

void InitSedMap(MAPSIZE *Map, SEDPIX ***SedMap);

void InitShadeMap(OPTIONSTRUCT *Options, int NDaySteps, int NY, int NX,
		  unsigned char ****ShadowMap, float ***SkyViewMap);

void InitPrecipMap(MAPSIZE *Map, PRECIPPIX ***PrecipMap, VEGPIX **VegMap,
		   LAYER *Veg, TOPOPIX **TopoMap);

void InitRadar(LISTPTR Input, MAPSIZE *Map, TIMESTRUCT *Time,
	       INPUTFILES *InFiles, MAPSIZE *Radar);

void InitRadarMap(MAPSIZE *Radar, RADARPIX ***RadarMap);

void InitRadMap(MAPSIZE *Map, RADCLASSPIX ***RadMap);

void InitSatVaporTable(void);

void InitSnowMap(MAPSIZE *Map, SNOWPIX ***SnowMap);

void InitSoilMap(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		 LAYER *Soil, TOPOPIX **TopoMap, SOILPIX ***SoilMap);

int InitSoilTable(SOILTABLE **SType, LISTPTR Input, LAYER *Soil);

void InitSnowTable(SNOWTABLE **SnowAlbedo, int Dt);

void InitStateDump(LISTPTR Input, int NStates, DATE **DState);

void InitGraphicsDump(LISTPTR Input, int NGraphics, int ***which_graphics);

void InitStations(LISTPTR Input, MAPSIZE *Map, int NDaySteps,
		  OPTIONSTRUCT *Options, int *NStats, METLOCATION **Stat);

void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options,
		SOILTABLE **SType, LAYER *Soil, VEGTABLE **VType,
		LAYER *Veg, SNOWTABLE **SnowAlbedo);

void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		     LAYER *Soil, TOPOPIX ***TopoMap, SOILPIX ***SoilMap,
		     VEGPIX ***VegMap);

void InitTopoMap(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		 TOPOPIX ***TopoMap);

void InitUnitHydrograph(LISTPTR Input, MAPSIZE *Map, TOPOPIX **TopoMap,
			UNITHYDR ***UnitHydrograph, float **Hydrograph,
			UNITHYDRINFO *HydrographInfo);

void InitVegMap( OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map, VEGPIX ***VegMap);

int InitVegTable(VEGTABLE **VType, LISTPTR Input, OPTIONSTRUCT *Options,
		 LAYER *Veg);

float evalexpint(int n, float x);

void InitWindModel(LISTPTR Input, INPUTFILES *InFiles, int NStats,
		   METLOCATION *Stat);

void InitWindModelMaps(char *WindPath, int NY, int NX, float ****WindModel);

uchar IsStationLocation(COORD *Loc, int NStats, METLOCATION *Station,
			int *WhichStation);

void InitXGraphics(int argc, char **argv,
		   int ny, int nx, int nd, MET_MAP_PIX ***MetMap);

float LapsePrecip(float Precip, float FromElev, float ToElev,
		  float PrecipLapse);

float LapseT(float Temp, float FromElev, float ToElev, float LapseRate);

void MainMWM(SEDPIX **SedMap, FINEPIX ***FineMap, VEGTABLE *VType, SEDTABLE *SedType,
	    CHANNEL *ChannelData, char *DumpPath, SOILPIX **SoilMap, TIMESTRUCT *Time,
	    MAPSIZE *Map, TOPOPIX **TopoMap, SOILTABLE *SType, VEGPIX **VegMap, 
	     int MaxStreamID, SNOWPIX **SnowMap);
 
PIXMET MakeLocalMetData(int y, int x, MAPSIZE *Map, int DayStep,
			OPTIONSTRUCT *Options, int NStats,
			METLOCATION *Stat, uchar *MetWeights,
			float LocalElev, RADCLASSPIX *RadMap,
			PRECIPPIX *PrecipMap, MAPSIZE *Radar,
			RADARPIX **RadarMap, float **PrismMap,
			SNOWPIX *LocalSnow, SNOWTABLE *SnowAlbedo,
			float ***MM5Input, float ***WindModel,
			float **PrecipLapseMap, MET_MAP_PIX ***MetMap,
			int NGraphics, int Month, float skyview,
			unsigned char shadow, float SunMax,
			float SineSolarAltitude);

void MassBalance(DATE *Current, FILES *Out, FILES *SedOut, AGGREGATED *Total,
		 WATERBALANCE *Mass, OPTIONSTRUCT * Options);

void MassEnergyBalance(int y, int x, float SineSolarAltitude, float DX, 
		       float DY, int Dt, int HeatFluxOption, 
		       int CanopyRadAttOption, int RoadRouteOption,
		       int InfiltOption,int MaxVegLayers,  PIXMET *LocalMet,
		       ROADSTRUCT *LocalNetwork, PRECIPPIX *LocalPrecip,
		       VEGTABLE *VType, VEGPIX *LocalVeg, SOILTABLE *SType,
		       SOILPIX *LocalSoil, SNOWPIX *LocalSnow,
		       EVAPPIX *LocalEvap, PIXRAD *TotalRad,
		       CHANNEL *ChannelData);

float MaxRoadInfiltration(ChannelMapPtr **map, int col, int row);

void OutputChannelSediment(Channel * Head, TIMESTRUCT Time, DUMPSTRUCT *Dump);

double pow (double a, double b);

void quick(ITEM *OrderedCells, int count);

void qs(ITEM *OrderedCells, int left, int right);

void ReadChannelState(char *Path, DATE *Current, Channel *Head);

void ReadMetRecord(OPTIONSTRUCT *Options, DATE *Current, int NSoilLayers,
		   FILES *InFile, unsigned char IsWindModelLocation,
		   MET *MetRecord);

void ReadRadarMap(DATE *Current, DATE *StartRadar, int Dt, MAPSIZE *Radar,
		  RADARPIX **RadarMap, char *HDFFileName);

void ReadPRISMMap(DATE *Current, DATE *StartRadar, int Dt, MAPSIZE *Radar,
		  RADARPIX **RadarMap, char *HDFFileName);

void ResetAggregate(LAYER *Soil, LAYER *Veg, AGGREGATED *Total,
                    OPTIONSTRUCT *Options);

void ResetValues(MAPSIZE *Map, SOILPIX **SoilMap);

int Round(double x);

void RouteChannelSediment(Channel * Head, TIMESTRUCT Time, 
			  DUMPSTRUCT *Dump, AGGREGATED * Total, float *SedDiams);

void RouteCulvertSediment(CHANNEL * ChannelData, MAPSIZE * Map,
			  TOPOPIX ** TopoMap, SEDPIX ** SedMap, 
			  AGGREGATED * Total, float *SedDiams);

void RouteDebrisFlow(float *SedimentToChannel, int prevy, 
		     int prevx, float SlopeAspect, CHANNEL *ChannelData, 
		     MAPSIZE *Map);

void RouteRoad(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap, 
	       SOILPIX ** SoilMap, ROADSTRUCT ** Network, SOILTABLE * SType, 
	       CHANNEL * ChannelData, PRECIPPIX ** PrecipMap, SEDPIX **SedMap,
	       float Tair, float Rh, float *SedDiams); 

void RouteSubSurface(int Dt, MAPSIZE *Map, TOPOPIX **TopoMap,
		     VEGTABLE *VType, VEGPIX **VegMap,
		     ROADSTRUCT **Network, SOILTABLE *SType,
		     SOILPIX **SoilMap, CHANNEL *ChannelData, 
		     TIMESTRUCT *Time, OPTIONSTRUCT *Options, 
		     char *DumpPath, SEDPIX **SedMap, FINEPIX ***FineMap,
		     SEDTABLE *SedType, int MaxStreamID, SNOWPIX **SnowMap);

void RouteSurface(MAPSIZE *Map, TIMESTRUCT *Time, TOPOPIX **TopoMap,
		  SOILPIX **SoilMap, OPTIONSTRUCT *Options,
		  UNITHYDR **UnitHydrograph,
		  UNITHYDRINFO *HydrographInfo, float *Hydrograph,
		  DUMPSTRUCT *Dump, VEGPIX **VegMap, VEGTABLE *VType,
		  SOILTABLE *SType, CHANNEL *ChannelData, SEDPIX **SedMap,
		  PRECIPPIX **PrecipMap, SEDTABLE *SedType, float Tair, 
		  float Rh, float *SedDiams);

float SatVaporPressure(float Temperature);

int SaveChannelSedInflow(Channel * Head, AGGREGATED * Total);

int ScanInts(FILE *FilePtr, int *X, int N);

int ScanDoubles(FILE *FilePtr, double *X, int N);

int ScanFloats(FILE *FilePtr, float *X, int N);

uchar ScanUChars(FILE *FilePtr, uchar *X, int N);

void SedimentFlag(OPTIONSTRUCT *Options,  TIMESTRUCT * Time);

void SkipHeader(FILES *InFile, int NLines);

void SkipLines(FILES *InFile, int NLines);

void StoreChannelState(char *Path, DATE *Current, Channel *Head);

void StoreModelState(char *Path, DATE * Current, MAPSIZE * Map,
		     OPTIONSTRUCT * Options, TOPOPIX ** TopoMap,
		     PRECIPPIX ** PrecipMap, SNOWPIX ** SnowMap,
		     MET_MAP_PIX ** MetMap, RADCLASSPIX ** RadMap,
		     VEGPIX ** VegMap, LAYER * Veg, SOILPIX ** SoilMap,
		     LAYER * Soil, ROADSTRUCT ** Network, 
		     UNITHYDRINFO * HydrographInfo, float *Hydrograph,
		     CHANNEL * ChannelData);

float viscosity(float Tair, float Rh);

#endif
