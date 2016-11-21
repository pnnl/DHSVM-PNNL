/*
 * SUMMARY:      functions.h - header file for number of DHSVM functions 
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD: Fri Jul 25 13:28:12 1997 by Bart Nijssen <nijssen@u.washington.edu>
 * DESCRIPTION:  header file for large number of DHSVM functions 
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 */

/* 	$Id: functions.h,v 1.13 1997/04/18 00:43:07 nijssen Exp $	 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "data.h"
#include "DHSVMChannel.h"

void Aggregate(MAPSIZE Map, OPTIONSTRUCT *Options, TOPOPIX **TopoMap, 
	       LAYER Soil, LAYER Veg, VEGPIX **VegMap, EVAPPIX **Evap, 
	       PRECIPPIX **Precip, RADCLASSPIX **RadMap, SNOWPIX **Snow, 
	       SOILPIX **SoilMap, AGGREGATED *Total, VEGTABLE *VType, 
	       ROADSTRUCT **Network);

double CalcDistance(COORD LocA, COORD LocB);

float CalcSnowAlbedo(float TSurf, unsigned short Last, SNOWTABLE *SnowAlbedo);


void DumpMap(MAPSIZE Map, TIMESTRUCT Time, MAPDUMP DMap, EVAPPIX **EvapMap, 
	     PRECIPPIX **PrecipMap, RADCLASSPIX **RadMap, SNOWPIX **SnowMap, 
	     SOILPIX **SoilMap, LAYER Soil, VEGPIX **VegMap, LAYER Veg);

void DumpPix(DATE Current, FILES OutFile, EVAPPIX EvapMap, 
	     PRECIPPIX PrecipMap, RADCLASSPIX RadMap, SNOWPIX SnowMap, 
	     SOILPIX SoilMap, int NSoil, int NVeg);

void ExecDump(MAPSIZE Map, TIMESTRUCT Time, OPTIONSTRUCT *Options, 
              DUMPSTRUCT Dump, TOPOPIX **TopoMap, EVAPPIX **EvapMap, 
              PRECIPPIX **PrecipMap, RADCLASSPIX **RadMap, SNOWPIX **SnowMap, 
              VEGPIX **VegMap, LAYER Veg, SOILPIX **SoilMap, LAYER Soil, 
              AGGREGATED Total, UNITHYDRINFO *HydrographInfo,
	      float *Hydrograph);

void FinalMassBalance(FILES Out, AGGREGATED Total, WATERBALANCE Mass);

void GenerateScales(MAPSIZE Map, int NumberType, void **XScale, void **YScale);

uchar InArea(MAPSIZE Map, COORD Loc);

uchar InBasin(uchar MaskValue);

void IncreaseTime(TIMESTRUCT *Time);

void InitAggregated(int MaxVegLayers, int MaxSoilLayers, AGGREGATED *Total);

void InitConstants(char *InFileName, OPTIONSTRUCT *Options, MAPSIZE *Map, 
		   SOLARGEOMETRY *SolarGeo, TIMESTRUCT *Time);

void InitDump(char *InFileName, OPTIONSTRUCT *Options, MAPSIZE *Map, 
	      int MaxSoilLayers, int MaxVegLayers, TIMESTRUCT *Time, 
	      TOPOPIX **TopoMap, DUMPSTRUCT *Dump);

void InitEvapMap(MAPSIZE Map, EVAPPIX ***EvapMap, SOILPIX **SoilMap, 
		 LAYER Soil, VEGPIX **VegMap, LAYER Veg, TOPOPIX **TopoMap);

void InitImageDump(char *InFileName, TIMESTRUCT *Time, int MaxSoilLayers, 
		   int MaxVegLayers, char *Path, int NMaps, int NImages, 
		   MAPDUMP **DMap);

void InitInFiles(INPUTFILES *InFiles);

void InitMapDump(char *InFileName, TIMESTRUCT *Time, int MaxSoilLayers, 
		 int MaxVegLayers, char *Path, int TotalMapImages,
		 int NMaps, MAPDUMP **DMap);

void InitMetMaps(DATE Start, int NDaySteps, MAPSIZE Map, MAPSIZE Radar, 
                 OPTIONSTRUCT *Options, char *WindPath, EVAPPIX ***EvapMap, 
		 PRECIPPIX ***PrecipMap, RADARPIX ***RadarMap, 
		 RADCLASSPIX ***RadMap, SOILPIX **SoilMap, LAYER Soil, 
		 VEGPIX **VegMap, LAYER Veg, TOPOPIX **TopoMap, 
		 float ****MM5Input, float ****WindModel);

void InitMM5(char *InFileName, int NSoilLayers, TIMESTRUCT *Time, 
	     INPUTFILES *InFiles);

void InitMM5Maps(int NSoilLayers, int NY, int NX, float ****MM5Input,
		 RADCLASSPIX ***RadMap);

void InitWindModelMaps(char *WindPath, int NY, int NX, float ****WindModel);

void InitModelState(TIMESTRUCT *Time, MAPSIZE *Map, OPTIONSTRUCT *Options,
		    PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, 
		    SOILPIX **SoilMap, LAYER Soil, SOILTABLE *SType, 
		    VEGPIX **VegMap, LAYER Veg, VEGTABLE *VType, char *Path, 
		    SNOWTABLE *SnowAlbedo, TOPOPIX **TopoMap, 
		    ROADSTRUCT **Network, UNITHYDRINFO *HydrographInfo,
		    float *Hydrograph);

void InitNetwork(int HasNetwork, int NY, int NX, float DX, TOPOPIX **TopoMap,
                 SOILPIX **SoilMap, VEGPIX **VegMap, VEGTABLE *VType,
                 ROADSTRUCT ***Network, CHANNEL *ChannelData);

void InitNewDay(int DayOfYear, SOLARGEOMETRY *SolarGeo);


void InitPixDump(char *InFileName, MAPSIZE *Map, uchar **BasinMask, char *Path, 
		 int NPix, PIXDUMP **Pix);

void InitPrecipMap(MAPSIZE Map, PRECIPPIX ***PrecipMap, VEGPIX **VegMap, 
		   LAYER Veg, TOPOPIX **TopoMap);

void InitRadar(char *InFileName, MAPSIZE *Map, TIMESTRUCT *Time, 
	       INPUTFILES *InFiles, MAPSIZE *Radar);

void InitRadarMap(MAPSIZE Radar, RADARPIX ***RadarMap);

void InitRadMap(unsigned char RadType, DATE Start, int NDaySteps,
                MAPSIZE Map, TOPOPIX **TopoMap, RADCLASSPIX ***RadMap);

void InitSatVaporTable(void);

void InitSnowMap(MAPSIZE Map, SNOWPIX ***SnowMap);

void InitSoilMap(char *InFileName, MAPSIZE *Map, LAYER *Soil, 
		 TOPOPIX **TopoMap, SOILPIX ***SoilMap);

int InitSoilTable(SOILTABLE **SType, char *InFileName, LAYER *Soil);

void InitSnowTable(SNOWTABLE **SnowAlbedo, int Dt);

void InitStateDump(char *InFileName, TIMESTRUCT *Time, int NStates, 
		   DATE **DState);

void InitStations(char *InFileName, MAPSIZE *Map, int NDaySteps, 
                  unsigned char RadType, int *NStats, METLOCATION **Stat);

void InitTables(int Dt, char *InFileName, SOILTABLE **SType, LAYER *Soil, 
		VEGTABLE **VType, LAYER *Veg, SNOWTABLE **SnowAlbedo);

void InitTerrainMaps(char *InFileName, OPTIONSTRUCT *Options, MAPSIZE *Map, 
		     LAYER *Soil, TOPOPIX ***TopoMap, SOILPIX ***SoilMap, 
		     VEGPIX ***VegMap);

void InitTopoMap(char *InFileName, OPTIONSTRUCT *Options, MAPSIZE *Map, 
		 TOPOPIX ***TopoMap);

void InitUnitHydrograph(char *InFileName, MAPSIZE *Map, TOPOPIX **TopoMap,
			UNITHYDR ***UnitHydrograph, float **Hydrograph, 
			UNITHYDRINFO *HydrographInfo);

void InitVegMap(char *InFileName, MAPSIZE *Map, VEGPIX ***VegMap);

int InitVegTable(VEGTABLE **VType, char *InFileName, LAYER *Veg);

void InitWindModel(char *InFileName, INPUTFILES *InFiles, int NStats,
		   METLOCATION *Stat);

uchar IsMultiLayer(int ID, int MaxSoilLayers, int MaxVegLayers, int
                   *MaxLayers); 

uchar IsStationLocation(COORD Loc, int NStats, METLOCATION *Station, 
			int *WhichStation);

uchar IsValidDumpID(int ID);

float LapsePrecip(float Precip, float FromElev, float ToElev, float PrecipLapse);

float LapseT(float Temp, float FromElev, float ToElev, float LapseRate);

PIXMET MakeLocalMetData(int y, int x, MAPSIZE Map, int DayStep, 
                        OPTIONSTRUCT *Options, int NStats, 
                        METLOCATION *Stat, uchar *MetWeights, 
                        float LocalElev, RADCLASSPIX *RadMap, 
                        PRECIPPIX *PrecipMap, MAPSIZE Radar, 
                        RADARPIX **RadarMap, SNOWPIX *LocalSnow,
                        SNOWTABLE *SnowAlbedo, float ***MM5Input,
                        float ***WindModel);

void MakeVarAttr(int ID, char *Var, char *FileLabel, int *NumberType, 
		 int Layer);

void MassBalance(DATE Current, FILES Out, AGGREGATED Total, WATERBALANCE *Mass);


void MassEnergyBalance(int y, int x, float DX, float Dt, int HeatFluxOption,
                       int MaxVegLayers, PIXMET LocalMet, 
                       ROADSTRUCT LocalNetwork, PRECIPPIX *LocalPrecip, 
                       VEGTABLE *VType, VEGPIX *LocalVeg, SOILTABLE SType, 
                       SOILPIX *LocalSoil, SNOWPIX *LocalSnow, 
                       EVAPPIX *LocalEvap, PIXRAD *TotalRad);

void ReadMetRecord(OPTIONSTRUCT *Options, DATE Current, int NSoilLayers, 
                   FILES InFile, unsigned char IsWindModelLocation, 
                   MET *MetRecord);

void ReadNetwork(void);

void ReadRadarMap(DATE Current, DATE StartRadar, int Dt, MAPSIZE Radar, 
		  RADARPIX **RadarMap, char *HDFFileName);

void ReadRadMap(MAPSIZE Map, int EndStep, int NDaySteps, int NStats, 
                char *FileName, TOPOPIX **TopoMap, RADCLASSPIX **RadMap, 
                METLOCATION *Stat);
 
int ReadStatInfo(FILES *InFile, uchar PrecipType, uchar RadType, 
                 int NDaySteps, METLOCATION **Stat, MAPSIZE Map); 

void ResetAggregate(LAYER Soil, LAYER Veg, AGGREGATED *Total);

void ResetValues(MAPSIZE Map, SOILPIX **SoilMap);

void RouteSubSurface(float Dt, MAPSIZE Map, TOPOPIX **TopoMap, 
		     VEGTABLE *VType, VEGPIX ** VegMap, 
                     ROADSTRUCT **Network, SOILTABLE *SType, 
                     SOILPIX **SoilMap, CHANNEL *ChannelData);

void RouteSurface(MAPSIZE Map, TIMESTRUCT TIME, TOPOPIX **TopoMap, 
                  SOILPIX **SoilMap, int HasNetwork, 
                  UNITHYDR **UnitHydrograph, 
                  UNITHYDRINFO *HydrographInfo, float *Hydrograph, 
                  FILES StreamFile);

int ScanInts(FILE *FilePtr, int *X, int N);

int ScanDoubles(FILE *FilePtr, double *X, int N);

int ScanFloats(FILE *FilePtr, float *X, int N);

uchar ScanUChars(FILE *FilePtr, uchar *X, int N);

void SkipHeader(FILES *InFile, int NLines);

void SkipLines(FILES *InFile, int NLines);

void StoreModelState(char *Path, TIMESTRUCT *Time, MAPSIZE *Map, 
		     OPTIONSTRUCT *Options, TOPOPIX **TopoMap, 
		     PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, 
		     VEGPIX **VegMap, LAYER Veg, SOILPIX **SoilMap, 
		     LAYER Soil, UNITHYDRINFO *HydrographInfo,
		     float *Hydrograph);
#endif


