/*
 * SUMMARY:      MainDHSVM.c - Distributed Hydrology-Soil-Vegetation Model
 * USAGE:        DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Main routine to drive DHSVM, the Distributed 
 *               Hydrology-Soil-Vegetation Model  
 * DESCRIP-END.cd
 * FUNCTIONS:    main()
 * COMMENTS:
 * $Id: MainDHSVM.c,v 1.42 2006/10/12 20:38:11 nathalie Exp $
 */

/******************************************************************************/
/*				    INCLUDES                                  */
/******************************************************************************/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "fileio.h"
#include "getinit.h"
#include "DHSVMChannel.h"
#include "channel.h"

/******************************************************************************/
/*				GLOBAL VARIABLES                              */
/******************************************************************************/

/* global function pointers */
void (*CreateMapFile) (char *FileName, ...);
int (*Read2DMatrix) (char *FileName, void *Matrix, int NumberType, int NY, int NX, int NDataSet, ...);
int (*Write2DMatrix) (char *FileName, void *Matrix, int NumberType, int NY, int NX, ...);

/* global strings */
char *version = "Version 3.1.1 ";               /* store version string */
char commandline[BUFSIZE + 1] = "";		/* store command line */
char fileext[BUFSIZ + 1] = "";			/* file extension */
char errorstr[BUFSIZ + 1] = "";			/* error message */
/******************************************************************************/
/*				      MAIN                                    */
/******************************************************************************/
int main(int argc, char **argv)
{
  float *Hydrograph = NULL;
  float ***MM5Input = NULL;
  float **PrecipLapseMap = NULL;
  float **PrismMap = NULL;
  unsigned char ***ShadowMap = NULL;
  float **SkyViewMap = NULL;
  float ***WindModel = NULL;
  int MaxStreamID, MaxRoadID;
  float SedDiams[NSEDSIZES];     /* Sediment particle diameters (mm) */
  float roadarea;
  time_t tloc;
  int flag;
  int i;
  int j;
  int x;						/* row counter */
  int y;						/* column counter */
  int shade_offset;				/* a fast way of handling arraay position given the number of mm5 input options */
  int NStats;					/* Number of meteorological stations */
  uchar ***MetWeights = NULL;	/* 3D array with weights for interpolating meteorological variables between the stations */

  int NGraphics;				/* number of graphics for X11 */
  int *which_graphics;			/* which graphics for X11 */
  char buffer[32];

  AGGREGATED Total = {			/* Total or average value of a  variable over the entire basin */

    {0.0, NULL, NULL, NULL, NULL, 0.0},												/* EVAPPIX */
    {0.0, 0.0, 0.0, 0.0, NULL, NULL, 0.0, 0, 0.0},									/* PRECIPPIX */
    {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 0.0, 0.0, 0.0},							/* PIXRAD */
    {0.0, 0.0},																		/* RADCLASSPIX */
    {0.0, 0.0, 0, NULL, NULL, 0.0, 0, 0.0, 0.0, 0.0, 0.0, NULL, 
	NULL, NULL, NULL, NULL, NULL, 0.0},												/* ROADSTRUCT*/
    {0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},		/* SNOWPIX */
    {0, 0.0, NULL, NULL, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},			/*SOILPIX */
    { 0.0, 0.0, 0.0, 0.0, 0.0},														/*SEDPIX */
    { 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},							/*FINEPIX */
    0.0, 0.0, 0.0, 0.0, 0.0, 0l, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  };
  CHANNEL ChannelData = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
  DUMPSTRUCT Dump;
  EVAPPIX **EvapMap = NULL;
  INPUTFILES InFiles;
  LAYER Soil;
  LAYER Veg;
  LISTPTR Input = NULL;			/* Linked list with input strings */
  MAPSIZE Map;					/* Size and location of model area */
  MAPSIZE Radar;				/* Size and location of area covered by precipitation radar */
  MAPSIZE MM5Map;				/* Size and location of area covered by MM5 input files */
  METLOCATION *Stat = NULL;
  OPTIONSTRUCT Options;			/* Structure with information which program options to follow */
  PIXMET LocalMet;				/* Meteorological conditions for current pixel */
  FINEPIX ***FineMap	= NULL;
  PRECIPPIX **PrecipMap = NULL;
  RADARPIX **RadarMap	= NULL;
  RADCLASSPIX **RadMap	= NULL;
  ROADSTRUCT **Network	= NULL;	/* 2D Array with channel information for each pixel */
  SNOWPIX **SnowMap		= NULL;
  MET_MAP_PIX **MetMap	= NULL;
  SNOWTABLE *SnowAlbedo = NULL;
  SOILPIX **SoilMap		= NULL;
  SEDPIX **SedMap		= NULL;
  SOILTABLE *SType	    = NULL;
  SEDTABLE *SedType		= NULL;
  SOLARGEOMETRY SolarGeo;		/* Geometry of Sun-Earth system (needed for INLINE radiation calculations */
  TIMESTRUCT Time;
  TOPOPIX **TopoMap = NULL;
  UNITHYDR **UnitHydrograph = NULL;
  UNITHYDRINFO HydrographInfo;	/* Information about unit hydrograph */
  VEGPIX **VegMap = NULL;
  VEGTABLE *VType = NULL;
  WATERBALANCE Mass =			/* parameter for mass balance calculations */
    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

/*****************************************************************************
  Initialization Procedures 
*****************************************************************************/
  if (argc != 2) {
    fprintf(stderr, "\nUsage: %s inputfile\n\n", argv[0]);
    fprintf(stderr, "DHSVM uses two output streams: \n");
    fprintf(stderr, "Standard Out, for the majority of output \n");
    fprintf(stderr, "Standard Error, for the final mass balance \n");
    fprintf(stderr, "\nTo pipe output correctly to files: \n");
    fprintf(stderr, "(cmd > f1) >& f2 \n");
    fprintf(stderr, "where f1 is stdout_file and f2 is stderror_file\n");
    exit(EXIT_FAILURE);
  }

  sprintf(commandline, "%s %s", argv[0], argv[1]);
  printf("%s \n", commandline);
  fprintf(stderr, "%s \n", commandline);
  strcpy(InFiles.Const, argv[1]);

  printf("\nRunning DHSVM %s\n", version);
  printf("\nSTARTING INITIALIZATION PROCEDURES\n\n");

  ReadInitFile(InFiles.Const, &Input);

  InitConstants(Input, &Options, &Map, &SolarGeo, &Time);

  InitFileIO(Options.FileFormat);
  InitTables(Time.NDaySteps, Input, &Options, &SType, &Soil, &VType, &Veg,
	     &SnowAlbedo);

  InitTerrainMaps(Input, &Options, &Map, &Soil, &TopoMap, &SoilMap, &VegMap);

  CheckOut(Options.CanopyRadAtt, Veg, Soil, VType, SType, &Map, TopoMap, 
	   VegMap, SoilMap);

  if (Options.HasNetwork)
    InitChannel(Input, &Map, Time.Dt, &ChannelData, SoilMap, &MaxStreamID, &MaxRoadID, &Options);
  else if (Options.Extent != POINT)
    InitUnitHydrograph(Input, &Map, TopoMap, &UnitHydrograph,
		       &Hydrograph, &HydrographInfo);
 
  InitNetwork(Map.NY, Map.NX, Map.DX, Map.DY, TopoMap, SoilMap, 
	      VegMap, VType, &Network, &ChannelData, Veg, &Options);

  InitMetSources(Input, &Options, &Map, Soil.MaxLayers, &Time,
		 &InFiles, &NStats, &Stat, &Radar, &MM5Map);

  /* the following piece of code is for the UW PRISM project */
  /* for real-time verification of SWE at Snotel sites */
  /* Other users, set OPTION.SNOTEL to FALSE, or use TRUE with caution */

  if (Options.Snotel == TRUE && Options.Outside == FALSE) {
    printf
      ("Warning: All met stations locations are being set to the vegetation class GLACIER\n");
    printf
      ("Warning: This requires that you have such a vegetation class in your vegetation table\n");
    printf("To disable this feature set Snotel OPTION to FALSE\n");
    for (i = 0; i < NStats; i++) {
      printf("veg type for station %d is %d ", i,
	     VegMap[Stat[i].Loc.N][Stat[i].Loc.E].Veg);
      for (j = 0; j < Veg.NTypes; j++) {
	if (VType[j].Index == GLACIER) {
	  VegMap[Stat[i].Loc.N][Stat[i].Loc.E].Veg = j;
	  break;
	}
      }
      if (j == Veg.NTypes) {	/* glacier class not found */
	ReportError("MainDHSVM", 62);
      }
      printf("setting to glacier type (assumed bare class): %d\n", j);
    }
  }

  InitMetMaps(Time.NDaySteps, &Map, &Radar, &Options, InFiles.WindMapPath,
	      InFiles.PrecipLapseFile, &PrecipLapseMap, &PrismMap,
	      &ShadowMap, &SkyViewMap, &EvapMap, &PrecipMap,
	      &RadarMap, &RadMap, SoilMap, &Soil, VegMap, &Veg, TopoMap,
	      &MM5Input, &WindModel);

  InitInterpolationWeights(&Map, &Options, TopoMap, &MetWeights, Stat, NStats);

  InitDump(Input, &Options, &Map, Soil.MaxLayers, Veg.MaxLayers, Time.Dt,
	   TopoMap, &Dump, &NGraphics, &which_graphics);

  if (Options.HasNetwork == TRUE) {
    InitChannelDump(&ChannelData, Dump.Path);
    ReadChannelState(Dump.InitStatePath, &(Time.Start), ChannelData.streams);
  }

  InitSnowMap(&Map, &SnowMap);
  InitAggregated(Veg.MaxLayers, Soil.MaxLayers, &Total);

  InitModelState(&(Time.Start), &Map, &Options, PrecipMap, SnowMap, SoilMap,
		 Soil, SType, VegMap, Veg, VType, Dump.InitStatePath,
		 SnowAlbedo, TopoMap, Network, &HydrographInfo, Hydrograph);

  InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, ShadowMap,
	       RadMap, &InFiles, Veg.NTypes, VType, NStats, Stat, 
	       Dump.InitStatePath);

  InitNewDay(Time.Current.JDay, &SolarGeo);

  if (NGraphics > 0) {
    printf("Initialzing X11 display and graphics \n");
    InitXGraphics(argc, argv, Map.NY, Map.NX, NGraphics, &MetMap);
  }

  shade_offset = FALSE;
  if (Options.Shading == TRUE)
    shade_offset = TRUE;

  /* Done with initialization, delete the list with input strings */
  DeleteList(Input);

  /*****************************************************************************
  Sediment Initialization Procedures 
  *****************************************************************************/
  if(Options.Sediment) {
     time (&tloc);
     srand (tloc);
  /* Randomize Random Generator */
 
  /* Commenting the line above and uncommenting the line below 
     allows for the comparison of scenarios. */
  /*  srand48 (0);  */
 

    printf("\nSTARTING SEDIMENT INITIALIZATION PROCEDURES\n\n");

    ReadInitFile(Options.SedFile, &Input);

    InitParameters(Input, &Options, &Map, &Network, &ChannelData, TopoMap,
		   &Time, SedDiams);

    InitSedimentTables(Time.NDaySteps, Input, &SedType, &SType, &VType, &Soil, &Veg);

    InitFineMaps(Input, &Options, &Map, &Soil, &TopoMap, &SoilMap, 
		  &FineMap);

    if (Options.HasNetwork){ 
      printf("Initializing channel sediment\n\n");
      InitChannelSedimentDump(&ChannelData, Dump.Path, Options.ChannelRouting); 
      InitChannelSediment(ChannelData.streams, &Total);
      InitChannelSediment(ChannelData.roads, &Total);
    }

    InitSedMap( &Map, &SedMap);

    /* Done with initialization, delete the list with input strings */
    DeleteList(Input);
  }

  /* setup for mass balance calculations */
  Aggregate(&Map, &Options, TopoMap, &Soil, &Veg, VegMap, EvapMap, PrecipMap,
	      RadMap, SnowMap, SoilMap, &Total, VType, Network, SedMap, FineMap,
	      &ChannelData, &roadarea);

  Mass.StartWaterStorage =
    Total.Soil.IExcess + Total.CanopyWater + Total.SoilWater + Total.Snow.Swq +
    Total.Soil.SatFlow;
  Mass.OldWaterStorage = Mass.StartWaterStorage;

  if (Options.Sediment) {
    Mass.StartChannelSedimentStorage = Total.ChannelSedimentStorage;
    Mass.LastChannelSedimentStorage = Mass.StartChannelSedimentStorage;
  }

/*****************************************************************************
  Perform Calculations 
*****************************************************************************/

  while (Before(&(Time.Current), &(Time.End)) ||
	 IsEqualTime(&(Time.Current), &(Time.End))) {
    ResetAggregate(&Soil, &Veg, &Total, &Options);

    if (IsNewMonth(&(Time.Current), Time.Dt))
      InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, ShadowMap,
		   RadMap, &InFiles, Veg.NTypes, VType, NStats, Stat, 
		   Dump.InitStatePath);

    if (IsNewDay(Time.DayStep)) {
      InitNewDay(Time.Current.JDay, &SolarGeo);
      PrintDate(&(Time.Current), stdout);
      printf("\n");
    }

/*     PrintDate(&(Time.Current),stdout); */
/*     printf("\n"); */
/*     uncomment the above lines to print the time at every step */

    /* determine surface erosion and routing scheme */
    SedimentFlag(&Options, &Time); 

    InitNewStep(&InFiles, &Map, &Time, Soil.MaxLayers, &Options, NStats, Stat,
		InFiles.RadarFile, &Radar, RadarMap, &SolarGeo, TopoMap, RadMap,
                SoilMap, MM5Input, WindModel, &MM5Map);

    /* initialize channel/road networks for time step */

    if (Options.HasNetwork) {
      channel_step_initialize_network(ChannelData.streams);
      channel_step_initialize_network(ChannelData.roads);
    }

    for (y = 0; y < Map.NY; y++) {
      for (x = 0; x < Map.NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {

	  if (Options.Shading)
	    LocalMet =
	      MakeLocalMetData(y, x, &Map, Time.DayStep, &Options, NStats,
			       Stat, MetWeights[y][x], TopoMap[y][x].Dem,
			       &(RadMap[y][x]), &(PrecipMap[y][x]), &Radar,
			       RadarMap, PrismMap, &(SnowMap[y][x]),
			       SnowAlbedo, MM5Input, WindModel, PrecipLapseMap,
			       &MetMap, NGraphics, Time.Current.Month,
			       SkyViewMap[y][x], ShadowMap[Time.DayStep][y][x],
			       SolarGeo.SunMax, SolarGeo.SineSolarAltitude);
	  else
	    LocalMet =
	      MakeLocalMetData(y, x, &Map, Time.DayStep, &Options, NStats,
			       Stat, MetWeights[y][x], TopoMap[y][x].Dem,
			       &(RadMap[y][x]), &(PrecipMap[y][x]), &Radar,
			       RadarMap, PrismMap, &(SnowMap[y][x]),
			       SnowAlbedo, MM5Input, WindModel, PrecipLapseMap,
			       &MetMap, NGraphics, Time.Current.Month, 0.0,
			       0.0, SolarGeo.SunMax,
			       SolarGeo.SineSolarAltitude);

	  for (i = 0; i < Soil.MaxLayers; i++) {
	    if (Options.HeatFlux == TRUE) {
	      if (Options.MM5 == TRUE)
		SoilMap[y][x].Temp[i] =
		  MM5Input[shade_offset + i + N_MM5_MAPS][y][x];
	      else
		SoilMap[y][x].Temp[i] = Stat[0].Data.Tsoil[i];
	    }
	    else
	      SoilMap[y][x].Temp[i] = LocalMet.Tair;
	  }
	  
	  MassEnergyBalance(y, x, SolarGeo.SineSolarAltitude, Map.DX, Map.DY, 
			    Time.Dt, Options.HeatFlux, Options.CanopyRadAtt, 
			    Options.RoadRouting, Options.Infiltration,
			    Veg.MaxLayers, &LocalMet, 
			    &(Network[y][x]), &(PrecipMap[y][x]), 
			    &(VType[VegMap[y][x].Veg-1]), &(VegMap[y][x]),
			    &(SType[SoilMap[y][x].Soil-1]), &(SoilMap[y][x]), 
			    &(SnowMap[y][x]), &(EvapMap[y][x]), &(Total.Rad),
			    &ChannelData);
	}
      }
    }

#ifndef SNOW_ONLY

    /* set sediment inflows to zero - they are incremented elsewhere */
  if ((Options.HasNetwork) && (Options.Sediment)){ 
    InitChannelSedInflow(ChannelData.streams);
    InitChannelSedInflow(ChannelData.roads);
  }
    
    RouteSubSurface(Time.Dt, &Map, TopoMap, VType, VegMap, Network,
		    SType, SoilMap, &ChannelData, &Time, &Options, Dump.Path,
		    SedMap, FineMap, SedType, MaxStreamID, SnowMap);

    if (Options.HasNetwork)
      RouteChannel(&ChannelData, &Time, &Map, TopoMap, SoilMap, &Total, 
		   &Options, Network, SType, PrecipMap, SedMap,
		   LocalMet.Tair, LocalMet.Rh, SedDiams);

    /* Sediment Routing in Channel and output to sediment files */
    if ((Options.HasNetwork) && (Options.Sediment)){
      SPrintDate(&(Time.Current), buffer);
      flag = IsEqualTime(&(Time.Current), &(Time.Start));
      if(Options.ChannelRouting){
	if (ChannelData.roads != NULL) {
	  RouteChannelSediment(ChannelData.roads, Time, &Dump, &Total, SedDiams);
	  channel_save_sed_outflow_text(buffer, ChannelData.roads,
					ChannelData.sedroadout,
					ChannelData.sedroadflowout, flag);
	  RouteCulvertSediment(&ChannelData, &Map, TopoMap, SedMap, 
			       &Total, SedDiams);
	}
	RouteChannelSediment(ChannelData.streams, Time, &Dump, &Total, SedDiams);
	channel_save_sed_outflow_text(buffer, ChannelData.streams,
				      ChannelData.sedstreamout,
				      ChannelData.sedstreamflowout, flag);
      }
      else{
	if (ChannelData.roads != NULL) {
	  channel_save_sed_inflow_text(buffer, ChannelData.roads,
				     ChannelData.sedroadinflow, SedDiams,
				     flag);
	}
	channel_save_sed_inflow_text(buffer, ChannelData.streams,
				     ChannelData.sedstreaminflow, SedDiams,
				     flag);
      }
      SaveChannelSedInflow(ChannelData.roads, &Total);
      SaveChannelSedInflow(ChannelData.streams, &Total);
    }
    
    if (Options.Extent == BASIN)
      RouteSurface(&Map, &Time, TopoMap, SoilMap, &Options,
		   UnitHydrograph, &HydrographInfo, Hydrograph,
		   &Dump, VegMap, VType, SType, &ChannelData, SedMap,
		   PrecipMap, SedType, LocalMet.Tair, LocalMet.Rh, SedDiams);

#endif

    if (NGraphics > 0)
      draw(&(Time.Current), IsEqualTime(&(Time.Current), &(Time.Start)),
	   Time.DayStep, &Map, NGraphics, which_graphics, VType,
	   SType, SnowMap, SoilMap, SedMap, FineMap, VegMap, TopoMap, PrecipMap,
	   PrismMap, SkyViewMap, ShadowMap, EvapMap, RadMap, MetMap, Network,
	   &Options);
    
    Aggregate(&Map, &Options, TopoMap, &Soil, &Veg, VegMap, EvapMap, PrecipMap,
	      RadMap, SnowMap, SoilMap, &Total, VType, Network, SedMap, FineMap,
	      &ChannelData, &roadarea);
    
    MassBalance(&(Time.Current), &(Dump.Balance), &(Dump.SedBalance), &Total, 
		&Mass, &Options);
    
    ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
	     EvapMap, PrecipMap, RadMap, SnowMap, MetMap, VegMap, &Veg, SoilMap,
	     SedMap, Network, &ChannelData, FineMap, &Soil, &Total, &HydrographInfo,
	     Hydrograph);
    
    IncreaseTime(&Time);

  }

  ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
	   EvapMap, PrecipMap, RadMap, SnowMap, MetMap, VegMap, &Veg, SoilMap,
	   SedMap, Network, &ChannelData, FineMap, &Soil, &Total, &HydrographInfo,
	   Hydrograph);

  FinalMassBalance(&(Dump.FinalBalance), &Total, &Mass, &Options, roadarea);

/*****************************************************************************
  Cleanup
*****************************************************************************/

  printf("\nSTARTING CLEANUP\n\n");

  printf("\nEND OF MODEL RUN\n\n");

  return EXIT_SUCCESS;

}
