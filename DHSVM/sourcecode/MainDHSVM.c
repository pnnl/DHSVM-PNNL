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
 * $Id: MainDHSVM.c,v 3.2 2018/2/22 16:58 ning Exp $
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

/* global strings */
char *version = "Version 3.2";        /* store version string */
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
  float **PptMultiplierMap = NULL;                                  
  int MaxStreamID, MaxRoadID;
  clock_t start, finish1;
  double runtime = 0.0;
  int t = 0;
  float roadarea;
  int i;
  int j;
  int x;						/* row counter */
  int y;						/* column counter */
  int shade_offset;				/* a fast way of handling arraay position given the number of mm5 input options */
  int NStats;					/* Number of meteorological stations */
  uchar ***MetWeights = NULL;	/* 3D array with weights for interpolating meteorological variables between the stations */

  int NGraphics;				/* number of graphics for X11 */
  int *which_graphics;			/* which graphics for X11 */

  AGGREGATED Total = {			/* Total or average value of a  variable over the entire basin */
    {0.0, NULL, NULL, NULL, NULL, 0.0},												/* EVAPPIX */
    {0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, 0.0, 0, 0.0},								/* PRECIPPIX */
    {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 0.0, {0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                                                                                /* PIXRAD */
    {0.0, 0.0, 0, NULL, NULL, 0.0, 0, 0.0, 0.0, 0.0, 0.0, NULL, NULL},				  /* ROADSTRUCT*/
	  {0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},     /* SNOWPIX */ 
    {0, 0.0, NULL, NULL, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL},			                /* SOILPIX */
    {0, 0, 0.0, 0.0, 0.0, NULL, NULL, NULL, NULL, 0.0},                             /* VEGPIX */
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0l, 0.0, 0.0
  };
  CHANNEL ChannelData = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  DUMPSTRUCT Dump;
  EVAPPIX **EvapMap = NULL;
  INPUTFILES InFiles;
  LAYER Soil;
  LAYER Veg;
  LISTPTR Input = NULL;			/* Linked list with input strings */
  MAPSIZE Map;					/* Size and location of model area */
  MAPSIZE Radar;				/* Size and location of area covered by precipitation radar */
  MAPSIZE MM5Map;				/* Size and location of area covered by MM5 input files */
  GRID Grid;
  METLOCATION *Stat = NULL;
  OPTIONSTRUCT Options;			/* Structure with information which program options to follow */
  PIXMET LocalMet;				/* Meteorological conditions for current pixel */
  PRECIPPIX **PrecipMap = NULL;
  RADARPIX **RadarMap	= NULL;
  PIXRAD **RadiationMap = NULL;
  ROADSTRUCT **Network	= NULL;	/* 2D Array with channel information for each pixel */
  SNOWPIX **SnowMap		= NULL;
  MET_MAP_PIX **MetMap	= NULL;
  SOILPIX **SoilMap		= NULL;
  SOILTABLE *SType	    = NULL;
  SOLARGEOMETRY SolarGeo;		/* Geometry of Sun-Earth system (needed for INLINE radiation calculations */
  TIMESTRUCT Time;
  TOPOPIX **TopoMap = NULL;
  UNITHYDR **UnitHydrograph = NULL;
  UNITHYDRINFO HydrographInfo;	/* Information about unit hydrograph */
  VEGPIX **VegMap = NULL;
  VEGTABLE *VType = NULL;
  WATERBALANCE Mass =			/* parameter for mass balance calculations */
    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };


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
#ifdef SNOW_ONLY
  printf("----------------------------------\n");
  printf("WARNING: USING SNOW ONLY MODULES (prescribed in makefile)!\n");
  printf("----------------------------------\n");
#endif
  printf("\nSTARTING INITIALIZATION PROCEDURES\n\n");

  /* Start recording time */
  start = clock();

  /* initiate input/output format */

  ReadInitFile(InFiles.Const, &Input);
  InitConstants(Input, &Options, &Map, &SolarGeo, &Time);

  InitFileIO(Options.FileFormat);
  InitTables(Time.NDaySteps, Input, &Options, &Map, &SType, &Soil, &VType, &Veg); 

  InitTerrainMaps(Input, &Options, &Map, &Soil, &Veg, &TopoMap, SType, &SoilMap, VType, &VegMap);

  InitSnowMap(&Map, &SnowMap, &Time);

  InitMappedConstants(Input, &Options, &Map, &SnowMap);

  CheckOut(&Options, Veg, Soil, VType, SType, &Map, TopoMap, VegMap, SoilMap);

#ifdef TOPO_DUMP
  DumpTopo(&Map, TopoMap);
#endif
  
  if (Options.HasNetwork)
    InitChannel(Input, &Map, Time.Dt, &ChannelData, SoilMap, &MaxStreamID, &MaxRoadID, &Options);
  else if (Options.Extent != POINT)
    InitUnitHydrograph(Input, &Map, TopoMap, &UnitHydrograph,
		       &Hydrograph, &HydrographInfo);
 
  InitNetwork(Map.NY, Map.NX, Map.DX, Map.DY, TopoMap, SoilMap, 
	      VegMap, VType, &Network, &ChannelData, Veg, &Options);

  InitMetSources(Input, &Options, &Map, TopoMap, Soil.MaxLayers, &Time,
		 &InFiles, &NStats, &Stat, &Radar, &MM5Map, &Grid);

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

  InitMetMaps(Input, Time.NDaySteps, &Map, &Radar, &Options, InFiles.WindMapPath,
	      InFiles.PrecipLapseFile, &PrecipLapseMap, &PrismMap,
	      &ShadowMap, &SkyViewMap, &EvapMap, &PrecipMap, &PptMultiplierMap,
	      &RadarMap, &RadiationMap, SoilMap, &Soil, VegMap, &Veg, TopoMap,
	      &MM5Input, &WindModel);

  InitInterpolationWeights(&Map, &Options, TopoMap, &MetWeights, Stat, NStats);

  InitDump(Input, &Options, &Map, Soil.MaxLayers, Veg.MaxLayers, Time.Dt,
	   TopoMap, &Dump, &NGraphics, &which_graphics);

#ifndef SNOW_ONLY
  if (Options.HasNetwork == TRUE) {
    InitChannelDump(&Options, &ChannelData, Dump.Path);
    ReadChannelState(Dump.InitStatePath, &(Time.Start), ChannelData.streams);
	if (Options.StreamTemp && Options.CanopyShading)
	  InitChannelRVeg(&Time, ChannelData.streams);
  }
#endif

  InitAggregated(&Options, Veg.MaxLayers, Soil.MaxLayers, &Total);

  InitModelState(&(Time.Start), Time.NDaySteps, &Map, &Options, PrecipMap, SnowMap, SoilMap,
		 Soil, SType, VegMap, Veg, VType, Dump.InitStatePath,
		 TopoMap, Network, &HydrographInfo, Hydrograph);

  InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, ShadowMap,
	       &InFiles, Veg.NTypes, VType, NStats, Stat, Dump.InitStatePath, &VegMap);

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

  /* setup for mass balance calculations */
  Aggregate(&Map, &Options, TopoMap, &Soil, &Veg, VegMap, EvapMap, PrecipMap,
	      RadiationMap, SnowMap, SoilMap, &Total, VType, Network, &ChannelData, &roadarea, Time.Dt);

  Mass.StartWaterStorage =
    Total.Soil.IExcess + Total.CanopyWater + Total.SoilWater + Total.Snow.Swq +
    Total.Soil.SatFlow;
  Mass.OldWaterStorage = Mass.StartWaterStorage;

  /* computes the number of grid cell contributing to one segment */
  if (Options.StreamTemp) 
	Init_segment_ncell(TopoMap, ChannelData.stream_map, Map.NY, Map.NX, ChannelData.streams);

/*****************************************************************************
  Perform Calculations 
*****************************************************************************/
  while (Before(&(Time.Current), &(Time.End)) ||
	 IsEqualTime(&(Time.Current), &(Time.End))) {

    /* reset aggregated variables */
    ResetAggregate(&Soil, &Veg, &Total, &Options);
    
    /* redistribute snow based on snow surface slope etc */
    if (Options.SnowSlide)
	    Avalanche(&Map, TopoMap, &Time, &Options, SnowMap);
    
    if (IsNewWaterYear(&(Time.Current)))
      InitNewWaterYear(&Time, &Options, &Map, TopoMap, SnowMap);

    if (IsNewMonth(&(Time.Current), Time.Dt))
      InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, ShadowMap,
		   &InFiles, Veg.NTypes, VType, NStats, Stat, Dump.InitStatePath, &VegMap);

    if (IsNewDay(Time.DayStep)) {
      InitNewDay(Time.Current.JDay, &SolarGeo);
      PrintDate(&(Time.Current), stdout);
      printf("\n");
    }

    InitNewStep(&InFiles, &Map, &Time, Soil.MaxLayers, &Options, NStats, Stat,
		InFiles.RadarFile, &Radar, RadarMap, &SolarGeo, TopoMap, 
                SoilMap, MM5Input, PrecipLapseMap, WindModel, &MM5Map);

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
	        MakeLocalMetData(y, x, &Map, Time.DayStep, Time.NDaySteps, &Options, NStats,
			       Stat, MetWeights[y][x], TopoMap[y][x].Dem,
			       &(RadiationMap[y][x]), &(PrecipMap[y][x]), &Radar,
			       RadarMap, PrismMap, &(SnowMap[y][x]),
			       &(VegMap[y][x].Type), &(VegMap[y][x]), 
             MM5Input, WindModel, PrecipLapseMap,
			       &MetMap, PptMultiplierMap[y][x], NGraphics, Time.Current.Month,
			       SkyViewMap[y][x], ShadowMap[Time.DayStep][y][x],
			       SolarGeo.SunMax, SolarGeo.SineSolarAltitude);
		  else
	        LocalMet =
	        MakeLocalMetData(y, x, &Map, Time.DayStep, Time.NDaySteps, &Options, NStats,
			       Stat, MetWeights[y][x], TopoMap[y][x].Dem,
			       &(RadiationMap[y][x]), &(PrecipMap[y][x]), &Radar,
			       RadarMap, PrismMap, &(SnowMap[y][x]),
			       &(VegMap[y][x].Type), &(VegMap[y][x]), 
             MM5Input, WindModel, PrecipLapseMap,
			       &MetMap, PptMultiplierMap[y][x],NGraphics, Time.Current.Month, 0.0,
			       0.0, SolarGeo.SunMax,
			       SolarGeo.SineSolarAltitude);

		  /* get surface tempeature of each soil layer */
		  for (i = 0; i < Soil.MaxLayers; i++) {
	        if (Options.HeatFlux == TRUE) {
	          if (Options.MM5 == TRUE)
		        SoilMap[y][x].Temp[i] =
				MM5Input[shade_offset + i + N_MM5_MAPS][y][x];

              /* read tempeature of each soil layer from met station input */
			  else
		        SoilMap[y][x].Temp[i] = Stat[0].Data.Tsoil[i];
			}
            /* if heat flux option is turned off, soil temperature of all 3 layers 
            is taken equal to air tempeature */
	        else
	          SoilMap[y][x].Temp[i] = LocalMet.Tair;
		  }
		  
          MassEnergyBalance(&Options, y, x, SolarGeo.SineSolarAltitude, Map.DX, Map.DY,
            Time.Dt, Options.HeatFlux, Options.CanopyRadAtt, Options.Infiltration, Soil.MaxLayers,
            Veg.MaxLayers, &LocalMet, &(Network[y][x]), &(PrecipMap[y][x]),
            &(VType[VegMap[y][x].Veg - 1]), &(VegMap[y][x]), &(SType[SoilMap[y][x].Soil - 1]),
            &(SoilMap[y][x]), &(SnowMap[y][x]), &(RadiationMap[y][x]), &(EvapMap[y][x]),
            &(Total.Rad), &ChannelData, SkyViewMap);
	 
		  PrecipMap[y][x].SumPrecip += PrecipMap[y][x].Precip;
		}
	  }
    }

	/* Average all RBM inputs over each segment */
	if (Options.StreamTemp) {
	  channel_grid_avg(ChannelData.streams);
      if (Options.CanopyShading)
	    CalcCanopyShading(&Time, ChannelData.streams, &SolarGeo);
	}

 #ifndef SNOW_ONLY
    
    RouteSubSurface(Time.Dt, &Map, TopoMap, VType, VegMap, Network,
		    SType, SoilMap, &ChannelData, &Time, &Options, Dump.Path,
		    MaxStreamID, SnowMap);

    if (Options.HasNetwork)
      RouteChannel(&ChannelData, &Time, &Map, TopoMap, SoilMap, &Total, 
		   &Options, Network, SType, PrecipMap, LocalMet.Tair, LocalMet.Rh, SnowMap);

    if (Options.Extent == BASIN)
      RouteSurface(&Map, &Time, TopoMap, SoilMap, &Options,
        UnitHydrograph, &HydrographInfo, Hydrograph,
        &Dump, VegMap, VType, &ChannelData);


#endif

    if (NGraphics > 0)
      draw(&(Time.Current), IsEqualTime(&(Time.Current), &(Time.Start)),
	   Time.DayStep, &Map, NGraphics, which_graphics, VType,
	   SType, SnowMap, SoilMap, VegMap, TopoMap, PrecipMap,
	   PrismMap, SkyViewMap, ShadowMap, EvapMap, RadiationMap, 
	   MetMap, Network, &Options);
    
    Aggregate(&Map, &Options, TopoMap, &Soil, &Veg, VegMap, EvapMap, PrecipMap,
	      RadiationMap, SnowMap, SoilMap, &Total, VType, Network, &ChannelData, &roadarea, Time.Dt);
    
    if (Options.SnowStats)
      SnowStats(&(Time.Current), &Map, &Options, TopoMap, SnowMap, Time.Dt);
    
    MassBalance(&(Time.Current), &(Time.Start), &(Dump.Balance), &Total, &Mass);

    ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
	     EvapMap, RadiationMap, PrecipMap, SnowMap, MetMap, VegMap, &Veg, 
		 SoilMap, Network, &ChannelData, &Soil, &Total, &HydrographInfo,Hydrograph);
	
    IncreaseTime(&Time);
	t += 1;
  }

  ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
	   EvapMap, RadiationMap, PrecipMap, SnowMap, MetMap, VegMap, &Veg, SoilMap,
	   Network, &ChannelData, &Soil, &Total, &HydrographInfo, Hydrograph);

#ifndef SNOW_ONLY
  FinalMassBalance(&(Dump.FinalBalance), &Total, &Mass);
#endif

  printf("\nEND OF MODEL RUN\n\n");

  /* record the run time at the end of each time loop */
  finish1 = clock ();
  runtime = (finish1-start)/CLOCKS_PER_SEC;
  printf("***********************************************************************************");
  printf("\nRuntime Summary:\n");
  printf("%6.2f hours elapsed for the simulation period of %d hours (%.1f days) \n", 
	  runtime/3600, t*Time.Dt/3600, (float)t*Time.Dt/3600/24);

  return EXIT_SUCCESS;
}
/*****************************************************************************
  Cleanup
*****************************************************************************/
void cleanup(DUMPSTRUCT *Dump, CHANNEL *ChannelData, OPTIONSTRUCT *Options)
{
	if (Dump->Aggregate.FilePtr != NULL) 
	  fclose(Dump->Aggregate.FilePtr);
	if (Dump->Balance.FilePtr != NULL) 
	  fclose(Dump->Balance.FilePtr);
	if (Dump->FinalBalance.FilePtr != NULL) 
	  fclose(Dump->FinalBalance.FilePtr);
	if (ChannelData->streamflowout != NULL)
	  fclose(ChannelData->streamflowout);
	if (ChannelData->streamout != NULL)
	  fclose(ChannelData->streamout);
	if (ChannelData->roadflowout != NULL)
	  fclose(ChannelData->roadflowout );
	if (ChannelData->roadout != NULL)
	  fclose(ChannelData->roadout);

	if (Options->StreamTemp) {
	  if (ChannelData->streaminflow != NULL) 
      fclose(ChannelData->streaminflow);
	  if (ChannelData->streamoutflow != NULL) 
      fclose(ChannelData->streamoutflow);
    if (ChannelData->streamMelt != NULL)
      fclose(ChannelData->streamMelt);                                    
	  if (ChannelData->streamNSW != NULL) 
      fclose(ChannelData->streamNSW);
	  if (ChannelData->streamNLW!= NULL) 
      fclose(ChannelData->streamNLW);								  
	  if (ChannelData->streamVP!= NULL) 
      fclose(ChannelData->streamVP);	
	  if (ChannelData->streamWND!= NULL) 
      fclose(ChannelData->streamWND);	
	  if (ChannelData->streamATP!= NULL) 
      fclose(ChannelData->streamATP);
	}
}
