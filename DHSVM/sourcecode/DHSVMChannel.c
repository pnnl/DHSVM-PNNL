/* -------------------------------------------------------------
   file: DHSVMChannel.c
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created August 30, 1996 by  William A Perkins
   $Id: DHSVMChannel.c, v3.1.2  2013/12/20   Ning Exp $
   ------------------------------------------------------------- */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "getinit.h"
#include "DHSVMChannel.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "settings.h"
#include "errorhandler.h"
#include "fileio.h"
#include "ParallelDHSVM.h"
#include "ParallelChannel.h"

/* -----------------------------------------------------------------------------
   InitChannel
   Reads stream and road files and builds the networks.
   -------------------------------------------------------------------------- */
void
InitChannel(LISTPTR Input, MAPSIZE *Map, int deltat, CHANNEL *channel,
	    SOILPIX ** SoilMap, int *MaxStreamID, int *MaxRoadID, OPTIONSTRUCT *Options)
{
  int i;
  STRINIENTRY StrEnv[] = {
    {"ROUTING", "STREAM NETWORK FILE", "", ""},
    {"ROUTING", "STREAM MAP FILE", "", ""},
    {"ROUTING", "STREAM CLASS FILE", "", ""},
    {"ROUTING", "RIPARIAN VEG FILE", "", ""},
    {"ROUTING", "ROAD NETWORK FILE", "", "none"},
    {"ROUTING", "ROAD MAP FILE", "", "none"},
    {"ROUTING", "ROAD CLASS FILE", "", "none"},
    {NULL, NULL, "", NULL}
  };

  printf("\nInitializing Road/Stream Networks\n");

  /* Read the key-entry pairs from the ROUTING section in the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  channel->stream_class = NULL;
  channel->road_class = NULL;
  channel->streams = NULL;
  channel->roads = NULL;
  channel->stream_map = NULL;
  channel->road_map = NULL;

  channel_init();
  channel_grid_init(Map->NX, Map->NY);

  if (strncmp(StrEnv[stream_class].VarStr, "none", 4)) {

    printf("\tReading Stream data\n");

    if ((channel->stream_class =
	 channel_read_classes(StrEnv[stream_class].VarStr, stream_class)) == NULL) {
      ReportError(StrEnv[stream_class].VarStr, 5);
    }
    if ((channel->streams =
	 channel_read_network(StrEnv[stream_network].VarStr,
			      channel->stream_class, MaxStreamID)) == NULL) {
      ReportError(StrEnv[stream_network].VarStr, 5);
    }
    if ((channel->stream_map =
	 channel_grid_read_map(Map, channel->streams,
			       StrEnv[stream_map].VarStr, SoilMap)) == NULL) {
      ReportError(StrEnv[stream_map].VarStr, 5);
    }
    error_handler(ERRHDL_STATUS,
		  "InitChannel: computing stream network routing coefficients");
    channel_routing_parameters(channel->streams, (double) deltat);
  }

  if (Options->StreamTemp) {
    if (strncmp(StrEnv[riparian_veg].VarStr, "none", 4)) {
      printf("\tReading channel riparian vegetation params\n");
      channel_read_rveg_param(channel->streams, StrEnv[riparian_veg].VarStr, MaxStreamID);
    }
  }

  if (strncmp(StrEnv[road_class].VarStr, "none", 4)) {

    printf("\tReading Road data\n");

    if ((channel->road_class =
	 channel_read_classes(StrEnv[road_class].VarStr, road_class) == NULL)) {
      ReportError(StrEnv[road_class].VarStr, 5);
    }
    if ((channel->roads =
	 channel_read_network(StrEnv[road_network].VarStr,
			      channel->road_class, MaxRoadID)) == NULL) {
      ReportError(StrEnv[road_network].VarStr, 5);
    }
    if ((channel->road_map =
	 channel_grid_read_map(Map, channel->roads,
			       StrEnv[road_map].VarStr, SoilMap)) == NULL) {
      ReportError(StrEnv[road_map].VarStr, 5);
    }
    error_handler(ERRHDL_STATUS,
		  "InitChannel: computing road network routing coefficients");
    channel_routing_parameters(channel->roads, (double) deltat);
  }

  ParallelBarrier();

  if (channel->streams != NULL) {
    channel->stream_state_ga = ChannelStateGA(channel->streams);
  }
  if (channel->roads != NULL) {
    channel->road_state_ga = ChannelStateGA(channel->roads);
  }

  ParallelBarrier();
}

/* -------------------------------------------------------------
   InitChannelDump
   ------------------------------------------------------------- */
void InitChannelDump(OPTIONSTRUCT *Options, CHANNEL * channel, 
                     char *DumpPath)
{
  char buffer[NAMESIZE];

  if (ParallelRank() == 0) {

    if (channel->streams != NULL) {
      sprintf(buffer, "%sStream.Flow", DumpPath);
      OpenFile(&(channel->streamout), buffer, "w", TRUE);
      sprintf(buffer, "%sStreamflow.Only", DumpPath);
      OpenFile(&(channel->streamflowout), buffer, "w", TRUE);
      /* output files for John's RBM model */
      if (Options->StreamTemp) {
        //inflow to segment
        sprintf(buffer, "%sInflow.Only", DumpPath);
        OpenFile(&(channel->streaminflow), buffer, "w", TRUE);
        // outflow ( redundant but it's a check
        sprintf(buffer, "%sOutflow.Only", DumpPath);
        OpenFile(&(channel->streamoutflow), buffer, "w", TRUE);
        // total incoming short wave
        sprintf(buffer, "%sISW.Only", DumpPath);
        OpenFile(&(channel->streamISW), buffer, "w", TRUE);
        //net incoming short wave
        sprintf(buffer, "%sNSW.Only", DumpPath);
        OpenFile(&(channel->streamNSW), buffer, "w", TRUE);
        // total incoming long wave
        sprintf(buffer, "%sILW.Only", DumpPath);
        OpenFile(&(channel->streamILW), buffer, "w", TRUE);
        // net incoming long wave
        sprintf(buffer, "%sNLW.Only", DumpPath);
        OpenFile(&(channel->streamNLW), buffer, "w", TRUE);
        //Vapor pressure
        sprintf(buffer, "%sVP.Only", DumpPath);
        OpenFile(&(channel->streamVP), buffer, "w", TRUE);
        //wind speed
        sprintf(buffer, "%sWND.Only", DumpPath);
        OpenFile(&(channel->streamWND), buffer, "w", TRUE);
        //air temperature
        sprintf(buffer, "%sATP.Only", DumpPath);
        OpenFile(&(channel->streamATP), buffer, "w", TRUE);
        //beam radiation
        sprintf(buffer, "%sBeam.Only", DumpPath);
        OpenFile(&(channel->streamBeam), buffer, "w", TRUE);
        //diffuse radiation
        sprintf(buffer, "%sDiffuse.Only", DumpPath);
        OpenFile(&(channel->streamDiffuse), buffer, "w", TRUE);
        //skyview
        sprintf(buffer, "%sSkyview.Only", DumpPath);
        OpenFile(&(channel->streamSkyView), buffer, "w", TRUE);
      }
    }
    if (channel->roads != NULL) {
      sprintf(buffer, "%sRoad.Flow", DumpPath);
      OpenFile(&(channel->roadout), buffer, "w", TRUE);
      sprintf(buffer, "%sRoadflow.Only", DumpPath);
      OpenFile(&(channel->roadflowout), buffer, "w", TRUE);

    }
  }
}

/* -------------------------------------------------------------
   ChannelCulvertFlow    
   computes outflow of channel/road network to a grid cell, if it
   contains a sink
   ------------------------------------------------------------- */
double ChannelCulvertFlow(int y, int x, CHANNEL * ChannelData)
{
  if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
    return channel_grid_outflow(ChannelData->road_map, x, y);
  }
  else {
    return 0;
  }
}

/* -------------------------------------------------------------
   RouteChannel
   ------------------------------------------------------------- */
void
RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
             TOPOPIX **TopoMap, SOILPIX **SoilMap, AGGREGATED *Total, 
	     OPTIONSTRUCT *Options, ROADSTRUCT **Network, SOILTABLE *SType, 
             PRECIPPIX **PrecipMap, float Tair, float Rh)
{
  int x, y;
  int flag;
  char buffer[32];
  float CulvertFlow;

  /* set flag to true if it's time to output channel network results */
  SPrintDate(&(Time->Current), buffer);
  flag = IsEqualTime(&(Time->Current), &(Time->Start));

  ParallelBarrier();

  if (ChannelData->roads != NULL) {

    /* give any surface water to roads w/o sinks */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          if (channel_grid_has_channel(ChannelData->road_map, x, y) && 
              !channel_grid_has_sink(ChannelData->road_map, x, y)) {	/* road w/o sink */
            SoilMap[y][x].RoadInt += SoilMap[y][x].IExcess; 
            channel_grid_inc_inflow(ChannelData->road_map, x, y, SoilMap[y][x].IExcess * Map->DX * Map->DY);
            SoilMap[y][x].IExcess = 0.0f;
          }
        }
      }
    }

    ParallelBarrier();

    ChannelGatherLateralInflow(ChannelData->roads, ChannelData->road_state_ga);

    if (ParallelRank() == 0) {
      /* route the road network and save results */
      channel_route_network(ChannelData->roads, Time->Dt);
      channel_save_outflow_text(buffer, ChannelData->roads,
                                ChannelData->roadout, ChannelData->roadflowout, flag);
    }

    ParallelBarrier();

    ChannelDistributeState(ChannelData->roads, ChannelData->road_state_ga);
  
    /* add culvert outflow to surface water */
    Total->CulvertReturnFlow = 0.0;
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          CulvertFlow = ChannelCulvertFlow(y, x, ChannelData);
          CulvertFlow /= Map->DX * Map->DY;
		
          /* CulvertFlow = (CulvertFlow > 0.0) ? CulvertFlow : 0.0; */
          if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
            channel_grid_inc_inflow(ChannelData->stream_map, x, y,
                                    (SoilMap[y][x].IExcess + CulvertFlow) * Map->DX * Map->DY);
            SoilMap[y][x].ChannelInt += SoilMap[y][x].IExcess;
            Total->CulvertToChannel += CulvertFlow;
            SoilMap[y][x].IExcess = 0.0f;
          }
          else {
            SoilMap[y][x].IExcess += CulvertFlow;
            Total->CulvertReturnFlow += CulvertFlow;
          }
        }
      }
    }
  }

  /* route stream channels */
  if (ChannelData->streams != NULL) {
    
    ChannelGatherLateralInflow(ChannelData->streams, ChannelData->stream_state_ga);

    if (ParallelRank() == 0) {
      channel_route_network(ChannelData->streams, Time->Dt);
      channel_save_outflow_text(buffer, ChannelData->streams,
                                ChannelData->streamout,
                                ChannelData->streamflowout, flag);
      /* save parameters for John's RBM model */
      if (Options->StreamTemp)
        channel_save_outflow_text_cplmt(Time, buffer,ChannelData->streams,ChannelData, flag);
    }

    ParallelBarrier();

    ChannelDistributeState(ChannelData->streams, ChannelData->stream_state_ga);
  }
  ParallelBarrier();
}

/* -------------------------------------------------------------
   ChannelCut
   computes necessary parameters for cell storage adjustment from
   channel/road dimensions
   ------------------------------------------------------------- */
void ChannelCut(int y, int x, CHANNEL * ChannelData, ROADSTRUCT * Network)
{
  float bank_height = 0.0;
  float cut_area = 0.0;

  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->stream_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->stream_map, x, y) * 
      channel_grid_cell_length(ChannelData->stream_map, x, y);
  }
  else if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->road_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->road_map, x, y) * 
      channel_grid_cell_length(ChannelData->road_map, x, y);
  }
  Network->Area = cut_area;
  Network->BankHeight = bank_height;
}

/* -------------------------------------------------------------
   ChannelFraction
   This computes the (sub)surface flow fraction for a road
   ------------------------------------------------------------- */
uchar ChannelFraction(TOPOPIX * topo, ChannelMapRec * rds)
{
  float effective_width = 0;
  float total_width;
  float sine, cosine;
  ChannelMapRec *r;
  float fract = 0.0;

  if (rds == NULL) {
    return 0;
  }
  cosine = cos(topo->Aspect);
  sine = sin(topo->Aspect);
  total_width = topo->FlowGrad / topo->Slope;
  effective_width = 0.0;

  for (r = rds; r != NULL; r = r->next) {
    effective_width += r->length * sin(fabs(topo->Aspect - r->aspect));
  }
  fract = effective_width / total_width * 255.0;
  fract = (fract > 255.0 ? 255.0 : floor(fract + 0.5));

  return (uchar) fract;
}

/* -------------------------------------------------------------
   DestroyChannel
   This completely destroys channel network data.
   ------------------------------------------------------------- */
void
DestroyChannel(OPTIONSTRUCT *Options, MAPSIZE *Map, CHANNEL *channel)
{
  ParallelBarrier();

  if (channel->streams != NULL) {
    channel_free_classes(channel->stream_class);
    channel_free_network(channel->streams);
    channel_grid_free_map(Map, channel->stream_map);
    /* GA_Destroy(channel->stream_state_ga); */
    if (ParallelRank() == 0) {
      fclose(channel->streamout);
      fclose(channel->streamflowout);
    }    
  }
  if (channel->roads != NULL) {
    channel_free_classes(channel->road_class);
    channel_free_network(channel->roads);
    channel_grid_free_map(Map, channel->road_map);
    /* GA_Destroy(channel->road_state_ga); */
    if (ParallelRank() == 0) {
      fclose(channel->roadout);
      fclose(channel->roadflowout);
    }    
  }

}
