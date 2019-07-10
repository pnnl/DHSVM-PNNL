/* -------------------------------------------------------------
   file: DHSVMChannel.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created August 30, 1996 by  William A Perkins
   $Id: DHSVMChannel.h,v3.1.2 2013/10/03  Ning Exp $
   ------------------------------------------------------------- */

#ifndef _DHSVMChannel_h_
#define _DHSVMChannel_h_

#include "settings.h"		/* for data.h */
#include "data.h"
#include "getinit.h"
#include "channel.h"
#include "channel_grid.h"

/* -------------------------------------------------------------
   struct CHANNEL
   ------------------------------------------------------------- */
typedef struct {
  ChannelClass *stream_class;
  ChannelClass *road_class;
  Channel *streams;
  Channel *roads;
  ChannelMapPtr **stream_map;
  ChannelMapPtr **road_map;
  int stream_state_ga;
  int road_state_ga;
  FILE *streamout;
  FILE *roadout;
  FILE *streamflowout;
  FILE *roadflowout;
  /* new output files for John's RBM model */
  FILE *streaminflow;
  FILE *streamoutflow;
  FILE *streamNSW;
  FILE *streamNLW;
  FILE *streamVP;
  FILE *streamWND;
  FILE *streamATP;
} CHANNEL;

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void InitChannel(LISTPTR Input, MAPSIZE *Map, int deltat, CHANNEL *channel,
		 SOILPIX **SoilMap, int *MaxStreamID, int *MaxRoadID, OPTIONSTRUCT *Options);
void InitChannelDump(OPTIONSTRUCT *Options, CHANNEL *channel, char *DumpPath);
double ChannelCulvertFlow(int y, int x, CHANNEL *ChannelData);
void RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
		  TOPOPIX **TopoMap, SOILPIX **SoilMap, AGGREGATED *Total, 
		  OPTIONSTRUCT *Options, ROADSTRUCT **Network, SOILTABLE *SType, 
		  PRECIPPIX **PrecipMap, float Tair, float Rh);
void ChannelCut(int y, int x, CHANNEL *ChannelData, ROADSTRUCT *Network);
uchar ChannelFraction(TOPOPIX *topo, ChannelMapRec *rds);
void DestroyChannel(OPTIONSTRUCT *Options, MAPSIZE *Map, CHANNEL *channel);

#endif
