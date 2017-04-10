/*
 * SUMMARY:      InitNetwork.c - Initialize road/channel work
 * USAGE:
 *
 * AUTHOR:       DHSVM Project (Bart Nijssen)
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    27-Aug-1996 at 18:34:01
 * DESCRIPTION:  Initialize road/channel work.  Memory is allocated, and the
 *               necessary adjustments for the soil profile are calculated
 * DESCRIP-END.
 * FUNCTIONS:    InitNetwork()
 * COMMENTS:
 * $Id: InitNetwork.c,v 1.8 2004/05/03 03:28:45 colleen Exp $
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "settings.h"
#include "soilmoisture.h"
#include "DHSVMChannel.h"
#include "ParallelDHSVM.h"

 /*****************************************************************************
   Function name: InitNetwork()

   Purpose      : Initialize road/channel work.  Memory is allocated, and the
                  necessary adjustments for the soil profile are calculated

   Comments     :
 *****************************************************************************/
void InitNetwork(MAPSIZE *Map, TOPOPIX **TopoMap,
  SOILPIX **SoilMap, VEGPIX **VegMap, VEGTABLE *VType,
  ROADSTRUCT ***Network, CHANNEL *ChannelData,
  LAYER Veg, OPTIONSTRUCT *Options)
{
  const char *Routine = "InitNetwork";
  int i;			/* counter */
  int x;			/* column counter */
  int y;			/* row counter */
  int sx, sy;
  int minx, miny;
  int doimpervious;
  int numroads;          /* Counter of number of pixels
                with a road and channel */
  int numroadschan;      /* Counter of number of pixels
                with a road */
  FILE *inputfile;
  /* Allocate memory for network structure */

  if (!(*Network = (ROADSTRUCT **)calloc(Map->NY, sizeof(ROADSTRUCT *))))
    ReportError((char *)Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*Network)[y] = (ROADSTRUCT *)calloc(Map->NX, sizeof(ROADSTRUCT))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        if (!((*Network)[y][x].Adjust =
          (float *)calloc((VType[VegMap[y][x].Veg - 1].NSoilLayers + 1), sizeof(float))))
          ReportError((char *)Routine, 1);

        if (!((*Network)[y][x].PercArea =
          (float *)calloc((VType[VegMap[y][x].Veg - 1].NSoilLayers + 1), sizeof(float))))
          ReportError((char *)Routine, 1);
      }
    }
  }

  numroadschan = 0;
  numroads = 0;

  /* If a road/channel Network is imposed on the area, read the Network
     information, and calculate the storage adjustment factors */
  if (Options->HasNetwork) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          ChannelCut(y, x, ChannelData, &((*Network)[y][x]));
          AdjustStorage(VType[VegMap[y][x].Veg - 1].NSoilLayers,
            SoilMap[y][x].Depth,
            VType[VegMap[y][x].Veg - 1].RootDepth,
            (*Network)[y][x].Area, Map->DX, Map->DY,
            (*Network)[y][x].BankHeight,
            (*Network)[y][x].PercArea,
            (*Network)[y][x].Adjust,
            &((*Network)[y][x].CutBankZone));
          (*Network)[y][x].IExcess = 0.;
          if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
            numroads++;
            if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
              numroadschan++;
            }
            (*Network)[y][x].fraction =
              ChannelFraction(&(TopoMap[y][x]), ChannelData->road_map[x][y]);
            (*Network)[y][x].MaxInfiltrationRate =
              MaxRoadInfiltration(ChannelData->road_map, x, y);
            (*Network)[y][x].RoadClass =
              channel_grid_class(ChannelData->road_map, x, y);
            (*Network)[y][x].FlowSlope =
              channel_grid_flowslope(ChannelData->road_map, x, y);
            (*Network)[y][x].FlowLength =
              channel_grid_flowlength(ChannelData->road_map, x, y, (*Network)[y][x].FlowSlope);
            (*Network)[y][x].RoadArea = channel_grid_cell_width(ChannelData->road_map, x, y) * channel_grid_cell_length(ChannelData->road_map, x, y);
          }
          else {
            (*Network)[y][x].MaxInfiltrationRate = DHSVM_HUGE;
            (*Network)[y][x].FlowSlope = 0.;
            (*Network)[y][x].FlowLength = 0.;
            (*Network)[y][x].RoadArea = 0.;
            (*Network)[y][x].RoadClass = NULL;
            (*Network)[y][x].IExcess = 0.;
          }
        }
      }
    }
  }
  /* if no road/channel Network is imposed, set the adjustment factors to the
     values they have in the absence of an imposed network */
  else {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          for (i = 0; i <= VType[VegMap[y][x].Veg - 1].NSoilLayers; i++) {
            (*Network)[y][x].Adjust[i] = 1.0;
            (*Network)[y][x].PercArea[i] = 1.0;
            (*Network)[y][x].CutBankZone = NO_CUT;
            (*Network)[y][x].MaxInfiltrationRate = 0.;
          }
          (*Network)[y][x].FlowSlope = 0.;
          (*Network)[y][x].FlowLength = 0.;
          (*Network)[y][x].RoadArea = 0.;
          (*Network)[y][x].RoadClass = NULL;
          (*Network)[y][x].IExcess = 0.;
        }
      }
    }
  }

  if (numroads > 0) {
    printf("There are %d pixels with a road and %d with a road and a channel.\n",
      numroads, numroadschan);
  }

  /* this all pertains to the impervious surface */

  doimpervious = 0;
  for (i = 0; i < Veg.NTypes; i++)
    if (VType[i].ImpervFrac > 0.0)
      doimpervious = 1;

  if (doimpervious) {
    if (!(inputfile = fopen(Options->ImperviousFilePath, "rt"))) {
      fprintf(stderr,
        "User has specified a percentage impervious area \n");
      fprintf(stderr,
        "To incorporate impervious area, DHSVM needs a file\n");
      fprintf(stderr,
        "identified by the key: \"IMPERVIOUS SURFACE ROUTING FILE\",\n");
      fprintf(stderr,
        "in the \"VEGETATION\" section.  This file is used to \n");
      fprintf(stderr,
        "determine the fate of surface runoff, i.e. where does it go \n");
      fprintf(stderr,
        "This file was not found: see InitNetwork.c \n");
      fprintf(stderr,
        "The code find_nearest_channel.c will make the file\n");
      ReportError(Options->ImperviousFilePath, 3);
    }

    /* when the impervious surface routing file is read by multiple
       processes, indexes cannot be checked as they are when the
       entire domain is on one process; so, run serial to check before
       using in parallel */

    /* Is this really the best approach for reading this information */

    if (ParallelSize() == 1) {
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN(TopoMap[y][x].Mask)) {
            if (fscanf(inputfile, "%d %d %d %d \n", &sy, &sx, &miny, &minx) !=
                EOF) {
              TopoMap[y][x].drains_x = minx;
              TopoMap[y][x].drains_y = miny;
            }
            else {
              ReportError(Options->ImperviousFilePath, 63);
            }
                                /* indexing check */
            if (sx != x || sy != y) {
              ReportError(Options->ImperviousFilePath, 64);
            }
          }
        }
      }
    } else {
                                /* there is supposed to be a line for
                                   each active cell in the domain (for
                                   some reason), so let's read that
                                   many lines */

      for (i = 0; i < Map->AllCells; ++i) {
        if (fscanf(inputfile, "%d %d %d %d \n", &sy, &sx, &miny, &minx) != EOF) {
          if (Global2Local(Map, sx, sy, &x, &y)) {
                                /* NOTE minx, miny are global indexes */
            TopoMap[y][x].drains_x = minx;
            TopoMap[y][x].drains_y = miny;
          }
        } else {
          ReportError(Options->ImperviousFilePath, 63);
        }
      }
    }
  }
}
