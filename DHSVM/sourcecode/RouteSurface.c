/*
 * SUMMARY:      RouteSurface.c - Route surface flow
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Route surface flow
 * DESCRIP-END.
 * FUNCTIONS:    RouteSurface()
 * Modification: Changes are made to exclude the impervious channel cell (with
 a non-zero impervious fraction) from surface routing. In the original
 code, some impervious channel cells are routed to themselves causing
 overestimated runoff in those cells (Ning, 2013).

 * $Id: RouteSurface.c, v3.1.2  2013/3/21   Ning Exp $
 */
#include <ga.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "ParallelDHSVM.h"
/*****************************************************************************
RouteSurface()
If the watertable calculated in WaterTableDepth() was negative, then water is
ponding on the surface.  At the moment no ponding of water is allowed in
DHSVM, and the "excess" water is routed to the outlet one pixel per time step
However, if the pixel contains an impervious fraction, then the surface water
is immediately routed to the nearest downslope pixel that contains a channel.
The net effect is that all pixels that have an impervious area are directly
connected (over the coarse of a single time step) to the channel network, this
assumption is likely to be true for small urban basins, and perhaps even for
large rural basins with some urban development
If Overland Routing = KINEMATIC, then "excess" water is routed to the outlet
using a infinite difference approximation to the kinematic wave solution of
the Saint-Venant equations.
*****************************************************************************/
void RouteSurface(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
                  SOILPIX ** SoilMap, OPTIONSTRUCT *Options,
                  UNITHYDR ** UnitHydrograph, UNITHYDRINFO * HydrographInfo, float *Hydrograph,
                  DUMPSTRUCT *Dump, VEGPIX ** VegMap, VEGTABLE * VType, CHANNEL *ChannelData)
{
  const char *Routine = "RouteSurface";
  static float one = 1.0;
  int Lag;			/* Lag time for hydrograph */
  int Step;
  float StreamFlow;
  int TravelTime;
  int WaveLength;
  int i, j, x, y, n;         /* Counters */
  int ga;
  float value;
  float theexcess;
  GA_Patch patch;


  /* Allocate memory for Runoff Matrix */
  if (Options->HasNetwork) {
    /* Option->Routing = false when routing = conventional */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          SoilMap[y][x].Runoff = SoilMap[y][x].IExcess;
          SoilMap[y][x].IExcess = 0;
          SoilMap[y][x].DetentionIn = 0;
        }
      }
    }

    /* ga = GA_Duplicate_type(Map->dist, "Subsurface Routing", C_FLOAT); */
    ga = Map->dist;
    value = 0.0;
    GA_Zero(ga);

    /* make a local array to store IExcess that covers the locally
       owned part of the domain plus ghost cells */
    GA_Alloc_patch_ghost(ga, Map, &patch);
    for (y = 0; y < patch.NY; ++y) {
      for (x = 0; x < patch.NX; ++x) {
        patch.patch[y][x] = 0.0;
      } 
    }

    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          if (!channel_grid_has_channel(ChannelData->stream_map, x, y)) {

            
            if (VType[VegMap[y][x].Veg - 1].ImpervFrac > 0.0) {

              theexcess = 0.0;

              /* Calculate the outflow from impervious portion of
                 urban cell straight to nearest channel cell */

              theexcess = (1 - VType[VegMap[y][x].Veg - 1].DetentionFrac) *
                VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].Runoff;

              /* Retained water in detention storage */

              SoilMap[y][x].DetentionIn = VType[VegMap[y][x].Veg - 1].DetentionFrac *
                VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].Runoff;

              /* Retained water in Detention storage routed to channel */

              SoilMap[y][x].DetentionStorage += SoilMap[y][x].DetentionIn;
              SoilMap[y][x].DetentionOut = 
                SoilMap[y][x].DetentionStorage * VType[VegMap[y][x].Veg - 1].DetentionDecay;

              theexcess += SoilMap[y][x].DetentionOut;

              SoilMap[y][x].DetentionStorage -= SoilMap[y][x].DetentionOut;
              if (SoilMap[y][x].DetentionStorage < 0.0)
                SoilMap[y][x].DetentionStorage = 0.0;

              /* there is no way to know what process owns the
                 (drains_x, drains_y) cell, so accumulate the GA value
                 for that cell. It shouldn't cost much if (drains_x,
                 drains_y) is local */

              GA_Acc_one_global(ga, Map, 
                         TopoMap[y][x].drains_x, TopoMap[y][x].drains_y, 
                         &theexcess, &one);

              /* Route the runoff from pervious portion of urban cell
                 to the neighboring cell, the local GA patch can be
                 used for this */
              
              for (n = 0; n < NDIRS; n++) {
                int xn = x + xdirection[n];
                int yn = y + ydirection[n];
                if (valid_cell(Map, xn, yn)) {
                  
                  xn += patch.ixoff;
                  yn += patch.iyoff;
                  patch.patch[yn][xn] += 
                    (1 - VType[VegMap[y][x].Veg - 1].ImpervFrac) * SoilMap[y][x].Runoff
                    *((float)TopoMap[y][x].Dir[n] / (float)TopoMap[y][x].TotalDir);
                }
              }
            }
            else {
              for (n = 0; n < NDIRS; n++) {
                int xn = x + xdirection[n];
                int yn = y + ydirection[n];
                if (valid_cell(Map, xn, yn)) {
                  /* SoilMap[yn][xn].IExcess += SoilMap[y][x].Runoff * */
                  /*   ((float)TopoMap[y][x].Dir[n] / (float)TopoMap[y][x].TotalDir); */
                  xn += patch.ixoff;
                  yn += patch.iyoff;
                  patch.patch[yn][xn] += SoilMap[y][x].Runoff *
                    ((float)TopoMap[y][x].Dir[n] / (float)TopoMap[y][x].TotalDir);
                }
              }
            }
          }
          else if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
            /* SoilMap[y][x].IExcess += SoilMap[y][x].Runoff; */
            patch.patch[y+patch.iyoff][x+patch.ixoff] += SoilMap[y][x].Runoff;
          }
        }
      }
    }

    GA_Acc_patch(ga, Map, &patch);
    ParallelBarrier();
    
    /* get the accumulated surface flow back from the GA (local array
       does not include ghosts) */

    GA_Get_patch(ga, Map, &patch);
    
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          SoilMap[y][x].IExcess = patch.patch[y+patch.iyoff][x+patch.ixoff];
        }
      }
    }

    GA_Free_patch(&patch);
    /* GA_Destroy(ga); */
    
  }/* end if Options->routing = conventional */

  /* MAKE SURE THIS WORKS WITH A TIMESTEP IN SECONDS */
  else {			/* No network, so use unit hydrograph method */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          TravelTime = (int)TopoMap[y][x].Travel;
          if (TravelTime != 0) {
            WaveLength = HydrographInfo->WaveLength[TravelTime - 1];
            for (Step = 0; Step < WaveLength; Step++) {
              Lag = UnitHydrograph[TravelTime - 1][Step].TimeStep;
              Hydrograph[Lag] += SoilMap[y][x].Runoff * UnitHydrograph[TravelTime - 1][Step].Fraction;

            }
            SoilMap[y][x].Runoff = 0.0;
          }
        }
      }
    }

    StreamFlow = 0.0;
    for (i = 0; i < Time->Dt; i++)
      StreamFlow += (Hydrograph[i] * Map->DX * Map->DY) / Time->Dt;

    /* Advance Hydrograph */
    for (i = 0; i < Time->Dt; i++) {
      for (j = 0; j < HydrographInfo->TotalWaveLength - 1; j++) {
        Hydrograph[j] = Hydrograph[j + 1];
      }

    }

    /* Set the last elements of the hydrograph to zero */
    for (i = 0; i < Time->Dt; i++)
      Hydrograph[HydrographInfo->TotalWaveLength - (i + 1)] = 0.0;

    PrintDate(&(Time->Current), Dump->Stream.FilePtr);
    fprintf(Dump->Stream.FilePtr, " %g\n", StreamFlow);
  }
}

