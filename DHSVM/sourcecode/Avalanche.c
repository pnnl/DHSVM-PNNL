/*
 * SUMMARY:      Avalanche.c - Downslope movement of snow on steep slopes
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Chris Frans
 * ORG:          University of Washington
 * E-MAIL:       chrisf2@uw.edu
 * ORIG-DATE:    Feb-15
 * DESCRIPTION:  Represent Snow Redistribution
 * DESCRIP-END.
 * FUNCTIONS:    Avalanche()
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "slopeaspect.h"
#include "array_alloc.h"
#include "ParallelDHSVM.h"

 /*****************************************************************************
   Avalanche()

   Sources:
   Bernhardt, M., and K. Schulz (2010), SnowSlide: A simple routine for calculating
   gravitational snow transport, Geophys. Res. Lett., 37, L11502,
   doi:10.1029/2010GL043086.

   This routine follows Bernhardt, M., and K. Schulz (2010) in calculating gravitational
   redistribution of snow. Routing algorithms are similar to those used in
   RouteSubsurface.c. The local gradient is taken to be equal to the slope of the land
   surface. The slope is not currently calculated based on surface elevation of
   the snowpack.

   Set the gradient with pixels that are outside tha basin to zero.  This
   ensures that there is no flux of water across the basin boundary.

   WORK IN PROGRESS:
   Calculate slope based on Ice and Snow on top of topography.
   Transfer Cold Content of Snowpack with Mass.
 *****************************************************************************/
void Avalanche(MAPSIZE *Map, TOPOPIX **TopoMap, TIMESTRUCT *Time, OPTIONSTRUCT *Options,
  SNOWPIX **Snow)
{

  float Shd;                     /*Snow Holding Depth of a cell(m) as a function slope*/
  unsigned char ***SubDir;       /* Fraction of flux moving in each direction*/
  unsigned int **SubTotalDir;    /* Sum of Dir array */
  float **slope_deg;             /* Surface Slope in Degrees */
  float **SubSnowGrad;           /* Snow Surface Slope*/
  const char *Routine = "Avalanche";
  int x;                         /* counter */
  int y;                         /* counter */
  int i, j, k;
  float Snowout;
  int ga;
  GA_Patch patch;

  /*****************************************************************************
     Allocate memory
    ****************************************************************************/

  if (!(SubSnowGrad = calloc_2D_float(Map->NY, Map->NX)))
    ReportError((char *)Routine, 1);
  if (!(slope_deg = calloc_2D_float(Map->NY, Map->NX)))
    ReportError((char *)Routine, 1);
  if (!((SubDir) = calloc_3D_uchar(Map->NY, Map->NX, NDIRS)))
    ReportError((char *)Routine, 1);
  if (!(SubTotalDir = calloc_2D_uint(Map->NY, Map->NX)))
    ReportError((char *)Routine, 1);

  /* calculate snow surface slope in the same approach as subflow direction */
  SnowSlopeAspect(Map, TopoMap, Snow, SubSnowGrad, SubDir, SubTotalDir);
  
  ga = Map->dist;
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
        /* convert slope from radian to degree */
        slope_deg[y][x] = atan(SubSnowGrad[y][x])*(180 / PI);

        /* snow holding depth as a function of slope and slide parameters */
        Shd = SNOWSLIDE1*exp(-slope_deg[y][x] * SNOWSLIDE2);

        /* only redistribute snow if Swq is above holding capacity */
        if (slope_deg[y][x] > 30. && Snow[y][x].Swq > Shd) {

          /*If avalanche occurs on glacier surface, Leave a 10mm of snow behind so that glacier 
          surface is not prematurely exposed */
          /* if (Snow[y][x].Iwq > 1.0) {
            Snowout = Snow[y][x].Swq - 0.01;
            Snow[y][x].Swq = 0.01;
          } */
          
          Snowout = Snow[y][x].Swq;
          Snow[y][x].Swq = 0.0;
          
          /* This seems premature. It probably needs to be done after
             redistribution */
          Snow[y][x].TSurf = 0.0;
          Snow[y][x].TPack = 0.0;
          Snow[y][x].PackWater = 0.0;
          Snow[y][x].SurfWater = 0.0;

          /* Assign the avalanched snow to appropriate surrounding pixels */
          if (SubTotalDir[y][x] > 0) {
            Snowout /= (float)SubTotalDir[y][x];

            for (k = 0; k < NDIRS; k++) {
              int nx = xdirection[k] + x;
              int ny = ydirection[k] + y;
              if (valid_cell(Map, nx, ny)) {
                nx += patch.ixoff;
                ny += patch.iyoff;
                patch.patch[ny][nx] += Snowout * SubDir[y][x][k];
                /* Snow[ny][nx].Swq += Snowout * SubDir[y][x][k]; */
              }
            }

          }
          else {
            /* (WAP) If we get here, there are apparently no neighbors
               to receive the slide, which seems really unlikely given
               the slope_deg condition above (the test is still
               necessary, though). This should mean that Snowout goes
               back into cell x,y Swq. It should not just disappear,
               as I found it:

               Snowout = 0.0; */
            int nx = x + patch.ixoff;
            int ny = y + patch.iyoff;
            patch.patch[x + patch.ixoff][y + patch.iyoff] += Snowout;
            /* Snow[y][x].Swq = Snowout; */
          }
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
        Snow[y][x].Swq = patch.patch[y+patch.iyoff][x+patch.ixoff];
      }
    }
  }

  GA_Free_patch(&patch);

  free_3D_uchar(SubDir);
  free_2D_uint(SubTotalDir);
  free_2D_float(slope_deg);
  free_2D_float(SubSnowGrad);
}
