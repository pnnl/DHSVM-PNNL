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
          

          Snow[y][x].TSurf = 0.0;
          Snow[y][x].TPack = 0.0;
          Snow[y][x].PackWater = 0.0;
          Snow[y][x].SurfWater = 0.0;
          
          /* Assign the avalanched snow to appropriate surrounding pixels */
          if (SubTotalDir[y][x] > 0) {
            Snowout /= (float)SubTotalDir[y][x];

          }
          else {
            Snowout = 0.0;
            Snow[y][x].Swq = Snowout;
          }
          for (k = 0; k < NDIRS; k++) {
            int nx = xdirection[k] + x;
            int ny = ydirection[k] + y;
            if (valid_cell(Map, nx, ny)) {
              Snow[ny][nx].Swq += Snowout * SubDir[y][x][k];
            }
          }
        }
      }
    }
  }
  free_3D_uchar(SubDir);
  free_2D_uint(SubTotalDir);
  free_2D_float(slope_deg);
  free_2D_float(SubSnowGrad);
}
