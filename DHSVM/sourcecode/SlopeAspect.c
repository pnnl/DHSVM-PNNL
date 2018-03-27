/*
 * SUMMARY:      SlopeAspect.c - Calculate slope and aspect of each pixel
 * USAGE:        Part of DHSVM/MWM
 *
 * AUTHOR:       William A Perkins
 * ORG:          Battelle Memorial Institute Pacific Northwest Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    21-May-96
 * DESCRIPTION:  This module contains two routines to compute "slope" and
 *               "aspect"  (direction of slope): one which uses only terrain
 *               elevations and another which uses water table elevations.
 * DESCRIP-END.
 * FUNCTIONS:    valid_cell()
 *               valid_cell_fine()
 *               slope_aspect()
 *               flow_fractions()
 *               ElevationSlopeAspect()
 *               HeadSlopeAspect()
 *               ElevationSlope()
 *               ElevationSlopeAspectfine()
 * COMMENTS:
                 This program is considerably changed to fix the problems including:
				 1) runoff from some basins cell is rounted to the neighnoring cells 
				    that are outside of basin boundary
				 2) unfilled sinks due to the D8 and D4 algorithm difference between 
				    ArcGIS and DHSVM.
				 Main changes are made to slope_aspect(), flow_fractions() &
				 ElevationSlopeAspect()
 * $Id: SlopeAspect.c, v 4.0  2013/1/2   Ning Exp $ 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <ga.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"
#include "ParallelDHSVM.h"

/* These indices are so neighbors can be looked up quickly */
int xdirection[NDIRS] = {
  0, 1, 0, -1
};
int ydirection[NDIRS] = {
  -1, 0, 1, 0
};

float temp_aspect[NNEIGHBORS] = {
  225., 180., 135., 90., 45., 0., 315., 270.
};

/* NNEIGHBORS used for redistribution of subsurface flow in topoindex
 * and for slope/aspect calculations, must equal 8. */
int xneighbor[NNEIGHBORS] = {
-1, 0, 1, 1, 1, 0, -1, -1
};

int yneighbor[NNEIGHBORS] = {
 1, 1, 1, 0, -1, -1, -1, 0
};

/* -------------------------------------------------------------
   valid_cell
   Checks to see if grid indices, x and y, are within the grid 
   defined by the specified Map
   ------------------------------------------------------------- */
int valid_cell(MAPSIZE *Map, int x, int y)
{
  int gx, gy;
  gx = Map->OffsetX + x;
  gy = Map->OffsetY + y;
  return (gx >= 0 && gy >= 0 && gx < Map->gNX && gy < Map->gNY);
}

/* -------------------------------------------------------------
   slope_aspect
   Calculation of slope and aspect given elevations of cell and neighbors
   ------------------------------------------------------------- */
static void slope_aspect(float dx, float dy, float celev, float
			 nelev[NNEIGHBORS], float *slope, float *aspect)
{
  int n;
  float dzdx, dzdy;
  float *dummyelev;
  /* this dummy varaible is added for calculation of elev difference,
  in which the elev of OUTSIDEBASIN cells (which is ZERO) is 
  replaced by the elev of the central cell */

  /* allocate memory */
  if (!(dummyelev = (float*) calloc(NNEIGHBORS, sizeof(float))))
	  ReportError("slope_aspect( )", 1);
  
  for (n = 0; n < NNEIGHBORS; n++) {
      if (nelev[n] == OUTSIDEBASIN) {
		  dummyelev[n] = celev;
      }
	  else
		  dummyelev[n] = nelev[n];
    }
  dzdx = ((dummyelev[0] + 2 * dummyelev[7] + dummyelev[6]) -
	    (dummyelev[2] + 2 * dummyelev[3] + dummyelev[4])) / (8 * dx);
  dzdy = ((dummyelev[0] + 2 * dummyelev[1] + dummyelev[2]) -
	    (dummyelev[4] + 2 * dummyelev[5] + dummyelev[6])) / (8 * dy);

  *slope = sqrt(dzdx * dzdx + dzdy * dzdy);
  if (fequal(dzdx, 0.0) && fequal(dzdy, 0.0)) {
    *aspect = 0.0;
  }
  else {
	  /* convert from radian to degree */
	  *aspect = atan2(dzdx, dzdy) ;
  }
  free(dummyelev);
  return;
}
/* -------------------------------------------------------------
   flow_fractions
   Computes subsurface flow fractions given the slope and aspect 

   Comment: this function is considerably modified to avoid any
   out flow to the cells outside of the basin mask (Ning, 2013)
------------------------------------------------------------- */
static void flow_fractions(float dx, float dy, float slope, float aspect,
			   float nelev[NDIRS], float *grad,
			   unsigned char dir[NDIRS], unsigned int *total_dir)
{
  float cosine = cos(aspect);
  float sine = sin(aspect);
  float total_width, effective_width;
  float *cos, *sin;
  int n;
 
  /* allocate memory */
  if (!(cos = (float*) calloc(NDIRS/2, sizeof(float))))
	  ReportError("slope_aspect( )", 1);
  if (!(sin = (float*) calloc(NDIRS/2, sizeof(float))))
	  ReportError("slope_aspect( )", 1);

 switch (NDIRS) {
  case 4:
    /* fudge any cells which flow outside the basin by just pointing the
       aspect in the opposite direction */
    if (cosine > 0 && nelev[5] == (float) OUTSIDEBASIN)
		cos[1] = -cosine;
	else cos[1] = cosine;
	if (cosine < 0 && nelev[1] == (float) OUTSIDEBASIN)
        cos[0] = -cosine;
	else cos[0] = cosine;
    if (sine > 0 && nelev[3] == (float) OUTSIDEBASIN) 
		sin[0] = -sine;
	else sin[0] = sine;
	if (sine < 0 && nelev[7] == (float) OUTSIDEBASIN)
        sin[1] = -sine;
	else sin[1] = sine;

    /* compute flow widths */
    total_width = fabs(sine) * dx + fabs(cosine) * dy;
    *grad = slope * total_width;
    *total_dir = 0;
    for (n = 0; n < NDIRS; n++) 
	{
      switch (n) {
      case 0:
		  effective_width = (cos[1] > 0 ? cos[1] * dx : 0.0);
		  break;
      case 2:
		  effective_width = (cos[0] < 0 ? -cos[0] * dx : 0.0);
		  break;
      case 1:
		  effective_width = (sin[0] > 0 ? sin[0] * dy : 0.0);
		  break;
      case 3:
		  effective_width = (sin[1] < 0 ? -sin[1] * dy : 0.0);
		  break;
      default:
		  ReportError("flow_fractions",65);
		  assert(0);		
    }
	dir[n] = (int) ((effective_width / total_width) * 255.0 + 0.5);
    *total_dir += dir[n];
	}
	break;
  case 8:
    ReportError("flow_fractions",65);
    assert(0);			
    break;
  default:
    ReportError("flow_fractions",65);
    assert(0);			/* other cases don't work either */
  }
  free(sin);
  free(cos);
  return;
}
/* -------------------------------------------------------------
   ElevationSlopeAspect
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap)
{
  const char *Routine = "ElevationSlopeAspect";
  int x;
  int y;
  int n;
  int k;
  int i;
  float neighbor_elev[NNEIGHBORS];
  int steepestdirection;
  float min;
  int xn, yn;
  int gaelev;
  float elev, outelev;
  GA_Patch patch;
  

  outelev = (float) OUTSIDEBASIN;

  /* Count the number of cells in the basin.  Need this to allocate
     memory for the new, smaller Elev[] and Coords[][].  */
  Map->NumCells = 0;
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	Map->NumCells++;
      }
    }
  }

  /* make a global array to hold elevation */

  gaelev = GA_Duplicate_type(Map->dist, "ElevationSlopeAspect", GA_Type(NC_FLOAT));
  GA_Fill(gaelev, &outelev);
  ParallelBarrier();

  /* make a local array in which to put local elevations ghost cells not needed */
  GA_Alloc_patch(gaelev, Map, &patch);
  GA_Get_patch(gaelev, Map, &patch);

  /* put elevations in the GA patch */

  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        elev = TopoMap[y][x].Dem;
      } else {
        elev = outelev;
      }
      patch.patch[y+patch.iyoff][x+patch.ixoff] = elev;
    }
  }
  GA_Put_patch(gaelev, Map, &patch);
  ParallelBarrier();
  GA_Free_patch(&patch);

  /* make a local array and fill it with elevations from GA that
     includes ghost cells */
  GA_Alloc_patch_ghost(gaelev, Map, &patch);
  GA_Get_patch(gaelev, Map, &patch);

  /* for each local, fill neighbor array with elevations from GA patch */
  
  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	for (n = 0; n < NNEIGHBORS; n++) {
	  xn = x + xneighbor[n];
	  yn = y + yneighbor[n];	  
	  if (valid_cell(Map, xn, yn)) {
            elev = patch.patch[yn+patch.iyoff][xn+patch.ixoff];
          } else {
            elev = outelev;
          }
          neighbor_elev[n] = elev;
        }

        slope_aspect(Map->DX, Map->DY, TopoMap[y][x].Dem, neighbor_elev,
		     &(TopoMap[y][x].Slope), &(TopoMap[y][x].Aspect));	

	/* fill Dirs in TopoMap too */
	flow_fractions(Map->DX, Map->DY, TopoMap[y][x].Slope,
		       TopoMap[y][x].Aspect,
		       neighbor_elev, &(TopoMap[y][x].FlowGrad),
		       TopoMap[y][x].Dir, &(TopoMap[y][x].TotalDir));
        
	/* If there is a sink, check again to see if there 
	   is a direction of steepest descent. Does not account 
	   for ties.*/
	if(TopoMap[y][x].TotalDir <= 0) {
	  steepestdirection = -99;
	  min = DHSVM_HUGE;	       
	  for (n = 0; n < NDIRS; n++) {
	    xn = x + xdirection[n];
	    yn = y + ydirection[n];	  
            if (valid_cell(Map, xn, yn)) {
              elev = patch.patch[yn+patch.iyoff][xn+patch.ixoff];
              /* elev = neighbor_elev[n]; */
              if (INBASIN((int)elev)) {
                if (elev < min) { 
                  min = elev;
                  steepestdirection = n;
                }
              }
            }
	  }
	  if(min < TopoMap[y][x].Dem) {
	    TopoMap[y][x].Dir[steepestdirection] = (int)(255.0 + 0.5);
	    TopoMap[y][x].TotalDir = (int)(255.0 + 0.5);
	  }
	  else {
	    /*  Last resort: set the Dir of the cell to the cell that is
		closest in elevation. This should only happen for the 
		basin outlet, unless the Dem wasn't filled. */	  
	    TopoMap[y][x].Dir[steepestdirection] = (int)(255.0 + 0.5);
	    TopoMap[y][x].TotalDir = (int)(255.0 + 0.5);	    
	  }
	}
      }
    }
  } 

  GA_Free_patch(&patch);  
  GA_Destroy(gaelev);
  ParallelBarrier();


  /* FIXME: Are the OrderedCells actually used anywhere? */
  /* Create a structure to hold elevations of only those cells
     within the basin and the y,x of those cells.*/
  if (!(Map->OrderedCells = (ITEM *) calloc(Map->NumCells, sizeof(ITEM))))
    ReportError((char *) Routine, 1);
  k = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      /* Save the elevation, y, and x in the ITEM structure. */
      if (INBASIN(TopoMap[y][x].Mask)) {
        Map->OrderedCells[k].Rank = TopoMap[y][x].Dem;
        Map->OrderedCells[k].y = y;
        Map->OrderedCells[k].x = x;
        k++;
      }
    }
  }
  /* Sort Elev in descending order-- Elev.x and Elev.y hold indices. */ 
  quick(Map->OrderedCells, Map->NumCells);

  /* End of modifications to create ordered cell coordinates.  SRW 10/02, LCB 03/03 */

  n = Map->NumCells;
  printf("%d: local NumCells = %d\n", ParallelRank(), n);

  i = n;
  k = n;
  GA_Igop(&n, 1, "+");
  GA_Igop(&i, 1, "max");
  GA_Igop(&k, 1, "min");
  if (ParallelRank() == 0) 
    printf("global NumCells = %d, local Maximum = %d, local Minimum = %d\n", n, i, k);
  Map->AllCells = n;

  return;
}

/* -------------------------------------------------------------
   QuickSort
   ------------------------------------------------------------- */

/**********************************************************************
        this subroutine starts the quick sort
**********************************************************************/

void quick(ITEM *OrderedCells, int count)
{
  qs(OrderedCells,0,count-1);
}

void qs(ITEM *item, int left, int right)
/**********************************************************************
        this is the quick sort subroutine - it returns the values in
        an array from high to low.
**********************************************************************/
{
  register int i,j;
  ITEM x,y;

  i=left;
  j=right;
  x=item[(left+right)/2];

  do {
    while(item[i].Rank<x.Rank && i<right) i++;
    while(x.Rank<item[j].Rank && j>left) j--;

    if (i<=j) {
      y=item[i];
      item[i]=item[j];
      item[j]=y;
      i++;
      j--;
    }
  } while (i<=j);
  if(left<j) qs(item,left,j);
  if(i<right) qs(item,i,right);
}
/* -------------------------------------------------------------
   HeadSlopeAspect
   This computes slope and aspect using the water table elevation. 

   Comment: rewritten to fill the sinks (Ning, 2013)
   ------------------------------------------------------------- */
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap,
		     float **FlowGrad, unsigned char ***Dir, unsigned int **TotalDir)
{
  int x;
  int y;
  int n;
  float neighbor_elev[NNEIGHBORS];
  float outelev, elev;
  int gaelev;

  /* make a global array to hold elevation */

  gaelev = GA_Duplicate_type(Map->dist, "WaterLevelSlopeAspect", GA_Type(NC_FLOAT));

  outelev = (float) OUTSIDEBASIN;

  GA_Fill(gaelev, &outelev);

  ParallelBarrier();

  /* put water levels in GA (presumably, all of these puts should be
     local, so it may be OK do put one value at a time; maybe try
     non-blocking too) */

  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        elev = SoilMap[y][x].WaterLevel;
        GA_Put_one(gaelev, Map, x, y, &elev);
      }
    }
  }
  GA_Update_ghosts(gaelev);


  /* let's assume for now that WaterLevel is the SOILPIX map is
     computed elsewhere */
  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        float slope, aspect;
        for (n = 0; n < NNEIGHBORS; n++) {
          int xn = x + xneighbor[n];
          int yn = y + yneighbor[n];			  
          if (valid_cell(Map, xn, yn)) {
            GA_Get_one(gaelev, Map, xn, yn, &elev);
            neighbor_elev[n] = elev;
          } 
        }
        slope_aspect(Map->DX, Map->DY, SoilMap[y][x].WaterLevel, neighbor_elev,
		     &slope, &aspect);
        flow_fractions(Map->DX, Map->DY, slope, aspect, neighbor_elev,
		       &(FlowGrad[y][x]), Dir[y][x], &(TotalDir[y][x])); 
      }
    }
  }
  GA_Destroy(gaelev);
  return;
}



/* -------------------------------------------------------------
SnowSlopeAspect
This computes slope and aspect using the SnowSurface Elevation.
------------------------------------------------------------- */
void SnowSlopeAspect(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow,
  float **SubSnowGrad, unsigned char ***Dir, unsigned int **TotalDir)
{
  int x;
  int y;
  int n;
  float neighbor_elev[NNEIGHBORS];

  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        float slope, aspect;
        for (n = 0; n < NNEIGHBORS; n++) {
          int xn = x + xneighbor[n];
          int yn = y + yneighbor[n];
          if (valid_cell(Map, xn, yn)) {
            /* snow elevation (swq+dem) of neighboring cells */
            neighbor_elev[n] =
              ((TopoMap[yn][xn].Mask) ? (TopoMap[yn][xn].Dem + Snow[yn][xn].Swq) : (float)OUTSIDEBASIN);
          }
          else {
            neighbor_elev[n] = (float)OUTSIDEBASIN;
          }
        }

        slope_aspect(Map->DX, Map->DY, (TopoMap[y][x].Dem + Snow[y][x].Swq), neighbor_elev,
          &slope, &aspect);
        flow_fractions(Map->DX, Map->DY, slope, aspect, neighbor_elev,
          &(SubSnowGrad[y][x]), Dir[y][x], &(TotalDir[y][x]));

        /* Reset SubSnowGrad to slope, don't want width in computation */
        SubSnowGrad[y][x] = slope;
      }
    }
  }
  return;
}




















