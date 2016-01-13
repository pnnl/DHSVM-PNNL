/*
 * SUMMARY:      SlopeAspect.c - Calculate slope and aspect of each pixel
 * USAGE:        Part of DHSVM/MWM
 *
 * AUTHOR:       William A Perkins
 * ORG:          Battelle Memorial Institute Pacific Northwest Laboratory
 * E-MAIL:       perk@clio.muse.pnl.gov
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
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"

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
int valid_cell(MAPSIZE * Map, int x, int y)
{
  return (x >= 0 && y >= 0 && x < Map->NX && y < Map->NY);
}
/******************************************************************************/
/*   valid_cell_fine                                                          */
/*   Checks to see if grid indices, x and y, are within the grid              */
/*   defined by the specified Map                                             */ 
/******************************************************************************/

int valid_cell_fine(MAPSIZE *Map, int x, int y) 
{
  return (x >= 0 && y >= 0 && x < Map->NXfine && y < Map->NYfine);
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
  float neighbor_elev[NNEIGHBORS];
  int steepestdirection;
  float min;
  int xn, yn;

  /* fill neighbor array */
  
  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	/* Count the number of cells in the basin.  
	   Need this to allocate memory for
	   the new, smaller Elev[] and Coords[][].  */
	Map->NumCells++;
	for (n = 0; n < NNEIGHBORS; n++) {
	  xn = x + xneighbor[n];
	  yn = y + yneighbor[n];	  
	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] = ((TopoMap[yn][xn].Mask) ? TopoMap[yn][xn].Dem : (float) OUTSIDEBASIN);
	  }
	  else {
	    neighbor_elev[n] = (float) OUTSIDEBASIN;
	  }
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
	if(TopoMap[y][x].TotalDir == 0) {
	  steepestdirection = -99;
	  min = DHSVM_HUGE;	       
	  for (n = 0; n < NDIRS; n++) {
	    xn = x + xdirection[n];
	    yn = y + ydirection[n];	  
	    if (valid_cell(Map, xn, yn)) {
	      if (INBASIN(TopoMap[yn][xn].Mask)) {
			  if(TopoMap[yn][xn].Dem < min) { 
				  min = TopoMap[yn][xn].Dem;
				  steepestdirection = n;}
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
	    xn = x + xdirection[steepestdirection];
	    yn = y + ydirection[steepestdirection];
	  }
	}
   } 
  }
 } 	
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
				  neighbor_elev[n] =
					  ((TopoMap[yn][xn].Mask) ? SoilMap[yn][xn].WaterLevel : (float) OUTSIDEBASIN);
			  }
			  else {
				  neighbor_elev[n] = (float) OUTSIDEBASIN;
			  }
		  }
		  slope_aspect(Map->DX, Map->DY, SoilMap[y][x].WaterLevel, neighbor_elev,
		     &slope, &aspect);
		  flow_fractions(Map->DX, Map->DY, slope, aspect, neighbor_elev,
		       &(FlowGrad[y][x]), Dir[y][x], &(TotalDir[y][x])); 
      }
    }
  }
  return;
}
/******************************************************************************/
/*			     ElevationSlope                            */
/* Part of MWM, should probably be merged w/ ElevationSlopeAspect function.   */
/******************************************************************************/

float ElevationSlope(MAPSIZE *Map, TOPOPIX **TopoMap, FINEPIX ***FineMap, int y, int x, int *nexty, 
		     int *nextx, int prevy, int prevx, float *Aspect) 
{
  int n, direction;
  float soil_elev[NNEIGHBORS];
  float bedrock_elev[NNEIGHBORS];
  float Slope;
  float temp_slope[NNEIGHBORS];
  double length_diagonal;
  float dx, dy, celev;
  int coarsej, coarsei;

  /* fill neighbor array */
  
  for (n = 0; n < NNEIGHBORS; n++) {

    int xn = x + xneighbor[n];
    int yn = y + yneighbor[n];
   
    // Initialize soil_elev and bedrock_elev
    soil_elev[n] = (float) OUTSIDEBASIN;
    bedrock_elev[n] = (float) OUTSIDEBASIN;

    // Check whether yn, xn are within FineMap array bounds
    if (valid_cell_fine(Map,xn,yn)){

      coarsej = floor(yn*Map->DMASS/Map->DY);
      coarsei = floor(xn*Map->DMASS/Map->DX);

      // Check whether FineMap element has been allocated for this cell
      // (equivalent to checking whether parent coarse grid cell is within coarse mask)
      if (INBASIN(TopoMap[coarsej][coarsei].Mask)) { 

	bedrock_elev[n] = (((*FineMap[yn][xn]).Mask) ? (*FineMap[yn][xn]).bedrock : (float) OUTSIDEBASIN);
	soil_elev[n] = (((*FineMap[yn][xn]).Mask) ? (*FineMap[yn][xn]).bedrock+(*FineMap[yn][xn]).sediment : (float) OUTSIDEBASIN);
	
      }
    }
    
  }       
  /*  Find bedrock slope in all directions. Negative slope = ascent, positive slope = descent.  */     
  dx = Map->DMASS;
  dy = Map->DMASS;
  celev = (*FineMap[y][x]).bedrock;


  length_diagonal = sqrt((pow((double)dx, (double)2)) + (pow((double)dy, (double)2))); 

  for (n = 0; n < NNEIGHBORS; n++) {
    if (bedrock_elev[n] == OUTSIDEBASIN) 
      bedrock_elev[n] = DHSVM_HUGE;
    
    if(n==0 || n==2 || n==4 || n==6)
      temp_slope[n] = (atan((celev - bedrock_elev[n]) / length_diagonal))
	* DEGPRAD;
    else if(n==1 || n==5)
      temp_slope[n] = (atan((celev - bedrock_elev[n]) / dy)) * DEGPRAD;
    else
      temp_slope[n] = (atan((celev - bedrock_elev[n]) / dx)) * DEGPRAD;
  }
    
 /* Find largest (positive) slope, this is the direction of failure along bedrock plain.  
     Backtracking isn't a problem if using the bedrock, but sinks may exist. */ 
   
  Slope = -999.;
  *Aspect = -99.;

  for (n = 0; n < NNEIGHBORS; n++){
    if(temp_slope[n] > Slope) {
      Slope = temp_slope[n];
      *Aspect = temp_aspect[n] * PI / 180.0;
      direction = n;
      *nexty = y + yneighbor[n];
      *nextx = x + xneighbor[n];
    }
  }

  /* If no positive slope found, a bedrock sink was encountered.  Assuming the 
     sink should be filled to the lowest "pour elevation", aspect should have 
     already been assigned correctly. */

  /* Find dynamic slope in direction of steepest descent. */
  
  celev = (*FineMap[y][x]).bedrock + (*FineMap[y][x]).sediment;
  if(direction==0 || direction==2 || direction==4 || direction==6)
    Slope = (atan((celev - soil_elev[direction]) / length_diagonal))
      * DEGPRAD;
  else if(direction==1 || direction==5)
    Slope = (atan((celev - soil_elev[direction]) / dy)) * DEGPRAD;
  else
    Slope = (atan((celev - soil_elev[direction]) / dx)) * DEGPRAD;

  /* It is possible that a "soil" sink could be encountered at this point.  
     This is not really an error, and is checked for in MainMWM. */
  // if(Slope < 0.0 ) {
  // fprintf(stderr, "Sink encountered in cell y= %d x= %d, all routes from here go up!\n", y,x);
  // }

  if(Slope == -999. || *Aspect == -99.) {
    fprintf(stderr, "Aspect not assigned, this shouldn't have happened.\n");
    exit(0);
  }

  return Slope;
}
























