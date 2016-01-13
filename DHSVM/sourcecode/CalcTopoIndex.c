/*
 * SUMMARY:      CalcTopoIndex.c - Calculate the topographic index
 * USAGE:        Part of MWM
 *
 * AUTHOR:       Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * DESCRIPTION:  Calculate the topographic index for the redistribution of 
                 soil moisture from the coarse grid to the fine grid. 
 * DESCRIP-END.
 * FUNCTIONS:    CalcTopoIndex()
 * COMMENTS:
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "DHSVMerror.h"
#include "constants.h"
#include "data.h"
#include "settings.h"
#include "slopeaspect.h"

#define VERTRES 1 /* vertical resolution of the dem */

/*****************************************************************************
  Function name: CalcTopoIndex()

  Purpose      : Calculate the topographic index for the redistribution of 
                 soil moisture from the coarse grid to the fine grid.  
                 
  Required     :
    MAPSIZE *Map      -  mass wasting resolution map data
    FINEPIX **FineMap -  mask and dem for mass wasting resoltuion map
   
  Returns      : float, topographic index

  Modifies     : none
   
  Comments     : Based on the topographic index, ln (a/tan beta) from TOPMODEL
                 (Beven & Kirkby 1979). Calculated according to Wolock 1995.

  The surrounding grid cells are numbered in the following way

                |-----| DX

          0-----1-----2  ---
          |\    |    /|   |
          | \   |   / |   |
          |  \  |  /  |   | DY
          |   \ | /   |   |
          |    \|/    |   |
          7-----*-----3  ---
          |    /|\    |
          |   / | \   |
          |  /  |  \  |
          | /   |   \ |
          |/    |    \|
          6-----5-----4

  For the current implementation it is assumed that the resolution is the 
  same in both the X and the Y direction.  

  Source       : Beven K.J. and M.J. Kirkby, 1979, A physically based, variable 
                 contributing area model of basin hydrology, Hydrol Sci Bull 24, 
                 43-69.
        
                 Wolock, David M. and Gregory J. McCabe, Jr., 1995, Comparison
                 of single and multiple flow direction algorithms for computing
                 topographic parameters in TOPMODEL, Water Resources Research,
                 31 (5), 1315-1324.
*****************************************************************************/

void CalcTopoIndex(MAPSIZE *Map, FINEPIX ***FineMap, TOPOPIX **TopoMap)
{
  FILE *fo;
  char topoindexmap[100];
  int printmap;
  int i, j, k, x, y, n, lower;  /* counters */
  int dx, dy;
  float celev;
  float neighbor_elev[NNEIGHBORS];
  float temp_slope[NNEIGHBORS];
  float delta_a[NNEIGHBORS];
  float length_diagonal;
  float **a;                    /* Area of hillslope per unit contour (m2) */
  float **tanbeta;
  float **contour_length;
  int coarsei, coarsej;
  ITEM *OrderedCellsfine;

  /* These indices are so neighbors can be looked up quickly */
  int xneighbor[NNEIGHBORS] = {
#if NNEIGHBORS == 4
    0, 1, 0, -1
#elif NNEIGHBORS == 8
    -1, 0, 1, 1, 1, 0, -1, -1
#endif
  };
  int yneighbor[NNEIGHBORS] = {
#if NNEIGHBORS == 4
    -1, 0, 1, 0
#elif NNEIGHBORS == 8
    1, 1, 1, 0, -1, -1, -1, 0
#endif
  };
  
  /**************************************************************************
   Allocate memory 
  ***************************************************************************/
  if (!(a = (float **)calloc(Map->NYfine, sizeof(float *))))
    ReportError("CalcTopoIndex", 1);
  for(i=0; i<Map->NYfine; i++) {
    if (!(a[i] = (float *)calloc(Map->NXfine, sizeof(float))))
      ReportError("CalcTopoIndex", 1);
  }
  
  if (!(tanbeta = (float **)calloc(Map->NYfine, sizeof(float *))))
    ReportError("CalcTopoIndex", 1);
  for(i=0; i<Map->NYfine; i++) {
    if (!(tanbeta[i] = (float *)calloc(Map->NXfine, sizeof(float))))
      ReportError("CalcTopoIndex", 1);
  }
  
  if (!(contour_length = (float **)calloc(Map->NYfine, sizeof(float *))))
    ReportError("CalcTopoIndex", 1);
  for(i=0; i<Map->NYfine; i++) {
    if (!(contour_length[i] = (float *)calloc(Map->NXfine, sizeof(float))))
      ReportError("CalcTopoIndex", 1);
  }
  
  if (!(OrderedCellsfine = (ITEM *) calloc(Map->NumCellsfine, sizeof(ITEM))))
    ReportError("CalcTopoIndex", 1);
  
  dx = Map->DMASS;
  dy = Map->DMASS;
  length_diagonal = sqrt((pow(dx, 2)) + (pow(dy, 2)));
  
  k = 0;
  for (y = 0; y < Map->NYfine; y++) {
    for (x = 0; x < Map->NXfine; x++) {
      
      coarsei = floor(y*Map->DMASS/Map->DY);
      coarsej = floor(x*Map->DMASS/Map->DX);
      
      /* Save the elevation, y, and x in the ITEM structure. */
      if (INBASIN(TopoMap[coarsei][coarsej].Mask)) {
	OrderedCellsfine[k].Rank = (*FineMap[y][x]).Dem;
	OrderedCellsfine[k].y = y;
	OrderedCellsfine[k].x = x;
	k++;
      }
    }
  }
 
  /* Sort Elev in descending order-- Elev.x and Elev.y hold indices. */
  
  quick(OrderedCellsfine, Map->NumCellsfine);

  for (k = (Map->NumCellsfine)-1; k >-1; k--) { 
    y = OrderedCellsfine[k].y;
    x = OrderedCellsfine[k].x;
    a[y][x]= dx * dy; /* Intialize cummulative area to cell area  */
  }
  
  /* Loop through all cells in descending order of elevation */
  for (k = (Map->NumCellsfine)-1; k >-1; k--) { 
    y = OrderedCellsfine[k].y;
    x = OrderedCellsfine[k].x;
    
    /* fill neighbor array */
    for (n = 0; n < NNEIGHBORS; n++) {

      int xn = x + xneighbor[n];
      int yn = y + yneighbor[n];
   
      // Initialize neighbor_elev
      neighbor_elev[n] = (float) OUTSIDEBASIN;

      // Check whether yn, xn are within FineMap array bounds
      if (valid_cell_fine(Map,xn,yn)){

        coarsei = floor(yn*Map->DMASS/Map->DY);
        coarsej = floor(xn*Map->DMASS/Map->DX);
      
        // Check whether FineMap element has been allocated for this cell
        // (equivalent to checking whether parent coarse grid cell is within coarse mask)
        if (INBASIN(TopoMap[coarsei][coarsej].Mask)){

	  // Solve for all grid cells within the coarse mask, not just the fine mask. 
	  neighbor_elev[n] = ((TopoMap[coarsei][coarsej].Mask) ? (*FineMap[yn][xn]).Dem : (float) OUTSIDEBASIN);

        }

      }

    }
    
    celev = (*FineMap[y][x]).Dem; 
    
    switch (NNEIGHBORS) { 
    case 8:
      lower = 0;
      for (n = 0; n < NNEIGHBORS; n++) {
	if(neighbor_elev[n] == OUTSIDEBASIN) {
	  neighbor_elev[n] = celev;
	}
	
	/* Calculating tanbeta as tanbeta * length of cell boundary between
	   the cell of interest and downsloping neighbor. */
	if(neighbor_elev[n] < celev){
	  if(n==0 || n==2 || n==4 || n==6){
	    temp_slope[n] = (celev - neighbor_elev[n]) /length_diagonal;
	    contour_length[y][x] += 0.4*Map->DMASS;
	    tanbeta[y][x] += temp_slope[n]*(0.4*Map->DMASS);
	    delta_a[n] = a[y][x] * temp_slope[n]*(0.4*Map->DMASS);
	  }
	  
	  else {
	    temp_slope[n] = (celev - neighbor_elev[n]) /Map->DMASS;
	    contour_length[y][x] += 0.6*Map->DMASS;
	    tanbeta[y][x] += temp_slope[n]*(0.6*Map->DMASS);
	    delta_a[n] = a[y][x] * temp_slope[n]*(0.6*Map->DMASS);

	  }
	}
	else
	  lower++;
      } /* end for (n = 0; n < NNEIGHBORS; n++) { */
      
      if (lower == 8){
	/* if this is a flat area then tanbeta = sum of (0.5 * vertical resolution of elevation 
	   data)/ honizontal distance between centers of neighboring grid cells */
	tanbeta[y][x] = ((NNEIGHBORS/2)*((0.5 * VERTRES)/length_diagonal)) +
	  ((NNEIGHBORS/2)*((0.5 * VERTRES)/(Map->DMASS)));
      }

      /* Distributing total upslope area to downslope neighboring cells */
      for (n = 0; n < NNEIGHBORS; n++) {
	if(neighbor_elev[n]<celev){  
	  switch (n) {
	  case 0:
	    a[y+1][x-1] += delta_a[n]/tanbeta[y][x];
	    break;
	  case 1:
	    a[y+1][x] += delta_a[n]/tanbeta[y][x];
	    break;
	  case 2:
	    a[y+1][x+1] += delta_a[n]/tanbeta[y][x];
	    break;
	  case 3:
	    a[y][x+1] += delta_a[n]/tanbeta[y][x];
	    break;
	  case 4:
	    a[y-1][x+1] += delta_a[n]/tanbeta[y][x];
	    break;
	  case 5:
	    a[y-1][x] += delta_a[n]/tanbeta[y][x];
	    break;
	  case 6:
	    a[y-1][x-1] += delta_a[n]/tanbeta[y][x];
	    break;
	  case 7:
	    a[y][x-1] += delta_a[n]/tanbeta[y][x];
	    break;
	  default:
	    ReportError("CalcTopoIndex", 65);
	    assert(0);
	  } /* end switch (n) {*/
	}/* end 	 if(neighbor_elev[n]<celev){ */	
      } /*  for (n = 0; n < NNEIGHBORS; n++) { */
      break; /*end case 8: */
		     
    case 4:
      ReportError("CalcTopoIndex", 65);
      assert(0);                       /* not set up to do this */
      break;
    default:
      ReportError("CalcTopoIndex", 65);
      assert(0);			/* other cases don't work either */
    } /* end  switch (NNEIGHBORS) {  */
  } /* end  for (k = 0; k < Map->NumCellsfine; k++) { */
  
  for (k = (Map->NumCellsfine)-1; k >-1; k--) { 
    y = OrderedCellsfine[k].y;
    x = OrderedCellsfine[k].x;
    
    (*FineMap[y][x]).TopoIndex = log(a[y][x]/tanbeta[y][x]);
  }
  
  /*************************************************************************/
  /* Create output files...currently hard-coded, chould be moved to dump   */
  /* functions for user specification. Creates the following file in the   */
  /* topoindex.asc - TopoIndex for the mass wasting resolution map        */
  /*************************************************************************/
   
  printmap = 0;

  if (printmap == 1){
    sprintf(topoindexmap, "topoindex.asc");
    
    if((fo=fopen(topoindexmap,"a")) == NULL)
      {
	printf("Cannot open TopoIndex output file.\n");
	exit(0);
      }
    
    fprintf(fo,"ncols %11d\n",Map->NXfine);
    fprintf(fo,"nrows %11d\n",Map->NYfine);
    fprintf(fo,"xllcorner %.1f\n",Map->Xorig);
    fprintf(fo,"yllcorner %.1f\n",Map->Yorig - Map->NY*Map->DY);
    fprintf(fo,"cellsize %.0f\n",Map->DMASS);
    fprintf(fo,"NODATA_value %d\n",0);
    
    for (y = 0; y < Map->NYfine; y++) {
      for (x  = 0; x < Map->NXfine; x++) {
	
	i = (int) floor(y/(Map->DY/Map->DMASS));
        j = (int) floor(x/(Map->DX/Map->DMASS));
	
        /* Check to make sure region is in the basin. */
        if (INBASIN(TopoMap[i][j].Mask)) 	
	   	fprintf(fo, "%2.3f ", (*FineMap[y][x]).TopoIndex);
	  /*   fprintf(fo, "%2.3f ", log(a[y][x])); */ 
	/*   fprintf(fo, "%2.3f ", log(1/tanbeta[y][x])); */
	else 
	  fprintf(fo, "0. "); 
	
      }
      fprintf(fo, "\n");
    }
    fclose(fo);
  }

  for(i=0; i<Map->NYfine; i++) { 
    free(a[i]);
    free(tanbeta[i]);
    free(contour_length[i]);
  }
  free(a);
  free(tanbeta);
  free(contour_length);
  free(OrderedCellsfine);
  
  return;
}
