/*
 * SUMMARY:      slopeaspect.h - header file for SlopeAspect.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for SlopeAspect.c
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: slopeaspect.h,v 1.10 2004/04/16 04:18:29 colleen Exp $     
 */

#ifndef SLOPEASPECT_H
#define SLOPEASPECT_H

#include "settings.h"
#include "data.h"

/* -------------------------------------------------------------
   available variables
   ------------------------------------------------------------- */
extern int xneighbor[NNEIGHBORS];
extern int yneighbor[NNEIGHBORS];
extern int xdirection[NDIRS];
extern int ydirection[NDIRS];

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap);
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap,
		     float **FlowGrad, unsigned char ***Dir, unsigned int **TotalDir);
int valid_cell(MAPSIZE * Map, int x, int y);
int valid_cell_fine(MAPSIZE *Map, int x, int y);
float ElevationSlope(MAPSIZE *Map, TOPOPIX ** TopoMap, FINEPIX ***FineMap, int y, int x, int *nexty, 
		     int *nextx, int prevy, int prevx, float *Aspect);
/* void ElevationSlopeAspectfine(MAPSIZE * Map, FINEPIX *** FineMap, TOPOPIX **TopoMap) ;*/
void quick(ITEM *OrderedCells, int count);
#endif

