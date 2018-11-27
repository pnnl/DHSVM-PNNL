/*
 * SUMMARY:      MapSize.h
 * USAGE:        Defines the MAPSIZE structure defining a 2D region and 
 *               (maybe) routines to operate on it.
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    October 2018
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-10-18 08:31:11 d3g096
 * COMMENTS:
 */
#ifndef _MapSize_h_
#define _MapSize_h_

#include "settings.h"

typedef struct {
  float Rank;
  int   x;
  int   y;
} ITEM;

typedef struct {
  char System[BUFSIZE + 1];     /* Coordinate system */
  double Xorig;                 /* X coordinate of Northwest corner */
  double Yorig;                 /* Y coordinate of Northwest corner */
  /* int X; */                        /* Current x position */
  /* int Y; */                       /* Current y position */
  int NX;                       /* Number of (local) pixels in x direction */
  int NY;                       /* Number of pixels in y direction */
  float DX;                     /* Pixel spacing in x-direction */
  float DY;                     /* Pixel spacing in y-direction */
  float DXY;                    /* Pixel spacing in diagonal */
  int OffsetX;                  /* Offset in x-direction compared to basemap */
  int OffsetY;                  /* Offset in y-direction compared to basemap */
  int NumCells;                 /* Number of active cells on this processor */
  int AllCells;                 /* Number of cells within the basin */
  ITEM *OrderedCells;           /* Structure array to hold the ranked elevations; NumCells in size */
  int dist;                     /* GA handle that provides parallel distribution */
  int gNX;                      /* Number of (global) pixels in x direction */
  int gNY;                      /* Number of (global) pixels in x direction */
} MAPSIZE;

#endif
