/*
 * SUMMARY:      InitSedMap() - Initialize sediment grid
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Nathalie Voisin
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Oct-11 2006
 * DESCRIPTION:   Initialize sediment grid when sediment options is on
 * DESCRIP-END.
 * FUNCTIONS:    InitSedMap()
 * COMMENTS:     
 */

#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "DHSVMerror.h"

/*****************************************************************************
  InitSedMap()
*****************************************************************************/
void InitSedMap(MAPSIZE *Map, SEDPIX *** SedMap )
{
  int   y;		/* Counters */
  
  if (!(*SedMap = (SEDPIX **) calloc(Map->NY, sizeof(SEDPIX *))))
     ReportError("InitSedMap", 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SedMap)[y] = (SEDPIX *) calloc(Map->NX, sizeof(SEDPIX))))
      ReportError("InitSedMap", 1);
    }
}

