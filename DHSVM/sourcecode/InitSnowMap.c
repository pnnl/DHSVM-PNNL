/*
 * SUMMARY:      InitSnowMap.c - Initialize snow coverage
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize snow coverage
 * DESCRIP-END.
 * FUNCTIONS:    InitSnowMap()
 * COMMENTS:
 * $Id: InitSnowMap.c,v 1.4 2003/07/01 21:26:17 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  Function name: InitSnowMap()

  Purpose      : Initialize the snow information for each pixel in the basin

  Required     :
    MAPSIZE Map        - Size and location of the model area
    SNOWPIX ***SnowMap - Address of array with snow information 

  Returns      : void

  Modifies     :
    Values stored at the locations pointed to by SnowMap

  Comments     :
*****************************************************************************/
void InitSnowMap(MAPSIZE * Map, SNOWPIX *** SnowMap)
{
  const char *Routine = "InitSnowMap";
  int y;			/* counter */

  printf("Initializing snow map\n");

  if (!(*SnowMap = (SNOWPIX **) calloc(Map->NY, sizeof(SNOWPIX *))))
    ReportError((char *) Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*SnowMap)[y] = (SNOWPIX *) calloc(Map->NX, sizeof(SNOWPIX))))
      ReportError((char *) Routine, 1);
  }
}
