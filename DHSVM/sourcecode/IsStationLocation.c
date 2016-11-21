/*
 * SUMMARY:      IsStationLocation.c - Check for station location
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Check whether the current pixel is a station location
 * DESCRIP-END.
 * FUNCTIONS:    IsStationLocation()
 * COMMENTS:
 * $Id: IsStationLocation.c,v 1.4 2003/07/01 21:26:19 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"

#define NOTSAME -1

/*****************************************************************************
function  : IsStationLocation()
input     : Current Location, address of structure with station data, and
            pointer to int indicating which station is at current location 
	    (if any)
output    : TRUE if the current location coincides with a station
            FALSE otherwise
modifies  : integer indicating which station is at current location
programmer: Bart Nijssen

This routine determines whether there is a station at the current location
*****************************************************************************/
uchar IsStationLocation(COORD * Loc, int NStats, METLOCATION * Station,
			int *WhichStation)
{
  int i;			/* Station Counter */

  for (i = 0, *WhichStation = NOTSAME; i < NStats &&
       *WhichStation == NOTSAME; i++) {
    if (Loc->N == Station[i].Loc.N && Loc->E == Station[i].Loc.E) {
      *WhichStation = i;
      return TRUE;
    }
  }

  return FALSE;
}
