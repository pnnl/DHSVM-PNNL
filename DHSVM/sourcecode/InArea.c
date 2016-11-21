/*
 * SUMMARY:      InArea.c - Check whether a pixel is in the area
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Check whether a pixel is in the area
 * DESCRIP-END.
 * FUNCTIONS:    InArea()
 * COMMENTS:
 * $Id $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "data.h"

/*****************************************************************************
  InArea()
*****************************************************************************/
uchar InArea(MAPSIZE * Map, COORD * Loc)
{
  if (Loc->N < 0 || Loc->N > (Map->NY - 1))
    return FALSE;
  else if (Loc->E < 0 || Loc->E > (Map->NX - 1))
    return FALSE;
  else
    return TRUE;
}
