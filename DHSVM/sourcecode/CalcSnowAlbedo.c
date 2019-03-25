/*
 * SUMMARY:      CalcSnowAlbedo.c - Calculate snow albedo
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the snow albedo as a function of snow age
 * DESCRIP-END.
 * FUNCTIONS:    CalcSnowAlbedo()
 * COMMENTS:
 * $Id: CalcSnowAlbedo.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "Calendar.h"
#include "functions.h"

/*****************************************************************************
  CalcSnowAlbedo()
*****************************************************************************/
float CalcSnowAlbedo(float TSurf, unsigned short Last, SNOWTABLE * SnowAlbedo)
{
  if (Last > (unsigned short) DAYPYEAR)
    Last = (unsigned short) DAYPYEAR;

  /* Accumulation season */
  if (TSurf < 0.0)
    return SnowAlbedo[Last].Freeze;

  /* Melt season */
  else
    return SnowAlbedo[Last].Thaw;
}
