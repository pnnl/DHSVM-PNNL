/*
 * SUMMARY:      Round.c - Round a double to the nearest int
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Round a double to the nearest int
 * DESCRIP-END.
 * FUNCTIONS:    Round()
 * COMMENTS:
 * $Id: Round.c,v 1.4 2003/07/01 21:26:23 olivier Exp $     
 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include "DHSVMerror.h"

/*****************************************************************************
  Round()
*****************************************************************************/
int Round(double x)
{
  int RoundX;

  if (x < INT_MIN || x > INT_MAX) {
    printf("%f\n", x);
    ReportError("Round()", 16);
  }
  if ((fabs(x) - floor(fabs(x))) < 0.5)
    RoundX = (int) ((x < 0) ? ceil(x) : floor(x));
  else
    RoundX = (int) ((x < 0) ? floor(x) : ceil(x));

  return RoundX;
}
