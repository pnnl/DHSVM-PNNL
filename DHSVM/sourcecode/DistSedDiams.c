/*
 * SUMMARY:      DistSedDiams.c - Distribute the sediment diameters
 * USAGE:        MWM
 *
 * AUTHOR:       Ed Maurer
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Sep-02
 * Last Change:  Thu Jun 19 09:27:02 2003 by Ed Maurer <edm@u.washington.edu>
 * DESCRIPTION:  For lateral sediment inflow, finds the particle diameters for each 
                  portion.  Assumes a lognormal distribution
 * DESCRIP-END.   
 * FUNCTIONS:    DistSedDiams()
 * COMMENTS:
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data.h"
#include "channel.h"
#include "constants.h"

/*****************************************************************************
  Function name: DistSedDiams()

  Purpose      : Calculate sediment diameters for NSEDSIZES bins
                 
  Required     :

  Returns      : Sediment diameters (mm)

  Modifies     : none
   
  Comments     : to be consistent with FindValue, use the Tukey (1960) 
                 approx to normal CDF -- y is probability, function returns Z(y)
*****************************************************************************/

void DistributeSedimentDiams(float SedDiams[NSEDSIZES])
{
  
#ifndef NORMALDIST
#define NORMALDIST(mean, stdev, y) (4.91 * stdev * (pow(y,.14) - pow(( 1 - y ),.14)) + mean )
#endif

  int i;
  float mn,std,z;
  float pctfiner;

  /* slope of the lognormal curve, log(d) on y-axis, Z on x-axis, is stdev */
  mn = log10(DEBRISd50);
  std = log10(DEBRISd90)-log10(DEBRISd50)/(NORMALDIST(0,1,0.9)-NORMALDIST(0,1,0.5));

  pctfiner = 1.0/(2.0*NSEDSIZES); /* midpoint of finest interval */

  for(i=0;i<NSEDSIZES;i++) {
    z = NORMALDIST(mn,std,pctfiner);
    SedDiams[i] = pow(10,mn+std*z);
    pctfiner += 1.0/NSEDSIZES;
  }
}

