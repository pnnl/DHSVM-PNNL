/*
 * SUMMARY:      AggregateRadiation.c - calculate basin-wide radiation
 * USAGE:        part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate basin-wide radiation
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: AggregateRadiation.c,v 1.4 2003/07/01 21:26:09 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "massenergy.h"

/*****************************************************************************
  AggregateRadiation()
  
  In the current implementation the local radiation elements are not stored 
  for the entire area.  Therefore these components are aggregated here.  They
  are averaged over the basin in Aggregate()
*****************************************************************************/
void AggregateRadiation(int MaxVegLayers, int NVegL, PIXRAD * Rad,
			PIXRAD * TotalRad)
{
  int i;			/* counter */

  /* aggregate radiation data */
  for (i = 0; i < NVegL; i++) {
    TotalRad->NetShort[i] += Rad->NetShort[i];
    TotalRad->LongIn[i] += Rad->LongIn[i];
    TotalRad->LongOut[i] += Rad->LongOut[i];
  }
  TotalRad->NetShort[MaxVegLayers] += Rad->NetShort[NVegL];
  TotalRad->LongIn[MaxVegLayers] += Rad->LongIn[NVegL];
  TotalRad->LongOut[MaxVegLayers] += Rad->LongOut[NVegL];
  TotalRad->PixelNetShort += Rad->PixelNetShort;
  TotalRad->PixelLongIn += Rad->PixelLongIn;
  TotalRad->PixelLongOut += Rad->PixelLongOut;

}
