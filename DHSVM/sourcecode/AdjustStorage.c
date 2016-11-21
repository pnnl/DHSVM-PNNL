/*
 * SUMMARY:      AdjustStorage.c - Adjust moisture storage for cut banks
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Batelle Pacific Northwest Laboratories
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Jul-96
 * DESCRIPTION:  Calculates corrections to adjust for the effects of 
 *               road cut-banks and channels in grid cells by calling the 
 *               function geometry().
 * DESCRIP-END.
 * FUNCTIONS:    AdjustStorage()
 * COMMENTS:
 * $Id: AdjustStorage.c,v 1.4 2003/07/01 21:26:08 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "soilmoisture.h"

/*****************************************************************************
  Function name: AdjustStorage()

  Purpose      : This subroutine calculates corrections to adjust for the
                 effects of road cut-banks and channels in grid cells by
                 calling the function CutBankGeometry(). 

  Required     :
    int NSoilLayers  - Number of soil layers
    float TotalDepth - Total thickness of the soil column (m)
    float *RootDepth - Array with thicknesses of root layers (m)
    float Area       - Area of channel or road surface (m)
    float DX         - Grid cell width (m)
    float DY         - Grid cell width (m)
    float BankHeight - Distance from ground surface to channel bed or bottom
                       of road-cut (m) 

  Returns      : void

  Modifies     :
    float *PercArea  - Array with area of the bottom of zone i, for perc
                       (m/m). This is number between 0 and 1
    float *Adjust    - Array with coefficients to correct for loss of soil
                       storage due to channel/road-cut for each soil layer.
                       Multiplied with RootDepth to give the zone thickness
                       for use in calculating soil moisture 
    int *CutBankZone - Number of the soil layer containing the bottom of the
                       cut-bank. Used in UnsaturatedFlow to check for surface
                       runoff.  If BankHeight = 0.; CutBankZone = NO_CUT. 

  Comments     :
*****************************************************************************/
void AdjustStorage(int NSoilLayers, float TotalDepth, float *RootDepth,
		   float Area, float DX, float DY, float BankHeight,
		   float *PercArea, float *Adjust, int *CutBankZone)
{
  float DeepLayerDepth;		/* depth of layer below deepest root zone
				   layer */
  float Depth;			/* depth below the ground surface (m) */
  int i;			/* counter */

  DeepLayerDepth = TotalDepth;
  Depth = 0.0;

  for (i = 0; i < NSoilLayers; i++) {
    CutBankGeometry(i, RootDepth[i], Depth, BankHeight, Area, DX, DY,
		    &PercArea[i], &Adjust[i], CutBankZone);
    Depth += RootDepth[i];
  }

  DeepLayerDepth = TotalDepth - Depth;

  CutBankGeometry(i, DeepLayerDepth, Depth, BankHeight, Area, DX, DY,
		  &PercArea[i], &Adjust[i], CutBankZone);
}
