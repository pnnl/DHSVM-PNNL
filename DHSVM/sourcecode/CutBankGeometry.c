/*
 * SUMMARY:      CutBankGeometry.c - Calculates effects of cut banks
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Batelle Pacific Northwest Laboratories
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Jul-96
 * DESCRIPTION:  This subroutine calculates corrections to adjust for the
 *               effects of road cut-banks and channels in grid cells.  It
 *               also add precip and updates the upper zone soil moisture.
 *               If the water table is below the road / channel bed precip is
 *               added to the coresponding zone. 
 * DESCRIP-END.
 * FUNCTIONS:    CutBankGeometry()
 * COMMENTS:
 * $Id: CutBankGeometry.c,v 1.4 2003/07/01 21:26:12 olivier Exp $     
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "settings.h"
#include "soilmoisture.h"

/*****************************************************************************
  Function name: CutBankGeometry()

  Purpose      : This subroutine calculates corrections to adjust for the
                 effects of road cut-banks and channels in grid cells.  It
                 also add precip and updates the upper zone soil
                 moisture.  If the water table is below the road /
                 channel bed precip is added to the coresponding zone.
	
  Required     :
    int i            - Number of the soil layer being processed.
    float RootDepth  - Depth of the soil layer being processed (m)
    float TopZone    - Distance from the ground surface to top of zone i (m) 
    float BankHeight - Distance from ground surface to channel bed or bottom
                       of road-cut (m) 
    float Area       - Area of channel or road surface (m)
    float DX         - Grid cell width (m)
    float DY         - Grid cell width (m)

  Returns      : void

  Modifies     :
    float *PercArea  - Area of the bottom of zone i, for perc (m/m)
    float *Adjust    - Corrects for loss of soil storage due to channel/
                       road-cut.  Multiplied with RootDepth to give the zone
                       thickness for use in calculating soil moisture
    int *CutBankZone - Number of the soil layer containing the bottom of the
                       cut-bank.  Used in UnsaturatedFlow to check for
                       surface runoff.  If BankHeight = 0.; CutBankZone =
                       NO_CUT 

  Comments     :
    Schematic: 
    
    <-----------------------------------DX---------------------------------->
    |====================|     -                     |======================|
    |      ^             |     |                     |   |                  |
    |      TopZone[0]    |<----|------Area---------->|   |- RootDepth[0]    |
    |                    |     |                     |   |                  |
    |<---PercArea*DX---->|     |                     |<--|--PercArea*DX---->|
    |                    |     |-BankHeight          |   V                  |
    |====================|     |                     |======================|
    |       ^            |     |                     |   |                  |
    |       TopZone[1]   |     |                     |   |-RootDepth[1]     |
    |                    |     V                     |   |                  |
    |                    |---------------------------|   |                  |
    |                                                    |                  |
    |                             CutBankZone            V                  |
    |=======================================================================|
    |                                                                       |
    |                                                                       |
    |=======================================================================|    

    Note: if the cut bank is cause by a road instead of a channel, then the 
          schematic would look as follows:

          |=================|
          |                 |			  
          |=================|			  
          |                 |_______________
          |                                 |
          |=================================|
          |                                 |
          |=================================|
*****************************************************************************/
void CutBankGeometry(int i, float RootDepth, float TopZone, float BankHeight,
		     float Area, float DX, float DY, float *PercArea, 
		     float *Adjust, int *CutBankZone)
{
  *PercArea = 1.0;
  *Adjust = 1.0;

  if (BankHeight > 0.0) {
    if (BankHeight <= TopZone) {
      /* below cut depth - full area */
      *PercArea = 1.0;
      *Adjust = 1.0;
    }
    else {
      if (BankHeight <= (TopZone + RootDepth)) {
	/* cut depth in this zone partial area */
	*PercArea = 1.0;
	*Adjust = 1.0 - (Area * (BankHeight - TopZone) / (RootDepth * DX * DY));
	*CutBankZone = i;
      }
      else {
	/* above cut depth - less than full area  */
	assert(DX*DY >= Area);
	*PercArea = 1 - Area / (DX * DY);
	*Adjust = *PercArea;
      }
    }
  }
}
