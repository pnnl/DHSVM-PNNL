/*
 * SUMMARY:      CalcBagnold.c - Calculate the total sediment transport capacity
 * USAGE:        Called by MainMWM.c
 *
 * AUTHOR:       Ed Maurer
 * ORG:          University of Washington, Department of Civil Engineering
 * DESCRIPTION:  Calculate the total sediment transport capacity
 * DESCRIP-END.
 * FUNCTIONS:    CalcBagnold()
 * COMMENTS:
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "DHSVMerror.h"
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMChannel.h"

#define DEPTHTHRESHOLD 0.0001   /* Min. depth below which no transport occurs
				   To avoid divide by zero */

/*****************************************************************************
  Function name: CalcBagnold()

  Purpose      : Calculate the total sediment transport capacity
 
                 
  Required     :


  Returns      : sediment transport capacity in kg (dry mass) per second

  Modifies     : none
   
  Comments     : Most equations: Hydralics of Sediment Transport, Graf(1971).
                 Analytical approximations made here for Figs 9.3 and 9.4.
                 This uses the simplification for fully turbulent conditions.
                 A check of Reynold's number should eventually be added.
*****************************************************************************/
float CalcBagnold(float DS,TIMESTRUCT * Time, float outflow, float width, float n, float slope)
{

  float Q, V, flowdepth;
  float visc, settling;
  float streampower;
  float tau0, taustar, tanalpha,tanalphamax, eb;
  float A,B;
  float TotalLoad;
 
  /* settling velocity uses Rubey's formula -- result in m/s */
  visc=VISCOSITY/1000000.0; /* convert mm2/s to m2/s to use SI units*/
  /* note -- this differs from the solution used in RouteSurface */
  settling = sqrt(36*visc*visc/(DS*DS)+0.667*(PARTDENSITY-WATER_DENSITY)*G*DS/WATER_DENSITY)-6*visc/DS;

  /* flow depth by Manning's equation; flow velocity */
  Q = outflow/Time->Dt;

  flowdepth = pow(Q*n/(width*sqrt(slope)),(0.6));
  if(flowdepth < DEPTHTHRESHOLD) TotalLoad = 0;
  else {
    V = Q/(flowdepth*width);

    /*printf ( "DS= %.5f   Vs= %.5f  Dt= %d  n= %.3f flow= %.3f",DS,settling,Time->Dt ,n,outflow);
      printf(" depth= %.3f V= %.3f\n",flowdepth,V); */

    /* streampower per area in J/s/m2 see eq. 9.10 in Graf (1971)*/
    streampower= WATER_DENSITY*G*flowdepth*V*slope;
    
    /* average shear stress tau0 and dimensionless taustar */
    tau0 = WATER_DENSITY*G*flowdepth*slope;
    taustar = tau0/(DS*(PARTDENSITY-WATER_DENSITY)*G);
    
    /* Now use approximations for Figs 9.3 and 9.4 for eb and tanalpha */
    A = -0.00125-0.0132*DS/MMTOM;
    B = 0.147-0.0132*DS/MMTOM;
    eb = A*log10(V*3.28)+B; /* original chart had V in ft/s */

    if(DS/MMTOM <= 0.6) {
      A = 0.142-0.71*DS/MMTOM;
      B = 0.808+0.11*DS/MMTOM;
      tanalphamax = 0.75;
      tanalpha = A*log10(taustar) + B;
      tanalpha = ( tanalpha > tanalphamax ) ? tanalphamax : tanalpha;
    }
    else if (DS/MMTOM > 0.6 && DS/MMTOM <= 2.0 ) {
      A = -0.46+0.23*DS/MMTOM;
      B = 1.12 - 0.44*DS/MMTOM;
      tanalphamax = (0.85-0.29*DS/MMTOM > 0.75) ? 0.75 : 0.85-0.29*DS/MMTOM;
      tanalpha = A*log10(taustar) + B;
      tanalpha = ( tanalpha > tanalphamax ) ? tanalphamax : tanalpha;
    }
    else {
      tanalpha = 0.375;
    }
    if(tanalpha < 0.375) tanalpha = 0.375;
    
    /* printf("tanalpha= %.4f eb= %.4f\n",tanalpha,eb);*/
    /* Calculate the total load (rate) as immersed weight per unit width */
    TotalLoad = streampower*((eb/tanalpha)+0.01*V/settling);
    /* Convert to dry mass (kg) per unit width per second */
    TotalLoad /= ((1-WATER_DENSITY/PARTDENSITY)*G);
    /* Convert to dry mass transport rate kg/s */
    TotalLoad *= width;
    if(TotalLoad<0.0) TotalLoad=0.0;
  }
  return TotalLoad;
}



