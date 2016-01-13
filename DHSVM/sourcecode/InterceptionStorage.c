/*
 * SUMMARY:      InterceptionStorage.c - Calculate interception storage
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate interception storage
 * DESCRIP-END.
 * FUNCTIONS:    InterceptionStorage()
 * COMMENTS:
 * $Id: InterceptionStorage.c,v 1.5 2003/11/12 20:01:51 colleen Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  InterceptionStorage()
*****************************************************************************/
void InterceptionStorage(int NMax, int NAct, float *MaxInt, float *Fract,
			 float *Int, float *Precip, float *MomentSq, float *Height, 
			 unsigned char Understory, float Dt, float MS_Rainfall,
			 float LD_FallVelocity)
{
  float Available;		/* Available storage */
  float Intercepted;		/* Amount of water intercepted during this 
				   timestep */
  int i;			/* counter */
  float OriginalPrecip;

  OriginalPrecip = *Precip;

  
  /* The precipitation is multiplied by the fractional coverage, since if the 
     vegetation covers only 10% of the grid cell, only 10% can be intercepted as a 
     maximum */
  for (i = 0; i < NAct; i++) {
    Available = MaxInt[i] - Int[i];
    if (Available > *Precip * Fract[i])
      Intercepted = (*Precip) * Fract[i];
    else
      Intercepted = Available;
    *Precip -= Intercepted;
    Int[i] += Intercepted;
  }

  /* Find momentum squared of rainfall for use by the sediment model. */
  if(Understory) 
    /* Since the understory is assumed to cover the entire grid cell, all 
       momentum is associated with leaf drip, eq. 2, Wicks and Bathurst (1996) */
    *MomentSq = pow(LD_FallVelocity * WATER_DENSITY, 2) * PI/6 *
      pow(LEAF_DRIP_DIA, 3) * (*Precip)/Dt;
  else
    /* If no understory, part of the rainfall reaches the ground as direct throughfall. */
    *MomentSq = (pow(LD_FallVelocity * WATER_DENSITY, 2) * PI/6 *
		 pow(LEAF_DRIP_DIA, 3) * (*Precip)/Dt) + (1-Fract[0]) *
      MS_Rainfall;  
 
  /* WORK IN PROGESS */
  /* It should be checked whether the following statement can cause a "loss"
     of water.  When there is interception storage at timestep t, and snow
     will cover this vegetation layer at t = T+1, the amount of water in 
     storage will be lost */
/*   for (i = NAct; i < NMax; i++) */
/*     Int[i] = 0; */
}
