/*
 * SUMMARY:      Desorption.c - Calculate soil desorption
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the maximum amount of moisture the soil can
 *               deliver to the atmosphere in one time step
 * DESCRIP-END.
 * FUNCTIONS:    Desorption()
 * COMMENTS:
 * $Id: Desorption.c,v 1.4 2003/07/01 21:26:13 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  Desorption()
*****************************************************************************/
float Desorption(int Dt, float MoistContent, float Porosity, float Ks,
		 float Press, float m)
{
  float Sorptivity;		/* sorptivity */
  float DesorptionVolume;	/* total desorption volume during timestep */

  /* Eq. 46, Wigmosta et al [1994] */

  if (MoistContent > Porosity)
    MoistContent = Porosity;

/*   Sorptivity =  */
/*     pow((double) ((8 * Porosity * Ks * Press)/(3.0*(1 + 3*m) * (1 + 4*m))), */
/* 	(double) 0.5) *  */
/* 	  pow((double) (MoistContent/Porosity), (double) (1.0/(2.0 * m) + 2)); */

  Sorptivity = sqrt((double)
		    ((8 * Porosity * Ks * Press) /
		     (3.0 * (1 + 3 * m) * (1 + 4 * m)))) *
    pow((double) (MoistContent / Porosity), (double) (1.0 / (2.0 * m) + 2));

  /* Eq. 45, Wigmosta et al [1994] */

/*  DesorptionVolume = Sorptivity * pow((double) Dt * SECPHOUR, (double) 0.5); */
  DesorptionVolume = Sorptivity * sqrt((double) Dt);

  return DesorptionVolume;
}
