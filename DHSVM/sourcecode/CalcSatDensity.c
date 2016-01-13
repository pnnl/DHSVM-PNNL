/*
 * SUMMARY:      CalcSatDensity.c - Calculate saturated soil density
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       colleen@u.washington.edu
 * ORIG-DATE:    Sept-03
 * DESCRIPTION:  This function calculates the saturated soil density based
                 on the bulk density and particlc density.
 * DESCRIP-END.
 * FUNCTIONS:    CalcSatDensity) 
 * COMMENTS:
 * $Id: CalcSatDensity.c,v 1.1 2003/10/29 01:01:08 colleen Exp $     
 */

#include "settings.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  Function name: CalcSatDensity()

  Purpose      : This function calculates the daturated soil density. 

  Required     :
    float Density - Soil bulk density in kg/m3

  Returns      : SatDensity - Saturated Soil Density

  Modifies     : NA

  Comments     :

*****************************************************************************/
float CalcSatDensity(float Density)
{
  float SatDensity;		      /* Saturated Soil Density (kg/m3) */
  float Porosity;

  /* This could also be calculated using the porosity from the configuration file */
  Porosity = 1 - (Density/PARTDENSITY); 

  SatDensity = PARTDENSITY*(1-Porosity) + WATER_DENSITY*Porosity;

  return SatDensity;
}
