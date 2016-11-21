/*
 * SUMMARY:      CalcKhDry.c - Calculate thermal conductivity
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  This function calculates the thermal conductivity of a soil 
 *               under dry conditions, KhDry, based on the soil density.
 * DESCRIP-END.
 * FUNCTIONS:    CalcKhDry() 
 * COMMENTS:
 * $Id: CalcKhDry.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include "settings.h"
#include "functions.h"

/*****************************************************************************
  Function name: CalcKhDry()

  Purpose      : This function calculates the thermal conductivity of a soil 
                 under dry conditions, KhDry, based on the soil density. This
                 is an empirical relationship which supposedly is accurate to
                 within 20% [Farouki, 1986; section 3.3.2].

  Required     :
    float Density - Soil bulk density in kg/m3

  Returns      : KhDry - Soil thermal conductivity

  Modifies     : NA

  Comments     : Source: Farouki, O. T., 1986, Thermal properties of soils, 
                 Trans Tech Publications

*****************************************************************************/
float CalcKhDry(float Density)
{
  float KhDry;			/* Dry soil thermal conductivity (W/(m*K)) */

  KhDry = (0.135 * Density + 64.7) / (2700 - 0.947 * Density);

  return KhDry;
}
