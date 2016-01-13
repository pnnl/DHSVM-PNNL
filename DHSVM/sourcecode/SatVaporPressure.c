/*
 * SUMMARY:      SatVaporPressure.c - Calculate saturated vapor pressure
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:     4-Oct-1996 at 09:34:43
 * DESCRIPTION:  Calculates the saturated vapor pressure in Pa for a certain
 *               temperature 
 * DESCRIP-END.
 * FUNCTIONS:    SatVaporPressure()
 * COMMENTS:
 * $Id: SatVaporPressure.c,v 1.4 2003/07/01 21:26:23 olivier Exp $     
 */

#include <stdlib.h>
#include <math.h>
#include "lookuptable.h"

float CalcVaporPressure(float T);
static FLOATTABLE svp;		/* Table that contains saturated vapor 
				   pressures as a function of temperature 
				   in degrees C */

/*****************************************************************************
  Function name: InitSatVaporTable()

  Purpose      : Initialize lookup table for saturated vapor pressure as a 
                 function of temperature in degrees C.

  Required     : void

  Returns      : void

  Modifies     : none
  
  Comments     :  Table runs from -100 C to 100 C with an interval of 0.02 C
*****************************************************************************/
void InitSatVaporTable(void)
{
  InitFloatTable(10000L, -100., .02, CalcVaporPressure, &svp);
}

/*****************************************************************************
  Function name: CalcVaporPressure() 

  Purpose      : Calculates the saturated vapor pressure in Pa for a certain
                 temperature in degrees C 

  Required     : 
    float T    - Temperature (C)

  Returns      :
    float      - Saturated vapor pressure (Pa) 

  Modifies     : none
  
  Comments     :
    References: Shuttleworth, W.J., Evaporation,  In: Maidment, D. R. (ed.),
                  Handbook of hydrology,  1993, McGraw-Hill, New York, etc..
                Bras, R. A., Hydrology, an introduction to hydrologic
                  science, Addisson Wesley, Inc., Reading, etc., 1990.

*****************************************************************************/
float CalcVaporPressure(float T)
{
  float Pressure;

  Pressure = 610.78 * exp((double) ((17.269 * T) / (237.3 + T)));

  /* Calculate the saturated vapor pressure in the snow pack, 
     (Equation 3.32, Bras 1990) */

  if (T < 0.0)
    Pressure *= 1.0 + .00972 * T + .000042 * T * T;

  return Pressure;
}

/*****************************************************************************
  Function name: SatVaporPressure() - new version, using a lookup table

  Purpose      : Looks up the saturated vapor pressure in Pa for a certain
                 temperature in a table

  Required     : 
    float T    - Temperature (C)

  Returns      :
    float      - Saturated vapor pressure (Pa) 

  Modifies     : none
  
  Comments     : Uses lookup table
*****************************************************************************/
float SatVaporPressure(float T)
{
  return FloatLookup(T, &svp);
}

/*****************************************************************************
  Function name: SatVaporPressure() - old version, pre lookup table

  Purpose      : Calculates the saturated vapor pressure in Pa for a certain
                 temperature 

  Required     : 
    float T    - Temperature (C)

  Returns      :
    float      - Saturated vapor pressure (Pa) 

  Modifies     : none
  
  Comments     : No longer used
*****************************************************************************/
/* float SatVaporPressure(float T) */
/* { */
/*   return CalcVaporPressure(T); */
/* } */
