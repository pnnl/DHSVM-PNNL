/*
 * SUMMARY:      CalcKinViscosity - Calculate the kinematic viscosity
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Jordan Lanini
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    
 * DESCRIPTION:  estimate kinematic viscosity through interpolation
 *               based on dewpoint temperature
 * DESCRIP-END.
 * FUNCTIONS:    viscosity()
 * COMMENTS:
 * $Id: CalcKinViscosity.c,v 1.3 2004/08/18 01:01:26 colleen Exp $     
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "functions.h"

float viscosity(float Tair, float Rh)
{
  float knviscosity=0;              /* kinematic viscosity */
  float Tdew=0;                     /* Dew point temperature */
  float X;                          /* complement of relative humidity, fraction */
  
  /* calculate Dewpoint temperature Eq. 2-7 Linsley */
  X = 1.-Rh/100.;
  
  Tdew =- (14.55+.114*Tair)*X-pow(((2.5+0.007*Tair)*X),3)-
    (15.9+.117*Tair)*pow(X,14)+Tair;

  if (Tdew<0.)
    knviscosity=1.792;
  else if (Tdew<4. && Tdew>=0.)
    knviscosity=(1.792-(Tdew-0.)/4.*(1.792-1.567));
  else if  (Tdew>=4. && Tdew<10.)
    knviscosity=(1.567-(Tdew-4.)/6.*(1.567-1.371));
  else if (Tdew>=10. && Tdew<20.)
    knviscosity=(1.371-(Tdew-10.)/10.*(1.371-1.007));        
  else if  (Tdew>=20. && Tdew<25.)
    knviscosity=(1.007-(Tdew-20.)/5.*(1.007-.8963));
  else if  (Tdew>=25. && Tdew<30.)
    knviscosity=(.8963-(Tdew-25.)/5.*(.8963-.8042));
  else if  (Tdew>=30. && Tdew<40.)
    knviscosity=(.8042-(Tdew-30.)/10.*(.8042-.6611));
  else
    knviscosity=(.6611-(Tdew-40.)/10.*(.6611-.556));
  return knviscosity;
}

