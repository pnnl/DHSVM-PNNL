/*
 * SUMMARY:      CalcAerodynamic.c - Calculate the aerodynamic resistances
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu, pstorck@u.washington.edu
 * ORIG-DATE:    Thu Mar 27 18:00:10 1997
 * DESCRIPTION:  Calculate the aerodynamic resistances
 * DESCRIP-END.
 * FUNCTIONS:    CalcAerodynamic()
 * COMMENTS:
 * $Id: CalcAerodynamic.c,v 1.4 2003/07/01 21:26:09 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "DHSVMerror.h"
#include "settings.h"
#include "constants.h"
#include "functions.h"

/*****************************************************************************
  Function name: CalcAerodynamic()

  Purpose      : Calculate the aerodynamic resistance for each vegetation 
                 layer, and the wind 2m above the layer boundary.  In case of 
                 an overstory, also calculate the wind in the overstory.
                 The values are normalized based on a reference height wind 
                 speed, Uref, of 1 m/s.  To get wind speeds and aerodynamic 
                 resistances for other values of Uref, you need to multiply 
                 the here calculated wind speeds by Uref and divide the 
                 here calculated aerodynamic resistances by Uref
                 
  Required     :
    int NVegLayers - Number of vegetation layers
    char OverStory - flag for presence of overstory.  Only used if NVegLayers 
                     is equal to 1
    float Zref     - Reference height for windspeed
    float n        - Attenuation coefficient for wind in the overstory
    float *Height  - Height of the vegetation layers (top layer first)
    float Trunk    - Multiplier for Height[0] that indictaes the top of the 
                     trunk space
    float *U       - Vector of length 2, with wind for vegetation layers
                     If OverStory == TRUE the first value is the wind in
                     the overstory canopy, and the second value the wind 
                     2m above the lower boundary.  Otherwise the first 
                     value is the wind 2m above the lower boundary and 
                     the second value is not used.
    float *U2mSnow - Wind velocity 2m above the snow surface
    float *Ra      - Vector of length 2, with aerodynamic resistance values.  
                     If OverStory == TRUE the first value is the aerodynamic 
                     resistance for the the overstory canopy, and the second 
                     value the aerodynamic resistance for the lower boundary.  
                     Otherwise the first value is the aerodynamic resistance
                     for the lower boundary and the second value is not used.
    float *RaSnow  - Aerodynamic resistance for the snow surface.

  Returns      : void

  Modifies     :
    float *U
    float *U2mSnow
    float *Ra     
    float *RaSnow
   
  Comments     :
*****************************************************************************/
void CalcAerodynamic(int NVegLayers, unsigned char OverStory,
		     float n, float *Height, float Trunk, float *U,
		     float *U2mSnow, float *Ra, float *RaSnow)
{
  float d_Lower;
  float d_Upper;
  float K2;
  float Uh;
  float Ut;
  float Uw;
  float Z0_Lower;
  float Z0_Upper;
  float Zt;
  float Zw;

  K2 = VON_KARMAN * VON_KARMAN;

  /* No OverStory, thus maximum one soil layer */

  if (OverStory == FALSE) {

    if (NVegLayers == 0) {
      Z0_Lower = Z0_GROUND;
      d_Lower = 0;
    }
    else {
      Z0_Lower = Z0_MULTIPLIER * Height[0];
      d_Lower = D0_MULTIPLIER * Height[0];
    }

    /* No snow */
    U[0] = log((2. + Z0_Lower) / Z0_Lower) / log((Zref - d_Lower) / Z0_Lower);
    Ra[0] =
      log((2. + Z0_Lower) / Z0_Lower) * log((Zref - d_Lower) / Z0_Lower) / K2;

    /* Snow */
    *U2mSnow = log((2. + Z0_SNOW) / Z0_SNOW) / log(Zref / Z0_SNOW);
    *RaSnow = log((2. + Z0_SNOW) / Z0_SNOW) * log(Zref / Z0_SNOW) / K2;
  }

  /* Overstory present, one or two vegetation layers possible */
  else {
    Z0_Upper = Z0_MULTIPLIER * Height[0];
    d_Upper = D0_MULTIPLIER * Height[0];

    if (NVegLayers == 1) {
      Z0_Lower = Z0_GROUND;
      d_Lower = 0;
    }
    else {
      Z0_Lower = Z0_MULTIPLIER * Height[1];
      d_Lower = D0_MULTIPLIER * Height[1];
    }

    Zw = 1.5 * Height[0] - 0.5 * d_Upper;
    Zt = Trunk * Height[0];
    if (Zt < (Z0_Lower + d_Lower))
      ReportError("Trunk space height below \"center\" of lower boundary", 48);

    /* Resistance for overstory */
    Ra[0] = log((Zref - d_Upper) / Z0_Upper) / K2 *
      (Height[0] / (n * (Zw - d_Upper)) *
       (exp(n * (1 - (d_Upper + Z0_Upper) / Height[0])) - 1) + (Zw -
								Height[0]) /
       (Zw - d_Upper) + log((Zref - d_Upper) / (Zw - d_Upper)));

    /* Wind at different levels in the profile */
    Uw = log((Zw - d_Upper) / Z0_Upper) / log((Zref - d_Upper) / Z0_Upper);
    Uh =
      Uw - (1 -
	    (Height[0] - d_Upper) / (Zw - d_Upper)) / log((Zref -
							   d_Upper) / Z0_Upper);
    U[0] = Uh * exp(n * ((Z0_Upper + d_Upper) / Height[0] - 1.));
    Ut = Uh * exp(n * (Zt / Height[0] - 1.));

    /* resistance at the lower boundary */

    /* No snow */
    /* case 1: the wind profile to a height of 2m above the lower boundary is 
       entirely logarithmic */
    if (Zt > (2. + Z0_Lower + d_Lower)) {
      U[1] =
	Ut * log((2. + Z0_Lower) / Z0_Lower) / log((Zt - d_Lower) / Z0_Lower);
      Ra[1] =
	log((2. + Z0_Lower) / Z0_Lower) * log((Zt -
					       d_Lower) / Z0_Lower) / (K2 * Ut);
    }

    /* case 2: the wind profile to a height of 2m above the lower boundary 
       is part logarithmic and part exponential, but the top of the overstory 
       is more than 2 m above the lower boundary */
    else if (Height[0] > (2. + Z0_Lower + d_Lower)) {
      U[1] = Uh * exp(n * ((2. + Z0_Lower + d_Lower) / Height[0] - 1.));
      Ra[1] =
	log((Zt - d_Lower) / Z0_Lower) * log((Zt -
					      d_Lower) / Z0_Lower) / (K2 *
								      Ut) +
	Height[0] * log((Zref - d_Upper) / Z0_Upper) / (n * K2 *
							(Zw -
							 d_Upper)) * (exp(n *
									  (1 -
									   Zt
									   /
									   Height
									   [0]))
								      -
								      exp(n *
									  (1 -
									   (Z0_Lower
									    +
									    d_Lower
									    +
									    2.)
									   /
									   Height
									   [0])));
    }

    /* case 3: the top of the overstory is less than 2 m above the lower 
       boundary.  The wind profile above the lower boundary is part logarithmic
       and part exponential, but only extends to the top of the overstory */
    else {
      U[1] = Uh;
      Ra[1] =
	log((Zt - d_Lower) / Z0_Lower) * log((Zt -
					      d_Lower) / Z0_Lower) / (K2 *
								      Ut) +
	Height[0] * log((Zref - d_Upper) / Z0_Upper) / (n * K2 *
							(Zw -
							 d_Upper)) * (exp(n *
									  (1 -
									   Zt
									   /
									   Height
									   [0]))
								      - 1);
      fprintf(stderr,
	      "WARNING:  Top of overstory is less than 2 meters above the lower boundary\n");
    }

    /* Snow */
    /* case 1: the wind profile to a height of 2m above the lower boundary is 
       entirely logarithmic */
    if (Zt > (2. + Z0_SNOW)) {
      *U2mSnow = Ut * log((2. + Z0_SNOW) / Z0_SNOW) / log(Zt / Z0_SNOW);
      *RaSnow = log((2. + Z0_SNOW) / Z0_SNOW) * log(Zt / Z0_SNOW) / (K2 * Ut);
    }

    /* case 2: the wind profile to a height of 2m above the lower boundary 
       is part logarithmic and part exponential, but the top of the overstory 
       is more than 2 m above the lower boundary */
    else if (Height[0] > (2. + Z0_SNOW)) {
      *U2mSnow = Uh * exp(n * ((2. + Z0_SNOW) / Height[0] - 1.));
      *RaSnow = log(Zt / Z0_SNOW) * log(Zt / Z0_SNOW) /
	(K2 * Ut) +
	Height[0] * log((Zref - d_Upper) / Z0_Upper) / (n * K2 *
							(Zw -
							 d_Upper)) * (exp(n *
									  (1 -
									   Zt
									   /
									   Height
									   [0]))
								      -
								      exp(n *
									  (1 -
									   (Z0_SNOW
									    +
									    2.)
									   /
									   Height
									   [0])));
    }

    /* case 3: the top of the overstory is less than 2 m above the lower boundary.
       The wind profile above the lower boundary is part logarithmic and part 
       exponential, but only extends to the top of the overstory */
    else {
      *U2mSnow = Uh;
      *RaSnow = log(Zt / Z0_SNOW) * log(Zt / Z0_SNOW) /
	(K2 * Ut) +
	Height[0] * log((Zref - d_Upper) / Z0_Upper) / (n * K2 *
							(Zw -
							 d_Upper)) * (exp(n *
									  (1 -
									   Zt
									   /
									   Height
									   [0]))
								      - 1);
      fprintf(stderr,
	      "WARNING:  Top of overstory is less than 2 meters above the lower boundary\n");
    }
  }
}

/*****************************************************************************
  The following is a small main program that helps in testing the function
  CalcAerodynamic().

  To compile you also need the files ReportError.c and InitErrorMessage.c, and
  the appropriate header files.

  To compile the test program (using gcc):
  gcc -DTEST_MAIN -o test_aero ReportError.c InitErrorMessage.c CalcAerodynamic.c -lm

  then run the program by typing test_aero.
*****************************************************************************/
#ifdef TEST_MAIN

int main(void)
{
  float Ra[2] = { 0, 0 };
  int NVegLayers = 1;
  char OverStory = FALSE;
  float Zref = 80;
  float n = .5;
  float Height[2] = { .3, .2 };
  float Trunk = 0.3;
  float RaSnow;
  float UoverStory = 0.0;
  float U[2] = { 0, 0 };
  float U2mSnow = 0.0;

  InitErrorMessage();
  CalcAerodynamic(NVegLayers, OverStory, Zref, n, Height, Trunk, U, &U2mSnow,
		  Ra, &RaSnow);

  return EXIT_SUCCESS;
}

#endif
