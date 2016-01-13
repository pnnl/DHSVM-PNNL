/*
 * SUMMARY:      MassRelease.c - Calculates mass release of snow from canopy
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Brian Connelly and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     6-Oct-1996 at 15:42:13
 * DESCRIPTION:  Calculates mass release of snow from canopy
 * DESCRIP-END.
 * FUNCTIONS:    MassRelease()
 * COMMENTS:
 * $Id: MassRelease.c,v 1.4 2003/07/01 21:26:21 olivier Exp $     
 */

#include <stdarg.h>
#include "constants.h"
#include "settings.h"
#include "massenergy.h"
#include "snow.h"

/*****************************************************************************
  Function name: MassRelease()

  Purpose      : Calculates mass release of snow from canopy

  Required     :
    float *InterceptedSnow
    float *TempInterceptionStorage
    float *ReleasedMass
    float *Drip 

  Returns      : none

  Modifies     : see under required (i.e. all variables are modified)

  Comments     :
*****************************************************************************/
void MassRelease(float *InterceptedSnow, float *TempInterceptionStorage,
		 float *ReleasedMass, float *Drip, float MDRatio)
{
  float TempDrip;
  float TempReleasedMass;

  /* If the amount of snow in the canopy is greater than some minimum
     value, MIN_INTERCEPTION_STORAGE, then calculte mass release and Drip */

  if (*InterceptedSnow > MIN_INTERCEPTION_STORAGE) {

    if ((*TempInterceptionStorage) >= 0.0) {
      *Drip += *TempInterceptionStorage;
      *InterceptedSnow -= *TempInterceptionStorage;

      if (*InterceptedSnow < MIN_INTERCEPTION_STORAGE)
	TempReleasedMass = 0.0;
      else
	TempReleasedMass = MIN((*InterceptedSnow - MIN_INTERCEPTION_STORAGE),
			       *TempInterceptionStorage * MDRatio);
      *ReleasedMass += TempReleasedMass;
      *InterceptedSnow -= TempReleasedMass;
      *TempInterceptionStorage = 0;

    }

    else {
      TempDrip = MIN(*TempInterceptionStorage, *InterceptedSnow);
      *Drip += TempDrip;
      *InterceptedSnow -= TempDrip;
    }
  }

  /* (*InterceptedSnow < MIN_INTERCEPTION_STORAGE) If the amount of snow in
     the canopy is less than some minimum value, MIN_INTERCEPTION_STORAGE,
     then only melt can occur and there is no mass release. */

  else {
    TempDrip = MIN(*TempInterceptionStorage, *InterceptedSnow);
    *Drip += TempDrip;
    *InterceptedSnow -= TempDrip;
    *TempInterceptionStorage = 0.0;
  }
}
