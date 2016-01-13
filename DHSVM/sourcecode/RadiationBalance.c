/*
 * SUMMARY:      RadiationBalance.c - Calculate radiation balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate radiation balance at each pixel
 * DESCRIP-END.
 * FUNCTIONS:    RadiationBalance()
 *               LongwaveBalance()
 *               ShortwaveBalance()
 * COMMENTS:
 * $Id: RadiationBalance.c,v 1.4 2003/07/01 21:26:21 olivier Exp $     
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
  Function name: RadiationBalance()

  Purpose      : Calculate the radiation balance for the individual canopy 
                 layers

  Required     :
    int HeatFluxOption - TRUE if surface temperature is being calculated
    float Rs           - Incoming shortwave radiation (W/m2)
    float Ld           - Incoming longwave radiation (W/m2)
    float Tair         - Ambient air temperature (C)
    float Tcanopy      - Canopy temperature from the previous timestep (C)
    float Tsoil        - soil surface temperature from the previous timestep
                         (C)
    int NAct           - Number of vegetation layers above the snow
    VEGTABLE VType     - Information about number of veg layers
    SNOWPIX LocalSnow  - Information about snow conditions at current pixel
    PIXRAD *LocalRad   - Components of radiation balance for current pixel
    
  Returns      : void
  
  Modifies     :
    PIXRAD *LocalRad  - Components of radiation balance for current pixel

  Comments     :
    This routine is implemented according to Wigmosta et al. (1994), with a
    small change for the soil surface temperature if a sensible heat flux is
    calculated, and the canopy surface temperature.  
    The following assumptions are made:
      - No snow:
        soil temperature is TSurf from last timestep if HeatFluxOption ==
        TRUE, otherwise the soil temperature is the same as the air
        temperature
      - Snow:
        soil temperature is the snow surface temperature 
      - Canopy temperature is equal to the air temperature is there is no snow
        interception
      - There are at most two vegetation layers

  Reference:
    Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
    hydrology-vegetation model for complex terrain, Water Resour. Res.,
    30(6), 1665-1679, 1994.

  Reference:
    Nijssen and Lettenmaier, A simplified approach for predicting shortwave
    radiation transfer through boreal forest canopies, JGR, 1999.

*****************************************************************************/
void RadiationBalance(int HeatFluxOption, int CanopyRadAttOption, 
		      float SineSolarAltitude, float Rs,
		      float Rsd, float Rsb, float Ld, float Tair,
		      float Tcanopy, float Tsoil, float SoilAlbedo,
		      VEGTABLE *VType, SNOWPIX *LocalSnow, PIXRAD *LocalRad)
{
  float F;			/* Fraction of pixel covered by top canopy
				   layer [0-1] */
  float Albedo[2];		/* Albedo of each layer */
  float Tau;			/* Transmittance for overstory vegetation 
				   layer */
  float Taub, Taud;		/* Transmittance for overstory vegetation layer for
				   direct and diffuse radiation, respectively */
  float Tsurf;			/* Surface temperature (C) */

  F = VType->Fract[0];
  /*following added 08/13/2001 by Pascal Storck */
  if (CanopyRadAttOption == VARIABLE) {
    F = VType->HemiFract[0];
  }
  /*end of add 08/13/2001 */

  /* Determine Albedo */

  if (VType->OverStory == TRUE) {
    Albedo[0] = VType->Albedo[0];
    if (LocalSnow->HasSnow == TRUE)
      Albedo[1] = LocalSnow->Albedo;
    else if (VType->UnderStory == TRUE)
      Albedo[1] = VType->Albedo[1];
    else
      Albedo[1] = SoilAlbedo;
  }
  else if (LocalSnow->HasSnow == TRUE)
    Albedo[0] = LocalSnow->Albedo;
  else if (VType->UnderStory == TRUE)
    Albedo[0] = VType->Albedo[0];
  else
    Albedo[0] = SoilAlbedo;

  /* if the attenuation is fixed, calculate the canopy transmittance */

  if (CanopyRadAttOption == FIXED) {
    if (VType->OverStory == TRUE)
      Tau = exp(-VType->Atten*VType->LAI[0]);
    else
      Tau = 0.;
  }

  /* calculate canopy transmittance coefficient for overstory vegetation 
     layer */
  /* for the case where Bart Nijssen's simplified radiation scheme is used
     then k*LAI is assumed to be the effective Leaf Area Index (L in Nijssen
     and Lettenmaier, 2000) */

  else if (CanopyRadAttOption == VARIABLE) {
    if (VType->OverStory == TRUE) {
      Taub =
	exp(-VType->LAI[0] / VType->ClumpingFactor *
	    (VType->LeafAngleA / SineSolarAltitude + VType->LeafAngleB));
      /*formulation is typically based on the cos of the solar zenith angle,
	which is the sin of the solar altitude (SA = 90 - SZA) */
      Taud = VType->Taud;
      if (Rs > 0.0) {
	Tau = Taub * Rsb / Rs + Taud * Rsd / Rs;
	Tau = (float) pow((double) Tau, (double)(VType->Scat));
	/* VType->Scat can be specified as a scattering paramter 
	   or DHSVM will set it to 0.8 if not specified */
	Tau = Tau / (1 - Albedo[0] * Albedo[1]);
      }
      else
	Tau = 0.;
    }
    else
      Tau = 0.;
  }

  ShortwaveBalance(VType->OverStory, F, Rs, Tau, Albedo, LocalRad);

  if (LocalSnow->HasSnow == TRUE)
    Tsurf = LocalSnow->TSurf;
  else if (HeatFluxOption == TRUE)
    Tsurf = Tsoil;
  else
    Tsurf = Tair;

  LongwaveBalance(VType->OverStory, F, Ld, Tcanopy, Tsurf, LocalRad);

}

/*****************************************************************************
  Function name: LongwaveBalance() 
  
  Purpose      : Calculate the longwave radiation balance for the individual 
                 canopy layers

  Required     :
    unsigned char OverStory
                     - Flag to indicate presence/absence of overstory
    float F          - Fraction of pixel covered by top vegetation layer
    float Ld         - Downward (incoming) longwave radiation (W/m2)
    float Tcanopy    - Canopy temperature (C)
    float Tsurf      - Surface temperature (C)
    PIXRAD &LocalRad - Components of radiation balance for current pixel

  Returns      :
    void

  Modifies     :
    Elements of LocalRad

  Comments     : This function is used to update the longwave radiation
                 balance when new surface temperatures are calculated
*****************************************************************************/
void LongwaveBalance(unsigned char OverStory, float F, float Ld,
		     float Tcanopy, float Tsurf, PIXRAD * LocalRad)
{
  double Tmp;			/* temporary variable */

  /* calculate emitted longwave for each layer */

  if (OverStory == TRUE) {
    Tmp = Tcanopy + 273.15;
    LocalRad->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp);
    Tmp = Tsurf + 273.15;
    LocalRad->LongOut[1] = STEFAN * (Tmp * Tmp * Tmp * Tmp);
  }
  else {
    Tmp = Tsurf + 273.15;
    LocalRad->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp);
    LocalRad->LongOut[1] = 0.;
  }

  /* Calculate the incoming longwave for each layer */

  if (OverStory == TRUE) {
    LocalRad->LongIn[0] = (Ld + LocalRad->LongOut[1]) * F;
    LocalRad->LongIn[1] = Ld * (1 - F) + LocalRad->LongOut[0] * F;
  }
  else {
    LocalRad->LongIn[0] = Ld;
    LocalRad->LongIn[1] = 0.;
  }

  /* Calculate the radiative components for the entire pixel.  Use the
     snow/soil surface temperature as an estimate for the pixel temperature.
     LocalRad->PixelLongOut is calculated in the sensible heat flux routine,
     and is not needed otherwise.  Here it is initialized anyway, as if the
     surface temperature is already known */

  LocalRad->PixelLongIn = Ld;

  if (OverStory == TRUE) {
    LocalRad->PixelLongOut = LocalRad->LongOut[0] * F +
      LocalRad->LongOut[1] * (1 - F);
  }
  else {
    LocalRad->PixelLongOut = LocalRad->LongOut[0];
  }
}

/*****************************************************************************
  Function name: ShortwaveBalance()

  Purpose      : Calculate the shortwave radiation balance for the individual
                 pixels

  Required     :
    unsigned char OverStory
                     - Flag to indicate presence of overstory
    float F          - Fraction of pixel covered by top vegetation layer
    float Rs         - Incoming shortwave radiation (W/m2)
    float Tau        - Canopy transmittance coefficient for overstory
    float *Albedo    - Albedo of each layer
    PIXRAD *LocalRad - Components of radiation balance for current pixel

  Returns      : 
    void

  Modifies     :
    Elements of LocalRad

  Comments     : 
    This function needs to be called only once for each pixel for each 
    timestep, since the shortwave radiation balance is independent of
    the surface temperatures
*****************************************************************************/
void ShortwaveBalance(unsigned char OverStory, float F, float Rs, float Tau,
		      float *Albedo, PIXRAD * LocalRad)
{
  /* Calculate the net shortwave for each layer */
  /* Overstory present, i.e. two layers */
  if (OverStory == TRUE) {
    LocalRad->NetShort[0] = Rs * F * ((1 - Albedo[0]) - Tau * (1 - Albedo[1]));
    LocalRad->NetShort[1] = Rs * (1 - Albedo[1]) * ((1 - F) + (Tau * F));
  }
  else {
    LocalRad->NetShort[0] = Rs * (1 - Albedo[0]);
    LocalRad->NetShort[1] = 0.;
  }

  /* Calculate the net shortwave for the entire pixel */

  if (OverStory == TRUE) {
    LocalRad->PixelNetShort = Rs * (1 - Albedo[0] * F - Albedo[1] * (1 - F));
  }
  else {
    LocalRad->PixelNetShort = LocalRad->NetShort[0];
  }
}
