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
 * Reference:  

   Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed
   hydrology-vegetation model for complex terrain, Water Resour. Res.,
   30(6), 1665-1679, 1994.

   Nijssen and Lettenmaier, A simplified approach for predicting shortwave
   radiation transfer through boreal forest canopies, JGR, 1999.

   **references of improved radiation scheme **

   Reifsnyder, W. E., and H. W. Lull (1965), Radiant energy in relation to forests,
   USDA For. Serv. Tech. Bull., 1344, 111.

   Thyer M. et al. (2004), Diagnosing a distributed hydrologic model for two
   high-elevation forested catchments based on detailed stand- and basin-scale data.
   Water Resources Research. DOI:10.1029/2003WR002414.
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

  Purpose      : Calculate the radiation balance for the individual canopy layers

  Required     :
    int HeatFluxOption - TRUE if surface temperature is being calculated
    float Rs           - Incoming shortwave radiation adjusted by topo shading 
	                     if shading option is on (W/m2)
    float Rsb          - Incoming direct shortwave radiation separated from Rs 
    float Rsd          - Incoming diffuse shortwave radiation separated from Rs
	float VIC_Rs       - Input 'observed' shortwave radiation in the forcing file
    float Ld           - Incoming longwave radiation (W/m2)
    float Tair         - Ambient air temperature (C)
    float Tcanopy      - Canopy temperature from the previous timestep (C)
    float Tsoil        - soil surface temperature from the previous timestep (C)
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

*****************************************************************************/
void RadiationBalance(OPTIONSTRUCT *Options, int HeatFluxOption, 
              int CanopyRadAttOption, float SineSolarAltitude, 
              float VIC_Rs, float Rs, float Rsb, float Rsd, float Ld, 
              float Tair, float Tcanopy, float Tsoil, float SoilAlbedo,
		      VEGTABLE *VType, SNOWPIX *LocalSnow, PIXRAD *LocalRad)
{
  float F;			    /* Fraction of pixel covered by top canopy layer [0-1] */
  float h;              /* Canopy height (m) */
  float Albedo[2];		/* Albedo of each layer */
  float Tau;			/* Transmittance for overstory vegetation layer */
  float Taub, Taud;		/* Transmittance for overstory vegetation layer for
				           direct and diffuse radiation, respectively */
  float Tsurf;			/* Surface temperature (C) */

  F = VType->Fract[0];
  h = VType->Height[0];

  /*following added 08/13/2001 by Pascal Storck */
  if (CanopyRadAttOption == VARIABLE) {
    F = VType->HemiFract[0];
  }
  /*end of add 08/13/2001 */

  /* Determine Albedo */
  if (VType->OverStory == TRUE) {
    Albedo[0] = VType->Albedo[0];
    /* With snow, understory canopy albedo is set equal to snow albedo */
    if (LocalSnow->HasSnow == TRUE || LocalSnow->Iwq > 0.0)
      Albedo[1] = LocalSnow->Albedo;
    else if (VType->UnderStory == TRUE)
      Albedo[1] = VType->Albedo[1];
    else
      Albedo[1] = SoilAlbedo;
  }
  else if (LocalSnow->HasSnow == TRUE || LocalSnow->Iwq > 0.0)
    Albedo[0] = LocalSnow->Albedo;
  else if (VType->UnderStory == TRUE)
    Albedo[0] = VType->Albedo[0];
  else
    Albedo[0] = SoilAlbedo;
  
  /* calculate canopy transmittance coefficient for overstory layer */
  /* if the attenuation is fixed, calculate the canopy transmittance */
  if (CanopyRadAttOption == FIXED) {
    if (VType->OverStory == TRUE)
      Tau = exp(-VType->Atten * VType->LAI[0]);
    else
      Tau = 0.;
  }
  /* Nijssen's simplified radiation scheme as in Nijssen and Lettenmaier, 1999 */
  else if (CanopyRadAttOption == VARIABLE) {
    if (VType->OverStory == TRUE) {
      /* Calculate transmittance of overstory canopy for direct radiation:
         1) LAI * ClumpingFactor = Effective LAI 
         2) Formulation is typically based on the cos of the solar zenith angle,
         which is the sin of the solar altitude (SA = 90 - SZA) */
      Taub = exp(-VType->LAI[0] / VType->ClumpingFactor *
	    (VType->LeafAngleA / SineSolarAltitude + VType->LeafAngleB));
      
      /* transmittance for diffuse radiation (cacluated in CheckOut.c as a function of
         LeafAngleA and LeafAngleB and solar altitude) */
      Taud = VType->Taud;

      /* cacluate the total canopy transimittance for shortwave radiation (adjusted to 
         scattering and multiple reflection */
      if (Rs > 0.0) {
	    Tau = Taub * Rsb / Rs + Taud * Rsd / Rs;
        /* adjust Tau to scaterring parameter */
	    Tau = pow(Tau, (VType->Scat));
        /* adjust Tau to over- and under- story reflection */
	    Tau = Tau / (1 - Albedo[0] * Albedo[1]);
      }
      else
	    Tau = 0.;
    }
    else
      Tau = 0.;
  }

  /* A new version of radiation scheme -- by Ning Sun */
  if (Options->ImprovRadiation) {
    if (VType->OverStory == TRUE) {
      if (SineSolarAltitude > 0.)
        Tau = exp(-VType->ExtnCoeff * h * F / SineSolarAltitude);
      else
        Tau = 0.;
    }
    else
      Tau = 0.;
  }

  /* Calculate shortwave radiation balance */
  ShortwaveBalance(Options, VType->OverStory, F, Rs, Rsb, Rsd, Tau, Albedo, LocalRad);

  if (LocalSnow->HasSnow == TRUE || LocalSnow->Iwq > 0.0)
    Tsurf = LocalSnow->TSurf;
  else if (HeatFluxOption == TRUE)
    Tsurf = Tsoil;
  else
    Tsurf = Tair;
  
  /* Calculate longwave radiation balance */
  LongwaveBalance(Options, VType->OverStory, F, VType->Vf, Ld, Tcanopy, Tsurf, LocalRad);

  /* For John's RBM input */
  LocalRad->PixelLongIn = Ld;
  
  // copy the local Met variable value -- total shortwave radiation without any shading effect 
  LocalRad->ObsShortIn = VIC_Rs;

}

/************************************************************************************************
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
************************************************************************************************/
void LongwaveBalance(OPTIONSTRUCT *Options, unsigned char OverStory, 
            float F, float Vf, float Ld, float Tcanopy, float Tsurf, 
            PIXRAD *LocalRad)
{
  double Tmp;			

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
    /* In the improved radiation scheme, F is replaced by the canopy view 
    factor Vf that indicates the proportion of the sky view blocked by the 
    canopy. Here Vf is estimated by VfAdjust * F. IN the case of clear-cut, 
    F = Vf = 0. Details can be found in Thyer et al. (2004) and 
    Reifsnyder and Lull (1965) */
    if (Options->ImprovRadiation == TRUE) {
      LocalRad->LongIn[0] = (Ld + LocalRad->LongOut[1]) * Vf;
      LocalRad->LongIn[1] = Ld * (1 - Vf) + LocalRad->LongOut[0] * Vf;
    }
    else {
      LocalRad->LongIn[0] = (Ld + LocalRad->LongOut[1]) * F;
      LocalRad->LongIn[1] = Ld * (1 - F) + LocalRad->LongOut[0] * F;
    }
  }
  else {
    LocalRad->LongIn[0] = Ld;
    LocalRad->LongIn[1] = 0.;
  }

  /* Calculate the net longwave for the entire pixel */
  /* added by Ning */
  if (Options->StreamTemp) {
    if (OverStory == TRUE) {
      LocalRad->RBMNetLong = 
		Ld * (1 - F) + LocalRad->LongOut[0] * F;
    }
    else {
      LocalRad->RBMNetLong = Ld;
    }
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
void ShortwaveBalance(OPTIONSTRUCT *Options, unsigned char OverStory, 
              float F, float Rs, float Rsb, float Rsd, float Tau, 
              float *Albedo, PIXRAD *LocalRad)
{
  /* Calculate the net shortwave for each layer */
  /* Overstory present, i.e. two layers */
  if (OverStory == TRUE) {
    /* A new version of radiation scheme: 
    Beause F was factored in during the tau calculations, F is not used 
    repeatedly here */
    if (Options->ImprovRadiation == TRUE) {
      LocalRad->NetShort[0] = Rs * (1 - Albedo[0]) * (1 - Tau * (1 - Albedo[1]));
      LocalRad->NetShort[1] = Rs * (1 - Albedo[1]) * Tau;
    }
    else {
      LocalRad->NetShort[0] = Rs * F * ((1 - Albedo[0]) - Tau * (1 - Albedo[1]));
      LocalRad->NetShort[1] = Rs * (1 - Albedo[1]) * ((1 - F) + (Tau * F));
    }
  }
  else {
    LocalRad->NetShort[0] = Rs * (1 - Albedo[0]);
    LocalRad->NetShort[1] = 0.;
  }

  /* Calculate the net shortwave for the entire pixel */
  if (OverStory == TRUE) {
    if (Options->ImprovRadiation == TRUE) {
      LocalRad->PixelNetShort = Rs * (1 - Albedo[0] - Albedo[0] * (1-Albedo[1]));
    }
    else {
      LocalRad->PixelNetShort = Rs * (1 - Albedo[0] * F - Albedo[1] * (1 - F));
    }
  }
  else 
    LocalRad->PixelNetShort = LocalRad->NetShort[0];

  /* Calculate the incoming shortwave reaching the water surface */
  /* When the canopy shading option is off, the model still takes into account the shading
  created by the local vegetation defined by the input vegetation map */
  if (Options->StreamTemp && !Options->CanopyShading){
    if (OverStory == TRUE) {
      LocalRad->RBMNetShort = Rs * (1 - F) + Rs * Tau * F;
	  LocalRad->PixelBeam = Rsb * (1 - F) + Rsb * Tau * F;    // direct beam radiation 
      LocalRad->PixelDiffuse = LocalRad->RBMNetShort - LocalRad->PixelBeam; // diffuse radiation 
    }
    else {
      LocalRad->RBMNetShort = Rs;
	  LocalRad->PixelBeam = Rsb;
	  LocalRad->PixelDiffuse = Rsd;
    }
  }
  /* When turning on the canopy shading top, only riparian veg characterized by the height, width and 
  other parameters are used for calculations. The VEG type defined by the vegetation map is dismissed. */
  else if (Options->StreamTemp && Options->CanopyShading){ 
      LocalRad->RBMNetShort = Rs;
	  LocalRad->PixelBeam = Rsb;
	  LocalRad->PixelDiffuse = Rsd;
  }
}
