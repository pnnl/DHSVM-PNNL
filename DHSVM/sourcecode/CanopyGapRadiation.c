/*
* SUMMARY:      CanopyGapRadiation.c
* USAGE:        Part of DHSVM
*
* AUTHOR:       Ning Sun
* E-MAIL:       ning.sun@pnnl.gov
* ORIG-DATE:    Jul-16
* DESCRIPTION:  Calculate radiation balance under an idealized cylindrical canopy
*               gap/openning
* DESCRIP-END.
* FUNCTIONS:    CanopyGapRadiation()
* Reference:

C.R. Ellis, J.W. Pomeroy, and T.E. Link, Modeling increases in snowmelt
yield and desynchronization resulting from forest gap-thinning treatments
in a northern mountain headwater basin, Water Resour. Res., 49, 936-949, 2013.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "functions.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"

#define MAXIT 100

/********************************************************************************
Function Name: CanopyGapRadiation()

Purpose      : Calculate net radiation incident on canopy gap

Required     :
h        - canopy height
dm       - canopy opening diameter (assuming a cylindrical gap)
Rsb      - downward beam radiation (W/m^2)
Rsd      - downward diffuse radiation
Extn     - light extinction coeff.
SunAngle - sine of solar altitude

Returns      :
Comments     :
********************************************************************************/
void CanopyGapRadiation(CanopyGapStruct **Gap, float SunAngle, float Rs,
  float Rsb, float Rsd, float Ld, float TSurf, float Tcanopy, float SoilAlbedo,
  VEGTABLE *VType, SNOWPIX *LocalSnow, PIXRAD *LocalRad, float Gapping, VEGPIX *LocalVeg)
{

  /* ==== OPENING PORTION (no overstory if gap presents) ====*/
  (*Gap)[Opening].OverStory = FALSE;
  (*Gap)[Opening].UnderStory = VType->UnderStory;

  /* net shortwave received by the opening*/
  (*Gap)[Opening].NetShort[1] = CanopyGapShortRadiation((*Gap)[Opening].UnderStory,
    (*Gap)[Opening].GapView, VType->Height[0], Gapping, SunAngle, Rsb,
    Rsd, VType->ExtnCoeff, SoilAlbedo, VType, LocalSnow, LocalVeg->Fract[0]);
  (*Gap)[Opening].NetShort[0] = 0;

  /* net longwave received by the opening*/
  CanopyGapLongRadiation(Gap[Opening], VType->Height[0], Gapping, Ld, 
    TSurf, LocalVeg->Fract[0]);
  (*Gap)[Opening].LongIn[0] = 0.;
}

/********************************************************************************
Function Name: CanopyGapShortRadiation()
********************************************************************************/
float CanopyGapShortRadiation(int Understory, float GapView, float h, float dm, 
  float SunAngle, float Rsb, float Rsd, float Extn, float SoilAlbedo, 
  VEGTABLE *VType, SNOWPIX *LocalSnow, float Vf)
{

  float Dmax;       /* max length of shadow */
  float Lmax;       /* max attenuation length */
  float Albedo;	    /* Albedo of each layer */
  float R;          /* Gap radiaus */
  float I1 = 0.;    /* area receiving direct beam radiation*/
  float I2 = 0.;    /* area receiving incoming attenuated beam radiation */
  float Area;       /* gap area */
  float ls;         /* length along optical axis that defines shading area */
  float Rbg = 0.;   /* net beam radiation incident on canopy gap */
  float Rdg = 0.;   /* net diffuse radiation incident on canopy gap */
  float Rsg;         /* net shortwave radiation - beam + diffuse radiation */

  if (LocalSnow->HasSnow == TRUE)
    Albedo = LocalSnow->Albedo;
  else if (Understory == TRUE)
    Albedo = VType->Albedo[1]; //understory albedo
  else
    Albedo = SoilAlbedo;

  /* Gap radius */
  R = 0.5 * dm;
  Area = PI * R * R;

  /* Canopy view to sky is adjusted for canopy gap, and is estimated as a fraction
  of the overlying hemisphere opened to the base of the gap (Ellis et al., 2013) */
  Rdg = Rsd*(GapView + VType->Taud*(1-GapView));

  if (SunAngle > 0.) {
    /* Calculate max length of cast shadow in the canopy gap */
    Dmax = h / tan(SunAngle);
    /* Calcuate max attenuation distance to the canopy gap*/
    Lmax = h / sin(SunAngle);

    /* case A: Dmax >= Gap diameter -- attenuated light only */
    if (Dmax >= dm) {
      I2 = AreaIntegral(Extn, Lmax, SunAngle, R, R, 0, Rsb, Rdg, Albedo);
      Rbg = (2 * Rsb) / Area * I2;
    }
    /* case B: Dmax < Gap diameter -- both direct and attenuated light */
    else {
      ls = sqrt(dm*dm - Dmax*Dmax);
      I1 = Area - 0.5 * (dm*dm*asin(Dmax / dm) + Dmax*ls);
      /* debug - output gap-center shortwave radiation */
      //fprintf(stderr, "%f\n", (Rsb+Rdg));
      /* debug ends */
      I2 = 2 * (AreaIntegral(Extn, Lmax, SunAngle, R, R, 0.5*ls, Rsb, Rsd, Albedo) +
        (exp(-0.5*Lmax / Extn)) * Dmax * (ls*0.5));
      Rbg = Rsb / Area * (I1 + I2);
    }
  }
  else {
    Rbg = 0.;
    /* debug - output gap-center shortwave radiation */
    //fprintf(stderr, "%f\n", (Rsb + Rsd));
    /* debug ends */
  }

  /* net shortwave radiation received by the opening */
  Rsg = (Rdg + Rbg) * (1 - Albedo);

  return Rsg;
}

/********************************************************************************
Function Name: CanopyGapLongRadiation()
Purpose      : Calculate net long radiation incident on canopy gap
********************************************************************************/
void CanopyGapLongRadiation(CanopyGapStruct *Gap, float h, float dm, float Ld, 
  float TCanopy, float Vf)
{

  float Tmp;        /* surface temperature in K */
  float GapView;    /* canopy view to sky adjusted for canopy gap */
  
  GapView = Gap->GapView/Vf;

  if (GapView <= (1 - Vf)) {
    GapView = 1 - Vf;
  }
  if (GapView >= 1) {
    GapView = 0.99;
    //printf("Warning: GapView>=100% .. reset to 0.99 ..\n");
  }

  /* calculate emitted longwave for each layer */
  Tmp = TCanopy + 273.15;
  /* understory */
  Gap->LongOut[1] = STEFAN * (Tmp*Tmp*Tmp*Tmp);
  /* removed overstory */
  Gap->LongOut[0] = 0.;

  /* net longwave radiation */
  Gap->LongIn[1] = Ld*GapView + STEFAN*(Tmp*Tmp*Tmp*Tmp)*(1-GapView);
}

/********************************************************************************
Function Name: GapSurroundingLongRadiation()
********************************************************************************/
void GapSurroundingLongRadiation(CanopyGapStruct *Forest, float Ld, float Vf,
  float F, float Tcanopy, float Tsurf)
{

  double Tmp;

  Tmp = Tcanopy + 273.15;
  Forest->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp);
  Tmp = Tsurf + 273.15;
  Forest->LongOut[1] = STEFAN * (Tmp * Tmp * Tmp * Tmp);

  Forest->LongIn[0] = (Ld + Forest->LongOut[1]) * Vf;
  Forest->LongIn[1] = Ld * (1 - Vf) + Forest->LongOut[0] * Vf;

  Forest->PixelLongIn = Ld;
  Forest->PixelLongOut = Forest->LongOut[0] * F + Forest->LongOut[1] * (1 - F);
}

/********************************************************************************
Function Name: GapSurroundingShortRadiation()
********************************************************************************/
void GapSurroundingShortRadiation(CanopyGapStruct *Forest, VEGTABLE *VType, 
  SNOWPIX *LocalSnow, float SoilAlbedo, float SineSolarAltitude, float Rs, VEGPIX *LocalVeg)
{
  float F;
  float h;
  float Albedo[2];
  float Tau;

  F = LocalVeg->Fract[0];
  h = VType->Height[0];

  /* Determine Albedo */
  Albedo[0] = VType->Albedo[0];
  /* With snow, understory canopy albedo is set equal to snow albedo */
  if (Forest->HasSnow == TRUE)
    Albedo[1] = LocalSnow->Albedo;
  else if (VType->UnderStory == TRUE)
    Albedo[1] = VType->Albedo[1];
  else
    Albedo[1] = SoilAlbedo;

  /* Improved radiation scheme taking into account solar position */
  if (SineSolarAltitude > 0. && Rs > 0.) 
    Tau = exp(-VType->ExtnCoeff * h * F / SineSolarAltitude);
  else
    Tau = 0.;

  /* Calculate the net shortwave for each layer */
  Forest->NetShort[0] = Rs * (1 - Albedo[0]) * (1 - Tau * (1 - Albedo[1]));
  Forest->NetShort[1] = Rs * (1 - Albedo[1]) * Tau;
}

/********************************************************************************
Function Name: AreaIntegral()

Purpose      : Calculate the integral of an areal function

Required     :
               xmax - the upper bound of x, which is the distance between any point
                      along the optical axis to gap center
               xmin - the lower bound of x

Returns      :
Comments     :
********************************************************************************/
float AreaIntegral(float Extn, float Lmax, float SolarAltitude, float R,
  float xmax, float xmin, float Rsb, float Rdg, float Albedo)
{

  float x;      /* distance between any point along the optical axis to gap center */
  float u;      /* area defined by the function between two boundary limits*/
  float deltax; /* increment */
  float Rad;

  deltax = (xmax - xmin) / MAXIT;
 
  u = 0.;
  for (x = xmin; x <= xmax; x += deltax) {
    /* debug */
    Rad = (Rsb*exp(-Extn*(Lmax-sqrt(R*R-x*x)/cos(SolarAltitude))) + Rdg);
    /* if (xmin == 0 && x==0)
      fprintf(stderr, "%f \n", Rad); */
    /* debug ends*/
    
    u +=
      exp(-Extn*(Lmax - sqrt(R*R - x*x) / cos(SolarAltitude))) * sqrt(R*R - x*x) * deltax;
  }
  return u;
}


/********************************************************************************
Function Name: CalcGapView()
********************************************************************************/
float CalcGapView(float R, float H, float F) {

  float GapView;
  float SVF;
  float I;
  float delta_r;
  float r;
  float delta_alpha;
  float alpha;
  float junk;
  int iter = 20;

  delta_r = R / iter;
  delta_alpha = 2 * PI / iter;
  GapView = 0;

  for (r = 0; r <= R; r += delta_r) {
    SVF = 0;
    for (alpha = 0; alpha <= 2 * PI; alpha += delta_alpha) {
      I = sqrt(R*R - r*r*sin(alpha)*sin(alpha)) - r*cos(alpha);
      junk = atan2(I, H);
      SVF += 1 / (PI*PI) * atan2(I, H) * delta_alpha / F;
    }
    GapView += SVF * delta_r;
  }
  GapView /= R;

  return GapView;
}




