/*
 * SUMMARY:      CalcEffectiveKh.c - Calculate effective thermal conducitivity
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the effective thermal conductivity of a soil under
 *               dry conditions
 * DESCRIP-END.
 * FUNCTIONS:    CalcEffectiveKh()
 * COMMENTS:
 * $Id: CalcEffectiveKh.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "constants.h"
#include "DHSVMerror.h"
#include "functions.h"

/*****************************************************************************
  CalcEffectiveKh()

  Source: Farouki, O. T., 1986, Thermal properties of soils, 
                          Trans Tech Publications

  This function calculates the effective thermal conductivity of a soil 
  based on the thermal conductivity under dry conditions, KhDry, and the
  thermal conductivity under saturated conditions, KhSat.  The latter
  differs for frozen and unfrozen soils, and is function of the effective
  solids thermal conductivity, KhSol, and the saturates soil moisture or 
  ice content.  

  The method followed here is Johansen's method, section 7.11 [Farouk, 1986]
 
  First the effective conductivity is calculated for each layer, after which
  the total effective thermal conductivity is calculated for the specified 
  depth.
*****************************************************************************/
float CalcEffectiveKh(int NSoilLayers, float Top, float Bottom,
		      float *SoilDepth, float *KhDry, float *KhSol,
		      float *Moisture, float *Porosity, float *TSoil)
{
  char NoEndLayer;		/* flag to indicate whether an end layer
				   has been determined */
  char NoStartLayer;		/* flag to indicate whether a start layer 
				   has been determined */
  float Dz;			/* Depth from soil surface (m) */
  float Ke;			/* Kersten number (see reference) */
  float KhEff;			/* Effective soil thermal conductivity 
				   (W/(m*K)) */
  float KhSat;			/* Thermal conductivity for saturated soils
				   (W/(m*K)) */
  float *LayerDepth;		/* Depth of each layer above the specified
				   depth (m) */
  float *LayerKh;		/* Effective thermal conductivity for each 
				   soil layer above the specified depth 
				   (W/(m*K)) */
  float Sr;			/* degree of saturation */
  float TotalDepth;		/* Depth of soil column for which to 
				   calculate the effective thermal 
				   conductivity (m) */
  int i;			/* counter */
  int NLayers;			/* Number of soil layers in depth interval */
  int StartLayer = 0;		/* First layer below top */

  TotalDepth = Bottom - Top;

  Dz = 0.0;
  NLayers = 0;
  NoStartLayer = TRUE;
  NoEndLayer = TRUE;
  LayerDepth = NULL;

  for (i = 0; i < NSoilLayers && NoEndLayer; i++) {
    Dz += SoilDepth[i];

    if (Dz > Top) {
      NLayers++;
      if (!(LayerDepth = (float *) realloc(LayerDepth, NLayers *
					   sizeof(float))))
	ReportError("CalcEffectiveKh()", 1);

      if (NoStartLayer) {
	StartLayer = i;
	LayerDepth[i - StartLayer] = Dz - Top;
	NoStartLayer = FALSE;
      }
      else if (Dz > Bottom) {
	LayerDepth[i - StartLayer] = SoilDepth[i] + Bottom - Dz;
	NoEndLayer = FALSE;
      }
      else
	LayerDepth[i - StartLayer] = SoilDepth[i];
    }
  }

  /* If the soil column is thinner than Bottom, than assign the soil 
     properties of the bottom layer to the remainder of the soil profile */

  if (NoEndLayer)
    LayerDepth[NLayers - 1] += Bottom - Dz;

  if (!(LayerKh = (float *) calloc(NLayers, sizeof(float))))
    ReportError("CalcEffectiveKh()", 1);

  for (i = 0; i < NLayers; i++) {

    Sr = Moisture[i + StartLayer] / Porosity[i + StartLayer];

    /* Assume for now that either all the water is either frozen or unfrozen */

    /* frozen soil */

    if (TSoil[i + StartLayer] < 0) {
      Ke = Sr;
      KhSat = pow((double) KhSol[i + StartLayer],
		  (double) (1 - Porosity[i + StartLayer])) *
	pow((double) 2.2, (double) Porosity[i + StartLayer]);
    }

    /* unfrozen soil */

    else {
      if (Sr > 0.1)
	Ke = log10((double) Sr) + 1.0;
      else
	Ke = 0.0;
      KhSat = pow((double) KhSol[i + StartLayer],
		  (double) (1 - Porosity[i + StartLayer])) *
	pow((double) KhH2O, (double) Porosity[i + StartLayer]);
    }
    LayerKh[i] = (KhSat - KhDry[i + StartLayer]) * Ke + KhDry[i + StartLayer];
  }

  /* Now place the soil layers in series, and calculate the resulting effective
     thermal conductivity for the entire soil */

  KhEff = 0.0;
  for (i = 0; i < NLayers; i++)
    KhEff += LayerDepth[i] / TotalDepth / LayerKh[i];
  KhEff = 1.0 / KhEff;

  /* clean up */
  free(LayerKh);
  free(LayerDepth);

  return KhEff;
}
