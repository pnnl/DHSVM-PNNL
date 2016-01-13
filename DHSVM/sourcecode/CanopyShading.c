/*
 * SUMMARY:      CanopyShading.c - Calculate 
 * USAGE:        Part of DHSVM-RBM
 *               This function is used to output files for John's RBM model
 *
 * AUTHOR:       Ning Sun
 * ORG:          University of Washington, Department of Civil Engineering
 * ORIG-DATE:    June-17-2013
 * DESCRIPTION:  Calculate the shadow length resulting from riparian vegetation zone
 * DESCRIP-END.
 * FUNCTIONS:    
 * Modification 
 * $Id: CalcShadowLength.c, v3.1.2  2013/06/17 Ning Exp $    

 Reference:
    Chen et al., Stream temperature simulation of forested riparian
	area: I. Watershed-scale model development, Journal of Environmental 
	Engineering, 1998.

	Sridhar et al., Prediction of stream temperature in forested watersheds, 
	JAWRS, 2004. 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "Calendar.h"
#include "errorhandler.h"
#include "DHSVMerror.h"
#include "channel.h"
#include "channel_grid.h"
#include "functions.h"
#include "constants.h"
#include "tableio.h"
#include "settings.h"
  
/*****************************************************************************
  Function name: InitChannelRVeg()
  Purpose:	 This subroutine initiates the riparian extinction parameter for
             each channel segment
*****************************************************************************/
void InitChannelRVeg(TIMESTRUCT *Time, Channel *Head) 
{
	if (Head != NULL) { 
	  Channel *Current = NULL;
	  Current = Head;
	  while (Current) {
		Current->rveg.Extn = Current->rveg.ExtnCoeff[Time->Current.Month - 1];
		Current = Current->next;
	  }
	}
}
/*****************************************************************************
  Function name: CalcCanopyShading()
  Purpose:	 This subroutine calculates the shadow length
*****************************************************************************/
void CalcCanopyShading(TIMESTRUCT *Time, Channel *Channel, SOLARGEOMETRY *SolarGeo) 
{
	float Dx1 = 0., Dx2 = 0.;              /* shadow length (in meters) */
	float SolarAltitude = 0.;   /* in radians */
	float StreamAzimuth = 0.;   /* in radians */
	float HDEM = 0.;            /* average height of the canopy in reference to the 
				                elevation of stream surface */
	float Net_Shade_Fctr = 0.;  /* the effective shade density fro beam radiation */
	float SKOP = 0.;            /* sky openess ranging from 0 to 1 */
	int ShadeCase;
	float debug;

	/* compute solar altitude in radians */
	SolarAltitude = asin(SolarGeo->SineSolarAltitude);
	
	while (Channel) {
	  if (IsNewMonth(&(Time->Current), Time->Dt))
		 Channel->rveg.Extn = Channel->rveg.ExtnCoeff[Time->Current.Month - 1];

	  /* compute stream azimuth in radians */
	  StreamAzimuth = Channel->azimuth*PI/180.;

	  /* debug */
	  if (Channel->rveg.TREEHEIGHT < 0.|| Channel->rveg.BUFFERWIDTH < 0. ||
		  Channel->rveg.OvhCoeff < 0. || Channel->rveg.Extn < 0. ||
		  Channel->rveg.CanopyBankDist < 0. ) {
		ReportError("CalcCanopyShading()", 68);
	  }
	  /* debug ends */

	  /* compute the average height of the canopy */
	  HDEM = Channel->rveg.TREEHEIGHT;

	  /* if a vegetation polygon on the sunward bank is located very close 
	  to the stream, the overhanging canopy righ above the stream surface 
	  may exist and contribute shade.
	  The horizontal width of the protruding portion of the canopy toward
	  the stream is assumed to be a percentage of the tree height. */
	  Channel->rveg.BUFFERWIDTH += Channel->rveg.TREEHEIGHT * Channel->rveg.OvhCoeff;

	  /* compute the shadow length on the stream */
	  if (SolarAltitude > 0){
		/* Examine six (6) cases based on the shadow length, 
		canopy bank distance and buffer width */
	    Dx1 = HDEM * fabs(sin(SolarGeo->SolarAzimuth - StreamAzimuth)/tan(SolarAltitude))
			-(Channel->rveg.CanopyBankDist+Channel->rveg.StreamWidth);
		Dx2 = HDEM * fabs(sin(SolarGeo->SolarAzimuth - StreamAzimuth)/tan(SolarAltitude))
			-Channel->rveg.CanopyBankDist;
		
		/* Case 1 - No shade */
		if (Dx2 <= 0.0 || Channel->rveg.Extn == 0 || Channel->rveg.BUFFERWIDTH == 0)
          ShadeCase = 1;
	    /* Case 2 - Partial shade, sun above buffer */
	    else if (Dx1 <= 0.0 && Dx2 <= Channel->rveg.BUFFERWIDTH) 
          ShadeCase = 2;
	    /*Case 3 - Partial shade, sun below buffer */
	    else if (Dx1 <= 0.0 &&  Dx2 > Channel->rveg.BUFFERWIDTH) 
          ShadeCase = 3;
	    /* Case 4 - Full shade, sun above buffer */
	    else if (Dx1 > 0.0 && Dx2 <= Channel->rveg.BUFFERWIDTH) 
          ShadeCase = 4;
	    /* Case 5 - Full shade, sun partially below buffer */
	    else if (Dx1 > 0.0 && Dx1 <= Channel->rveg.BUFFERWIDTH && Dx2 > Channel->rveg.BUFFERWIDTH)
          ShadeCase = 5;
	    /* Case 6 - Full shade, sun entirely below buffer */
	    else if (Dx1 > Channel->rveg.BUFFERWIDTH && Dx2 > Channel->rveg.BUFFERWIDTH) 
          ShadeCase = 6;

		if (ShadeCase > 1) 
	      /* calculate the effective shade density */
		  Net_Shade_Fctr = CalcShadeDensity(ShadeCase, HDEM, Channel->rveg.StreamWidth,
					  SolarGeo->SolarAzimuth, StreamAzimuth, SolarAltitude, 
					  Channel->rveg.TREEHEIGHT, Channel->rveg.BUFFERWIDTH, Dx1, Dx2, Channel->rveg.Extn);
		else Net_Shade_Fctr = 0.;
		if (Net_Shade_Fctr > 1) {
		  printf("The shading density > 1! must be <=0\n");
		  exit(0);
		}

	    /* VEGSHD + OVHSHD is then divided by the stream surface width to 
	    approximate the fraction of stream surface covered by the composite shade. 
	    This ratio then is used to estimate the amount of incoming direct beam 
	    radiation that actually reaches the water surface. */
	    Channel->Beam *= (1 - Net_Shade_Fctr);
	    if (Channel->Beam < 0)
		  Channel->Beam = 0.;
	  }

	  /* compute shading effect on diffusive radiation */
	  if (HDEM > 0 && Channel->rveg.Extn != 0 && Channel->rveg.BUFFERWIDTH != 0) {
		SKOP = CalcCanopySkyView(HDEM, Channel->rveg.CanopyBankDist);
		Channel->Diffuse *= MIN(Channel->skyview, SKOP);
	    /* compute long-wave radiation */
	  }
	  else {
		SKOP = 1;
		Channel->Diffuse *= Channel->skyview;
	  }

	  /* compute the net shortwave raidation adjusted by canopy shading */
	  Channel->NSW = Channel->Diffuse + Channel->Beam;
	  debug = Channel->NLW;
	  Channel->NLW = Channel->NLW * MIN(Channel->skyview, SKOP) + 
		   0.96*(1-MIN(Channel->skyview, SKOP))*0.96*STEFAN*pow((double)(Channel->ATP+273.15),4);

	  Channel = Channel->next;
	}
}
/*****************************************************************************
  Function name: CalcShadeDensity()
  Purpose:	 This subroutine calculates the effective shade density
  Note: The current program assumes same vegetation type with same height
        are plants along the riparaian zone. 
        The program does not handle multiple vegeatation buffers with 
        varying vegetation height. If needed, the code can be easily modified to
		account for this scenario.
******************************************************************************/
float CalcShadeDensity(int ShadeCase, float HDEM, float WStream, 
					  float SunAzimuth, float StreamAzim, float SunAltitude,
					  float TREEHEIGTH, float BUFFERWIDTH, float Dx1, float Dx2, 
					  float Ext_Ceoff) 
{
	double Pavg;   // the average path length of a sunbeam through the buffer
	double SHDDEN;  // the effective shade density
	float Shaded, Not_Shaded;
	float Net_Shade_Fctr;
	
	/* Initialize the parameters */
	Shaded = WStream;
	Not_Shaded = 0;

	switch (ShadeCase) {
	  case 2:
		Not_Shaded = WStream -Dx2;
        Shaded = Dx2;
        Pavg = 0.5 * Dx2 /cos(SunAltitude)/fabs(sin(SunAzimuth - StreamAzim));
		break;
	  case 3:
		Not_Shaded = WStream - Dx2;
        Shaded = Dx2;
        Pavg = 0.5 * BUFFERWIDTH /cos(SunAltitude)/fabs(sin(SunAzimuth - StreamAzim));
		break;
	  case 4:
        Pavg = 0.5 * (Dx1 + Dx2) /cos(SunAltitude)/fabs(sin(SunAzimuth - StreamAzim));
		break;
	  case 5:
        Pavg = 0.5 * (Dx1 + BUFFERWIDTH) /cos(SunAltitude)/fabs(sin(SunAzimuth - StreamAzim));
		break;
	  case 6:
        Pavg = BUFFERWIDTH /cos(SunAltitude)/fabs(sin(SunAzimuth - StreamAzim));
		break;
	}

	/* compute the effective shade density */
	SHDDEN = 1 - exp((float)(-Ext_Ceoff * Pavg));
	Net_Shade_Fctr = SHDDEN * Shaded /WStream;

	return Net_Shade_Fctr;
}
/*****************************************************************************
  Function name: CalcCanopySkyView()
  Purpose:	 This subroutine calculates the skyview factor above the riparaian 
             vegetation
******************************************************************************/
float CalcCanopySkyView(float HDEM, float dist) 
{
	float VSA;  /* vegetation shade angle (degrees) */
	float SKOP; /* sky openess (ranges from 0 ~ 1) */

	/* compute the vegeation shade angle in degrees */
	VSA = atan2(HDEM, dist)*180./PI;

	/* assume the vegetation zone on both side of the bank have
	the same height and distance to the bank */
	SKOP = (180 - 2*VSA)/180;

	return SKOP;
}
