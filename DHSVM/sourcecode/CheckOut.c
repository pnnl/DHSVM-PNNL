/*
 * SUMMARY:      CheckOut.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    May-2000
 * DESCRIPTION:  Check stuff out for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:   CheckOut()
 * COMMENTS:
 * $Id: CheckOut.c,v3.1.2 2013/11/13 Ning Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "DHSVMerror.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "constants.h"

void CheckOut(int CanopyRadAttOption, LAYER Veg, LAYER Soil, 
	      VEGTABLE *VType, SOILTABLE *SType, MAPSIZE *Map, 
	      TOPOPIX **TopoMap, VEGPIX **VegMap, SOILPIX **SoilMap)
{

  int y, x, i, j;
  int *count = NULL, *scount = NULL;
  float a, b, l, Taud, Taub20, Taub40, Taub60, Taub80;

  int npixels;

  if (!(count = calloc(Veg.NTypes, sizeof(int)))) {
    ReportError("Checkout", 1);
  }
  if (!(scount = calloc(Soil.NTypes, sizeof(int)))) {
    ReportError("Checkout", 1);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	if (VegMap[y][x].Veg < 1 || VegMap[y][x].Veg > Veg.NTypes) {
	  printf("veg value %d out of range \n", VegMap[y][x].Veg);
	  exit(-1);
	}

	count[VegMap[y][x].Veg - 1]++;

	if (SoilMap[y][x].Soil < 1 || SoilMap[y][x].Soil > Soil.NTypes) {
	  printf("soil value %d out of range \n", SoilMap[y][x].Soil);
	  exit(-1);
	}
	scount[SoilMap[y][x].Soil - 1]++;

      }
    }
  }

  i = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	i = i + 1;
      }
    }
  }
  printf("\nBasin has %d active pixels \n", i);
  npixels = i;

  printf("\nThe following VEG types are in the current basin \n");

  for (i = 0; i < Veg.NTypes; i++) {
    if (count[i] > 0)
      printf
	("Class # %d of Type: %s has fraction basin area: %5.3f\n",
	 i + 1, VType[i].Desc, (float) count[i] / (float) npixels);
    VType[i].TotalDepth = 0.0;
    for (y = 0; y < VType[i].NSoilLayers; y++) {
      VType[i].TotalDepth += VType[i].RootDepth[y];
    }
  }

  printf("\nThe following SOIL types are in the current basin \n");
  for (i = 0; i < Soil.NTypes; i++)
    if (scount[i] > 0)
      printf
	("Class # %d of Type: %s has fraction basin area: %5.3f\n",
	 i + 1, SType[i].Desc, (float) scount[i] / (float) npixels);

  printf("\nSome estimates for current vegetation specification\n");
  for (i = 0; i < Veg.NTypes; i++) {
    if (count[i] > 0) {

      printf("\nVegetation Type: %s\n", VType[i].Desc);
      printf("2meter    wind speed fraction of ref level %1.3f\n", 
	     VType[i].USnow);
      if (VType[i].OverStory) {
	for (j = 0; j < 12; j++) {
	  if (fequal(VType[i].LAIMonthly[0][j], 0.0)) {
	    printf("Overstory LAI must be > 0\n");
	    exit(-1);
	  }
	}
	/* printf("Overstory LAI July %2.3f Effective LAI July %2.3f\n", VType[i].LAIMonthly[0][6]);*/
	if (CanopyRadAttOption == VARIABLE) {
	  a = VType[i].LeafAngleA;
	  b = VType[i].LeafAngleB;
	  l = VType[i].LAIMonthly[0][6] / VType[i].ClumpingFactor;
	  if (l == 0)
	    Taud = 1.0;
	  else
	    Taud =
	      exp(-b * l) * ((1 - a * l) * exp(-a * l) +
			     (a * l) * (a * l) * evalexpint(1, a * l));
	  Taub20 = exp(-l * (VType[i].LeafAngleA / 
			     0.342 + VType[i].LeafAngleB));
	  Taub40 = exp(-l * (VType[i].LeafAngleA / 
			     0.642 + VType[i].LeafAngleB));
	  Taub60 = exp(-l * (VType[i].LeafAngleA / 
			     0.866 + VType[i].LeafAngleB));
	  Taub80 = exp(-l * (VType[i].LeafAngleA / 
			     0.984 + VType[i].LeafAngleB));
	  printf("Solar Altitude 20 deg Tbeam %f Tdiff %f\n", Taub20, Taud);
	  printf("Solar Altitude 40 deg Tbeam %f Tdiff %f\n", Taub40, Taud);
	  printf("Solar Altitude 60 deg Tbeam %f Tdiff %f\n", Taub60, Taud);
	  printf("Solar Altitude 80 deg Tbeam %f Tdiff %f\n", Taub80, Taud);
	}
      }
    }
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	if (SoilMap[y][x].Depth <= VType[VegMap[y][x].Veg - 1].TotalDepth) {
	  printf("Error for class %d of Type %s  \n", VegMap[y][x].Veg,
		 VType[VegMap[y][x].Veg - 1].Desc);
	  printf("%d %d Soil depth is %f, Root depth is %f \n", y,x,SoilMap[y][x].Depth,
		 VType[VegMap[y][x].Veg - 1].TotalDepth);
	  exit(-1);
	}
      }
    }
  }
  if (count) {
    free(count);
  }
  if (scount) {
    free(scount);
  }
}
