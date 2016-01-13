
/*
 * SUMMARY:      Draw.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    2000
 * DESCRIPTION:  X11 routines for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:   
 * COMMENTS:
 * $Id: Draw.c,v 1.12 2006/10/03 22:50:22 nathalie Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "snow.h"
#include "Calendar.h"

#ifdef HAVE_X11
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

extern Display *display;
extern Window window;
extern GC gc;
extern XColor my_color[50];
extern float **temp_array;
extern long black, white;
extern int e, ndx;
#endif

void draw(DATE * Day, int first, int DayStep, MAPSIZE *Map, int NGraphics,
	  int *which_graphics, VEGTABLE * VType, SOILTABLE * SType,
	  SNOWPIX ** SnowMap, SOILPIX ** SoilMap, SEDPIX ** SedMap, FINEPIX *** FineMap,
	  VEGPIX ** VegMap, TOPOPIX ** TopoMap, PRECIPPIX ** PrecipMap, float **PrismMap,
	  float **SkyViewMap, unsigned char ***ShadowMap, EVAPPIX ** EvapMap,
	  RADCLASSPIX ** RadMap, MET_MAP_PIX ** MetMap, ROADSTRUCT **Network,
	  OPTIONSTRUCT * Options)
{				/*begin */
  int i, j, k, ie, je, ir, jr;
  int ii, jj, yy, xx;
  int PX, PY;
  int MapNumber;
  float min, max, scale;
  float temp, surf_swe, pack_swe;
  int index, skip_it;
  char *text;
  char text2[6];
  char text3[20];
  int length;
  float max_temp;
  float re;
  int buf = 50;
  int sample = 1;		/* if sample =0, then average if compression needed, otherwise if =1 then sample */
  /*  obviously, for really large domains, the 1 option is much faster */
  int expand;
  int draw_static_colorbar;
#ifdef HAVE_X11
  XWindowAttributes windowattr;

  expand = e;
  draw_static_colorbar = 1;
  if (XGetWindowAttributes(display, window, &windowattr) == 0) {
    printf("failed to get window attributes in draw \n");
    exit(-1);
  }

  /* windowatt.map_state = 0 if DHSVM realtime display is set to an icon */
  /* windowatt.map_state = 2 if DHSVM realtime is active */
  /* if the user iconizes DHSVM display then */
  /* turn the graphics off and let DHSVM fly (or at least try to fly) */

  if (windowattr.map_state > 0) {

    XSetForeground(display, gc, black);
    SPrintDate(Day, text3);
    XClearArea(display, window, 10, 0, 100, 20, False);
    XDrawString(display, window, gc, 10, 20, text3, 19);

    /* if the user changes window size below 300 by 300 then */
    /* turn the graphics off and let DHSVM fly (or at least try to fly) */
    /* but at least print the date and time to the display */

    if (windowattr.width > 300 && windowattr.height > 300) {

      if (first == 1 || Day->Hour == 0)
	draw_static_colorbar = 1;

      for (k = 0; k < NGraphics; k++) {
	/* this is the beginning of the master loop which tries to draw */
	/* all the graphic variables */
	/* however we override the static fields, 3, 4, 5 and 6 such that */
	/* they are only drawn on the first call or on each new day */
	/* the same limitation is given to drawing all the color bars */
	MapNumber = which_graphics[k];
	if (MapNumber < 3 || MapNumber > 6 || draw_static_colorbar == 1) {
	  PY = k / ndx;
	  PX = k - ndx * PY;

	  if (expand > 0) {
	    PX = PX * (Map->NX * expand + buf) + 10;
	    PY = PY * (Map->NY * expand + buf) + 20;	/*top 20 pixels reserved for date stamp */
	  }
	  else {
	    PX = PX * (Map->NX * (1.0 / ((float) (-expand))) + buf) + 10;
	    PY = PY * (Map->NY * (1.0 / ((float) (-expand))) + buf) + 20;
	  }
	  max = -1000000.;
	  min = 1000000.;

	  if (MapNumber == 1) {
	    text = "SWE (mm)";
	    length = 8;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SnowMap[j][i].Swq * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 2) {
	    text = "Water Table Depth (mm)";
	    length = 22;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].TableDepth * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 3) {
	    text = "Digital Elevation Model (m)";
	    length = 27;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = TopoMap[j][i].Dem;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 4) {
	    text = "Vegetation Class";
	    length = 16;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = VegMap[j][i].Veg;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 5) {
	    text = "Soil Class";
	    length = 10;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].Soil;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 6) {
	    text = "Soil Depth (mm)";
	    length = 15;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].Depth * 1000.;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 7) {
	    text = "Precipitation (mm)";
	    length = 18;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = PrecipMap[j][i].Precip * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 8) {
	    text = "Incoming Shortwave (W/sqm)";
	    length = 26;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = RadMap[j][i].Beam + RadMap[j][i].Diffuse;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 9) {
	    text = "Intercepted Snow (mm)";
	    length = 21;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)
		    && VType[VegMap[j][i].Veg - 1].OverStory == 1) {
		  temp = PrecipMap[j][i].IntSnow[0] * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		if (VType[VegMap[j][i].Veg - 1].OverStory == 1 && temp > 0.0)
		  temp_array[j][i] = temp;
		else
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 10) {
	    text = "Snow Surface Temp (C)";
	    length = 21;
	    max = 0.0;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SnowMap[j][i].TSurf;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}

		temp_array[j][i] = temp;
		if (fequal(SnowMap[j][i].Swq, 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 11) {
	    text = "Cold Content (kJ)";
	    length = 17;
	    max = 0.0;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  if (SnowMap[j][i].Swq > MAX_SURFACE_SWE) {
		    pack_swe = SnowMap[j][i].Swq - MAX_SURFACE_SWE;
		    surf_swe = SnowMap[j][i].Swq - pack_swe;
		    temp =
		      2.10e3 * (SnowMap[j][i].TSurf * surf_swe +
				SnowMap[j][i].TPack * pack_swe);
		  }
		  else {
		    temp = 2.10e3 * SnowMap[j][i].Swq * SnowMap[j][i].TSurf;
		  }
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(SnowMap[j][i].Swq, 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 12) {
	    text = "Snow Melt (mm)";
	    length = 14;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SnowMap[j][i].Melt * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 13) {
	    text = "Snow Pack Outflow (mm)";
	    length = 22;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SnowMap[j][i].Outflow * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 14) {
	    text = "Sat. Subsurf Flow (mm) 0=white";
	    length = 30;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].SatFlow * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 15) {
	    text = "Overland Flow (mm)";
	    length = 18;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].Runoff * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 16) {
	    text = "Total EvapoTranspiration (mm)";
	    length = 29;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = EvapMap[j][i].ETot * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 17) {
	    text = "Snow Pack Vapor Flux (mm)";
	    length = 25;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SnowMap[j][i].VaporMassFlux * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 18) {
	    text = "Int Snow Vapor Flux (mm)";
	    length = 24;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SnowMap[j][i].CanopyVaporMassFlux * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 19) {
	    text = "Soil Moist L1 (% Sat)";
	    length = 21;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  /*  temp=(SoilMap[j][i].Moist[0]-SType[SoilMap[j][i].Soil-1].FCap[0])/
		     (SType[SoilMap[j][i].Soil-1].Porosity[0]-
		     SType[SoilMap[j][i].Soil-1].FCap[0])*100.0; */
		  temp =
		    SoilMap[j][i].Moist[0] / SType[SoilMap[j][i].Soil -
						   1].Porosity[0] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 20) {
	    text = "Soil Moist L2 (% Sat)";
	    length = 21;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  /*      temp=(SoilMap[j][i].Moist[1]-SType[SoilMap[j][i].Soil-1].FCap[1])/
		     (SType[SoilMap[j][i].Soil-1].Porosity[1]-
		     SType[SoilMap[j][i].Soil-1].FCap[1])*100.0; */
		  temp =
		    SoilMap[j][i].Moist[1] / SType[SoilMap[j][i].Soil -
						   1].Porosity[1] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 21) {
	    text = "Soil Moist L3 (% Sat)";
	    length = 21;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  /*  temp=(SoilMap[j][i].Moist[2]-SType[SoilMap[j][i].Soil-1].FCap[2])/
		     (SType[SoilMap[j][i].Soil-1].Porosity[2]-
		     SType[SoilMap[j][i].Soil-1].FCap[2])*100.0; */
		  temp = SoilMap[j][i].Moist[2] /
		    SType[SoilMap[j][i].Soil - 1].Porosity[2] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 22) {
	    text = "Accumulated Precip (mm)";
	    length = 23;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = MetMap[j][i].accum_precip * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 23) {
	    text = "Air Temp (C) 0=white";
	    length = 20;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = MetMap[j][i].air_temp;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (temp_array[j][i] > -0.5 && temp_array[j][i] < 0.0)
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 24) {
	    text = "Wind Speed (m/s)";
	    length = 16;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = MetMap[j][i].wind_speed;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 25) {
	    text = "RH";
	    length = 2;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = MetMap[j][i].humidity;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}

		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 26) {
	    text = "Prism Precip (mm)";
	    length = 17;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = PrismMap[j][i] / 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 27) {
	    text = "Deep Layer Storage (% Sat)";
	    length = 26;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {

		  temp =
		    SoilMap[j][i].Moist[3] / SType[SoilMap[j][i].Soil -
						   1].Porosity[2] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 28) {
	    text = "Surface runoff from HOF and Return Flow (mm)";
	    length = 22;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].IExcess * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 29 && Options->Infiltration == DYNAMIC) {
	    text = "Infiltration Accumulation (mm)";
	    length = 22;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].TableDepth * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 30 && Options->MassWaste) {
	    text = "Sediment to Channel (m)";
	    length = 23;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {

		  // FineMap quantities must be aggregated to coarse grid
		  temp = 0.0;
		  for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
		    for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		      yy = (int) j*Map->DY/Map->DMASS + ii;
		      xx = (int) i*Map->DX/Map->DMASS + jj;
		      temp += (*FineMap[yy][xx]).SedimentToChannel;
		    }
		  }
		  // Normalize by # FineMap cells in a pixel
		  temp /= Map->DMASS*Map->DMASS;

		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 31) {
	    text = "Overstory Trans (mm)";
	    length = 20;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = EvapMap[j][i].EAct[0] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 32) {
	    text = "Understory Trans (mm)";
	    length = 21;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = EvapMap[j][i].EAct[1] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 33) {
	    text = "Soil Evaporation (mm)";
	    length = 21;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = EvapMap[j][i].EvapSoil * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 34) {
	    text = "Overstory Int Evap (mm)";
	    length = 23;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = EvapMap[j][i].EInt[0] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 35) {
	    text = "Understory Int Evap (mm)";
	    length = 24;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = EvapMap[j][i].EInt[1] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 36 && Options->MassWaste) {
	    text = "Fine Map Elevation (m)";
	    length = 20;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {

		  // FineMap quantities must be aggregated to coarse grid
		  temp = 0.0;
		  for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
		    for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		      yy = (int) j*Map->DY/Map->DMASS + ii;
		      xx = (int) i*Map->DX/Map->DMASS + jj;
		      temp += (*FineMap[yy][xx]).Dem;
		    }
		  }
		  // Normalize by # FineMap cells in a pixel
		  temp /= Map->DMASS*Map->DMASS;

		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  /* if (MapNumber == 37 && Options->Sediment) { */
/* 	    text = "Slope"; */
/* 	    length = 5; */
/* 	    for (i = 0; i < Map->NX; i++) { */
/* 	      for (j = 0; j < Map->NY; j++) { */

/* 		if (INBASIN(TopoMap[j][i].Mask)) { */

		  // FineMap quantities must be aggregated to coarse grid
/* 		  temp = 0.0; */
/* 		  for (ii=0; ii< Map->DY/Map->DMASS; ii++) { */
/* 		    for (jj=0; jj< Map->DX/Map->DMASS; jj++) { */
/* 		      yy = (int) j*Map->DY/Map->DMASS + ii; */
/* 		      xx = (int) i*Map->DX/Map->DMASS + jj; */
/* 		      temp += (*FineMap[yy][xx]).Slope; */
/* 		    } */
/* 		  } */
		  // Normalize by # FineMap cells in a pixel
/* 		  temp /= Map->DMASS*Map->DMASS; */

/* 		  if (temp > max) */
/* 		    max = temp; */
/* 		  if (temp < min) */
/* 		    min = temp; */
/* 		} */
/* 		temp_array[j][i] = temp; */
/* 	      } */
/* 	    } */
/* 	  } */

	  if (MapNumber == 37 && Options->MassWaste) {
	    text = "Water Table Thickness (m)";
	    length = 25;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {

		  // FineMap quantities must be aggregated to coarse grid
		  temp = 0.0;
		  for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
		    for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		      yy = (int) j*Map->DY/Map->DMASS + ii;
		      xx = (int) i*Map->DX/Map->DMASS + jj;
		      temp += (*FineMap[yy][xx]).SatThickness;
		    }
		  }
		  // Normalize by # FineMap cells in a pixel
		  temp /= Map->DMASS*Map->DMASS;

		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 38 && Options->MassWaste) {
	    text = "Change in Sediment Depth (m)";
	    length = 28;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {

		  // FineMap quantities must be aggregated to coarse grid
		  temp = 0.0;
		  for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
		    for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		      yy = (int) j*Map->DY/Map->DMASS + ii;
		      xx = (int) i*Map->DX/Map->DMASS + jj;
		      temp += (*FineMap[yy][xx]).DeltaDepth;
		    }
		  }
		  // Normalize by # FineMap cells in a pixel
		  temp /= Map->DMASS*Map->DMASS;

		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 39 && Options->MassWaste) {
	    text = "Failure Probability";
	    length = 19;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {

		  // FineMap quantities must be aggregated to coarse grid
		  temp = 0.0;
		  for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
		    for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
		      yy = (int) j*Map->DY/Map->DMASS + ii;
		      xx = (int) i*Map->DX/Map->DMASS + jj;
		      temp += (*FineMap[yy][xx]).Probability;
		    }
		  }
		  // Normalize by # FineMap cells in a pixel
		  temp /= Map->DMASS*Map->DMASS;

		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 41) {
	    text = "Sky View Factor (%)";
	    length = 19;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SkyViewMap[j][i] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 42) {
	    text = "Shade Map  (%)";
	    length = 14;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = (float) ShadowMap[DayStep][j][i] / 0.2223191;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (temp_array[j][i] < 0.0)
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 43) {
	    text = "Direct Shortwave (W/sqm)";
	    length = 24;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = RadMap[j][i].Beam;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 44) {
	    text = "Diffuse Shortwave (W/sqm)";
	    length = 25;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = RadMap[j][i].Diffuse;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 45) {
	    text = "Aspect (degrees)";
	    length = 16;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = TopoMap[j][i].Aspect * 57.2957;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 46) {
	    text = "Slope (percent)";
	    length = 15;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = TopoMap[j][i].Slope * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 47 && Options->InitSedFlag) {
	    text = "Total Sediment (m3)";
	    length = 22;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SedMap[j][i].SedFluxOut;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 48 && Options->InitSedFlag) {
	    text = "Erosion (mm)";
	    length = 12;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SedMap[j][i].Erosion;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 49 && Options->RoadRouting) {
	    text = "Road Erosion (mm)";
	    length = 17;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {

		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = Network[j][i].Erosion * 1000.;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }


	  if (MapNumber == 50) {
	    text = "Channel Sub Surf Int (mm)";
	    length = 25;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].ChannelInt * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 51) {
	    text = "Road Sub Surf Inter (mm)";
	    length = 24;
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (INBASIN(TopoMap[j][i].Mask)) {
		  temp = SoilMap[j][i].RoadInt * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (fequal(max, min))
	    scale = 0.0;
	  else
	    scale = 50 / (max - min);

	  /* draw the raster image for the current data set */
	  /* all values set to -9999.0 will be drawn as white */
	  /* each image is left and bottom justified in its drawing area */
	  /* i.e. buf pixels are available on the top and right for text */
	  /* and the color bar */

	  if (expand > 0) {
	    for (i = 0; i < Map->NX; i++) {
	      for (j = 0; j < Map->NY; j++) {
		if (!fequal(temp_array[j][i], -9999.0) && 
		    INBASIN(TopoMap[j][i].Mask)) {
		  index = (int) (scale * (temp_array[j][i] - min));
		  if (index > 49)
		    index = 49;
		  XSetForeground(display, gc, my_color[index].pixel);
		}
		else {
		  XSetForeground(display, gc, white);
		}
		for (ie = PX + i * expand; ie < PX + i * expand + expand; ie++)
		  for (je = PY + j * expand; je < PY + j * expand + expand;
		       je++) { {
		      XDrawPoint(display, window, gc, ie, je + buf);
		  }
		  }
	      }
	    }
	  }
	  else {		/* expand < 0 need to average image */

	    for (i = 0; i < Map->NX / (-expand); i++) {
	      for (j = 0; j < Map->NY / (-expand); j++) {
		jr = j * (-expand);
		ir = i * (-expand);
		temp = 0.0;
		skip_it = 0;
		/* for map numbers less than 26 just get the average or a sample */

		if (MapNumber < 50 && sample == 0) {
		  for (ie = 0; ie < (-expand); ie++) {
		    for (je = 0; je < (-expand); je++) {
		      if (temp_array[je + jr][ie + ir] != -9999.0
			  && (INBASIN(TopoMap[je + jr][ie + ir].Mask)))
			temp = temp + temp_array[je + jr][ie + ir];
		      else
			skip_it = 1;
		    }
		  }
		  temp = temp / ((float) (expand * expand));
		}

		if (MapNumber < 50 && sample == 1) {

		  if (temp_array[jr][ir] != -9999.0
		      && (INBASIN(TopoMap[jr][ir].Mask)))
		    temp = temp_array[jr][ir];
		  else
		    skip_it = 1;

		}

		/* for map numbers of 50 or larger, which are the channel and runoff */
		/* subsurface interception, get the max for each aggregated pixel */
		if (MapNumber > 49) {
		  max_temp = -10000.0;
		  for (ie = 0; ie < (-expand); ie++) {
		    for (je = 0; je < (-expand); je++) {
		      if (INBASIN(TopoMap[je + jr][ie + ir].Mask)) {
			if (temp_array[je + jr][ie + ir] > max_temp)
			  max_temp = temp_array[je + jr][ie + ir];
		      }
		      else {
			skip_it = 1;
		      }
		    }
		  }
		  temp = max_temp;
		  if (temp == -9999.0)
		    skip_it = 1;
		}

		index = (int) (scale * (temp - min));
		if (index > 49)
		  index = 49;

		if (skip_it == 0) {
		  XSetForeground(display, gc, my_color[index].pixel);
		}
		else {
		  XSetForeground(display, gc, white);
		}

		XDrawPoint(display, window, gc, i + PX, j + PY + buf);

	      }
	    }
	  }

	  if (expand > 0)
	    re = (float) expand;
	  else
	    re = 1.0 / ((float) -expand);

	  if (draw_static_colorbar == 1) {
	    /* write the title */
	    XSetForeground(display, gc, black);
	    XSetBackground(display, gc, white);
	    XDrawString(display, window, gc, PX, PY + 40, text, length);

	    /* draw the color bar */
	    for (j = 0; j < Map->NY * re; j++) {
	      XSetForeground(display, gc,
			     my_color[(int) (50 * j / (Map->NY * re))].pixel);
	      /*    if((int)((float)(j*50/(Map->NY*re))/scale+min)==0) XSetForeground(display,gc,white); */
	      XDrawLine(display, window, gc, (int) (PX + Map->NX * re + 10),
			(int) (PY + Map->NY * re - j + buf),
			(int) (PX + Map->NX * re + 20),
			(int) (PY + Map->NY * re - j + buf));
	    }
	  }
	  /* label the color bar */
	  sprintf(text2, "%6f", max);
	  XSetForeground(display, gc, black);
	  XClearArea(display, window, (int) (PX + Map->NX * re),
		     (int) (PY - 20 + buf), 50, 20, False);
	  XDrawString(display, window, gc, (int) (PX + Map->NX * re),
		      (int) (PY - 10 + buf), text2, 6);
	  sprintf(text2, "%6.1f", min);
	  XClearArea(display, window, (int) (PX + Map->NX * re),
		     (int) (PY + Map->NY * re + buf), 50, 30, False);
	  XDrawString(display, window, gc, (int) (PX + Map->NX * re),
		      (int) (PY + Map->NY * re + 20 + buf), text2, 6);

	}
      }
    }
  }

#endif
}
