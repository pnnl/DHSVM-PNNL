/*
 * SUMMARY:      InitXGraphics.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    May-99
 * DESCRIPTION:  Initialize the X11 graphics for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    InitXGraphics()
 * COMMENTS:
 * $Id: InitXGraphics.c,v 1.4 2003/07/01 21:26:18 olivier Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"

#ifdef HAVE_X11
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

Display *display;
Window window;
GC gc;
XColor my_color[50];
float **temp_array;
long black, white;
int e, ndx;
#endif

void InitXGraphics(int argc, char **argv, int ny, int nx, int nd, 
		   MET_MAP_PIX *** MetMap)
{
  /* following is for the X11 libraries */

  int i, x, y, screen;		/* screen is an int. */
  int border_width;
  int c1, c2, c3;
  float re, best_re;
  long dy, dx;
  int best_e, best_ndx, best_ndy;
  int el, ndxl, ndyl;
  int buf = 50;

#ifdef HAVE_X11
  Colormap cmap;		/* map from pixel values to colors */
  char *window_name = "DHSVM Realtime Display", *display_name = NULL;
  XSizeHints size_hints;

  /* connect to Xserver */

  if ((display = XOpenDisplay(display_name)) == NULL) {
    fprintf(stderr, "InitXGraphics: cannot connect to X server %s\n",
	    XDisplayName(display_name));
    exit(-1);
  }

  /* get screen size */

  screen = DefaultScreen(display);

  dx = 0.95 * DisplayWidth(display, screen);
  dy = 0.95 * DisplayHeight(display, screen);

  border_width = 4;
  x = 0;
  y = 0;

  /* figure out the actual size of the window and stuff */

  /* now we need to figure out the expansion factor for the number of graphs */
  /* want to maximize the size of each graph while constrained by the
     number of graphs and the screen size */
  /* save top 20 pixels of the display for the date counter */
  /* save 40 pixels for the X11 window title bar- allocated by system */
  /* save left 10 pixels for border */
  dx = dx - 10;
  dy = dy - 60;

  best_e = -10;
  best_ndx = 1;

  for (ndxl = 1; ndxl <= nd; ndxl++) {
    for (el = -10; el <= 10; el++) {

      if (el != 0) {
	if (el < 0)
	  re = 1 / ((float) (-el));
	if (el > 0)
	  re = (float) el;
	ndyl = nd / (ndxl);
	if (ndxl * ndyl < nd)
	  ndyl = ndyl + 1;
	c1 = nd * (nx * re + buf) * (ny * re + buf);
	c2 = ndxl * (nx * re + buf);
	c3 = ndyl * (ny * re + buf);
	if (c1 <= dx * dy && c2 <= dx && c3 <= dy) {
	  /*      printf("ndx %d ndy %d e %d \n",ndxl,ndy,el);
	     printf("nx ny buf re %d %d %d %f\n",nx,ny,buf,re);
	     printf("c1 c2 c3 rh1 rh2 rh3 %d %d %d %d %d %d \n",c1,c2,c3,dx*dy,dx,dy); */
	  if (el > best_e) {
	    best_e = el;
	    best_ndx = ndxl;
	    best_re = re;
	  }
	}
      }
    }
  }

  printf("best use of display for %d images: \n", nd);
  printf("Expand images by factor %f with %d columns\n", best_re, best_ndx);

  e = best_e;
  ndx = best_ndx;

  best_ndy = nd / best_ndx;
  if (best_ndy * best_ndx < nd)
    best_ndy = best_ndy + 1;

  /* now create the window given the new size classes */

  dx = best_ndx * (nx * best_re + buf) + 10;
  dy = best_ndy * (ny * best_re + buf) + 60;

  window = XCreateSimpleWindow(display, RootWindow(display, screen),
			       0, 0, dx, dy, border_width,
			       WhitePixel(display, screen),
			       WhitePixel(display, screen));
  size_hints.flags = PPosition | PSize | PMinSize | PMaxSize;
  size_hints.x = 0;
  size_hints.y = 0;
  size_hints.width = dx;
  size_hints.height = dy;
  size_hints.min_width = 200;
  size_hints.min_height = 200;
  size_hints.max_width = DisplayWidth(display, screen);
  size_hints.max_height = DisplayHeight(display, screen);

  XSetStandardProperties(display, window, window_name,
			 window_name, None, NULL, 0, &size_hints);

  black = XBlackPixel(display, screen);
  white = XWhitePixel(display, screen);

  gc = XCreateGC(display, window, 0, 0);

  XMapWindow(display, window);

  cmap = XDefaultColormap(display, screen);
  /* black to blue */
  for (i = 0; i < 10; i++) {
    my_color[i].red = 0;
    my_color[i].green = 0;
    my_color[i].blue = 65535 * i / 9;
  }
  /* blue to cyan */
  for (i = 10; i < 20; i++) {
    my_color[i].red = 0;
    my_color[i].green = i * 65535 / 19;
    my_color[i].blue = 65535;
  }
  /*cyan to green */
  for (i = 20; i < 25; i++) {
    my_color[i].red = 0;
    my_color[i].green = 65535;
    my_color[i].blue = 65535 - (65535 * (i - 20) / 5);
  }
  /*green to yellow */
  for (i = 25; i < 30; i++) {
    my_color[i].red = (i - 25) * 65535 / 5;
    my_color[i].green = 65535;
    my_color[i].blue = 0;
  }
  /*yellow to magenta */
  for (i = 30; i < 40; i++) {
    my_color[i].red = 65535;
    my_color[i].green = 65535 - (65535 * (i - 30) / 9);
    my_color[i].blue = 65535 * (i - 30) / 9;
  }
  /*magenta to red */
  for (i = 40; i < 50; i++) {
    my_color[i].red = 65335;
    my_color[i].green = 0;
    my_color[i].blue = 65535 - (65535 * (i - 40) / 9);
  }

  for (i = 0; i < 50; i++) {
    if (XAllocColor(display, cmap, &my_color[i]) == 0) {
      printf("Can't do DHSVM colors\n");
      exit(3);
    }
  }

  /* done initializing the X11 Display, available for drawing */

  /* initialize the memory used solely by the drawing functions */

  if (!((*MetMap) = (MET_MAP_PIX **) calloc(ny, sizeof(MET_MAP_PIX *))))
    ReportError("InitXGraphics", 1);

  for (y = 0; y < ny; y++) {
    if (!((*MetMap)[y] = (MET_MAP_PIX *) calloc(nx, sizeof(MET_MAP_PIX))))
      ReportError("InitXGraphics", 1);
  }

  if ((temp_array = (float **) malloc(ny * sizeof(float *))) == NULL) {
    ReportError("draw.c", 1);
  }
  for (y = 0; y < ny; y++) {
    if ((temp_array[y] = (float *) malloc(nx * sizeof(float))) == NULL) {
      ReportError("draw.c", 1);
    }
  }

#endif
}
