/*
 * SUMMARY:      deg2tum.c - 
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Ning Sun
 * ORG:          PNNL
 * E-MAIL:       ning@pnnl.gov
 * ORIG-DATE:    May-2016
 * DESCRIPTION:  Convert latitude and longitude to UTM coordinates
 * DESCRIP-END.
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "functions.h"
#include "constants.h"

void deg2utm(float la, float lo, float *x, float *y)
{
   float a, c, v, e2, e2cuadrada, S, deltaS, epsilon, a1, a2, j2, j4, j6;
   float alfa, beta, gama, Bm, xx, yy, ta;
   float lat, lon, nu, Huso;
   float sa = 6378137.000000 ; 
   float sb = 6356752.314245;
   char Letra;
         
   e2 = pow((pow(sa, 2) - pow(sb, 2)), 0.5) / sb;
   e2cuadrada = pow(e2, 2);
   c = pow(sa, 2) / sb;
   
   lat = la * ( PI / 180 );
   lon = lo * ( PI / 180 );

   Huso = floor ( ( lo / 6 ) + 31);
   S = 6 * Huso  - 183;
   deltaS = lon - ( S * ( PI / 180 ) );


   a = cos(lat) * sin(deltaS);
   epsilon = 0.5 * log((1. + a) / (1. - a));
   nu = atan(tan(lat) / cos(deltaS)) - lat;
   v = c / pow((1. + e2cuadrada * pow(cos(lat), 2.) ), 0.5) * 0.9996;
   ta = (e2cuadrada / 2.) * pow(epsilon, 2.) * pow(cos(lat), 2.);
   a1 = sin(2. * lat);
   a2 = a1 * pow(cos(lat), 2.);
   j2 = lat + ( a1 / 2. );
   j4 = (3. * j2  + a2 ) / 4.;
   j6 = (5. * j4  + (a2 * pow(cos(lat), 2.)) ) / 3.;
   alfa = 3. / 4.  * e2cuadrada;
   beta = 5. / 3.  * pow(alfa, 2.);
   gama = 35. / 27. * pow(alfa, 3.);
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
   xx = epsilon * v * ( 1. + ( ta / 3. ) ) + 500000;
   yy = nu * v * ( 1. + ta ) + Bm;

   if (yy < 0)
       yy += 9999999;

   *x = xx;
   *y = yy;
}