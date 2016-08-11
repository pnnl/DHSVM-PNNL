#ifdef HAVE_GLACIER

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "glacier.h"
#include "settings.h"
#include "constants.h"
#include "data.h"

void main_spinup(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, GLPIX **GlacierMap, double dt_year, double year_min, double year_max,DATE *Current, DUMPSTRUCT *Dump, OPTIONSTRUCT * Options)
{
 
  FILE   *f_out1 ;
  double  *b, *s_init, *s_out,  *b_dot; 
  double  dt_yr;
  double  yr_min;
  double  yr_max;
  float gl_cov;
  float h_max;
  int  k;
  int x;
  int y;
  char file_out1[80];
  float *Array1;
  int n_mask;
  /* Use Glen parameter with Annual units */
  double  A_GLEN  = 7.5738e-17; // 7.5738e-17 Cuffey & Paterson (4th ed) Glen's law parameter in Pa^{-3} yr^{-1} units (same as A_GLEN=2.4e-24 Pa^{-3} s^{-1})

  extern double  RHO;     
  extern double  n_GLEN;
  extern double  C_SLIDE;
  extern double  m_SLIDE;
  extern double  A_tilde;
  extern double  C_tilde;
  extern double g;
  extern double dx;
  extern int N;
  extern double  nm_half;
  extern double  np1;
  extern double  mm_half;
  extern double  m1; 

  
  
  // h      = malloc(N*sizeof(double));         // ice thickness (vectorized)
  b      = malloc(N*sizeof(double));         // bed surface elevation (vectorized)
  s_init = malloc(N*sizeof(double));         // initial ice surface elevation (vectorized)
  s_out  = malloc(N*sizeof(double));         // output ice surface elevation (vectorized)
        
  b_dot      = malloc(N*sizeof(double));
  
  // Array  = malloc(N*sizeof(float));
  // if (!(Array1 = (float *) calloc(Map->NX * Map->NY, sizeof(float))));
  
  dt_yr = dt_year;
 
  yr_min = year_min;
  yr_max =  year_max;
  k = 0;
  nm_half = (n_GLEN-1)/2;
  np1     = n_GLEN+1;
  mm_half = (m_SLIDE-1)/2;
  m1      = m_SLIDE;
  n_mask = 0;
  gl_cov = 0;
  h_max = 0;
  
  //  printf("dt_yr= %f yr_min= %f yr_max= %f", dt_yr, yr_min, yr_max);ic_jp[k] = ny*ic[i] + jp[j];for (i=0; i<nx; i++)
   for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
       k    = x*Map->NY + y;
       if (INBASIN(TopoMap[y][x].Mask)){
	 b_dot[k] = GlacierMap[y][x].Mbal;
	 b[k] = GlacierMap[y][x].b;
	 s_init[k]= b[k];
	 n_mask +=1;
       }
       else {
	 /*The following line forces mass balance to be very negative outside of glacier mask*/
        //b_dot[k] = -20;
	 b[k] = GlacierMap[y][x].b;
	 s_init[k] =  b[k];
       }
	    
     }
   }

   A_tilde = 2*A_GLEN*pow(RHO*g,n_GLEN)/((n_GLEN+2)*pow(dx,2));
   C_tilde = C_SLIDE*pow(RHO*g,m_SLIDE)/pow(dx,2);
   
   SetupIndexArrays();

   RunGlacier(b, s_init, s_out, yr_min, yr_max, dt_yr, b_dot, Options);
  
   gl_cov = 0;
   h_max = 0;
   for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
       if (INBASIN(TopoMap[y][x].Mask)){
	 k    = x*Map->NY + y;
	 GlacierMap[y][x].s_init = s_out[k];
	 GlacierMap[y][x].s_out = s_out[k];
	 GlacierMap[y][x].h = GlacierMap[y][x].s_out - GlacierMap[y][x].b;
	 Snow[y][x].Iwq = GlacierMap[y][x].h * (900./1000.);
	 Snow[y][x].iweold = Snow[y][x].Iwq;
	 GlacierMap[y][x].totmbal = 0.0;
	 if( GlacierMap[y][x].s_out<=GlacierMap[y][x].b)
	   GlacierMap[y][x].s_out  =  GlacierMap[y][x].b;
	 else
	   {
	     if (GlacierMap[y][x].s_out- GlacierMap[y][x].b > h_max)
	       h_max =  GlacierMap[y][x].s_out - GlacierMap[y][x].b;
	   }
	 
	 if (GlacierMap[y][x].s_out > GlacierMap[y][x].b)
	  gl_cov += 1;
       }
       //else{
	 //b_dot[k] = -20;
	 //b[k] = GlacierMap[y][x].b;
	 //s_init[k] =  b[k];
	 
      // } 
       
     }
   }

   sprintf(file_out1, "%sh_spinup.bin", Dump->Path);

   if((f_out1=fopen(file_out1,"wb"))==NULL)
     {
       printf("main(): Cannot open output file %s\n", file_out1);
       exit(1);
     }
   
   for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++){
       k    = x*Map->NY + y;
       //((float *) Array)[y * Map->NX + x] = b_dot[k];
       ((float *) Array1)[y * Map->NX + x] = GlacierMap[y][x].h;
     }
   }
   fwrite(Array1, sizeof(Array1), Map->NY * Map->NX, f_out1);

   free(Array1);
   free(b);
   free(s_out);
   free(s_init);
   
   free(b_dot); 
   fclose(f_out1);
   printf("ALL DONE: %.2f yr integration\n", yr_max);
}
#endif
