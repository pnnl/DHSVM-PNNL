#ifdef HAVE_GLACIER

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "snow.h"
#include "glacier.h"
#include "settings.h"
#include "constants.h"
#include "data.h"



void main_gl(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, GLPIX **GlacierMap, double dt_year, double year_min, double year_max,DATE *Current, DUMPSTRUCT *Dump,OPTIONSTRUCT * Options)
{
 
  FILE   *f_out1, *f_out2 ;
  double  *b, *s_init, *s_out, *swecc, *iwecc, *h_old,*hcc, *b_dot, *s_initcc; 
  double  dt_yr;
  double  yr_min;
  double  yr_max;
  float gl_cov,wsh_gl_cov;
  float sn_cov, sn_cov_gl;
  float  wsh_sn_cov,wsh_sn_cov_gl,wsh_mask;
  float h_max;
  int  k;
  int x;
  int y;
  char file_out1[150];
  char file_out2[150];
  float *Array1;
  int n_mask;
  float gl_vol;


    
  extern double  RHO;     
  extern double  A_GLEN ;
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
  swecc      = malloc(N*sizeof(double));
  iwecc      = malloc(N*sizeof(double));
  h_old      = malloc(N*sizeof(double));
  hcc      = malloc(N*sizeof(double));
  s_initcc = malloc(N*sizeof(double));
  // Array  = malloc(N*sizeof(float));
  if (!(Array1 = (float *) calloc(Map->NX * Map->NY, sizeof(float))));
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
  sn_cov=sn_cov_gl= 0;
  h_max = 0;
  gl_cov = 0;
  gl_vol = 0;
  wsh_gl_cov = 0;
  sn_cov= 0;
  h_max = 0;
  wsh_sn_cov= wsh_sn_cov_gl=0;
  wsh_mask = 0;
   for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
       k    = x*Map->NY + y;
       if (INBASIN(TopoMap[y][x].Mask)){

	 swecc[k] = (double)Snow[y][x].Swq - (double)Snow[y][x].sweold;  //change in swe in previous month
	 Snow[y][x].sweold = Snow[y][x].Swq;
	 iwecc[k] = (double)Snow[y][x].Iwq - (double)Snow[y][x].iweold;  // change in iwe in previous month
	 

	 b_dot[k] = iwecc[k];    //mass gain or loss of glacier layer used in dynamics
	 GlacierMap[y][x].Mbal = iwecc[k]+swecc[k];
	 GlacierMap[y][x].totmbal =  GlacierMap[y][x].Mbal + GlacierMap[y][x].totmbal; 
	 b[k] = GlacierMap[y][x].b;
	 s_initcc[k] = GlacierMap[y][x].s_init + b_dot[k]; // change in surface topography due to mass balance
	 if(s_initcc[k] <= b[k])
	   s_initcc[k] = b[k];
	 h_old[k] = s_initcc[k] - b[k];
	 s_init[k] =  GlacierMap[y][x].s_init;
	 if (Options->Glacier == GLSTATIC)
	   GlacierMap[y][x].s_out = GlacierMap[y][x].s_init;
	 
	 
       }
       else {
	 //b_dot[k] = -20;
	 b[k] = GlacierMap[y][x].b;
	 h_old[k]= 0.0;
	 s_init[k] =  b[k];
	 //printf("h_old= %f", h_old[k]);
       }
	    
     }
   }
   if (Options->Glacier == GLDYNAMIC)
     {
       printf("Glacier Model monthly run (Dynamic)\n");
       A_tilde = 2*A_GLEN*pow(RHO*g,n_GLEN)/((n_GLEN+2)*pow(dx,2));
       C_tilde = C_SLIDE*pow(RHO*g,m_SLIDE)/pow(dx,2);
       SetupIndexArrays();
       RunGlacier(b, s_init, s_out, yr_min, yr_max, dt_yr, b_dot, Options); 
     }
   
   
   for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
       if (INBASIN(TopoMap[y][x].Mask)){ 
	 if (Options->Glacier == GLDYNAMIC)
     {
	 k    = x*Map->NY + y;
	 GlacierMap[y][x].s_init = s_out[k];
	 GlacierMap[y][x].s_out = s_out[k];
	 GlacierMap[y][x].h = GlacierMap[y][x].s_out - GlacierMap[y][x].b;
	 hcc[k] = GlacierMap[y][x].h - h_old[k]; //change in ice thickness due to glacier movement
	 /*Update Glacier Ice Base on Ice Flow Only (Surface Balance Change Already Accounted for)*/
	 Snow[y][x].Iwq =  Snow[y][x].Iwq + (hcc[k] * (900./1000.)); 
	 }
	 if(GlacierMap[y][x].GlMask>=1){
	   Snow[y][x].Iwq = Snow[y][x].Iwq;
	   Snow[y][x].iweold = Snow[y][x].Iwq;
	   
	 }
	 	 
	 else
	   {
     /*The following was added to not allow glaciers to exist outside of Glacier Mask*/
     /*Do to unavoidable inaccuracies in model inputs, in some cases small glaciers */
     /*sometime grow outside of the historical glacier footprint, since they do not contribute to runoff*/
     /*but would erroneously contribute in a future warmer climate they are deleted to avoid error*/
     /*in future glacier contribution. The amount removed is tracked with the IceRemoved variable*/
	     Snow[y][x].IceRemoved += Snow[y][x].Iwq;
             Snow[y][x].Iwq = 0;
	     Snow[y][x].iweold = 0.0;
             b[k] = GlacierMap[y][x].b;
             h_old[k]= 0.0;
             s_init[k] =  b[k];


	   }
	 
	 if(Snow[y][x].Iwq < 0.00){
	   Snow[y][x].Iwq = 0;
	   Snow[y][x].iweold = 0.0;
	 }
	/*Only Inluding Pixels with more than 1 meter I.W.E. in extent calculation*/ 
         if(GlacierMap[y][x].GlMask>=1){
             n_mask += 1;
	     gl_vol += Snow[y][x].Iwq * dx * dx;
	   if(Snow[y][x].Iwq > 1.0)
	     gl_cov +=1;    
	 }
	 if(Snow[y][x].Iwq > 1.0 && GlacierMap[y][x].WshMask==1) 
	   wsh_gl_cov +=1; 
	 if(Snow[y][x].Swq > 0){  
	   sn_cov +=1;
	   if(Snow[y][x].Iwq > 1.0) 
	     sn_cov_gl +=1;
	 }
	 if(GlacierMap[y][x].WshMask==1 && Snow[y][x].Swq > 0){
	   wsh_sn_cov +=1;
	   if(Snow[y][x].Iwq > 1.0) 
	     wsh_sn_cov_gl +=1;
	 }
	 if(GlacierMap[y][x].WshMask==1)
	   wsh_mask +=1;
       }
       else {
	 
	 b[k] = GlacierMap[y][x].b;
	 h_old[k]= 0.0;
	 s_init[k] =  b[k];
         Snow[y][x].Iwq = 0;
	 Snow[y][x].iweold = 0.0;
       }
     }
   }

     
   

   /********************************OUTPUTS*****************************/
   sprintf(file_out2, "%sgl_sn_cov.txt", Dump->Path);
   if((f_out2=fopen(file_out2,"ab"))==NULL)
     {
       printf("main(): Cannot open output file %s\n", file_out2);
       exit(1);
     }

   
   
   //gl_cov = 100*gl_cov/(float)n_mask;
   //sn_cov = 100*sn_cov/(float)n_mask;
   //wsh_gl_cov = 100*wsh_gl_cov/wsh_mask;
   //wsh_sn_cov = 100*wsh_sn_cov/wsh_mask;
   //printf("Ice Extent (percent of initial) = %.3f percent\n", 100*gl_cov/(float)n_mask);
   printf("Watershed Snow-covered area      = %.3f percent\n", 100*wsh_sn_cov/wsh_mask);
   printf("Watershed glacier-covered area   = %.3f percent\n", 100*wsh_gl_cov/wsh_mask);

   fprintf(f_out2,"Ice-covered area (km^2)                                   = %04d %02d %02d %.3f\n",Current->Year,Current->Month,Current->Day, (gl_cov*dx*dx)/1000000);
   fprintf(f_out2, "Snow-covered area (km^2)                                = %04d %02d %02d %.3f\n", Current->Year,Current->Month,Current->Day,(sn_cov*dx*dx)/1000000);
   fprintf(f_out2, "Snow-covered area on glacier surfaces (km^2)            = %04d %02d %02d %.3f\n", Current->Year,Current->Month,Current->Day,(sn_cov_gl*dx*dx)/1000000);
   fprintf(f_out2, "Watershed Snow-covered (km^2)                           = %04d %02d %02d %.3f\n", Current->Year,Current->Month,Current->Day,(wsh_sn_cov*dx*dx)/1000000);
   fprintf(f_out2, "Watershed Snow-covered area on glacier surface (km^2)   = %04d %02d %02d %.3f\n", Current->Year,Current->Month,Current->Day,(wsh_sn_cov_gl*dx*dx)/1000000);
   fprintf(f_out2, "Watershed glacier-covered area (km^2)                   = %04d %02d %02d %.3f\n", Current->Year,Current->Month,Current->Day,(wsh_gl_cov*dx*dx)/1000000);

    fprintf(f_out2, "Glacier Volume                                         = %04d %02d %02d %.3f\n", Current->Year,Current->Month,Current->Day,gl_vol);

   
   sprintf(file_out1, "%sbalance_sum.bin", Dump->Path);


     if((f_out1=fopen(file_out1,"wb"))==NULL)
     {
       printf("main(): Cannot open output file %s\n", file_out1);
       exit(1);
     }
   
   for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++){
       k    = x*Map->NY + y;
       //((float *) Array)[y * Map->NX + x] = b_dot[k];
       ((float *) Array1)[y * Map->NX + x] = GlacierMap[y][x].totmbal;
     }
   }
   fwrite(Array1, sizeof(Array1), Map->NY * Map->NX, f_out1);

   free(Array1);
   free(b);
   free(s_out);
   free(s_init);
   
   free(b_dot); 
   free(swecc);
   free(iwecc);
   free(hcc);
   free(h_old);
   free(s_initcc);
   fclose(f_out1);
   fclose(f_out2);
   printf("ALL DONE: %.2f yr integration\n", yr_max);
}
#endif
