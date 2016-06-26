#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "snow.h"
#include "glacier.h"
#include "settings.h"
#include "constants.h"
#include "data.h"

void gl_massbalance(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, GLPIX **GlacierMap, double dt_year, double year_min, double year_max,DATE *Current, DUMPSTRUCT *Dump)
{
 
  FILE   *f_out, *f_out2 ;
  float gl_cov_peyto,gl_cov_bow;
  int  k;
  int x;
  int y;
  char file_out[150];
  char file_out2[150];
  int n_mask;
  float avg_6, avg_7;
  float sum_6,sum_7;
  float count_p,count_b;
  float gl_vol_b,gl_vol_p;
  extern double dx;    
  int mth;
  int yr;

  // h      = malloc(N*sizeof(double));         // ice thickness (vectorized)
 
  k = 0;
 
  n_mask = 0;
  gl_cov_peyto = 0;
  gl_cov_bow = 0;
  sprintf(file_out, "%sbalance_glac6_7.txt", Dump->Path);
   if((f_out=fopen(file_out,"ab"))==NULL)
     {
       printf("main(): Cannot open output file %s\n", file_out);
       exit(1);
     }
   //for(i=0;i<=10;i++){
   /* Identify Individual Glacier in Glacier Mask with Unique Integers
 * To Calculate the cumulative mass balance for each, this is currently
 * set up for glacier identidied with a 6 and 7 in the Glacier Mask */
    sum_6 = sum_7=0.0;
    count_p= count_b=0;
     for (x = 0; x < Map->NX; x++) {
       for (y = 0; y < Map->NY; y++) {
	 if (INBASIN(TopoMap[y][x].Mask)){ 
	   if(GlacierMap[y][x].GlMask==6 && Snow[y][x].Iwq > 0.01){
	     //if(TopoMap[y][x].Dem <= uplim && TopoMap[y][x].Dem >= llim){
	       sum_6 = sum_6 + GlacierMap[y][x].Mbal;
	       count_p +=1;
	   }
	   
	   if(GlacierMap[y][x].GlMask==7 && Snow[y][x].Iwq > 0.01){
	     // if(TopoMap[y][x].Dem <= uplim && TopoMap[y][x].Dem >= llim){
	       sum_7 = sum_7 + GlacierMap[y][x].Mbal;
	      
	       count_b +=1;
	       //}
	   }
	 }
       }
     }
     
     avg_6 = sum_6/count_p;
     avg_7 = sum_7/count_b;

/* Mass Balance estimate represents previous month */
/* Create variable so that the proper month is output */
if(Current->Month>1)
mth = Current->Month-1;
     else
mth = 12;

if(Current->Month==1)
yr = Current->Year-1;
else
yr = Current->Year;

     
     printf("Area average mass balance Glacier 6= %.3f, Glacier 7= %.3f\n",avg_6,avg_7);
     fprintf(f_out, " %04d %02d %02d Glacier 6 = %.3f Glacier 7 = %.3f\n",yr,mth,Current->Day,avg_6,avg_7);
  // printf("Area average mass balance = %.3f,swe= %.3f and iwe=%.3f between %.3f annd %.3f\n",avg,avgswe,avgiwe,uplim,llim);
     
       //fprintf(f_out, " %04d %02d %02d %.3f %.3f %.3f %.3f %.3f\n",Current->Year,Current->Month,Current->Day,avg,avgswe,avgiwe,uplim,llim);

     //uplim = uplim - 100;
     //llim = llim - 100;
     

     //}
     gl_vol_b =  gl_vol_p = 0.0; 
   for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
       if (INBASIN(TopoMap[y][x].Mask)){ 
	 if(Snow[y][x].Iwq > 0.01){
	   if(GlacierMap[y][x].GlMask==6){
	     gl_cov_peyto +=1; 
	     gl_vol_p = gl_vol_p+Snow[y][x].Iwq;
	   }
	   if(GlacierMap[y][x].GlMask==7){
	     gl_cov_bow +=1; 
	     gl_vol_b = gl_vol_b+Snow[y][x].Iwq;
	   }
	 }
       }
       else {
	 Snow[y][x].Iwq = 0;
	 Snow[y][x].iweold = 0.0;
       }
     }
   }
   sprintf(file_out2, "%sgl_cov_glac_6_7.txt", Dump->Path);
   if((f_out2=fopen(file_out2,"ab"))==NULL)
     {
       printf("main(): Cannot open output file %s\n", file_out2);
       exit(1);
     }

   
   printf("Ice-covered area and volume Glacier 6 = %.3f %.3f\n",gl_cov_peyto*dx*dx/1000000,gl_vol_p*dx*dx);
   printf("Ice-covered area and volume Glacier 7   = %.3f %.3f\n",gl_cov_bow*dx*dx/1000000,gl_vol_b*dx*dx);
   
   fprintf(f_out2,"Ice-covered area and volume = %04d %02d %02d Glacier 6 = %.3f %.3f Glacier 7= %.3f %.3f\n",Current->Year,mth,Current->Day, (gl_cov_peyto*dx*dx)/1000000,gl_vol_p*dx*dx,(gl_cov_bow*dx*dx)/1000000,gl_vol_b*dx*dx);
   
      
   fclose(f_out);
   fclose(f_out2);
}
