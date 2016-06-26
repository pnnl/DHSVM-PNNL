#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "snow.h"
#include "glacier.h"
#include "settings.h"
#include "constants.h"
#include "data.h"

 
 int calc_ela(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, GLPIX **GlacierMap, DATE *Current, DUMPSTRUCT *Dump)
 {
 FILE *f_out3;
 float ELA_all;
 float ELA_all_yr;
 float cnt2;
 float ELA_6;
 float ELA_6_yr;
 float cnt6;
 float ELA_7;
 float ELA_7_yr;
 float cnt7;


 float deltaswe;
 float deltaice;
 char file_out3[100];
 int x;
 int y;
 float mbal;

 ELA_all = 0.0;
 ELA_all_yr = 0.0;
 cnt2 = 0.0;

 ELA_6 = 0.0;
 ELA_6_yr = 0.0;
 cnt6 = 0.0;

 ELA_7 = 0.0;
 ELA_7_yr = 0.0;
 cnt7 = 0.0;
 


 for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
if (INBASIN(TopoMap[y][x].Mask)){
/*Calculate ELA for glacier 6*/
 if (GlacierMap[y][x].GlMask == 6){
mbal = 0.0;
deltaice = 0.0;
deltaswe = 0.0;

 deltaswe = (double)Snow[y][x].Swq - (double)Snow[y][x].sweann;
 deltaice = (double)Snow[y][x].Iwq - (double)Snow[y][x].iweann;
 mbal = deltaswe + deltaice;
  
 
 if(mbal > -0.1 && mbal < 0.1 && Snow[y][x].Iwq > 10.0) {
	 ELA_6 += TopoMap[y][x].Dem;
	 cnt6 += 1;
	 }
 
 
 
}
}
}
}
/*Calculate ELA for glacier 7*/
 for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
if (INBASIN(TopoMap[y][x].Mask)){

 if (GlacierMap[y][x].GlMask == 7){
mbal = 0.0;
deltaice = 0.0;
deltaswe = 0.0;

 deltaswe = (double)Snow[y][x].Swq - (double)Snow[y][x].sweann;
 deltaice = (double)Snow[y][x].Iwq - (double)Snow[y][x].iweann;
 mbal = deltaswe + deltaice;
 
/* A larger range of M.B. centered on zero may be required to capture ELA */
/* when using coarse model resolutions */   
 if(mbal > -0.1 && mbal < 0.1 && Snow[y][x].Iwq > 10.0) {
	 ELA_7 += TopoMap[y][x].Dem;
	 cnt7 += 1;
	 }
 
 
 
 
}
}
}
}
/*Calculate ELA for all glaciers    */ 
for (x = 0; x < Map->NX; x++) {
     for (y = 0; y < Map->NY; y++) {
if (INBASIN(TopoMap[y][x].Mask)){
   
mbal = 0.0;
deltaice = 0.0;
deltaswe = 0.0;

 deltaswe = (double)Snow[y][x].Swq - (double)Snow[y][x].sweann;
 deltaice = (double)Snow[y][x].Iwq - (double)Snow[y][x].iweann;
 mbal = deltaswe + deltaice;
 
 Snow[y][x].sweann = Snow[y][x].Swq;
 Snow[y][x].iweann = Snow[y][x].Iwq;
 
 
 if(mbal > -0.1 && mbal < 0.1 && Snow[y][x].Iwq > 10.0) {
	 ELA_all += TopoMap[y][x].Dem;
	 cnt2 += 1;
	 }
 
 
 
 

}

}
}
ELA_all_yr = ELA_all/cnt2;
ELA_6_yr = ELA_6/cnt6;
ELA_7_yr = ELA_7/cnt7;


   sprintf(file_out3, "%sELA_all_6_7.txt", Dump->Path);
   if((f_out3=fopen(file_out3,"ab"))==NULL)
     {
       printf("main(): Cannot open output file %s\n", file_out3);
       exit(1);
     }
         fprintf(f_out3," %04d %.3f %.3f %.3f meters\n", Current->Year, ELA_all_yr, ELA_6_yr, ELA_7_yr);
         fclose(f_out3);
		 printf("ELA = %.2f \n", ELA_all_yr);
}
