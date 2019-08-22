/*
 * SUMMARY:      SnowStats.c - calculate snow stats
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Zhuoran Duan
 * ORG:          Pacific Northwest National Laboratory
 * E-MAIL:       zhuoran.duan@pnnl.gov
 * ORIG-DATE:    Jul-2019
 * DESCRIPTION:  Calculate the values for the snow 
 *               state variables over the basin.
 * DESCRIP-END.
 * FUNCTIONS:    SnowStats()
 * COMMENTS:
 * 
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  SnowStats()
  
  Calculate the statistics for SWE analysis (peak, peak date, melt out date).  
 
  The intial values are set to zero in the function InitNewYear on very first 
  timestep of a new water year.

  Dates were converted to unsigned int in format of YYYYMMDD.
*****************************************************************************/
void SnowStats(DATE *Now, MAPSIZE *Map, OPTIONSTRUCT *Options, 
        TOPOPIX **TopoMap, SNOWPIX **Snow, int Dt)
{
  int x;
  int y;
  int DNum; 
  //printf("updating SWE stats map\n");
 
  DNum = Now->Year * 10000 + Now->Month * 100 + Now->Day;
  // printf("currnet year %d \n", Now->Year);
  // printf("currnet month is %d \n", Now->Month);
  // printf("currnet day is %d \n", Now->Day);
  // printf("currnet DNum is %d \n", DNum);
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
           //printf("currnet SWE is %f \n", Snow[y][x].Swq);
          // Update Peak SWE and Peak SWE date
          if ( Snow[y][x].Swq > Snow[y][x].MaxSwe){
            Snow[y][x].MaxSwe = Snow[y][x].Swq;
            Snow[y][x].MaxSweDate = DNum; 
            /* When the MaxSwe is updated, reset the melt out date to 0 so that it
            overwrites previous in-corret dates*/
            Snow[y][x].MeltOutDate = 0; 
          }

          // Update Peak SWE Date
          /* Criteria :
            1. If snow < 5mm
            2. First date past the peak SWE date
            3. And Preceding 7/15 day has snow  //for now this was not implimented
          */
         if ((Snow[y][x].Swq < MIN_SWE) && (DNum > Snow[y][x].MaxSweDate) && (Snow[y][x].MeltOutDate == 0)){
            Snow[y][x].MeltOutDate = DNum;    
            printf("SWE Melt out date is %d \n", Snow[y][x].MeltOutDate);
            }
      }
    }
  }
}