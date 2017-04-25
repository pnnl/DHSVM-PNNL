/*
 * SUMMARY:      InitInterpolationWeights.c - Initialize interpolation
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize the interpolation weights
 * DESCRIP-END.
 * FUNCTIONS:    InitInterpolationWeights()
 * COMMENTS:
 * $Id: InitInterpolationWeights.c,v 1.5 2003/10/28 20:02:43 colleen Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "ParallelDHSVM.h"

 /*****************************************************************************
   InitInterpolationWeights()
 *****************************************************************************/
void InitInterpolationWeights(MAPSIZE *Map, OPTIONSTRUCT *Options,
  TOPOPIX **TopoMap, uchar ****MetWeights, METLOCATION *Stats, int NStats)
{
  const char *Routine = "InitInterpolationWeights";
  uchar **BasinMask;
  int x;			/* counter */
  int y;			/* counter */
  int i;

  if (Options->GRIDMET) {
    for (i = 0; i < NStats; i++) {
      Stats[i].Elev = 0.0;
      if (Global2Local(Map, Stats[i].Loc.E, Stats[i].Loc.N, &x, &y)) {
        Stats[i].Elev = TopoMap[y][x].Dem;
      }
      GA_Fgop(&Stats[i].Elev, 1, "+");
    }
  }

  if (Options->MM5 == TRUE && Options->QPF == FALSE) {
    if (!((*MetWeights) = (uchar ***)calloc(Map->NY, sizeof(uchar **))))
      ReportError("CalcWeights()", 1);

    for (y = 0; y < Map->NY; y++)
      if (!((*MetWeights)[y] = (uchar **)calloc(Map->NX, sizeof(uchar *))))
        ReportError("CalcWeights()", 1);

    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
        (*MetWeights)[y][x] = NULL;
  }
  else {
    if (!(BasinMask = (uchar **)calloc(Map->NY, sizeof(uchar *))))
      ReportError((char *)Routine, 1);
    for (y = 0; y < Map->NY; y++) {
      if (!(BasinMask[y] = (uchar *)calloc(Map->NX, sizeof(uchar))))
        ReportError((char *)Routine, 1);
    }

    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
        BasinMask[y][x] = TopoMap[y][x].Mask;

    CalcWeights(Stats, NStats, Map, BasinMask, MetWeights, Options);

    ParallelBarrier();
    ParallelBarrier();
    if (ParallelRank() == 0) {
      printf("\nSummary info on met stations used for current model run \n");
      printf("        Name\t\tY\tX\tIn Mask\tDefined Elev\tActual Elev\n");
    }
    for (i = 0; i < NStats; i++) {
      if ((Stats[i].Loc.N > Map->gNY || Stats[i].Loc.N < 0 ||
           Stats[i].Loc.E > Map->gNX || Stats[i].Loc.E < 0)) {
        printf("%20s\t%d\t%d\t%5s\t%5.1f\t\t%5s\n",
               Stats[i].Name, Stats[i].Loc.N, Stats[i].Loc.E,
               "NA", Stats[i].Elev, "NA");
      } else {
        if (Global2Local(Map, Stats[i].Loc.E, Stats[i].Loc.N, &x, &y)) {
          printf("%20s\t%d\t%d\t%d\t%5.1f\t\t%5.1f\n",
                 Stats[i].Name, Stats[i].Loc.N, Stats[i].Loc.E,
                 BasinMask[y][x], Stats[i].Elev, TopoMap[y][x].Dem);
        } 
      }
    }
    printf("\n");
    ParallelBarrier();

    for (y = 0; y < Map->NY; y++)
      free(BasinMask[y]);

    free(BasinMask);
  }
}
