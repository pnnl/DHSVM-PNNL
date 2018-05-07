/*
 * SUMMARY:      patch_test.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    February 2017
 * DESCRIPTION:  Tests of GA patch functions
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-05-07 12:47:25 d3g096
 * COMMENTS:
 */

#undef GPTL_TIMING

#include <stdio.h>
#include <ga.h>
#include "ParallelDHSVM.h"

/* -------------------------------------------------------------
   Main Program
   ------------------------------------------------------------- */
int
main(int argc, char **argv)
{
  int me;
  MAPSIZE gMap, Map;
  GA_Patch patch;
  float value;
  int ga;
  int x, y;
  
  ParallelInitialize(&argc, &argv);
  me = ParallelRank();

  gMap.NX = 40;
  gMap.NY = 5;
  gMap.gNX = gMap.NX;
  gMap.gNY = gMap.NY;
  gMap.DX = 10;
  gMap.DY = 10;
  gMap.DXY = 0;
  gMap.OffsetX = 0;
  gMap.OffsetY = 0;
  SimpleDomainDecomposition(&gMap, &Map);
  DomainSummary(&gMap, &Map);

  ga = GA_Duplicate_type(Map.dist, "Patch Test", C_FLOAT);
  value = 0.0;
  GA_Fill(ga, &value);
  
  GA_Alloc_patch(ga, &Map, &patch);

  for (y = 0; y < patch.NY; ++y) {
    for (x = 0; x < patch.NX; ++x) {
      patch.patch[y][x] = (float)me + 1;
    }
  }
  GA_Put_patch(ga, &Map, &patch);
  GA_Free_patch(&patch);

  ParallelBarrier();
  GA_Print(ga);

  value = 0.0;
  GA_Fill(ga, &value);

  GA_Alloc_patch_ghost(ga, &Map, &patch);
  for (y = 0; y < patch.NY; ++y) {
    for (x = 0; x < patch.NX; ++x) {
      patch.patch[y][x] = 1.0;
    }
  }

  GA_Acc_patch(ga, &Map, &patch);
  GA_Free_patch(&patch);

  ParallelBarrier();
  GA_Print(ga);

  ParallelFinalize();
  return 0;
}


