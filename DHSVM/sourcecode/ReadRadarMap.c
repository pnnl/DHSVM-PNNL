/*
 * SUMMARY:      ReadRadarMap.c - Read radar precipitation map from file
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Read radar precipitation map from file
 * DESCRIP-END.
 * FUNCTIONS:    ReadRadarMap()
 * COMMENTS:
 * $Id: ReadRadarMap.c,v 1.4 2003/07/01 21:26:22 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"

/*****************************************************************************
  ReadRadarMap()
*****************************************************************************/

void ReadRadarMap(DATE * Current, DATE * StartRadar, int Dt, MAPSIZE * Radar,
		  RADARPIX ** RadarMap, char *FileName)
{
  const char *Routine = "ReadRadarMap";
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int RadarStep;		/* Location of current timestep in radarfile */
  int NumberType;		/* number type */
  void *Array;

  if (DEBUG)
    printf("Reading precipitation radar data from file: %s\n", FileName);

  NumberType = NC_FLOAT;

  if (!(Array = (float *) calloc(Radar->NY * Radar->NX,
				 SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);

  RadarStep = NumberOfSteps(StartRadar, Current, Dt);

  /* Read the precipitation */

  Read2DMatrix(FileName, Array, NumberType, Radar->NY, Radar->NX, RadarStep);

  for (y = 0, i = 0; y < Radar->NY; y++)
    for (x = 0; x < Radar->NX; x++, i++)
      RadarMap[y][x].Precip = ((float *) Array)[i];

  free(Array);
}
