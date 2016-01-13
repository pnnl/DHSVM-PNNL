/*
 * SUMMARY:      InitUnitHydrograph.c - Initialize unit hydrograph
 * USAGE:        DHSVM
 *
 * AUTHOR:       Bart Nijssen based on code by Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Dec-11-96
 * DESCRIPTION:  Initialize the unit hydrograph components
 *               
 * DESCRIP-END.
 * FUNCTIONS:    InitUnitHydrograph()
 * COMMENTS:
 * $Id: InitUnitHydrograph.c,v 1.4 2003/07/01 21:26:18 olivier Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "fileio.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "varid.h"

enum KEY { travel_file, hydrograph_file };

/*****************************************************************************
  Function name: InitUnitHydrograph()

  Purpose      :

  Required     :

  Returns      : void

  Modifies     :

  Comments     :
*****************************************************************************/
void InitUnitHydrograph(LISTPTR Input, MAPSIZE * Map, TOPOPIX ** TopoMap,
			UNITHYDR *** UnitHydrograph, float **Hydrograph,
			UNITHYDRINFO * HydrographInfo)
{
  const char *Routine = "InitUnitHydrograph()";
  char VarName[BUFSIZE + 1];	/* Variable name */
  FILE *HydrographFile;
  int MaxTravelTime;
  int NumberType;
  int TravelTimeStep;
  int WaveLength;
  int i;
  int j;
  int x;
  int y;
  STRINIENTRY StrEnv[] = {
    {"ROUTING", "TRAVEL TIME FILE", "", NULL},
    {"ROUTING", "UNIT HYDROGRAPH FILE", "", NULL},
    {NULL, NULL, "", NULL}
  };
  unsigned short *Travel;

  printf("Initializing unit hydrograph\n");

  /* Read the key-entry pairs from the [ROUTING] section of the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
    if (!StrEnv[i].VarStr)
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the travel times */
  GetVarName(006, 0, VarName);
  GetVarNumberType(006, &NumberType);
  if (!(Travel = (unsigned short *) calloc(Map->NX * Map->NY,
					   SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);

  Read2DMatrix(StrEnv[travel_file].VarStr, Travel, NumberType, Map->NY,
	       Map->NX, 0, VarName);

  /* Assign the travel times to the correct pixels */
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      TopoMap[y][x].Travel = Travel[i];
  free(Travel);

  /* Read the unit hydrograph file */
  OpenFile(&HydrographFile, StrEnv[hydrograph_file].VarStr, "r", FALSE);
  fscanf(HydrographFile, "%d", &MaxTravelTime);
  HydrographInfo->MaxTravelTime = MaxTravelTime;

  if (!(HydrographInfo->WaveLength =
	(int *) calloc(MaxTravelTime, sizeof(int))))
    ReportError((char *) Routine, 1);
  if (!(*UnitHydrograph =
	(UNITHYDR **) calloc(MaxTravelTime, sizeof(UNITHYDR *))))
    ReportError((char *) Routine, 1);

  for (i = 0; i < MaxTravelTime; i++) {

    fscanf(HydrographFile, "%d %d", &TravelTimeStep, &WaveLength);
    if (i != TravelTimeStep - 1)
      ReportError(StrEnv[hydrograph_file].VarStr, 46);
    HydrographInfo->WaveLength[i] = WaveLength;

    if (!((*UnitHydrograph)[i] =
	  (UNITHYDR *) calloc(WaveLength, sizeof(UNITHYDR))))
      ReportError((char *) Routine, 1);

    for (j = 0; j < WaveLength; j++) {
      fscanf(HydrographFile, "%d %f",
	     &((*UnitHydrograph)[i][j].TimeStep),
	     &((*UnitHydrograph)[i][j].Fraction));
    }
  }

  HydrographInfo->TotalWaveLength =
    (*UnitHydrograph)[MaxTravelTime - 1][WaveLength - 1].TimeStep + 1;

  if (!(*Hydrograph = (float *) calloc(HydrographInfo->TotalWaveLength,
				       sizeof(float))))
    ReportError((char *) Routine, 1);

  fclose(HydrographFile);
}
