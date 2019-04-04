
/*
* SUMMARY:      InitParameterMaps.c - Initialize spatial input of parameters
* USAGE:        Part of DHSVM
*
* AUTHOR:       Ning Sun
* ORG:          PNNL
* E-MAIL:       ning.sun@pnnl.gov
* ORIG-DATE:    Mar-2019
* DESCRIPTION:  
* DESCRIP-END.
* FUNCTIONS:    InitParameterMaps()
* COMMENTS:
* $Id: InitParameterMaps.c,v3.2 ning Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "slopeaspect.h"
#include "varid.h"

/*****************************************************************************
InitParameterMaps()
*****************************************************************************/
void InitParameterMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, int Id, 
  char *FileName, SNOWPIX ***SnowMap, int ParamType, float temp) {


  const char *Routine = "InitParameterMaps";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  int flag;         /* either or not reverse the matrix */
  int NumberType;		/* Number type of data set */
  float *Array;

  /* Read the map */

  if (ParamType == MAP) {
	GetVarName(Id, 0, VarName);
	GetVarNumberType(Id, &NumberType);
	if (!(Array = (float *)calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType))))
	  ReportError((char *)Routine, 1);

	flag = Read2DMatrix(FileName, Array, NumberType, Map, 0, VarName, 0);

	/* Assign the attributes to the map pixel */
	/* Reverse the matrix is flag = 1 & netcdf option is selected */
	if ((Options->FileFormat == NETCDF && flag == 0) || (Options->FileFormat == BIN)) {
	  for (y = 0, i = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++, i++) {
		  switch (Id) {
		  case 800:
			(*SnowMap)[y][x].Ts = Array[i];
			break;
		  case 801:
			(*SnowMap)[y][x].Tr = Array[i];
			break;
		  case 802:
			(*SnowMap)[y][x].amax = Array[i];
			break;
		  case 803:
			(*SnowMap)[y][x].LamdaAcc = Array[i];
			break;
		  case 804:
			(*SnowMap)[y][x].LamdaMelt = Array[i];
			break;
		  case 805:
			(*SnowMap)[y][x].AccMin = Array[i];
			break;
		  case 806:
			(*SnowMap)[y][x].MeltMin = Array[i];
			break;
		  default:
			printf("%s: Map ID %d not found", Routine, Id);
			exit(74);
		  } /* end switch (n) */
		}
	  }
	}
	else if (Options->FileFormat == NETCDF && flag == 1) {
	  for (y = Map->NY - 1, i = 0; y >= 0; y--) {
		for (x = 0; x < Map->NX; x++, i++) {
		  switch (Id) {
		  case 800:
			(*SnowMap)[y][x].Ts = Array[i];
			break;
		  case 801:
			(*SnowMap)[y][x].Tr = Array[i];
			break;
		  case 802:
			(*SnowMap)[y][x].amax = Array[i];
			break;
		  case 803:
			(*SnowMap)[y][x].LamdaAcc = Array[i];
			break;
		  case 804:
			(*SnowMap)[y][x].LamdaMelt = Array[i];
			break;
		  case 805:
			(*SnowMap)[y][x].AccMin = Array[i];
			break;
		  case 806:
			(*SnowMap)[y][x].MeltMin = Array[i];
			break;
		  default:
			printf("%s: Map ID %d not found", Routine, Id);
			exit(74);
		  } /* end switch (n) */
		}
	  }
	}
	else
	  ReportError((char *)Routine, 57);

	free(Array);
  }
  /* assign a constant parameter to all model grids*/
  else if (ParamType == CONSTANT) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
		switch (Id) {
		case 800:
		  (*SnowMap)[y][x].Ts = temp;
		  break;
		case 801:
		  (*SnowMap)[y][x].Tr = temp;
		  break;
		case 802:
		  (*SnowMap)[y][x].amax = temp;
		  break;
		case 803:
		  (*SnowMap)[y][x].LamdaAcc = temp;
		  break;
		case 804:
		  (*SnowMap)[y][x].LamdaMelt = temp;
		  break;
		case 805:
		  (*SnowMap)[y][x].AccMin = temp;
		  break;
		case 806:
		  (*SnowMap)[y][x].MeltMin = temp;
		  break;
		default:
		  printf("%s: Map ID %d not found", Routine, Id);
		  exit(0);
		}
	  }
	}
  }
  else {
	printf("%s: Parameter type %d not found", Routine, ParamType);
	exit(75);
  }
}