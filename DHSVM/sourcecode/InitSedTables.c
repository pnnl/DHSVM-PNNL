/*
 * SUMMARY:      InitSedTables.c - Initialize lookup tables
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD: Sat Feb  7 17:23:51 1998 by Bart Nijssen <nijssen@u.washington.edu>
 * DESCRIPTION:  Initialize lookup tables
 * DESCRIP-END.
 * FUNCTIONS:    InitSedimentTables() 
 *               InitSedTable() 
 *               InitVegStats() 
 * COMMENTS:     
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "DHSVMerror.h"
#include "Calendar.h"
#include "data.h"
#include "constants.h"
#include "fileio.h"
#include "getinit.h"

int InitSedTable(SEDTABLE **SedType, LISTPTR Input, SOILTABLE **SType);
int InitVegStats(VEGTABLE **VType, LISTPTR Input);
float CalcSatDensity(float AveDensity);

/*******************************************************************************/
/*		      InitSedimentTables()                                 */
/*******************************************************************************/
void InitSedimentTables(int StepsPerDay, LISTPTR Input, SEDTABLE **SedType, SOILTABLE **SType, 
			VEGTABLE **VType, LAYER *Soil, LAYER *Veg)
{
  int NSedimentTypes, NVegTypes;

  printf("Initializing sediment tables\n");

  if((NSedimentTypes = InitSedTable(SedType, Input, SType)) == 0)
    ReportError("Input Sediment File", 8);

  if(Soil->NTypes != NSedimentTypes)
     ReportError("Input Sediment File", 2);

  if((NVegTypes = InitVegStats(VType, Input)) == 0)
    ReportError("Input Vegetation File", 8); 

  if(Veg->NTypes != NVegTypes)
     ReportError("Input Vegetation File", 2);
}

/********************************************************************************
  Function Name: InitSedTable()

  Purpose      : Initialize the sediment lookup table
                 Processes most of the following section in InFileName:
		 [SEDIMENT]

  Required     :
    SEDTABLE **SedType - Pointer to lookup table
    LISTPTR Input     - Pointer to linked list with input info
    LAYER *Soil       - Pointer to structure with soil layer information

  Returns      : Number of soil layers

  Modifies     : SedTable

  Comments     :
********************************************************************************/
int InitSedTable(SEDTABLE **SedType, LISTPTR Input, SOILTABLE **SType)
{
  const char *Routine = "InitSedTable";
  int i;			       /* counter */
  int j;                        /* counter */
  int NSoils;			/* Number of soil types */
  char KeyName[fa_mode+1][BUFSIZE+1];
  char *KeyStr[] = {
    "SOIL DESCRIPTION",
    "KINDEX",
    "D50", 
    "SOIL COHESION DISTRIBUTION",
    "SC MIN",
    "SC MAX",
    "SC MEAN",
    "SC DEV",
    "SC MODE",
    "ANGLE OF INTERNAL FRICTION DISTRIBUTION",
    "AIF MIN",
    "AIF MAX",
    "AIF MEAN",
    "AIF DEV",
    "AIF MODE"
  };
  char SectionName[] = "SEDIMENT";
  char VarStr[fa_mode+1][BUFSIZE+1];
  float AveDensity;        /* Average soil density (km/m^3) */
 
  /* Get the number of different soil types */
  GetInitString(SectionName, "NUMBER OF SOIL TYPES", "", VarStr[0], 
		(unsigned long) BUFSIZE, Input);
  if (!CopyInt(&NSoils, VarStr[0], 1))
    ReportError("NUMBER OF SOIL TYPES", 51);

  if (NSoils == 0)
    return NSoils;

  if (!(*SedType = (SEDTABLE *) calloc(NSoils, sizeof(SEDTABLE))))
    ReportError((char *) Routine, 1);
  
  /********** Read information and allocate memory for each soil type *********/

  for (i = 0; i < NSoils; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= fa_mode; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i+1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j], 
		    (unsigned long) BUFSIZE, Input);
    }
    
    /* Assign the entries to the appropriate variables */
    if (IsEmptyStr(VarStr[sed_description]))
      ReportError(KeyName[sed_description], 51);
    strcpy((*SedType)[i].Desc, VarStr[sed_description]);

    if (!CopyFloat(&((*SedType)[i].KIndex), VarStr[kindex], 1))
      ReportError(KeyName[kindex], 51);

    if (!CopyFloat(&((*SedType)[i].d50), VarStr[dfifty], 1))
      ReportError(KeyName[dfifty], 51);

    if (IsEmptyStr(VarStr[cohesion]))
      ReportError(KeyName[cohesion], 51);
    strcpy((*SedType)[i].Cohesion.Distribution, VarStr[cohesion]);

    if(strcmp((*SedType)[i].Cohesion.Distribution, "NORMAL")==0) {

      if (!CopyFloat(&((*SedType)[i].Cohesion.mean), VarStr[coh_mean], 1))
	ReportError(KeyName[coh_mean], 51);

      if (!CopyFloat(&((*SedType)[i].Cohesion.stdev), VarStr[coh_dev], 1))
	ReportError(KeyName[coh_dev], 51);
    }
    else {
      
      if (!CopyFloat(&((*SedType)[i].Cohesion.min), VarStr[coh_min], 1))
	ReportError(KeyName[coh_min], 51);

      if (!CopyFloat(&((*SedType)[i].Cohesion.max), VarStr[coh_max], 1))
	ReportError(KeyName[coh_max], 51);
    }

    if(strcmp((*SedType)[i].Cohesion.Distribution, "TRIANGULAR")==0) {

      if (!CopyFloat(&((*SedType)[i].Cohesion.mode), VarStr[coh_mode], 1))
	ReportError(KeyName[coh_mode], 51);
    }

    if (IsEmptyStr(VarStr[friction_angle]))
      ReportError(KeyName[friction_angle], 51);
    strcpy((*SedType)[i].Friction.Distribution, VarStr[friction_angle]);

    if(strcmp((*SedType)[i].Friction.Distribution, "NORMAL")==0) {

      if (!CopyFloat(&((*SedType)[i].Friction.mean), VarStr[fa_mean], 1))
	ReportError(KeyName[fa_mean], 51);
      
      if (!CopyFloat(&((*SedType)[i].Friction.stdev), VarStr[fa_dev], 1))
	ReportError(KeyName[fa_dev], 51);
    }
    else {
      
      if (!CopyFloat(&((*SedType)[i].Friction.min), VarStr[fa_min], 1))
	ReportError(KeyName[fa_min], 51);

      if (!CopyFloat(&((*SedType)[i].Friction.max), VarStr[fa_max], 1))
	ReportError(KeyName[fa_max], 51);
    }

    if(strcmp((*SedType)[i].Friction.Distribution, "TRIANGULAR")==0) {

      if (!CopyFloat(&((*SedType)[i].Friction.mode), VarStr[fa_mode], 1))
	ReportError(KeyName[fa_mode], 51);
    }
    
    /* Calculating the saturated soil density based on average density - 
       this is not a weighted average */
    AveDensity = 0.0;
    for (j = 0; j < (*SType)[i].NLayers; j++){
      AveDensity += (*SType)[i].Dens[j];
    }
    AveDensity /= j;
    (*SedType)[i].SatDensity = CalcSatDensity(AveDensity);

  }

  return NSoils;
}

/********************************************************************************
  Function Name: InitVegTable()

  Purpose      : Initialize the vegetation lookup table
                 Processes most of the following section in the input file:
		 [VEGETATION]

  Required     :
    VEGTABLE **VType - Pointer to lookup table
    LISTPTR Input    - Pointer to linked list with input info
    LAYER *Veg       - Pointer to structure with veg layer information

  Returns      : Number of vegetation types

  Modifies     : VegTable and Veg

  Comments     :
********************************************************************************/
int InitVegStats(VEGTABLE **VType, LISTPTR Input)
{
  const char *Routine = "InitVegStats";
  int i;			/* Counter */
  int j;			/* Counter */
  int NVegs;			/* Number of vegetation types */
  char KeyName[vs_mode+1][BUFSIZE+1];
  char *KeyStr[] = {
    "ROOT COHESION DISTRIBUTION",
    "RC MIN",
    "RC MAX",
    "RC MEAN",
    "RC DEV",
    "RC MODE",
    "VEGETATION SURCHARGE DISTRIBUTION",
    "VS MIN",
    "VS MAX",
    "VS MEAN",
    "VS DEV",
    "VS MODE",
  };
  char SectionName[] = "VEGETATION";
  char VarStr[vs_mode+1][BUFSIZE+1];
  

  /* Get the number of different vegetation types */
  GetInitString(SectionName, "NUMBER OF VEGETATION TYPES", "", VarStr[0], 
		(unsigned long) BUFSIZE, Input);
  if (!CopyInt(&NVegs, VarStr[0], 1))
    ReportError("NUMBER OF VEGETATION TYPES", 51);

  if (NVegs == 0)
    return NVegs;
  
  /******* Read information for each vegetation type ******/

  for (i = 0; i < NVegs; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= vs_mode; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i+1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j], 
		    (unsigned long) BUFSIZE, Input);
    }
    
    /* Assign the entries to the appropriate variables */
      
    if (IsEmptyStr(VarStr[root_cohesion]))
      ReportError(KeyName[root_cohesion], 51);
    strcpy((*VType)[i].RootCoh.Distribution, VarStr[root_cohesion]);

    if(strcmp((*VType)[i].RootCoh.Distribution, "NORMAL")==0) {

      if (!CopyFloat(&((*VType)[i].RootCoh.mean), VarStr[rc_mean], 1))
	ReportError(KeyName[rc_mean], 51);

      if (!CopyFloat(&((*VType)[i].RootCoh.stdev), VarStr[rc_dev], 1))
	ReportError(KeyName[rc_dev], 51);
    }
    else {
      
      if (!CopyFloat(&((*VType)[i].RootCoh.min), VarStr[rc_min], 1))
	ReportError(KeyName[rc_min], 51);

      if (!CopyFloat(&((*VType)[i].RootCoh.max), VarStr[rc_max], 1))
	ReportError(KeyName[rc_max], 51);
    }

    if(strcmp((*VType)[i].RootCoh.Distribution, "TRIANGULAR")==0) {

      if (!CopyFloat(&((*VType)[i].RootCoh.mode), VarStr[rc_mode], 1))
	ReportError(KeyName[rc_mode], 51);
    }

    if (IsEmptyStr(VarStr[veg_surcharge]))
      ReportError(KeyName[veg_surcharge], 51);
    strcpy((*VType)[i].VegSurcharge.Distribution, VarStr[veg_surcharge]);

    if(strcmp((*VType)[i].VegSurcharge.Distribution, "NORMAL")==0) {

      if (!CopyFloat(&((*VType)[i].VegSurcharge.mean), VarStr[vs_mean], 1))
	ReportError(KeyName[vs_mean], 51);
      
      if (!CopyFloat(&((*VType)[i].VegSurcharge.stdev), VarStr[vs_dev], 1))
	ReportError(KeyName[vs_dev], 51);
    }
    else {

      if (!CopyFloat(&((*VType)[i].VegSurcharge.min), VarStr[vs_min], 1))
	ReportError(KeyName[vs_min], 51);

      if (!CopyFloat(&((*VType)[i].VegSurcharge.max), VarStr[vs_max], 1))
	ReportError(KeyName[vs_max], 51);
    }

    if(strcmp((*VType)[i].VegSurcharge.Distribution, "TRIANGULAR")==0) {

      if (!CopyFloat(&((*VType)[i].VegSurcharge.mode), VarStr[vs_mode], 1))
	ReportError(KeyName[vs_mode], 51);
    }
  }
  
  return NVegs;
}

