/*
 * SUMMARY:      InitParameters.c - Initialize constants for sediment tranport model
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Laura C. Bowling/Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Apr-96 
 * DESCRIPTION:  Initialize constants for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    InitParameters()
 * COMMENTS:
 * $Id: InitParameters.c,v 1.11 2004/08/10 20:07:47 tbohn Exp $     
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "Calendar.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "getinit.h"
#include "constants.h"
#include "rad.h"

/*****************************************************************************
  Function name: InitParameterss()

  Purpose      : Initialize constants and settings for DHSVM run
                 Processes the following sections in InFile:
                 [PARAMETERS]

  Required     :
    LISTPTR Input          - Linked list with input strings
    OPTIONSTRUCT *Options   - Structure with different program options
    MAPSIZE *Map            - Coverage and resolution of model area
4
  Returns      : void

  Modifies     : (see list of required above)

  Comments     :
*****************************************************************************/
void InitParameters(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
		    ROADSTRUCT ***Network, CHANNEL *ChannelData, TOPOPIX **TopoMap,
		    TIMESTRUCT * Time, float *SedDiams)
{
  int i, x, y;			/* counter */
  const char *Routine = "InitParameters";
  STRINIENTRY StrEnv[] = {
    {"SEDOPTIONS", "MASS WASTING", "", ""},
    {"SEDOPTIONS", "SURFACE EROSION", "", ""},
    {"SEDOPTIONS", "ROAD EROSION", "", ""},
    {"SEDOPTIONS", "CHANNEL ROUTING", "", ""},
    {"PARAMETERS", "MASS WASTING SPACING", "", ""},
    {"PARAMETERS", "MAXIMUM ITERATIONS", "", ""},
// Channel Parent parameters not currently used
//    {"PARAMETERS", "CHANNEL PARENT D50", "", ""},
//    {"PARAMETERS", "CHANNEL PARENT D90", "", ""},
    {"PARAMETERS", "DEBRIS FLOW D50", "", ""},
    {"PARAMETERS", "DEBRIS FLOW D90", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default, 
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
  }

  /**************** Determine model options ****************/

 /* Determine whether mass wasting model should be run */
  if (strncmp(StrEnv[mass_wasting].VarStr, "TRUE", 4) == 0) {
    printf("Sediment Mass Wasting component will be run\n");
    Options->MassWaste = TRUE;
    InitMassWaste(Input, Time);
  }
  else if (strncmp(StrEnv[mass_wasting].VarStr, "FALSE", 5) == 0)
    Options->MassWaste = FALSE;
  else
    ReportError(StrEnv[mass_wasting].KeyName, 51);

  /* Determine whether surface erosion model should be run */
  if (strncmp(StrEnv[surface_erosion].VarStr, "TRUE", 4) == 0){
    Options->ErosionPeriod = TRUE;
    Options->SurfaceErosion = TRUE;
    printf("Sediment Surface Erosion component will be run\n");
  }
  else if (strncmp(StrEnv[surface_erosion].VarStr, "FALSE", 5) == 0){
    Options->ErosionPeriod = FALSE;
    Options->SurfaceErosion = FALSE;
  }
  else
    ReportError(StrEnv[surface_erosion].KeyName, 51);

 /* Determine whether road erosion model should be run */
  if (strncmp(StrEnv[road_erosion].VarStr, "TRUE", 4) == 0){
    Options->RoadRouting = TRUE;
    if ((ChannelData->roads) == NULL) {
      printf("Cannot route the road network without the network files!\n");
      Options->RoadRouting = FALSE;
    }
   }
  else if (strncmp(StrEnv[road_erosion].VarStr, "FALSE", 5) == 0)
    Options->RoadRouting = FALSE;
  else
    ReportError(StrEnv[road_erosion].KeyName, 51);

  if(Options->RoadRouting){ 
  printf("Sediment Road Erosion component will be run\n");
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  if ((*Network)[y][x].RoadArea > 0) {
	    if (!((*Network)[y][x].h = 
		  (float *) calloc(CELLFACTOR, sizeof(float))))
	      ReportError((char *) Routine, 1);
	    if (!((*Network)[y][x].startRunoff =
		  (float *) calloc(CELLFACTOR, sizeof(float))))
	      ReportError((char *) Routine, 1);
	    if (!((*Network)[y][x].startRunon =
		  (float *) calloc(CELLFACTOR, sizeof(float))))
	      ReportError((char *) Routine, 1);
	    if (!((*Network)[y][x].OldSedIn = 
		  (float *) calloc(CELLFACTOR, sizeof(float))))
	      ReportError((char *) Routine, 1);
	    if (!((*Network)[y][x].OldSedOut =
		  (float *) calloc(CELLFACTOR, sizeof(float))))
	      ReportError((char *) Routine, 1);
	  }
	}
      }
    }
  }
  
  /* Determine whether channel routing  model should be run */
  if (strncmp(StrEnv[channel_routing].VarStr, "TRUE", 4) == 0){
    Options->ChannelRouting = TRUE;
    printf("Sediment Channel Routing component will be run\n");
  }
  else if (strncmp(StrEnv[channel_routing].VarStr, "FALSE", 5) == 0){
    Options->ChannelRouting = FALSE;
  }
  else
    ReportError(StrEnv[channel_routing].KeyName, 51);
  
  if (!CopyFloat(&(Map->DMASS), StrEnv[mass_spacing].VarStr, 1))
    ReportError(StrEnv[mass_spacing].KeyName, 51);
  
  Map->NYfine = Map->NY * (Map->DY/Map->DMASS);
  Map->NXfine = Map->NX * (Map->DY/Map->DMASS);
  Map->NumCellsfine = 0;

  if (!CopyFloat(&MASSITER, StrEnv[max_iterations].VarStr, 1))
  ReportError(StrEnv[max_iterations].KeyName, 51);

// Channel Parent parameters not currently used
//  if (!CopyFloat(&CHANNELd50, StrEnv[channeld50].VarStr, 1))
//    ReportError(StrEnv[channeld50].KeyName, 51);
//
//  if (!CopyFloat(&CHANNELd90, StrEnv[channeld90].VarStr, 1))
//    ReportError(StrEnv[channeld90].KeyName, 51);

  if (!CopyFloat(&DEBRISd50, StrEnv[debrisd50].VarStr, 1))
    ReportError(StrEnv[debrisd50].KeyName, 51);
  
  if (!CopyFloat(&DEBRISd90, StrEnv[debrisd90].VarStr, 1))
    ReportError(StrEnv[debrisd90].KeyName, 51);

  DistributeSedimentDiams(SedDiams);  /* find diameter for each portion */

  /* Determine surface erosion period */
  if (Options->SurfaceErosion == TRUE){
    InitSurfaceSed(Input, Time);
  }
  /* Store initial sediment routing option for dumping */
  Options->InitSedFlag=Options->SurfaceErosion;
}

/*******************************************************************************
  Function name: InitMassWaste()

  Purpose      : Read in the dates to run the mass wasting model. This information
                 is in the [SEDTIME] section of the sediment input file

  Required     : 
    LISTPTR Input         - Linked list with input strings

  Returns      : void

  Modifies     : Members of Time

  Comments     : 
*******************************************************************************/
void InitMassWaste(LISTPTR Input, TIMESTRUCT *Time)
{
  const char *Routine = "InitMassWaste";
  int i,j;			/* counter */
  char KeyName[mass_wasting_date+1][BUFSIZE + 1];
  char *KeyStr[] = {
    "MASS WASTING DATE", 
   };

  char SectionName[] = "SEDTIME";
  char VarStr[mass_wasting_date+1][BUFSIZE+1];

  /* Get the number of calculation periods */
  GetInitString(SectionName, "MWM TIME STEPS", "", VarStr[0], 
		(unsigned long) BUFSIZE, Input);

  if (!CopyInt(&(Time->NMWMTotalSteps), VarStr[0], 1))
    ReportError("MWM TIME STEPS", 51);

  if ((Time->NMWMTotalSteps) < 0)
    ReportError("MWM TIME STEPS", 51);
 
  if(Time->NMWMTotalSteps > 0){

    if (!((*Time).MWM = (DATE *) calloc(Time->NMWMTotalSteps, sizeof(DATE)))) 
      ReportError((char *)Routine, 1);  
    /* Read the key-entry pairs from the input file */
    for (i = 0; i < Time->NMWMTotalSteps; i++) {
      
      /* Read the key-entry pairs from the input file */
      for (j = 0; j <= mass_wasting_date; j++) {
	sprintf(KeyName[j], "%s %d",KeyStr[j], i+1);
	GetInitString(SectionName, KeyName[j], "", VarStr[j],
		      (unsigned long) BUFSIZE, Input);
      }
      
      if (!SScanDate(VarStr[mass_wasting_date], &((*Time).MWM[i])))
	ReportError(KeyName[mass_wasting_date], 51);
    }
    
    /* Initialize first run date */
    if ((*Time).Start.Julian > (*Time).MWM[0].Julian)
      ReportError("First Mass Wasting Date is before the beginning of the model run", 51);
    
    (*Time).MWMnext = (*Time).MWM[0];
  }
}

/*******************************************************************************
  Function name: InitSurfaceSed()

  Purpose      : Initialize the image dumps.  This information is in the 
		 [OUTPUT] section of the input file

  Required     : 
    LISTPTR Input         - Linked list with input strings
     int NSteps          - Number of images to dump 
 

  Returns      : void

  Modifies     : Members of Time

  Comments     : 
*******************************************************************************/
void InitSurfaceSed(LISTPTR Input, TIMESTRUCT *Time)
{
  const char *Routine = "InitSurfaceSed";
  int i,j;			/* counter */
  char KeyName[erosion_end+1][BUFSIZE + 1];
  char *KeyStr[] = {
    "EROSION START", 
    "EROSION END",
  };

  char SectionName[] = "SEDTIME";
  char VarStr[erosion_end+1][BUFSIZE+1];

  /* Get the number of calculation periods */
  GetInitString(SectionName, "SE TIME STEPS", "", VarStr[0], 
		(unsigned long) BUFSIZE, Input);

  if (!CopyInt(&(Time->NSETotalSteps), VarStr[0], 1))
    ReportError("SE TIME STEPS", 51);

  if ((Time->NSETotalSteps) < 0)
    ReportError("SE TIME STEPS", 51);

  if (!((*Time).StartSed = (DATE *) calloc(Time->NSETotalSteps, sizeof(DATE)))) 
    ReportError((char *)Routine, 1);  

  if (!((*Time).EndSed = (DATE *) calloc(Time->NSETotalSteps, sizeof(DATE)))) 
      ReportError((char *)Routine, 1);

  /* Read the key-entry pairs from the input file */
  for (i = 0; i < Time->NSETotalSteps ; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= erosion_end; j++) {
      sprintf(KeyName[j], "%s %d",KeyStr[j], i+1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
		    (unsigned long) BUFSIZE, Input);
    }
    
    if (!SScanDate(VarStr[erosion_start], &((*Time).StartSed[i])))
      ReportError(KeyName[erosion_start], 51);
        
    if (!SScanDate(VarStr[erosion_end], &((*Time).EndSed[i])))
      ReportError(KeyName[erosion_end], 51);

    /* Ensure that end times are berfore start tiems */
    if (After(&((*Time).StartSed[i]),&((*Time).EndSed[i])))
       ReportError(SectionName, 23);

   }

}
