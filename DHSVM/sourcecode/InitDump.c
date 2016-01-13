/*
 * SUMMARY:      InitDump.c - Initialize output settings
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Functions to initialize which variables need to be output when
 * DESCRIP-END.
 * FUNCTIONS:    InitDump()
 *               InitStateDump()
 *               InitImageDump()
 *               InitMapDump()
 *               InitPixDump()
 * COMMENTS:
 * $Id: InitDump.c,v 1.11 2004/08/18 01:01:29 colleen Exp $     
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
#include "varid.h"

/*******************************************************************************
  Function name: InitDump()

  Purpose      : Read the model output information from the options file, and 
                 organize what to output when.  This information is in the 
		 [OUTPUT] section

  Required     : 
    LISTPTR Input         - Linked list with input strings
    OPTIONSTRUCT *Options - Mode options
    MAPSIZE *Map          - Information about basin area
    int MaxSoilLayers     - Maximum number of soil layers
    int MaxVegLayers      - Maximum number of vegetation layers
    int Dt                - Time step in seconds
    TOPOPIX **TopoMap     - Information about terrain characteristics
    DUMPSTRUCT *Dump      - Information on what to output when

  Returns      : void

  Modifies     : Members of Dump

  Comments     :
*****************************************************************************/
void InitDump(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
	      int MaxSoilLayers, int MaxVegLayers, int Dt,
	      TOPOPIX ** TopoMap, DUMPSTRUCT * Dump, int *NGraphics,
	      int **which_graphics)
{
  char *Routine = "InitDump";
  int i;
  int x;			/* counter */
  int y;			/* counter */
  int NImageVars;		/* Number of different variables for which to 
				   dump images */
  int NMapVars;			/* Number of different variables for which to 
				   dump maps */
  int temp_count;
  uchar **BasinMask;
  char sumoutfile[100];
  

  STRINIENTRY StrEnv[] = {
    {"OUTPUT", "OUTPUT DIRECTORY", "", ""},
    {"OUTPUT", "INITIAL STATE DIRECTORY", "", ""},
    {"OUTPUT", "NUMBER OF OUTPUT PIXELS", "", ""},
    {"OUTPUT", "NUMBER OF MODEL STATES", "", ""},
    {"OUTPUT", "NUMBER OF MAP VARIABLES", "", ""},
    {"OUTPUT", "NUMBER OF IMAGE VARIABLES", "", ""},
    {"OUTPUT", "NUMBER OF GRAPHICS", "", ""},
    {NULL, NULL, "", NULL},
  };

  printf("Initializing dump procedures\n");

  /* Get the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);

  /* Assign the entries to the variables */
  if (IsEmptyStr(StrEnv[output_path].VarStr))
    ReportError(StrEnv[output_path].KeyName, 51);
  strcpy(Dump->Path, StrEnv[output_path].VarStr);

  // delete any previous failure_summary.txt file
  sprintf(sumoutfile, "%sfailure_summary.txt", Dump->Path);
  if (remove(sumoutfile) != -1)
    printf(" - removed old version of failure_summary.txt\n");

  if (IsEmptyStr(StrEnv[initial_state_path].VarStr))
    strcpy(Dump->InitStatePath, Dump->Path);
  strcpy(Dump->InitStatePath, StrEnv[initial_state_path].VarStr);

  if (IsEmptyStr(StrEnv[npixels].VarStr))
    Dump->NPix = 0;
  else if (!CopyInt(&(Dump->NPix), StrEnv[npixels].VarStr, 1) || Dump->NPix < 0)
    ReportError(StrEnv[npixels].KeyName, 51);

  if (IsEmptyStr(StrEnv[nstates].VarStr))
    Dump->NStates = 0;
  else if (!CopyInt(&(Dump->NStates), StrEnv[nstates].VarStr, 1))
    ReportError(StrEnv[nstates].KeyName, 51);

  if (IsEmptyStr(StrEnv[nmapvars].VarStr))
    NMapVars = 0;
  else if (!CopyInt(&NMapVars, StrEnv[nmapvars].VarStr, 1) || NMapVars < 0)
    ReportError(StrEnv[nmapvars].KeyName, 51);

  if (IsEmptyStr(StrEnv[nimagevars].VarStr))
    NImageVars = 0;
  else if (!CopyInt(&NImageVars, StrEnv[nimagevars].VarStr, 1) ||
	   NImageVars < 0)
    ReportError(StrEnv[nimagevars].KeyName, 51);

  if (IsEmptyStr(StrEnv[ngraphics].VarStr))
    *NGraphics = 0;
  else if (!CopyInt(NGraphics, StrEnv[ngraphics].VarStr, 1) || *NGraphics < 0)
    ReportError(StrEnv[ngraphics].KeyName, 51);

  if (Options->Extent == POINT)
    *NGraphics = 0;

  Dump->NMaps = NMapVars + NImageVars;

  // Open file for recording aggregated values for entire basin
  sprintf(Dump->Aggregate.FileName, "%sAggregated.Values", Dump->Path);
  OpenFile(&(Dump->Aggregate.FilePtr), Dump->Aggregate.FileName, "w", TRUE);

  // If specified, open file for recording aggregated sediment values for entire basin
  if (Options->Sediment) {
    sprintf(Dump->AggregateSediment.FileName, "%sAggregatedSediment.Values", Dump->Path);
    OpenFile(&(Dump->AggregateSediment.FilePtr), Dump->AggregateSediment.FileName, "w", TRUE);

    // Open file for recording mass balance for entire basin
    sprintf(Dump->SedBalance.FileName, "%sMassSediment.Balance", Dump->Path);
    OpenFile(&(Dump->SedBalance.FilePtr), Dump->SedBalance.FileName, "w", TRUE);
  }

  // Open file for recording mass balance for entire basin
  sprintf(Dump->Balance.FileName, "%sMass.Balance", Dump->Path);
  OpenFile(&(Dump->Balance.FilePtr), Dump->Balance.FileName, "w", TRUE);
  
  sprintf(Dump->FinalBalance.FileName, "%sMass.Final.Balance", Dump->Path);
  OpenFile(&(Dump->FinalBalance.FilePtr), Dump->FinalBalance.FileName, "w", TRUE);

  if (Options->Extent != POINT) {

    /* Read remaining information from dump info file */

    if (Dump->NStates > 0)
      InitStateDump(Input, Dump->NStates, &(Dump->DState));

    /* if Dump->NStates < 0, the state will be dumped every time step, this is
       done directly in ExecDump */

    if (!(BasinMask = (uchar **) calloc(Map->NY, sizeof(uchar *))))
		ReportError(Routine, 1);
    for (y = 0; y < Map->NY; y++)
      if (!(BasinMask[y] = (uchar *) calloc(Map->NX, sizeof(uchar))))
		  ReportError(Routine, 1);

    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
	BasinMask[y][x] = TopoMap[y][x].Mask;

    if (Dump->NPix > 0) {
      temp_count = InitPixDump(Input, Map, BasinMask, Dump->Path, Dump->NPix,
			       &(Dump->Pix), Options);

      if (temp_count == 0) {
	Dump->NPix = 0;
	printf("no candidate dump pixels accepted \n");
      }
      else {
	Dump->NPix = temp_count;
	printf("total number of accepted dump pixels %d \n", Dump->NPix);
      }
    }
    for (y = 0; y < Map->NY; y++)
      free(BasinMask[y]);
    free(BasinMask);

    if (Dump->NMaps > 0)
      InitMapDump(Input, Map, MaxSoilLayers, MaxVegLayers, Dump->Path,
		  Dump->NMaps, NMapVars, &(Dump->DMap));
    if (NImageVars > 0)
      InitImageDump(Input, Dt, Map, MaxSoilLayers, MaxVegLayers, Dump->Path,
		    Dump->NMaps, NImageVars, &(Dump->DMap));

    if (*NGraphics > 0)
      InitGraphicsDump(Input, *NGraphics, &which_graphics);

    /* if no network open unit hydrograph file */
    if (!(Options->HasNetwork)) {
      sprintf(Dump->Stream.FileName, "%sStream.Flow", Dump->Path);
      OpenFile(&(Dump->Stream.FilePtr), Dump->Stream.FileName, "w", TRUE);
    }
  }
}

/*******************************************************************************
  Function name: InitGraphicsDump()

  Purpose      : Initialize the model state dumps.  This information is in the 
		 [OUTPUT] section of the input file

  Required     : 
    LISTPTR Input         - Linked list with input strings
    int NStates           - Number of graphics to display
    int *which_graphics        - Array with graphic id's

  Returns      : void

  Modifies     : which_graphics, NStates

  Comments     :
*****************************************************************************/
void InitGraphicsDump(LISTPTR Input, int NGraphics, int ***which_graphics)
{
  char *Routine = "InitGraphicsDump";
  int i;			/* counter */
  char KeyName[BUFSIZE + 1];
  char *KeyStr = "GRAPHICS ID";
  char *SectionName = "OUTPUT";
  char VarStr[BUFSIZE + 1];

  if (((**which_graphics) = (int *) malloc(NGraphics * sizeof(int))) == NULL)
    ReportError(Routine, 1);

  for (i = 0; i < NGraphics; i++) {
    sprintf(KeyName, "%s %d", KeyStr, i + 1);
    GetInitString(SectionName, KeyName, "", VarStr,
		  (unsigned long) BUFSIZE, Input);

    if (!CopyInt(&(**which_graphics)[i], VarStr, 1))
      ReportError("GRAPHICS ID", 51);

  }
}

/*******************************************************************************
  Function name: InitStateDump()

  Purpose      : Initialize the model state dumps.  This information is in the 
		 [OUTPUT] section of the input file

  Required     : 
    LISTPTR Input         - Linked list with input strings
    int NStates           - Number of model states to dump
    DATE **DState         - Array with dump dates

  Returns      : void

  Modifies     : DState and its members

  Comments     :
*****************************************************************************/
void InitStateDump(LISTPTR Input, int NStates, DATE ** DState)
{
  char *Routine = "InitStateDump";
  int i;			/* counter */
  char KeyName[BUFSIZE + 1];
  char *KeyStr = "STATE DATE";
  char *SectionName = "OUTPUT";
  char VarStr[BUFSIZE + 1];

  if (!(*DState = (DATE *) calloc(NStates, sizeof(DATE))))
    ReportError(Routine, 1);

  for (i = 0; i < NStates; i++) {
    sprintf(KeyName, "%s %d", KeyStr, i + 1);
    GetInitString(SectionName, KeyName, "", VarStr,
		  (unsigned long) BUFSIZE, Input);
    if (!SScanDate(VarStr, &((*DState)[i])))
      ReportError(KeyName, 51);
  }
}

/*******************************************************************************
  Function name: InitImageDump()

  Purpose      : Initialize the image dumps.  This information is in the 
		 [OUTPUT] section of the input file

  Required     : 
    LISTPTR Input         - Linked list with input strings
    int Dt                - Model timestep in seconds
    MAPSIZE *MapDump      - Information about areal extent
    int MaxSoilLayers     - Maximum number of soil layers
    int MaxVegLayers      - Maximum number of vegetation layers
    char *Path            - Directory to write output to
    int NMaps             - Number of maps to dump
    int NImages           - Number of images to dump 
    MAPDUMP **DMap        - Array of maps and images to dump

  Returns      : void

  Modifies     : Members of DMap

  Comments     : InitImageDump must be preceded by a call to InitMapDump, since
                 the necessary memory is allocated there
*******************************************************************************/
void InitImageDump(LISTPTR Input, int Dt, MAPSIZE * Map, int MaxSoilLayers,
		   int MaxVegLayers, char *Path, int NMaps, int NImages,
		   MAPDUMP ** DMap)
{
  char *Routine = "InitImageDump";
  DATE End;			/* End of low resolution map dump period */
  DATE Start;			/* Start of low resolution map dump period */
  int i;			/* counter */
  int j;			/* counter */
  int Interval;			/* Interval between low resolution map dumps */
  int MaxLayers;		/* Maximum number of layers allowed for this 
				   variable */
  char KeyName[image_lower + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "IMAGE VARIABLE",
    "IMAGE LAYER",
    "IMAGE START",
    "IMAGE END",
    "IMAGE INTERVAL",
    "IMAGE UPPER LIMIT",
    "IMAGE LOWER LIMIT"
  };
  char *SectionName = "OUTPUT";
  char VarStr[image_lower + 1][BUFSIZE + 1];
  float tmpInterval;

  for (i = NMaps - NImages; i < NMaps; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= image_lower; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i - (NMaps - NImages) + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
		    (unsigned long) BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (!CopyInt(&((*DMap)[i].ID), VarStr[image_variable], 1))
      ReportError(KeyName[image_variable], 51);

    if (!IsValidID((*DMap)[i].ID))
      ReportError("Input Options File", 19);

    if (IsMultiLayer((*DMap)[i].ID)) {
      MaxLayers = GetVarNLayers((*DMap)[i].ID, MaxSoilLayers, MaxVegLayers);
      if (!CopyInt(&((*DMap)[i].Layer), VarStr[image_layer], 1))
	ReportError(KeyName[image_layer], 51);
      if ((*DMap)[i].Layer < 1 || (*DMap)[i].Layer > MaxLayers)
	ReportError("Input Options File", 20);
    }
    else
      (*DMap)[i].Layer = 1;

    (*DMap)[i].Resolution = IMAGE_OUTPUT;

    strncpy((*DMap)[i].FileName, Path, BUFSIZE);
    GetVarAttr(&((*DMap)[i]));
    (*DMap)[i].NumberType = NC_BYTE;
    strcpy((*DMap)[i].Format, "%d");

    CreateMapFile((*DMap)[i].FileName, (*DMap)[i].FileLabel, Map);

    if (!SScanDate(VarStr[image_start], &Start))
      ReportError(KeyName[image_start], 51);

    if (!SScanDate(VarStr[image_end], &End))
      ReportError(KeyName[image_end], 51);

    if (!CopyFloat(&tmpInterval, VarStr[image_interval], 1))
      ReportError(KeyName[image_interval], 51);
    Interval = SECPHOUR * tmpInterval;

    if (Interval % Dt != 0 || Interval <= 0)
      ReportError("Input Options File", 24);

    if (((*DMap)[i].N = NumberOfSteps(&Start, &End, Interval)) < 1)
      ReportError("Input Options File", 25);

    if (!((*DMap)[i].DumpDate = (DATE *) calloc((*DMap)[i].N, sizeof(DATE))))
      ReportError(Routine, 1);

    CopyDate(&((*DMap)[i].DumpDate[0]), &Start);

    for (j = 1; j < (*DMap)[i].N; j++)
      (*DMap)[i].DumpDate[j] =
	NextDate(&((*DMap)[i].DumpDate[j - 1]), Interval);

    if (!CopyFloat(&((*DMap)[i].MaxVal), VarStr[image_upper], 1))
      ReportError(KeyName[image_upper], 51);

    if (!CopyFloat(&((*DMap)[i].MinVal), VarStr[image_lower], 1))
      ReportError(KeyName[image_lower], 51);
  }
}

/*******************************************************************************
  Function name: InitMapDump()

  Purpose      : Initialize the map dumps.  This information is in the 
		 [OUTPUT] section of the input file

  Required     : 
    LISTPTR Input         - Linked list with input strings
    MAPSIZE *MapDump      - Information about areal extent
    int MaxSoilLayers     - Maximum number of soil layers
    int MaxVegLayers      - Maximum number of vegetation layers
    char *Path            - Directory to write output to
    int NTotalMapImages   - Total number of maps and images to dump
    int NMaps             - Number of maps to dump 
    MAPDUMP **DMap        - Array of maps and images to dump

  Returns      : void

  Modifies     : DMap and its members

  Comments     :
*******************************************************************************/
void InitMapDump(LISTPTR Input, MAPSIZE * Map, int MaxSoilLayers,
		 int MaxVegLayers, char *Path, int TotalMapImages, int NMaps,
		 MAPDUMP ** DMap)
{
  char *Routine = "InitMapDump";
  int i;			/* counter */
  int j;			/* counter */
  int MaxLayers;		/* Maximum number of layers allowed for this 
				   variable */
  char KeyName[map_date + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "MAP VARIABLE",
    "MAP LAYER",
    "NUMBER OF MAPS",
    "MAP DATE",
  };
  char *SectionName = "OUTPUT";
  char VarStr[map_date + 1][BUFSIZE + 1];

  if (!(*DMap = (MAPDUMP *) calloc(TotalMapImages, sizeof(MAPDUMP))))
    ReportError(Routine, 1);

  for (i = 0; i < NMaps; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= nmaps; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
		    (unsigned long) BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (!CopyInt(&((*DMap)[i].ID), VarStr[map_variable], 1))
      ReportError(KeyName[map_variable], 51);

    if (!IsValidID((*DMap)[i].ID))
      ReportError("Input Options File", 19);

    if (IsMultiLayer((*DMap)[i].ID)) {
      MaxLayers = GetVarNLayers((*DMap)[i].ID, MaxSoilLayers, MaxVegLayers);
      if (!CopyInt(&((*DMap)[i].Layer), VarStr[map_layer], 1))
	ReportError(KeyName[map_layer], 51);
      if ((*DMap)[i].Layer < 1 || (*DMap)[i].Layer > MaxLayers)
	ReportError("Input Options File", 20);
    }
    else
      (*DMap)[i].Layer = 1;

    (*DMap)[i].Resolution = MAP_OUTPUT;

    strncpy((*DMap)[i].FileName, Path, BUFSIZE);
    GetVarAttr(&((*DMap)[i]));

    CreateMapFile((*DMap)[i].FileName, (*DMap)[i].FileLabel, Map);

    if (!CopyInt(&((*DMap)[i].N), VarStr[nmaps], 1))
      ReportError(KeyName[nmaps], 51);

    if ((*DMap)[i].N < 1)
      ReportError("Input Options File", 22);

    if (!((*DMap)[i].DumpDate = (DATE *) calloc((*DMap)[i].N, sizeof(DATE))))
      ReportError(Routine, 1);

    for (j = 0; j < (*DMap)[i].N; j++) {
      sprintf(KeyName[map_date], "%s %d %d", KeyStr[map_date], j + 1, i + 1);
      GetInitString(SectionName, KeyName[map_date], "", VarStr[map_date],
		    (unsigned long) BUFSIZE, Input);
      if (!SScanDate(VarStr[map_date], &((*DMap)[i].DumpDate[j])))
	ReportError(KeyName[map_date], 51);
    }

    (*DMap)[i].MinVal = 0.0;
    (*DMap)[i].MaxVal = 0.0;
  }
}

/*******************************************************************************
  Function name: InitPixDump()

  Purpose      : Initialize the pixel dumps.  This information is in the 
		 [OUTPUT] section of the input file

  Required     : 
    LISTPTR Input         - Linked list with input strings
    MAPSIZE *Map          - Information about basin extent
    uchar **BasinMask     - Basin mask
    char *Path            - Directory to write output to
    int NPix              - Number of pixels to dump 
    PIXDUMP **Pix         - Array of pixels to dump
    OPTIONSTRUCT *Options - Mode options; affects whether sediment files are initialized

  Returns      : number of accepted dump pixels (i.e. in the mask, etc)

  Modifies     : NPix and its members

  Comments     :
*******************************************************************************/
int InitPixDump(LISTPTR Input, MAPSIZE * Map, uchar ** BasinMask, char *Path,
		int NPix, PIXDUMP ** Pix, OPTIONSTRUCT *Options)
{
  char *Routine = "InitPixDump";
  char Str[BUFSIZE + 1];
  double North;
  double East;
  int i;			/* counter */
  int j;
  int ok;
  char temp_name[BUFSIZE + 1];
  char KeyName[name + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "NORTH COORDINATE",
    "EAST COORDINATE",
    "NAME"
  };
  char *SectionName = "OUTPUT";
  char VarStr[name + 1][BUFSIZE + 1];

  ok = 0;

  if (!(*Pix = (PIXDUMP *) calloc(NPix, sizeof(PIXDUMP))))
    ReportError(Routine, 1);

  for (i = 0; i < NPix; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= name; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
		    (unsigned long) BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (!CopyDouble(&North, VarStr[north], 1))
      ReportError(KeyName[north], 51);

    if (!CopyDouble(&East, VarStr[east], 1))
      ReportError(KeyName[east], 51);

    if (IsEmptyStr(VarStr[name]))
      ReportError(KeyName[name], 51);
    strcpy(temp_name, VarStr[name]);

    /* Convert map coordinates to matrix coordinates */
    (*Pix)[i].Loc.N = Round(((Map->Yorig - 0.5 * Map->DY) - North) / Map->DY);
    (*Pix)[i].Loc.E = Round((East - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);

    if (!InArea(Map, &((*Pix)[i].Loc)) ||
	!INBASIN(BasinMask[(*Pix)[i].Loc.N][(*Pix)[i].Loc.E])) {
      printf("Ignoring dump command for pixel named %s \n", temp_name);
    }
    else {
      printf("Accepting dump command for pixel named %s \n", temp_name);
      sprintf(Str, "%s", temp_name);
      sprintf((*Pix)[ok].OutFile.FileName, "%sPixel.%s", Path, Str);
      if (Options->Sediment)
        sprintf((*Pix)[ok].OutFileSediment.FileName, "%sPixelSediment.%s", Path, Str);
      (*Pix)[ok].Loc.N = (*Pix)[i].Loc.N;
      (*Pix)[ok].Loc.E = (*Pix)[i].Loc.E;
      OpenFile(&((*Pix)[ok].OutFile.FilePtr), (*Pix)[ok].OutFile.FileName, "w", TRUE);
      if (Options->Sediment)
        OpenFile(&((*Pix)[ok].OutFileSediment.FilePtr), (*Pix)[ok].OutFileSediment.FileName, "w", TRUE);
      ok++;
    }
  }
  return ok;
}
