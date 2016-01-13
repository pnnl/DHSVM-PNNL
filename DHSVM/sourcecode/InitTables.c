/*
 * SUMMARY:      InitTables.c - Initialize lookup tables
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize lookup tables
 * DESCRIP-END.
 * FUNCTIONS:    InitTables() 
 *               InitSoilTable() 
 *               InitVegTable() 
 *               InitSnowTable()
 * COMMENTS:
 * $Id: InitTables.c,v3.1.2 2013/12/11 ning Exp $     
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "Calendar.h"
#include "constants.h"
#include "fileio.h"
#include "getinit.h"

/*******************************************************************************/
/*				  InitTables()                                 */
/*******************************************************************************/
void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options,
		SOILTABLE **SType, LAYER *Soil, VEGTABLE **VType,
		LAYER *Veg, SNOWTABLE **SnowAlbedo)
{
  printf("Initializing tables\n");

  if ((Soil->NTypes = InitSoilTable(Options, SType, Input, Soil,
	  Options->Infiltration)) == 0)
    ReportError("Input Options File", 8);

  if ((Veg->NTypes = InitVegTable(VType, Input, Options, Veg)) == 0)
    ReportError("Input Options File", 8);

  InitSnowTable(SnowAlbedo, StepsPerDay);
  InitSatVaporTable();
}

/********************************************************************************
  Function Name: InitSoilTable()

  Purpose      : Initialize the soil lookup table
                 Processes most of the following section in InFileName:
		 [SOILS]

  Required     :
    SOILTABLE **SType - Pointer to lookup table
    LISTPTR Input     - Pointer to linked list with input info
    LAYER *Soil       - Pointer to structure with soil layer information

  Returns      : Number of soil layers

  Modifies     : SoilTable and Soil

  Comments     :
********************************************************************************/
int InitSoilTable(OPTIONSTRUCT *Options, SOILTABLE ** SType, 
				  LISTPTR Input, LAYER * Soil, int InfiltOption)
{
  const char *Routine = "InitSoilTable";
  int i;			/* counter */
  int j;			/* counter */
  int NSoils;			/* Number of soil types */
  char KeyName[thermal_capacity + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "SOIL DESCRIPTION",
    "LATERAL CONDUCTIVITY",
    "EXPONENTIAL DECREASE",
    "DEPTH THRESHOLD",
    "MAXIMUM INFILTRATION",
    "CAPILLARY DRIVE",
    "SURFACE ALBEDO",
    "MANNINGS N",
    "NUMBER OF SOIL LAYERS",
    "POROSITY",
    "PORE SIZE DISTRIBUTION",
    "BUBBLING PRESSURE",
    "FIELD CAPACITY",
    "WILTING POINT",
    "BULK DENSITY",
    "VERTICAL CONDUCTIVITY",
    "THERMAL CONDUCTIVITY",
    "THERMAL CAPACITY"
  };
  char SectionName[] = "SOILS";
  char VarStr[thermal_capacity + 1][BUFSIZE + 1];

  /* Get the number of different soil types */
  GetInitString(SectionName, "NUMBER OF SOIL TYPES", "", VarStr[0],
		(unsigned long) BUFSIZE, Input);
  if (!CopyInt(&NSoils, VarStr[0], 1))
    ReportError("NUMBER OF SOIL TYPES", 51);

  if (NSoils == 0)
    return NSoils;

  if (!(Soil->NLayers = (int *) calloc(NSoils, sizeof(int))))
    ReportError((char *) Routine, 1);

  if (!(*SType = (SOILTABLE *) calloc(NSoils, sizeof(SOILTABLE))))
    ReportError((char *) Routine, 1);

  /********** Read information and allocate memory for each soil type *********/

  Soil->MaxLayers = 0;

  for (i = 0; i < NSoils; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= thermal_capacity; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
		    (unsigned long) BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (IsEmptyStr(VarStr[soil_description]))
      ReportError(KeyName[soil_description], 51);

    strcpy((*SType)[i].Desc, VarStr[soil_description]);
    (*SType)[i].Index = i;

    if (!CopyFloat(&((*SType)[i].KsLat), VarStr[lateral_ks], 1))
      ReportError(KeyName[lateral_ks], 51);

    if (!CopyFloat(&((*SType)[i].KsLatExp), VarStr[exponent], 1))
      ReportError(KeyName[exponent], 51);

    if (!CopyFloat(&((*SType)[i].DepthThresh), VarStr[depth_thresh], 1))
      ReportError(KeyName[depth_thresh], 51);

    if (!CopyFloat(&((*SType)[i].MaxInfiltrationRate), VarStr[max_infiltration], 1))
      ReportError(KeyName[max_infiltration], 51);

    if (InfiltOption == DYNAMIC) {
	  if (!CopyFloat(&((*SType)[i].G_Infilt), VarStr[capillary_drive], 1))
		ReportError(KeyName[capillary_drive], 51);
	}
    else (*SType)[i].G_Infilt = NOT_APPLICABLE;

    if (!CopyFloat(&((*SType)[i].Albedo), VarStr[soil_albedo], 1))
      ReportError(KeyName[soil_albedo], 51);

    if (!CopyInt(&(*SType)[i].NLayers, VarStr[number_of_layers], 1))
      ReportError(KeyName[number_of_layers], 51);
    Soil->NLayers[i] = (*SType)[i].NLayers;

	if (Options->Routing) {
	  if (!CopyFloat(&((*SType)[i].Manning), VarStr[manning], 1))
        ReportError(KeyName[manning], 51);
	}
	else if (!Options->Routing)
		(*SType)[i].Manning = NOT_APPLICABLE;
    
    if (Soil->NLayers[i] > Soil->MaxLayers)
      Soil->MaxLayers = Soil->NLayers[i];
    
    /* allocate memory for the soil layers */
    if (!((*SType)[i].Porosity = (float *) calloc((*SType)[i].NLayers,
						  sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].PoreDist = (float *) calloc((*SType)[i].NLayers,
						  sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].Press = (float *) calloc((*SType)[i].NLayers,
					       sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].FCap = (float *) calloc((*SType)[i].NLayers,
					      sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].WP = (float *) calloc((*SType)[i].NLayers,
					    sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].Dens = (float *) calloc((*SType)[i].NLayers,
					      sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].Ks = (float *) calloc((*SType)[i].NLayers,
					    sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].KhDry = (float *) calloc((*SType)[i].NLayers,
					       sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].KhSol = (float *) calloc((*SType)[i].NLayers,
					       sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*SType)[i].Ch = (float *) calloc((*SType)[i].NLayers,
					    sizeof(float))))
      ReportError((char *) Routine, 1);

    if (!CopyFloat((*SType)[i].Porosity, VarStr[porosity], (*SType)[i].NLayers))
      ReportError(KeyName[porosity], 51);

    if (!CopyFloat((*SType)[i].PoreDist, VarStr[pore_size],
		   (*SType)[i].NLayers))
      ReportError(KeyName[pore_size], 51);

    if (!CopyFloat((*SType)[i].Press, VarStr[bubbling_pressure],
		   (*SType)[i].NLayers))
      ReportError(KeyName[bubbling_pressure], 51);

    if (!CopyFloat((*SType)[i].FCap, VarStr[field_capacity],
		   (*SType)[i].NLayers))
      ReportError(KeyName[field_capacity], 51);

    if (!CopyFloat((*SType)[i].WP, VarStr[wilting_point], (*SType)[i].NLayers))
      ReportError(KeyName[wilting_point], 51);

    if (!CopyFloat((*SType)[i].Dens, VarStr[bulk_density], (*SType)[i].NLayers))
      ReportError(KeyName[bulk_density], 51);

    if (!CopyFloat((*SType)[i].Ks, VarStr[vertical_ks], (*SType)[i].NLayers))
      ReportError(KeyName[vertical_ks], 51);

    if (!CopyFloat((*SType)[i].KhSol, VarStr[solids_thermal],
		   (*SType)[i].NLayers))
      ReportError(KeyName[solids_thermal], 51);

    if (!CopyFloat((*SType)[i].Ch, VarStr[thermal_capacity],
		   (*SType)[i].NLayers))
      ReportError(KeyName[thermal_capacity], 51);
  }

  for (i = 0; i < NSoils; i++)
    for (j = 0; j < (*SType)[i].NLayers; j++) {
      (*SType)[i].KhDry[j] = CalcKhDry((*SType)[i].Dens[j]);
      if (((*SType)[i].Porosity[j] < (*SType)[i].FCap[j])
	  || ((*SType)[i].Porosity[j] < (*SType)[i].WP[j])
	  || ((*SType)[i].FCap[j] < (*SType)[i].WP[j]))
	ReportError((*SType)[i].Desc, 11);
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
int InitVegTable(VEGTABLE ** VType, LISTPTR Input, OPTIONSTRUCT * Options,
		 LAYER * Veg)
{
  const char *Routine = "InitVegTable";
  int i;			/* Counter */
  int j;			/* Counter */
  float impervious;		/* flag to check whether impervious layers are 
				   specified */
  int NVegs;			/* Number of vegetation types */
  char KeyName[understory_monalb + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "VEGETATION DESCRIPTION",
    "OVERSTORY PRESENT",
    "UNDERSTORY PRESENT",
    "FRACTIONAL COVERAGE",
    "HEMI FRACT COVERAGE",
    "TRUNK SPACE",
    "AERODYNAMIC ATTENUATION",
    "RADIATION ATTENUATION",
    "CLUMPING FACTOR",
    "LEAF ANGLE A",
    "LEAF ANGLE B",
    "SCATTERING PARAMETER",
    "MAX SNOW INT CAPACITY",
    "MASS RELEASE DRIP RATIO",
    "SNOW INTERCEPTION EFF",
    "IMPERVIOUS FRACTION",
	"DETENTION FRACTION",
    "DETENTION DECAY",
    "HEIGHT",
    "MAXIMUM RESISTANCE",
    "MINIMUM RESISTANCE",
    "MOISTURE THRESHOLD",
    "VAPOR PRESSURE DEFICIT",
    "RPC",
    "NUMBER OF ROOT ZONES",
    "ROOT ZONE DEPTHS",
    "OVERSTORY ROOT FRACTION",
    "UNDERSTORY ROOT FRACTION",
    "OVERSTORY MONTHLY LAI",
    "UNDERSTORY MONTHLY LAI",
    "OVERSTORY MONTHLY ALB",
    "UNDERSTORY MONTHLY ALB"
  };
  char SectionName[] = "VEGETATION";
  char VarStr[understory_monalb + 1][BUFSIZE + 1];

  /* Get the number of different vegetation types */
  GetInitString(SectionName, "NUMBER OF VEGETATION TYPES", "", VarStr[0],
		(unsigned long) BUFSIZE, Input);
  if (!CopyInt(&NVegs, VarStr[0], 1))
    ReportError("NUMBER OF VEGETATION TYPES", 51);

  if (NVegs == 0)
    return NVegs;

  if (!(Veg->NLayers = (int *) calloc(NVegs, sizeof(int))))
    ReportError((char *) Routine, 1);

  if (!(*VType = (VEGTABLE *) calloc(NVegs, sizeof(VEGTABLE))))
    ReportError((char *) Routine, 1);

  /******* Read information and allocate memory for each vegetation type ******/

  Veg->MaxLayers = 0;
  impervious = 0.0;
  for (i = 0; i < NVegs; i++) {

    /* Read the key-entry pairs from the input file */
    for (j = 0; j <= understory_monalb; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
		    (unsigned long) BUFSIZE, Input);
    }

    /* Assign the entries to the appropriate variables */
    if (IsEmptyStr(VarStr[veg_description]))
      ReportError(KeyName[veg_description], 51);
    strcpy((*VType)[i].Desc, VarStr[veg_description]);
    MakeKeyString(VarStr[veg_description]);	/* basically makes the string all
						   uppercase and removed spaces so
						   it is easier to compare */
    if (strncmp(VarStr[veg_description], "GLACIER", strlen("GLACIER")) == 0) {
      (*VType)[i].Index = GLACIER;
    }
    else
      (*VType)[i].Index = i;

    (*VType)[i].NVegLayers = 0;

    if (strncmp(VarStr[overstory], "TRUE", 4) == 0) {
      (*VType)[i].OverStory = TRUE;
      ((*VType)[i].NVegLayers)++;
    }
    else if (strncmp(VarStr[overstory], "FALSE", 5) == 0)
      (*VType)[i].OverStory = FALSE;
    else
      ReportError(KeyName[overstory], 51);

    if (strncmp(VarStr[understory], "TRUE", 4) == 0) {
      (*VType)[i].UnderStory = TRUE;
      ((*VType)[i].NVegLayers)++;
    }
    else if (strncmp(VarStr[understory], "FALSE", 5) == 0)
      (*VType)[i].UnderStory = FALSE;
    else
      ReportError(KeyName[understory], 51);

    Veg->NLayers[i] = (*VType)[i].NVegLayers;
    if ((*VType)[i].NVegLayers > Veg->MaxLayers)
      Veg->MaxLayers = (*VType)[i].NVegLayers;

    if (!CopyInt(&(*VType)[i].NSoilLayers, VarStr[number_of_root_zones], 1))
      ReportError(KeyName[number_of_root_zones], 51);

    if (!CopyFloat(&((*VType)[i].ImpervFrac), VarStr[imperv_frac], 1))
		ReportError(KeyName[imperv_frac], 51);
	impervious += (*VType)[i].ImpervFrac;

	if ((*VType)[i].ImpervFrac > 0) {
	  if (!CopyFloat(&((*VType)[i].DetentionFrac), VarStr[detention_frac], 1))
		ReportError(KeyName[detention_frac], 51);
	  if (!CopyFloat(&((*VType)[i].DetentionDecay), VarStr[detention_decay], 1))
		ReportError(KeyName[detention_decay], 51);
	}
	else {
	  (*VType)[i].DetentionFrac = 0.;
	  (*VType)[i].DetentionDecay = 0.;
	}

    /* allocate memory for the vegetation layers */

    if (!((*VType)[i].Fract = (float *) calloc((*VType)[i].NVegLayers,
					       sizeof(float))))
      ReportError((char *) Routine, 1);

    if (Options->CanopyRadAtt == VARIABLE) {
      if (!((*VType)[i].HemiFract = (float *) calloc((*VType)[i].NVegLayers,
						     sizeof(float))))
      ReportError((char *) Routine, 1);
    }
    else {
      (*VType)[i].HemiFract = NULL;
    }

    if (!((*VType)[i].Height = (float *) calloc((*VType)[i].NVegLayers,
						sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].RsMax = (float *) calloc((*VType)[i].NVegLayers,
					       sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].RsMin = (float *) calloc((*VType)[i].NVegLayers,
					       sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].MoistThres = (float *) calloc((*VType)[i].NVegLayers,
						    sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].VpdThres = (float *) calloc((*VType)[i].NVegLayers,
						  sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].Rpc = (float *) calloc((*VType)[i].NVegLayers,
					     sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].Albedo = (float *) calloc(((*VType)[i].NVegLayers + 1),
						sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].MaxInt = (float *) calloc((*VType)[i].NVegLayers,
						sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].LAI = (float *) calloc((*VType)[i].NVegLayers,
					     sizeof(float))))
      ReportError((char *) Routine, 1);
    if (!((*VType)[i].RootFract = (float **) calloc((*VType)[i].NVegLayers,
						    sizeof(float *))))
      ReportError((char *) Routine, 1);

    for (j = 0; j < (*VType)[i].NVegLayers; j++) {
      if (!((*VType)[i].RootFract[j] =
	    (float *) calloc((*VType)[i].NSoilLayers, sizeof(float))))
	ReportError((char *) Routine, 1);
    }
    if (!((*VType)[i].RootDepth = (float *) calloc((*VType)[i].NSoilLayers,
						   sizeof(float))))
      ReportError((char *) Routine, 1);

    if (!((*VType)[i].LAIMonthly = (float **) calloc((*VType)[i].NVegLayers,
						     sizeof(float *))))
      ReportError((char *) Routine, 1);
    for (j = 0; j < (*VType)[i].NVegLayers; j++) {
      if (!((*VType)[i].LAIMonthly[j] = (float *) calloc(12, sizeof(float))))
	ReportError((char *) Routine, 1);
    }
    
    if (!((*VType)[i].AlbedoMonthly = (float **) calloc((*VType)[i].NVegLayers,
							sizeof(float *))))
      ReportError((char *) Routine, 1);
    for (j = 0; j < (*VType)[i].NVegLayers; j++) {
      if (!((*VType)[i].AlbedoMonthly[j] = (float *) calloc(12, sizeof(float))))
	ReportError((char *) Routine, 1);
    }

    /* assign the entries to the appropriate variables */
    /* allocation of zero memory is not supported on some
       compilers */
    if ((*VType)[i].OverStory == TRUE) {
		if (!CopyFloat(&((*VType)[i].Fract[0]), VarStr[fraction], 1))
			ReportError(KeyName[fraction], 51);
		
		if (Options->CanopyRadAtt == VARIABLE) {
			if (!CopyFloat(&((*VType)[i].HemiFract[0]), VarStr[hemifraction], 1))
				ReportError(KeyName[hemifraction], 51);
			if (!CopyFloat(&((*VType)[i].ClumpingFactor),VarStr[clumping_factor], 1))
			   ReportError(KeyName[clumping_factor], 51);
			if (!CopyFloat(&((*VType)[i].LeafAngleA), VarStr[leaf_angle_a], 1))
				ReportError(KeyName[leaf_angle_a], 51);
			if (!CopyFloat(&((*VType)[i].LeafAngleB), VarStr[leaf_angle_b], 1))
				ReportError(KeyName[leaf_angle_b], 51);
			if (!CopyFloat(&((*VType)[i].Scat), VarStr[scat], 1))
				ReportError(KeyName[scat], 51);
			(*VType)[i].Atten = NOT_APPLICABLE;
		}
		else if (Options->CanopyRadAtt == FIXED) {
			if (!CopyFloat(&((*VType)[i].Atten), VarStr[radiation_att], 1))
				ReportError(KeyName[radiation_att], 51);
			(*VType)[i].ClumpingFactor = NOT_APPLICABLE;
			(*VType)[i].Scat = NOT_APPLICABLE;
			(*VType)[i].LeafAngleA = NOT_APPLICABLE;
			(*VType)[i].LeafAngleB = NOT_APPLICABLE;
		}
		
		if (!CopyFloat(&((*VType)[i].Trunk), VarStr[trunk_space], 1))
			ReportError(KeyName[trunk_space], 51);
		
		if (!CopyFloat(&((*VType)[i].Cn), VarStr[aerodynamic_att], 1))
			ReportError(KeyName[aerodynamic_att], 51);

      if (!CopyFloat(&((*VType)[i].MaxSnowInt), VarStr[snow_int_cap], 1))
		  ReportError(KeyName[snow_int_cap], 51);

      if (!CopyFloat(&((*VType)[i].MDRatio), VarStr[mass_drip_ratio], 1))
		  ReportError(KeyName[mass_drip_ratio], 51);

      if (!CopyFloat(&((*VType)[i].SnowIntEff), VarStr[snow_int_eff], 1))
	ReportError(KeyName[snow_int_eff], 51);

      if (!CopyFloat((*VType)[i].RootFract[0], VarStr[overstory_fraction],
		     (*VType)[i].NSoilLayers))
	ReportError(KeyName[overstory_fraction], 51);

      if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[overstory_monlai], 12))
	ReportError(KeyName[overstory_monlai], 51);

      if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[overstory_monalb],
		     12))
	ReportError(KeyName[overstory_monalb], 51);

      if ((*VType)[i].UnderStory == TRUE) {
	(*VType)[i].Fract[1] = 1.0;
	if (!CopyFloat((*VType)[i].RootFract[1], VarStr[understory_fraction],
		       (*VType)[i].NSoilLayers))
	  ReportError(KeyName[understory_fraction], 51);

	if (!CopyFloat((*VType)[i].LAIMonthly[1], VarStr[understory_monlai],
		       12))
	  ReportError(KeyName[understory_monlai], 51);

	if (!CopyFloat((*VType)[i].AlbedoMonthly[1], VarStr[understory_monalb],
		       12))
	  ReportError(KeyName[understory_monalb], 51);

      }
    }
    else {
      if ((*VType)[i].UnderStory == TRUE) {
	(*VType)[i].Fract[0] = 1.0;
	if (!CopyFloat((*VType)[i].RootFract[0], VarStr[understory_fraction],
		       (*VType)[i].NSoilLayers))
	  ReportError(KeyName[understory_fraction], 51);

	if (!CopyFloat((*VType)[i].LAIMonthly[0], VarStr[understory_monlai],
		       12))
	  ReportError(KeyName[understory_monlai], 51);

	if (!CopyFloat((*VType)[i].AlbedoMonthly[0], VarStr[understory_monalb],
		       12))
	  ReportError(KeyName[understory_monalb], 51);

      }
      (*VType)[i].Trunk = NOT_APPLICABLE;
      (*VType)[i].Cn = NOT_APPLICABLE;
      (*VType)[i].Atten = NOT_APPLICABLE;
      (*VType)[i].ClumpingFactor = NOT_APPLICABLE;
    }

    if (!CopyFloat((*VType)[i].Height, VarStr[height], (*VType)[i].NVegLayers))
      ReportError(KeyName[height], 51);

    if (!CopyFloat((*VType)[i].RsMax, VarStr[max_resistance],
		   (*VType)[i].NVegLayers))
      ReportError(KeyName[max_resistance], 51);

    if (!CopyFloat((*VType)[i].RsMin, VarStr[min_resistance],
		   (*VType)[i].NVegLayers))
      ReportError(KeyName[min_resistance], 51);

    if (!CopyFloat((*VType)[i].MoistThres, VarStr[moisture_threshold],
		   (*VType)[i].NVegLayers))
      ReportError(KeyName[moisture_threshold], 51);

    if (!CopyFloat((*VType)[i].VpdThres, VarStr[vpd], (*VType)[i].NVegLayers))
      ReportError(KeyName[vpd], 51);

    if (!CopyFloat((*VType)[i].Rpc, VarStr[rpc], (*VType)[i].NVegLayers))
      ReportError(KeyName[rpc], 51);

    if (!CopyFloat((*VType)[i].RootDepth, VarStr[root_zone_depth],
		   (*VType)[i].NSoilLayers))
      ReportError(KeyName[root_zone_depth], 51);

    /* Calculate the wind speed profiles and the aerodynamical resistances
       for each layer.  The values are normalized for a reference height wind
       speed of 1 m/s, and are adjusted each timestep using actual reference 
       height wind speeds */

    CalcAerodynamic((*VType)[i].NVegLayers, (*VType)[i].OverStory,
		    (*VType)[i].Cn, (*VType)[i].Height, (*VType)[i].Trunk,
		    (*VType)[i].U, &((*VType)[i].USnow), (*VType)[i].Ra,
		    &((*VType)[i].RaSnow));
  }

  if (impervious) {
    GetInitString(SectionName, "IMPERVIOUS SURFACE ROUTING FILE", "", VarStr[0],
		(unsigned long) BUFSIZE, Input);
    if (IsEmptyStr(VarStr[0]))
      ReportError("IMPERVIOUS SURFACE ROUTING FILE", 51);
    strcpy(Options->ImperviousFilePath, VarStr[veg_description]);
  }
  
  return NVegs;
}

/********************************************************************************
  InitSnowTable()

  Source:
  Laramie, R. L., and J. C. Schaake, Jr., Simulation of the continuous
      snowmelt process, Ralph M. Parsons Laboratory, Mass. Inst. of Technol.,
      1972

  Snow albedo is calculated as a function of the number of days since the
  last observed snow fall.  There are separete albedo curves for the freeze
  and thaw conditions.
********************************************************************************/
void InitSnowTable(SNOWTABLE ** SnowAlbedo, int StepsPerDay)
{
  const char *Routine = "InitSnowTable";
  int i;

  if (!
      (*SnowAlbedo =
       (SNOWTABLE *) calloc((int) ((DAYPYEAR + 1) * StepsPerDay),
			    sizeof(SNOWTABLE))))
    ReportError((char *) Routine, 1);

  /* Laramie and Schaake (1972) */
  /* Updated based on Storck (2000) */
  for (i = 0; i < ((DAYPYEAR + 1) * StepsPerDay); i++) {
    (*SnowAlbedo)[i].Freeze =
      0.85 * pow(0.92, pow(((float) i) / StepsPerDay, 0.58));
    if ((*SnowAlbedo)[i].Freeze < 0.4)
      (*SnowAlbedo)[i].Freeze = 0.4;
    (*SnowAlbedo)[i].Thaw =
      0.85 * pow(0.70, pow(((float) i) / StepsPerDay, 0.46));
    if ((*SnowAlbedo)[i].Thaw < 0.4)
      (*SnowAlbedo)[i].Thaw = 0.4;
  }
}
