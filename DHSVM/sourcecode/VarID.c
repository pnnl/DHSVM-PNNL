/*
 * SUMMARY:      VarID.c - Provide Info about variables
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * ORIG-DATE:    Tue Jan 26 18:26:15 1999
 * E-MAIL:       nijssen@u.washington.edu
 * DESCRIPTION:  Maintains a structure that acts as a database with info on each
 *               variable, and provides functions to query this database
 * DESCRIP-END.
 * FUNCTIONS:    MakeVarAttr()
 *               IsValidDumpID()
 *               IsMultiLayer()
 * COMMENTS:     If the number of IDs increases it might be worthwhile to use a
 *               better, faster search.  This is not done here, because in the 
 *               overall scheme of DHSVM it is not worth the programming effort
 *               right now.
 * $Id: VarID.c,v 1.7 2004/05/04 19:39:00 colleen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "sizeofnt.h"
#include "varid.h"

#ifdef TEST_VARID
char *fileext = ".test";
#else
extern char fileext[];
#endif

struct {
  int ID;
  char Name[BUFSIZE + 1];
  char LongName[BUFSIZE + 1];
  char Format[BUFSIZE + 1];
  char Units[BUFSIZE + 1];
  char FileLabel[BUFSIZE + 1];
  int NumberType;
  int IsMultiLayer;
  int IsVegLayer;
  int IsSoilLayer;
  int AddLayer;
} varinfo[] = {
  {
  001, "Basin.DEM",
      "DEM", "%.3f",
      "m", "Digital Elevation Model", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  002, "Basin.Mask",
      "Basin mask", "%d", "", "Basin mask", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  003, "Soil.Type",
      "Soil type", "%d", "", "Soil type", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  004, "Soil.Depth",
      "Soil depth", "%.3f",
      "m", "Total soil depth", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  005, "Veg.Type",
      "Vegetation type", "%d",
      "", "Vegetation type", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  006, "Travel.Time",
      "Travel time", "%d",
      "hours", "Travel time", NC_SHORT, FALSE, FALSE, FALSE, 0}, {
  101, "Evap.ETot",
      "Evapotranspiration (Total)", "%.4g",
      "m/timestep", "Total amount of evapotranspiration",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  102, "Evap.EPot",
      "Potential Evapotranspiration", "%.4g",
      "m/timestep", "Potential evaporation/transpiration",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  103, "Evap.EInt",
      "Interception Evaporation", "%.4g",
      "m/timestep", "Evaporation from interception",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  104, "Evap.ESoil",
      "Not implemented yet", "%.4g",
      "", "Not implemented yet", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  105, "Evap.EAct",
      "Evaporation", "%.4g",
      "m/timestep", "Actual evaporation/transpiration",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  201, "Precip",
      "Precipitation", "%.4g",
      "m/timestep", "Precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  202, "Precip.IntRain",
      "Interception Storage (liquid)", "%.4g",
      "m", "Interception storage (liquid)", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  203, "Precip.IntSnow",
      "Interception Storage (frozen)", "%.4g",
      "m", "Interception storage (frozen)", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  204, "Temp.Instor",
      "Temporary interception storage for top vegetation layer", "%.4g",
      "m", "Temporary interception storage for top vegetation layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  205, "PRISM.Precip",
      "PRISM Precipitation", "%.4g",
      "mm/month", "PRISM precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  206, "SumPrecip",
      "SumPrecipitation", "%.4g",
      "m", "Accumulated Precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  301, "Rad.ISW",
      "Incoming ShortWave Radiation", "%.4g",
      "W/m2", "Incoming ShortWave Radiation",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  302, "Rad.NSW",
      "Net Shortwave Radiation", "%.4g",
      "W/m2", "Net Shortwave solar radiation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  303, "Rad.Beam",
      "Net Beam Radiation", "%.4g",
      "W/m2", "Net Beam Radiation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  304, "Shade.Factor",
      "Shade Factor", "%d",
      "", "Shade Factor", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  305, "SkyView.Factor",
      "SkyView Factor", "%.4g",
      "-", "Skyview Factor", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  401, "Snow.HasSnow", "Snow Presence/Absence", "%1d", "", "Snow cover flag",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  402, "Snow.SnowCoverOver",
      "Overstory Snow Flag", "%1d", "", "Flag overstory can be covered",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  403, "Snow.LastSnow",
      "Last Snowfall", "%4d", "days", "Days since last snowfall",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  404, "Snow.Swq",
      "Snow Water Equivalent", "%.4g",
      "m", "Snow water equivalent", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  405, "Snow.Melt",
      "Snow Melt", "%.4g",
      "m/timestep", "Snow Melt", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  406, "Snow.PackWater",
      "Liquid Water Content (Deep Layer)", "%.4g",
      "m", "Liquid water content of snow pack",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  407, "Snow.TPack",
      "Snow Temperature (Deep Layer)", "%.4g",
      "C", "Temperature of snow pack", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  408, "Snow.SurfWater",
      "Liquid Water Content (Surface Layer)", "%.4g",
      "m", "Liquid water content of surface layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  409, "Snow.TSurf",
      "Snow Temperature (Surface Layer)", "%.4g",
      "C", "Temperature of snow pack surface layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  410, "Snow.ColdContent",
      "Snow Cold Content", "%.4g",
      "J", "Cold content of snow pack", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  501, "Soil.Moist",
      "Soil Moisture Content", "%.4g",
      "", "Soil moisture for layer %d", NC_FLOAT, TRUE, FALSE, TRUE, 0}, {
  502, "Soil.Perc",
      "Percolation", "%.4g",
      "m/timestep", "Percolation", NC_FLOAT, TRUE, FALSE, TRUE, 0}, {
  503, "Soil.TableDepth",
      "Water Table Depth", "%.4g",
      "m below surface", "Depth of water table",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  504, "Soil.NetFlux",
      "Net Water Flux", "%.4g",
      "m/timestep", "Net flux of water", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  505, "Soil.TSurf",
      "Surface Temperature", "%.4g",
      "C", "Soil surface temperature", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  506, "Soil.Qnet",
      "Net Radiation", "%.4g",
      "W/m2", "Net radiation exchange at surface",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  507, "Soil.Qs",
      "Sensible Heat Flux", "%.4g",
      "W/m2", "Sensible heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  508, "Soil.Qe",
      "Latent Heat Flux", "%.4g",
      "W/m2", "Latent heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  509, "Soil.Qg",
      "Ground Heat Flux", "%.4g",
      "W/m2", "Ground heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  510, "Soil.Qst",
      "Ground Heat Storage", "%.4g",
      "W/m2", "Ground heat storage", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  511, "Soil.Temp",
      "Soil Temperature", "%.4g",
      "C", "Soil Temperature", NC_FLOAT, TRUE, FALSE, TRUE}, {
  512, "Soil.Runoff",
      "Surface Ponding", "%.4g",
      "m", "Surface Ponding", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  513, "SoilMap.IExcess",
      "Surface runoff from HOF and Return Flow", "%.4g",
      "m", "Surface runoff from HOF and Return Flow",
       NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  514, "SoilMap.InfiltAcc",
      "Infiltration Accumulation", "%.4g",
      "m", "Accumulated water in top layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  601, "WindModel",
      "Wind Direction Multiplier", "%.5f",
      "", "Wind Direction Multiplier", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  602, "Precip.Lapse",
      "Precipitation Lapse Rate", "%.5f",
      "", "Precipitation Lapse Rate", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  605, "RadarMap.Precip",
      "Radar Precipitation", "%.4f",
      "m/timestep", "Radar precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  701, "MetMap.accum_precip",
      "Accumulated Precipitation", "%.5f",
      "m", "Accumulated Precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  702, "MetMap.air_temp",
      "Air Temperature", "%.2f",
      "C", "Air Temperature", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  703, "MetMap.windspeed",
      "Windspeed", "%.2f",
      "m/s", "Windspeed", NC_INT, FALSE, FALSE, FALSE, 0}, {
  704, "MetMap.humidity",
      "Humidity", "%.2f", "", "Humidity", NC_INT, FALSE, FALSE, FALSE, 0}, {
  801, "FineMap.Dem",
      "Elevation", "%.2f", "m", "Elevation",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  802, "FineMap.Slope",
      "Slope", "%.2f", "", "Land Surface Slope",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  803, "FineMap.SatThickness",
      "SatThickness", "%.2f", "m", "Water Table Thickness",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  804, "FineMap.DeltaDepth",
      "DeltaDepth", "%.2f", "m", "Cumulative Change in Sediment Depth",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  805, "FineMap.Probability",
      "Failure Probability", "%.2f", "", "Pixel Failure Probability",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  806, "FineMap.SedimentToChannel",
      "SedimentToChannel", "%.2f", "m", "Sediment from Mass-Wasting Deposited into Channels",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  901, "SedMap.SedFluxOut",
      "Sediment Flux Out", "%.4g",
      "(m3/m3)/timestep", "Total Sediment Flux from Grid Cell",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  902, "SedMap.Erosion",
      "Erosion", "%.4g",
      "m/timestep", "Change in Grid Cell Elevation due to Erosion",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  903, "Network.Erosion",
      "Road Erosion", "%.4g",
      "m/timestep", "Change in Road Elevation due to Erosion",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  ENDOFLIST, "", "", "", "", "",
      ENDOFLIST, ENDOFLIST, ENDOFLIST, ENDOFLIST, ENDOFLIST}
};

/*****************************************************************************
  GetVarAttr()
*****************************************************************************/
void GetVarAttr(MAPDUMP * DMap)
{
  GetVarName(DMap->ID, DMap->Layer, DMap->Name);
  GetVarLongName(DMap->ID, DMap->Layer, DMap->LongName);
  GetVarFormat(DMap->ID, DMap->Format);
  GetVarUnits(DMap->ID, DMap->Units);
  GetVarFileName(DMap->ID, DMap->Layer, DMap->Resolution, DMap->FileName);
  GetVarFileLabel(DMap->ID, DMap->FileLabel);
  GetVarNumberType(DMap->ID, &(DMap->NumberType));
}

/******************************************************************************/
/*				  GetVarName()                                */
/******************************************************************************/
void GetVarName(int ID, int Layer, char *Name)
{
  char *Routine = "GetVarName";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsMultiLayer == TRUE)
	sprintf(Name, "%d.%s", Layer, varinfo[i].Name);
      else
	strcpy(Name, varinfo[i].Name);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				GetVarLongName()                              */
/******************************************************************************/
void GetVarLongName(int ID, int Layer, char *LongName)
{
  char *Routine = "GetVarLongName";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsMultiLayer == TRUE)
	sprintf(LongName, "%s (Layer %d)", varinfo[i].LongName, Layer);
      else
	strcpy(LongName, varinfo[i].LongName);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				  GetVarFormat()                              */
/******************************************************************************/
void GetVarFormat(int ID, char *Format)
{
  char *Routine = "GetVarFormat";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(Format, varinfo[i].Format);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				  GetVarUnits()                               */
/******************************************************************************/
void GetVarUnits(int ID, char *Units)
{
  char *Routine = "GetVarUnits";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(Units, varinfo[i].Units);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				 GetVarFileName()                             */
/******************************************************************************/
void GetVarFileName(int ID, int Layer, unsigned char Resolution, char *FileName)
{
  char *Routine = "GetVarFileName";
  char Name[BUFSIZE + 1];
  char Str[BUFSIZE + 1];
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      GetVarName(ID, Layer, Name);
      if (Resolution == MAP_OUTPUT) {
	sprintf(Str, "%sMap.%s%s", FileName, Name, fileext);
      }
      else if (Resolution == IMAGE_OUTPUT) {
	sprintf(Str, "%sImage.%s%s", FileName, Name, fileext);
      }
      else
	ReportError((char *) Routine, 21);
      strncpy(FileName, Str, BUFSIZE);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				GetVarFileLabel()                             */
/******************************************************************************/
void GetVarFileLabel(int ID, char *FileLabel)
{
  char *Routine = "GetVarFileLabel";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(FileLabel, varinfo[i].FileLabel);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*			       GetVarNumberType()                             */
/******************************************************************************/
void GetVarNumberType(int ID, int *NumberType)
{
  char *Routine = "GetVarNumberType";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      *NumberType = varinfo[i].NumberType;
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				    IsValidID()                               */
/******************************************************************************/
unsigned char IsValidID(int ID)
{
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID)
      return TRUE;
    i++;
  }
  return FALSE;
}

/******************************************************************************/
/*				 IsMultiLayer()                               */
/******************************************************************************/
unsigned char IsMultiLayer(int ID)
{
  char *Routine = "IsMultiLayer";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID)
      return varinfo[i].IsMultiLayer;
    i++;
  }
  ReportError((char *) Routine, 26);
  return FALSE;
}

/******************************************************************************/
/*			 	 GetVarNLayers()                              */
/******************************************************************************/
int GetVarNLayers(int ID, int MaxSoilLayers, int MaxVegLayers)
{
  char *Routine = "GetVarNLayers";
  int NLayers = -1;
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsVegLayer == TRUE)
	NLayers = MaxVegLayers + varinfo[i].AddLayer;
      else if (varinfo[i].IsSoilLayer == TRUE)
	NLayers = MaxSoilLayers + varinfo[i].AddLayer;
      else
	NLayers = 1;
      return NLayers;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
  return NLayers;
}

/******************************************************************************/
/* Test main.  Compile by typing:                                             */
/* gcc -Wall -g -o test_varid -DTEST_VARID VarID.c ReportError.c              */
/* Then test by typing test_varid                                             */
/******************************************************************************/
#ifdef TEST_VARID

int main(void)
{
  MAPDUMP DMap;
  int i = 0;

  DMap.Layer = 2;
  while (varinfo[i].ID != ENDOFLIST) {
    strcpy(DMap.FileName, "<path>/");
    DMap.Resolution = i % 2;
    if (DMap.Resolution == 0)
      DMap.Resolution += 2;
    DMap.ID = varinfo[i].ID;
    if (IsValidID(DMap.ID)) {	/* only added to test IsvalidID */
      GetVarAttr(&DMap);
      printf("************************************************************\n");
      printf("ID        : %d\n", DMap.ID);
      printf("Name      : %s\n", DMap.Name);
      printf("LongName  : %s\n", DMap.LongName);
      printf("FileName  : %s\n", DMap.FileName);
      printf("FileLabel : %s\n", DMap.FileLabel);
      printf("Format    : %s\n", DMap.Format);
      printf("Units     : %s\n", DMap.Units);
      printf("NumberType: %d\n", DMap.NumberType);
      printf("NLayers   : %d\n", GetVarNLayers(DMap.ID, 2, 3));
      printf("FileName  : %s\n", DMap.FileName);
      printf("************************************************************\n");
      i++;
    }
  }
  DMap.ID = -1;
  if (IsValidID(DMap.ID)) {	/* only added to test IsvalidID */
    GetVarAttr(&DMap);
  }
  else
    return EXIT_SUCCESS;

  printf("Error: the test program should not have reached this line\n");

  return EXIT_FAILURE;
}

#endif
