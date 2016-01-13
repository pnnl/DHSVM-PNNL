/*
 * SUMMARY:      MakeModelStateBin.c - Create initial model state for DHSVM
 * USAGE:        MakeModelStateBin <infofile>
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD:     Mon Jul 28 20:59:20 1997 by Laura Bowling <lxb@u.washington.edu>
 * DESCRIPTION:  Create initial model state for DHSVM (binary I/O format)
 * DESCRIP-END.
 * FUNCTIONS:    ()
 * COMMENTS:
 */

static const char vcid[] = "$Id: MakeModelStateBin.c,v 3.1";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "sizeofnt.h"
#include "settings.h"
#include "fileio.h"
#include "fifobin.h"
#include "Calendar.h"
#include "DHSVMerror.h"

#define MAXLAYERS 10
#define MONTHSPYR 12

void StoreModelState(char *Path, DATE Current, int NY, int NX, 
                     int NVegLayers, float *RainInt, float *SnowInt,
                     float TempIntStorage, unsigned char SnowMask, 
                     unsigned short LastSnow, float Swq, float LwBottom, 
                     float TBottom, float LwTop, float TTop, float Cold, 
                     int NSoilLayers, float *Moist, float SoilTSurf,
                     float *Temp, float GroundHeat, float Runoff);

/*****************************************************************************
  MakeModelStateBin()

  Make a model state file to use as an initialization for DHSVM

  The state variables for DHSVM include the following variables:

    - Rain interception for each vegetation layer
    - Snow interception for top vegetation layer
    - Temporary interception (used for mass release algorithm)

    - Snow pack conditions:
      - presence/absence
      - number of days since last snowfall (used in albedo calculation)
      - snow water equivalent
      - for each layer of the snow pack:
        - liquid water content
        - temperature
      - cold content

    - Soil conditions:
      - for each soil layer:
        - soil moisture (also for the layer below the deepest root zone)
        - temperature
      - surface temperature
      - ground heat storage
*****************************************************************************/
int main(int argc, char *argv[])
{
  char InfoFileName[MAXSTRING+1];
  char Path[MAXSTRING+1];
  DATE Day;
  int NY;
  int NX;
  int NVegLayers;
  int i;
  float RainInt[MAXLAYERS];
  float SnowInt[MAXLAYERS];
  float TempIntStorage;
  unsigned char SnowMask;
  unsigned short LastSnow;
  float Swq;
  float LwBottom;
  float TBottom;
  float LwTop;
  float TTop;
  float Cold;
  int NSoilLayers;
  float Moist[MAXLAYERS];
  float SoilTSurf;
  float Temp[MAXLAYERS];
  float GroundHeat;
  float Runoff;
  int Junk;
  FILE *InfoFile;

  if (argc != 2) {
    fprintf(stderr, "Usage: MakeModelState <infofile>\n");
    fprintf(stderr, "The info file MUST contain the following information:\n");
    fprintf(stderr, " - path for output file\n");
    fprintf(stderr, " - date for the model state, in mm/dd/yyyy-hh\n");
    fprintf(stderr, " - number of rows (ny) and number of columns (nx)\n");
    fprintf(stderr, " - maximum number of vegetation layers\n");
    fprintf(stderr, " - rain interception in m for each vegetation layer\n");
    fprintf(stderr, " - snow interception in m for top vegetation layer\n");
    fprintf(stderr, " - snow cover mask\n");
    fprintf(stderr, " - number of days since last snow fall\n");
    fprintf(stderr, " - snow water equivalent in m\n");
    fprintf(stderr, " - liquid water content in m of bottom layer of snowpack\n");
    fprintf(stderr, " - temperature in C of bottom layer of snow pack\n");
    fprintf(stderr, " - liquid water content in m of top layer of snowpack\n");
    fprintf(stderr, " - temperature in C of top layer of snow pack\n");
    fprintf(stderr, " - cold content of snow pack\n");
    fprintf(stderr, " - maximum number of root zone layers\n");
    fprintf(stderr, " - volumetric soil moisture content for each layer\n");
    fprintf(stderr, "   (including the layer below the lowest root zone layer)\n");
    fprintf(stderr, " - temperature in C at soil surface\n");
    fprintf(stderr, " - soil temperature in C for each root zone layer\n");
    fprintf(stderr, " - ground heat storage\n");
    fprintf(stderr, " - runoff\n");
    exit(1);
  }
  strcpy(InfoFileName, argv[1]);
  
  if (!(InfoFile = fopen(InfoFileName, "r"))) {
    fprintf(stderr, "Canot open info file %s\n", InfoFileName);
    exit(1);
  }

  fscanf(InfoFile, "%s", Path);
  ScanDate(InfoFile, &Day);
  fscanf(InfoFile, "%d %d", &NY, &NX);
  fscanf(InfoFile, "%d", &NVegLayers);
  for (i = 0; i < NVegLayers; i++) 
    fscanf(InfoFile, "%f", &RainInt[i]);
  fscanf(InfoFile, "%f", &SnowInt[0]);
  for (i = 1; i < NVegLayers; i++) 
    SnowInt[i] = 0;
  TempIntStorage = 0.0;
  fscanf(InfoFile, "%d", &Junk);
  SnowMask = (unsigned char) Junk;
  fscanf(InfoFile, "%d", &Junk);
  LastSnow = (unsigned short) Junk;
  fscanf(InfoFile, "%f", &Swq);
  fscanf(InfoFile, "%f", &LwBottom);
  fscanf(InfoFile, "%f", &TBottom);
  fscanf(InfoFile, "%f", &LwTop);
  fscanf(InfoFile, "%f", &TTop);
  fscanf(InfoFile, "%f", &Cold);
  fscanf(InfoFile, "%d", &NSoilLayers);
  for (i = 0; i <= NSoilLayers; i++)
    fscanf(InfoFile, "%f", &Moist[i]);
  fscanf(InfoFile, "%f", &SoilTSurf);
  for (i = 0; i < NSoilLayers; i++)
    fscanf(InfoFile, "%f", &Temp[i]);
  fscanf(InfoFile, "%f", &GroundHeat);
  fscanf(InfoFile, "%f", &Runoff);

  //  InitErrorMessage();

  StoreModelState(Path, Day, NY, NX, NVegLayers, RainInt, SnowInt,
                  TempIntStorage, SnowMask,
                  LastSnow, Swq, LwBottom, TBottom, LwTop, TTop, Cold,
                  NSoilLayers, Moist, SoilTSurf, Temp, GroundHeat, Runoff);

  return 0;
}

/*****************************************************************************
  StoreModelState()

  Store the current state of the model.

  The state variables for DHSVM include the following variables:

    - Canopy interception for each vegetation layer

    - Snow pack conditions:
      - presence/absence
      - number of days since last snowfall (used in albedo calculation)
      - snow water equivalent
      - for each layer of the snow pack:
        - liquid water content
        - temperature
      - cold content

    - Soil conditions:
      - for each soil layer:
        - soil moisture (also for the layer below the deepest root zone)
        - temperature
      - surface temperature
      - ground heat storage
*****************************************************************************/

void StoreModelState(char *Path, DATE Current, int NY, int NX, 
                     int NVegLayers, float *RainInt, float *SnowInt,
                     float TempIntStorage, unsigned char SnowMask, 
                     unsigned short LastSnow, float Swq, float LwBottom, 
                     float TBottom, float LwTop, float TTop, float Cold, 
                     int NSoilLayers, float *Moist, float SoilTSurf,
                     float *Temp, float GroundHeat, float Runoff)
{
  const char *Routine = "StoreModelState";
  char Str[NAMESIZE+1];
  char DataLabel[MAXSTRING+1];
  char FileLabel[MAXSTRING+1];
  char FileName[NAMESIZE+1];
  char Units[MAXSTRING+1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  void *Array;

  /* print a message to stdout that state is being stored */

  printf("Storing model state\n");
 
  /* Store the canopy interception */

  sprintf(Str, "%02d.%02d.%04d.%02d.00.00", Current.Month, Current.Day, 
	  Current.Year, Current.Hour);
  MakeFileNameBin(Path, "Interception.State.", Str, FileName); 
  strcpy(FileLabel, "Interception storage for each vegetation layer");
  CreateFileBin(FileName, FileLabel);

  if (!(Array = (float *) calloc(NY * NX, sizeof(float)))) {
    perror(Routine);
    exit(1);
  }

  for (i = 0; i < NVegLayers; i++) {
    for (y = 0; y < NY; y++) 
      for (x = 0; x < NX; x++) 
        ((float *)Array)[y*NX + x] = RainInt[i];
    NumberType = NT_FLOAT32;
    sprintf(DataLabel, "Rain interception for vegetation Layer %d", i+1);
    strcpy(Units, "mm");
    Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                  FileName);
  }
  
  for (i = 0; i < NVegLayers; i++) {
    for (y = 0; y < NY; y++) 
      for (x = 0; x < NX; x++) 
        ((float *)Array)[y*NX + x] = SnowInt[i];
    NumberType = NT_FLOAT32;
    sprintf(DataLabel, "Snow interception for vegetation Layer %d", i+1);
    strcpy(Units, "mm");
    Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                  FileName);
  }

  for (y = 0; y < NY; y++) 
    for (x = 0; x < NX; x++) 
      ((float *)Array)[y*NX + x] = TempIntStorage;
  NumberType = NT_FLOAT32;
  sprintf(DataLabel, "Temporary snow interception for vegetation Layer %d", i+1);
  strcpy(Units, "mm");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                   FileName);
  
  free(Array);

  /* Store the snow pack conditions */
  
  MakeFileNameBin(Path, "Snow.State.", Str, FileName); 
  strcpy(FileLabel, "Snow pack moisture and temperature state");
  CreateFileBin(FileName, FileLabel);
  
  NumberType = NT_FLOAT32;
  if (!(Array = (float *) calloc(NY * NX, sizeof(SizeOfNumberType(NumberType))))) {
    perror(Routine);
    exit(1);
  }
  
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = (float)SnowMask;
  
  strcpy(DataLabel, "Snow Cover Mask");
  strcpy(Units, "");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);


  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = (float)LastSnow;

  strcpy(DataLabel, "Number of Days Since Last Snowfall");
  strcpy(Units, "Days");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = Swq;

  strcpy(DataLabel, "Snow Water Equivalent");
  strcpy(Units, "mm");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = LwBottom;

  strcpy(DataLabel, "Liquid Water Content of Bottom Layer");
  strcpy(Units, "mm");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);
    
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = TBottom;

  strcpy(DataLabel, "Temperature of Bottom Layer");
  strcpy(Units, "C");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = LwTop;

  strcpy(DataLabel, "Liquid Water Content of Surface Layer");
  strcpy(Units, "mm");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = TTop;

  strcpy(DataLabel, "Temperature of Surface Layer");
  strcpy(Units, "mm");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = Cold;;

  strcpy(DataLabel, "Cold Content of Snow Pack");
  strcpy(Units, "");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  free(Array);

  /* Store the soil conditions */ 
  
  MakeFileNameBin(Path, "Soil.State.", Str, FileName);
  strcpy(FileLabel, "Soil moisture and temperature state");
  CreateFileBin(FileName, FileLabel);
  
  if (!(Array = (float *) calloc(NY * NX, sizeof(float)))) {
    perror(Routine);
    exit(1);
  }
  
  for (i = 0; i < NSoilLayers+1; i++) {

    for (y = 0; y < NY; y++) 
      for (x = 0; x < NX; x++) 
        ((float *)Array)[y*NX + x] = Moist[i];
    
    NumberType = NT_FLOAT32;
    sprintf(DataLabel, "Soil Moisture Content of Layer %d", i);
    strcpy(Units, "");
    Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                  FileName);
  }
  
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = SoilTSurf;
  
  NumberType = NT_FLOAT32;
  strcpy(DataLabel, "Temperature of Soil Surface");
  strcpy(Units, "C");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  for (i = 0; i < NSoilLayers; i++) {

    for (y = 0; y < NY; y++) 
      for (x = 0; x < NX; x++) 
        ((float *)Array)[y*NX + x] = Temp[i];
    
    NumberType = NT_FLOAT32;
    sprintf(DataLabel, "Soil Temperature of Layer %d", i);
    strcpy(Units, "");
    Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                  FileName);
  }
  
  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = GroundHeat;
  
  NumberType = NT_FLOAT32;
  strcpy(DataLabel, "Ground Heat Storage");
  strcpy(Units, "");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  for (y = 0; y < NY; y++)
    for (x = 0; x < NX; x++)
      ((float *)Array)[y*NX + x] = Runoff;
  
  NumberType = NT_FLOAT32;
  strcpy(DataLabel, "Runoff");
  strcpy(Units, "");
  Write2DMatrixBin(NY, NX, NumberType, DataLabel, Units, Array, 
                FileName);

  free(Array);
}

