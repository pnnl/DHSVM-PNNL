/*
 * SUMMARY:      ReadMetRecord.c - Read station meteorological data
 * USAGE:        Part of DHSVM
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Read station meteorological data
 * DESCRIP-END.
 * FUNCTIONS:    ReadMetRecord()
 * COMMENTS:
 * $Id: ReadMetRecord.c,v 1.4 2003/07/01 21:26:22 olivier Exp $     
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

#define MAXMETVARS    21	/* Maximum Number of meteorological variables 
				   to read.  Hack to be replaced by something
				   better */

/*****************************************************************************
  ReadMetRecord()
*****************************************************************************/
void ReadMetRecord(OPTIONSTRUCT * Options, DATE * Current, int NSoilLayers,
		   FILES * InFile, unsigned char IsWindModelLocation,
		   MET * MetRecord)
{
  DATE MetDate;			/* Date of meteorological record */
  float Array[MAXMETVARS];	/* Temporary storage of met variables */
  int i;
  int NMetVars;			/* Number of meteorological variables to 
				   read */
  NMetVars = 5;
  /* these are - in order: 
     air temp,
     wind,
     humidity
     shortwave (total i.e. direct+diffuse)
     longwave */

  if (Options->HeatFlux == TRUE)
    NMetVars += NSoilLayers;
  /* expect to see temperature for each soil layer */
  if (Options->PrecipType == STATION)
    NMetVars++;
  if (Options->PrecipLapse == VARIABLE)
    NMetVars++;
  if (Options->TempLapse == VARIABLE)
    NMetVars++;
  if (IsWindModelLocation)
    NMetVars++;

  if (!ScanDate(InFile->FilePtr, &MetDate))
    ReportError(InFile->FileName, 23);

  while (!IsEqualTime(&MetDate, Current) && !feof(InFile->FilePtr)) {
    if (ScanFloats(InFile->FilePtr, Array, NMetVars) != NMetVars)
      ReportError(InFile->FileName, 5);
    if (!ScanDate(InFile->FilePtr, &MetDate))
      ReportError(InFile->FileName, 23);
  }

  if (!IsEqualTime(&MetDate, Current)) {
    if (DEBUG) {
      printf("Metfile: ");
      PrintDate(&MetDate, stdout);
      printf("Current: ");
      PrintDate(Current, stdout);
    }
    ReportError(InFile->FileName, 28);
  }

  if (ScanFloats(InFile->FilePtr, Array, NMetVars) != NMetVars)
    ReportError(InFile->FileName, 5);

  MetRecord->Tair = Array[0];
  MetRecord->Wind = Array[1];
  MetRecord->Rh = Array[2];
  if (MetRecord->Rh < 0.0 || MetRecord->Rh > 100.0) {
    printf("warning: RH out of bounds: %s\n", InFile->FileName);
    if (MetRecord->Rh < 0.0)
      MetRecord->Rh = 0.0;
    if (MetRecord->Rh > 100.0)
      MetRecord->Rh = 100.0;
  }
  MetRecord->Sin = Array[3];
  if (MetRecord->Sin > 1380.0) {
    printf("warning: Shortwave out of bounds: %s\n", InFile->FileName);
    MetRecord->Sin = 1380.0;
  }
  if (MetRecord->Sin < 0.0) {
    printf("Warning: Negative Shortwave, setting to zero: %s\n",
	   InFile->FileName);
    MetRecord->Sin = 0.0;
  }
  MetRecord->Lin = Array[4];
  if (MetRecord->Lin < 0.0 || MetRecord->Lin > 1800.0) {
    printf("warning: Longwave out of bounds: %s\n", InFile->FileName);
  }

  i = 0;
  if (Options->HeatFlux == TRUE)
    for (i = 0; i < NSoilLayers; i++)
      MetRecord->Tsoil[i] = Array[5 + i];

  if (Options->PrecipType == STATION) {
    MetRecord->Precip = Array[5 + i];
    if (MetRecord->Precip < 0) {
      printf("Warning: negative precip %s \n", InFile->FileName);
      MetRecord->Precip = 0.0;
    }
    i++;

  }
  else
    MetRecord->Precip = 0.0;

  if (Options->PrecipLapse == VARIABLE) {
    MetRecord->PrecipLapse = Array[5 + i];
    i++;
  }
  else
    MetRecord->PrecipLapse = PRECIPLAPSE;

  if (Options->TempLapse == VARIABLE) {
    MetRecord->TempLapse = Array[5 + i];
    i++;
  }
  else
    MetRecord->TempLapse = TEMPLAPSE;

  if (IsWindModelLocation) {
    MetRecord->WindDirection = (int) Array[5 + i];
    i++;
  }
  else
    MetRecord->WindDirection = NOT_APPLICABLE;

}
