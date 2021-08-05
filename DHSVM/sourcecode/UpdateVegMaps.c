/*
 * SUMMARY:      UpdateVegMaps() - Update Vegetation Maps at selected date
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Zhuoran Duan
 * ORG:          Pacific Northwest National Laboratory, Hydrology Group
 * E-MAIL:       zhuoran.duan@pnnl.gov
 * ORIG-DATE:    Mar-2020
 * DESCRIPTION:  Update vegetation map for type, fractional cover, LAI, and height
 * DESCRIP-END.
 * FUNCTIONS:    UpdateVegMaps()
 *               InitTopoMap()
 *               InitSoilMap()
 *               InitVegMap()
 * COMMENTS:
 * 
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
#include "Calendar.h"


/*******************************************************************************
  Function name: InitVegUpdate()

  Purpose      : Initialize the vegetation updates.  This information is in the
         [VEGETATION] section of the input file

  Required     :
    LISTPTR Input         - Linked list with input strings
    int NUpdate           - Number of vegetation layers to update
    DATE **DUpdate        - Array with update dates

  Returns      : void

  Modifies     : DUpdates and its members

  Comments     :
*****************************************************************************/
void InitVegUpdate(LISTPTR Input, int NUpdate, DATE ** DUpdate)
{
  char *Routine = "InitVegUpdate";
  int i;			/* counter */
  char KeyName[BUFSIZE + 1];
  char *KeyStr = "UPDATE DATE";
  char *SectionName = "VEGETATION";
  char VarStr[BUFSIZE + 1];

  if (!(*DUpdate = (DATE *)calloc(NUpdate, sizeof(DATE))))
    ReportError(Routine, 1);

  for (i = 0; i < NUpdate; i++) {  
    sprintf(KeyName, "%s %d", KeyStr, i + 1);
    GetInitString(SectionName, KeyName, "", VarStr,
      (unsigned long)BUFSIZE, Input);
    if (!SScanDate(VarStr, &((*DUpdate)[i])))
      ReportError(KeyName, 51);
  }
}

/*****************************************************************************
  IsVegDate()
*****************************************************************************/
uchar IsVegDate(DATE *Current, DYNAVEG *DVeg){
  //char *Routine = "IsVegDate";
  int i;			/* counter */

  for (i = 0; i < DVeg->NUpdate; i++) {
    if (IsEqualTime(Current, &(DVeg->DUpdate[i])))
    {
      return TRUE;
      break;
    }
  } 
  return FALSE;
}

/*****************************************************************************
  Function name: UpdateVegMap()

  Purpose      : Update Vegetation Maps at user defined date, this updates:
                  VEGPIX *** VegMap
  Required     :
    VEGPIX *** VegMap       - Updates vegetation information
    DYNAVEG *DVeg           - Dynamic Vegetaion, with input path and dates 

  Returns      : void

  Modifies     : (see list of required above)

  Comments     :
*****************************************************************************/
void UpdateVegMap(DATE *Current, OPTIONSTRUCT * Options, LISTPTR Input, MAPSIZE * Map,
                LAYER *Veg, VEGPIX *** VegMap, VEGTABLE *VType, DYNAVEG *DVeg)
{  
  const char *Routine = "UpdateVegMap";
  char VarName[BUFSIZE + 1];
  char Path[BUFSIZE + 1];
  char FileName[NAMESIZE + 1];
  char Str[NAMESIZE + 1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int j;			/* counter */
  int flag;
  int NumberType;		/* number type */
  unsigned char *Type;		/* Vegetation type */
  float *FC = NULL;		/* Vegetation Fractional Coverage */
  float *LAIMonthly= NULL; /* Vegetation Leaf Area Index, monthly */
  int NSet; /*Counter for LAI map month*/

  
  strcpy(Path, DVeg->DynaVegPath);

  printf("Updating vegetation maps\n");

  /*Vegetation Maps must be named as Veg.XXX.MM.DD.YYYY.HH.MM.SS.ext*/
  sprintf(Str, "%02d.%02d.%02d.%02d.%02d.%02d", Current->Month, Current->Day,
  Current->Year, Current->Hour, Current->Min, Current->Sec);
    
  /*Update Vegetation Type Map*/
  sprintf(FileName, "%sVegetation.Type.%s%s", Path, Str, fileext);
  printf("updating file %s\n",FileName);
  /* Read the vegetation type */
  GetVarName(005, 0, VarName);
  GetVarNumberType(005, &NumberType);
  if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(FileName, Type, NumberType, Map, 0, VarName, 0);
  
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Veg = Type[i];
        //printf("veg type at x %d y %d is %d \n",x,y,Type[i]);
        //(*VegMap)[y][x].Tcanopy = 0.0;
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Veg = Type[i];
        //(*VegMap)[y][x].Tcanopy = 0.0;
      }
    }
  }
  else ReportError((char *)Routine, 57);

  free(Type);

  /*Update Vegetation Fractional Cover Map*/
  sprintf(FileName, "%sVegetation.FC.%s%s", Path, Str, fileext);
  printf("updating file %s\n",FileName);
  /* Read the vegetation fractional coverage map */
  GetVarName(010, 0, VarName);
  GetVarNumberType(010, &NumberType);

  if (strncmp(FileName, "none", 4)) {
    printf("Spatial fractional cover map provided, reading FC from map\n");
    if (!(FC = (float *)calloc(Map->NX * Map->NY,
      SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(FileName, FC, NumberType, Map, 0, VarName, 0);

    if ((Options->FileFormat == NETCDF && flag == 0)
      || (Options->FileFormat == BIN))
    {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          /*Allocate Memory*/
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            if (FC[i] > 0.0)
              (*VegMap)[y][x].Fract[0] = FC[i];
            else
              (*VegMap)[y][x].Fract[0] = VType[(*VegMap)[y][x].Veg - 1].Fract[0];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[1] = 1.0;
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[0] = 1.0;
          }

          //printf("veg fc at x %d y %d is %f \n",x,y,(*VegMap)[y][x].Fract[0]);

        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
          /*Allocate memory*/
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {  
            if (FC[i] > 0.0)
              (*VegMap)[y][x].Fract[0] = FC[i];
            else
            /* If value from the fractional cover map is NaN, then read value from attribute table*/
              (*VegMap)[y][x].Fract[0] = VType[(*VegMap)[y][x].Veg - 1].Fract[0];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[1] = 1.0;
          }
          else{
            if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].Fract[0] = 1.0;	   
          }
        }
      }
    }
    else ReportError((char *)Routine, 57);
    free(FC);
  }
  else{
    //printf("Vegetation fractional coverage created from vegetation table\n");
    ReportError((char *)Routine, 57);
  }

  /*Calculate Vf */
  for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if ( VType[(*VegMap)[y][x].Veg - 1].NVegLayers >0) 
          (*VegMap)[y][x].Vf = (*VegMap)[y][x].Fract[0] * VType[(*VegMap)[y][x].Veg - 1].VfAdjust;
      }
  }


  /*Update Vegetation LAI Map*/
  sprintf(FileName, "%sVegetation.LAI.%s%s", Path, Str, fileext);
  printf("updating file %s\n",FileName);

  /* Read the vegetation LAI map */
  /*Instead of reading LAI data by month, read data all together as Bill suggested*/
    
  GetVarName(011, 0, VarName);
  GetVarNumberType(011, &NumberType);
 
  if (strncmp(FileName, "none", 4)) {
    printf("Spatial LAI provided, reading LAI from map\n");
  /*Read data monthy by month*/
  for (NSet = 0; NSet < 12; NSet++) {
    if (!(LAIMonthly = (float *)calloc(Map->NX * Map->NY,
      SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrix(FileName, LAIMonthly, NumberType, Map, NSet, VarName, 0);
    
    printf("begining month %d\n",NSet);
    
    if ((Options->FileFormat == NETCDF && flag == 0)
      || (Options->FileFormat == BIN))
    {
      for (y = 0, i = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++, i++) {
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
           if (LAIMonthly[i] > 0.0)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
            else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
           
            if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory  == TRUE )
              (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet];
          }
          else{
            if ( VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
          }

          //printf("veg lai at x %d y %d is %f \n",x,y,(*VegMap)[y][x].LAIMonthly[0][NSet]);

        }
      }
    }
    else if (Options->FileFormat == NETCDF && flag == 1) {
      for (y = Map->NY - 1, i = 0; y >= 0; y--) {
        for (x = 0; x < Map->NX; x++, i++) {
        
          if ( VType[(*VegMap)[y][x].Veg - 1].OverStory == TRUE) {
            if (LAIMonthly[i] > 0.0)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = LAIMonthly[i];
            else
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
            /*If understory exists, set default understory FC=1.0*/
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].LAIMonthly[1][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[1][NSet];
          }
          else{
            if (VType[(*VegMap)[y][x].Veg - 1].UnderStory == TRUE)
              (*VegMap)[y][x].LAIMonthly[0][NSet] = VType[(*VegMap)[y][x].Veg - 1].LAIMonthly[0][NSet];
          }
        }
      }
    }
    else ReportError((char *)Routine, 57);  

    free(LAIMonthly);     
  }   
  }
  else{
    printf("No spatial LAI provided, generating from vegetation table\n");
    ReportError((char *)Routine, 57);
  }
}