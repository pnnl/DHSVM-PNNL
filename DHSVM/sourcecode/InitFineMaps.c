/*
 * SUMMARY:      InitFineMaps() - Initialize fine resolution coverages
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Laura C. Bowling/Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       colleen@hydro.washington.edu
 * ORIG-DATE:    Oct-03
 * Last Change:  Mon Oct 27 14:40:44 2003 by Colleen O. Doten <colleen@hydro.washington.edu>
 * DESCRIPTION:  Initialize terrain coverages for fine resolution map
 * DESCRIP-END.
 * FUNCTIONS:    InitFineMaps()
 * COMMENTS:     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "constants.h"
#include "getinit.h"
#include "varid.h"
#include "sizeofnt.h"
#include "slopeaspect.h"

void CalcTopoIndex (MAPSIZE *Map, FINEPIX ***FineMap, TOPOPIX **TopoMap);

/*****************************************************************************
  InitFineMaps()
*****************************************************************************/
void InitFineMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, 
		     LAYER *Soil, TOPOPIX ***TopoMap, SOILPIX ***SoilMap, 
		     FINEPIX ****FineMap)
{
 
  const char *Routine = "InitFineMaps";
  char VarName[BUFSIZE+1];	/* Variable name */
  int i, k, x, y;		/* Counters */
  int ii, jj, xx, yy, xy;            /* Counters */
  int NumberType;		/* Number type of data set */
  float *Elev;                   /* Surface elevation */
  int MASKFLAG;
  unsigned char *Mask = NULL;          /* Fine resolution mask */

  STRINIENTRY StrEnv[] = {
    {"FINEDEM", "DEM FILE"        , ""  , ""},
    {"FINEDEM", "MASK FILE"        , ""  , ""},
    {NULL       , NULL            , ""  , NULL}
  };
  
  printf("Initializing mass wasting resolution maps\n");
  
  /* Process the [FINEDEM] section in the input file */
  
  /* Read the key-entry pair for the dem from the input file */
  GetInitString(StrEnv[0].SectionName, StrEnv[0].KeyName, StrEnv[0].Default,
		StrEnv[0].VarStr, (unsigned long) BUFSIZE, Input);
  if (IsEmptyStr(StrEnv[0].VarStr))
    ReportError(StrEnv[0].KeyName, 51);
  
  /* Read the elevation dataset */ 
  
  GetVarName(001, 0, VarName);
  GetVarNumberType(001, &NumberType);
  if (!(Elev = (float *) calloc(Map->NXfine * Map->NYfine, 
				SizeOfNumberType(NumberType)))) 
    ReportError((char *) Routine, 1);
  Read2DMatrix(StrEnv[demfile].VarStr, Elev, NumberType, Map->NYfine, Map->NXfine, 0,
	       VarName);
  
  /* Read the key-entry pair for the mask from the input file */
  GetInitString(StrEnv[1].SectionName, StrEnv[1].KeyName, StrEnv[1].Default,
		StrEnv[1].VarStr, (unsigned long) BUFSIZE, Input);
  if (IsEmptyStr(StrEnv[1].VarStr)) {
    printf("\nWARNING: Fine resolution mask not provided, will be set equal to \n");
    printf("coarse resolution mask.\n\n");
    MASKFLAG = FALSE;
  }
  else {
    printf("fine mask = %s\n",StrEnv[1].VarStr);
    /* Read the mask */
    GetVarName(002, 0, VarName);
    GetVarNumberType(002, &NumberType);
    if (!(Mask = (unsigned char *) calloc(Map->NXfine * Map->NYfine,
					  SizeOfNumberType(NumberType))))
      ReportError((char *) Routine, 1);
    Read2DMatrix(StrEnv[maskfile].VarStr, Mask, NumberType, Map->NYfine, Map->NXfine, 0,
		 VarName);
    MASKFLAG = TRUE;
    
  }
  
  /* Assign the attributes to the correct map pixel */
  if (!(*FineMap = (FINEPIX ***) calloc(Map->NYfine, sizeof(FINEPIX **))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NYfine; y++) {
    if (!((*FineMap)[y] = (FINEPIX **) calloc(Map->NXfine, sizeof(FINEPIX *))))
      ReportError((char *) Routine, 1);
  }
  // Only allocate a FINEPIX structure for a fine grid cell if that grid cell
  // is in the coarse grid mask
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) { 
      if (INBASIN((*TopoMap)[y][x].Mask)) {
	for (ii=0; ii< Map->DY/Map->DMASS; ii++) { 
 	  for (jj=0; jj< Map->DX/Map->DMASS; jj++) { 
	    yy = (int) y*Map->DY/Map->DMASS + ii; 
 	    xx = (int) x*Map->DX/Map->DMASS + jj; 
            if (!((*FineMap)[yy][xx] = (FINEPIX *) malloc(sizeof(FINEPIX)))) {
	      printf("error allocating FineMap[%d][%d]\n",yy,xx);
              ReportError((char *) Routine, 1);
            }
	  }
	}
      }
    }
  }
  
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) { 
      if (INBASIN((*TopoMap)[y][x].Mask)) {
	for (ii=0; ii< Map->DY/Map->DMASS; ii++) { 
 	  for (jj=0; jj< Map->DX/Map->DMASS; jj++) { 
	    yy = (int) y*Map->DY/Map->DMASS + ii; 
 	    xx = (int) x*Map->DX/Map->DMASS + jj; 
	    xy = (int) yy*Map->NXfine + xx;
	    (*(*FineMap)[yy][xx]).Dem  = Elev[xy]; 
	  }
	}
      }
    }
  }
  
  free(Elev);
  
  if(MASKFLAG == TRUE){
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) { 
	if (INBASIN((*TopoMap)[y][x].Mask)) {
	  for (ii=0; ii< Map->DY/Map->DMASS; ii++) { 
	    for (jj=0; jj< Map->DX/Map->DMASS; jj++) { 
	      yy = (int) y*Map->DY/Map->DMASS + ii; 
	      xx = (int) x*Map->DX/Map->DMASS + jj; 
	      xy = (int) yy*Map->NXfine + xx;
	      (*(*FineMap)[yy][xx]).Mask  = Mask[xy]; 
	    }
	  }
	}
      }
    }
  }

 free(Mask);  
  /* Create fine resolution mask, sediment and bedrock maps.
   Initialize other variables*/
  
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
	  if (INBASIN((*TopoMap)[y][x].Mask)) {
	    yy = (int) y*Map->DY/Map->DMASS + ii;
	    xx = (int) x*Map->DX/Map->DMASS + jj;
	    if(MASKFLAG==FALSE)
	      (*(*FineMap)[yy][xx]).Mask = (*TopoMap)[y][x].Mask;
	    /*   else { */
	      // Don't allow fine mask to extend beyond edges of coarse mask.
	      // This means that failures may still try to leave the basin if the coarse 
	      // mask isn't wide enough because some of the drainage area may be cropped.
	/*       if (!INBASIN((*TopoMap)[y][x].Mask)) */
/* 		(*(*FineMap)[yy][xx]).Mask = (*TopoMap)[y][x].Mask; */
	 /*    } */
	    (*(*FineMap)[yy][xx]).bedrock = (*(*FineMap)[yy][xx]).Dem - (*SoilMap)[y][x].Depth;
	    (*(*FineMap)[yy][xx]).sediment = (*SoilMap)[y][x].Depth;
	    (*(*FineMap)[yy][xx]).SatThickness = 0.;
	    (*(*FineMap)[yy][xx]).DeltaDepth = 0.;
	    (*(*FineMap)[yy][xx]).Probability = 0.;
	    (*(*FineMap)[yy][xx]).MassWasting = 0.;
	    (*(*FineMap)[yy][xx]).MassDeposition = 0.;
	    (*(*FineMap)[yy][xx]).SedimentToChannel = 0.;
	    (*(*FineMap)[yy][xx]).TopoIndex = 0.;
	  }
	}
      }
    }
  }
  
  Map->NumFineIn = (Map->DX/Map->DMASS) * (Map->DY/Map->DMASS);

 /* NumCellsfine is used in CalcTopoIndex.  The topo index is calculated for every 
     fine cell within the boundary of the coarse mask, so this number may exceed the 
     number of pixels within the fine resolution mask.  */

  Map->NumCellsfine = Map->NumCells*Map->NumFineIn;

  printf("Basin has %d active pixels in the mass wasting resolution map\n",
	 Map->NumCellsfine);
  
  /* Calculate the topographic index */
  CalcTopoIndex(Map, *FineMap, *TopoMap);
  
  for (y = 0; y < Map->NY; y++) {
    for (x  = 0; x < Map->NX; x++) {
      if (INBASIN((*TopoMap)[y][x].Mask)) {
	if (!((*TopoMap)[y][x].OrderedTopoIndex = (ITEM *) calloc(Map->NumFineIn, sizeof(ITEM))))
	  ReportError((char *) Routine, 1);
      }
    }
  }
  
  for (y = 0; y < Map->NY; y++) {
    for (x  = 0; x < Map->NX; x++) {
      if (INBASIN((*TopoMap)[y][x].Mask)) {
	k = 0;
	for(ii=0; ii< Map->DY/Map->DMASS; ii++) {
	  for(jj=0; jj< Map->DX/Map->DMASS; jj++) {
	    yy = (int) y*Map->DY/Map->DMASS + ii;
	    xx = (int) x*Map->DX/Map->DMASS + jj;
	    (*TopoMap)[y][x].OrderedTopoIndex[k].Rank = (*(*FineMap)[yy][xx]).TopoIndex;
	    (*TopoMap)[y][x].OrderedTopoIndex[k].y = yy;
	    (*TopoMap)[y][x].OrderedTopoIndex[k].x = xx;
	    k++;
	  }
	}
       	quick((*TopoMap)[y][x] .OrderedTopoIndex, Map->NumFineIn);
      }
    }
  }
}
   
