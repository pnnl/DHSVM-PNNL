/* $Id: InitGlacierMap.c,v 1.1.1.1 2002/09/24 04:58:50 nijssen Exp $     
 */
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "constants.h"
#include "glacier.h"

/*****************************************************************************
  Function name: InitGlacierMap()

  Purpose      : Initialize the glacier information for each pixel in the basin
  Author       : Bibi S. Naz
  Required     :
  MAPSIZE Map        - Size and location of the model area
      
  GLPIX ***GlacierMap - Address of array with glacier information 
  Returns      : void

  Modifies     :
    Values stored at the locations pointed to by GlacierMap

  Comments     :
*****************************************************************************/
void InitGlacierMap(MAPSIZE * Map, GLPIX *** GlacierMap)
{
 
  const char *Routine = "InitGlacierMap";
  int y;			/* counter */
  
  extern int  nx;
  extern int  ny;
  extern int  N;
  extern double  dx;
  extern double  dy;
  N  = Map->NX * Map->NY;
  nx = Map->NX;
  ny = Map->NY;
  dx = Map->DX;
  dy = Map->DY;
 
  printf("Initializing glacier map\n");

  if (!(*GlacierMap = (GLPIX **) calloc(Map->NY, sizeof(GLPIX *))))
    ReportError((char *) Routine, 1);

  for (y = 0; y < Map->NY; y++) {
    if (!((*GlacierMap)[y] = (GLPIX *) calloc(Map->NX, sizeof(GLPIX))))
      ReportError((char *) Routine, 1);
  }
  
  
  
}
