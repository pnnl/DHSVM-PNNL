/*
 * SUMMARY:      Map2D.h
 * USAGE:        Routines for dealing with 2D Map input and output
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    October 2018
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-10-19 07:58:44 d3g096
 * COMMENTS:
 */

#ifndef _Map2D_h_
#define _Map2D_h_

#include <stdio.h>
#include "settings.h"
#include "MapSize.h"
#include "fileio.h"

/* typedef void* InputMap2D; */
/* typedef void* OutputMap2D; */

/* InputMap2D *InputMap2DAlloc(const char* fname, const char* vname,  */
/*                             int NumberType, MAPSIZE *Map, int mirror); */
/* int InputMap2DOpen(InputMap2D *map2d); */
/* int InputMap2DRead(InputMap2D *map2d, int index, void *ldata); */
/* int InputMap2DClose(InputMap2D *map2d); */
/* void InputMap2DFree(InputMap2D *map2d); */

void Map2DInit(int FileFormat);


int Read2DMatrix(const char *FileName, void *Matrix, int NumberType, 
                 MAPSIZE *Map, int NDataSet, const char *VarName, int index);

int Read2DMatrixAll(const char *FileName, void *Matrix, int NumberType, 
                    MAPSIZE *Map, int NDataSet, const char *VarName, int index);


#endif 



