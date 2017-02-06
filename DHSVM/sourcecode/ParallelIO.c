/*
 * SUMMARY:      ParallelIO.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    January 2017
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2017-02-06 08:31:37 d3g096
 * COMMENTS:
 */

/******************************************************************************/
/*                                INCLUDES                                    */
/******************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <ga.h>

#include "data.h"
#include "sizeofnt.h"
#include "ParallelDHSVM.h"
#include "fileio.h"

/* global function pointers */
extern void (*CreateMapFileFmt) (char *FileName, ...);
extern int (*Read2DMatrixFmt) (char *FileName, void *Matrix, int NumberType, int NY, int NX, int NDataSet, ...);
extern int (*Write2DMatrixFmt) (char *FileName, void *Matrix, int NumberType, int NY, int NX, ...);

static const gaXdim = 1;
static const gaYdim = 0;

/******************************************************************************/
/*                            CreateMapFile                                   */
/******************************************************************************/
void
CreateMapFile(char *FileName, char *FileLabel, MAPSIZE *Map)
{
  int me = ParallelRank();
  if (me == 0) {
    CreateMapFileFmt(FileName, FileLabel, Map);
  }
  ParallelBarrier();
}

/******************************************************************************/
/*                           GlobalMapSize                                    */
/******************************************************************************/
void
GlobalMapSize(MAPSIZE *Map, int *globalNX, int *globalNY)
{
  int ndim;
  int dims[GA_MAX_DIM];
  int type;

  NGA_Inquire(Map->dist, &type, &ndim, &dims[0]);
  assert(ndim == 2);
  *globalNX = dims[gaXdim];
  *globalNY = dims[gaYdim];
  return;
}


/******************************************************************************/
/*                                GA_Type                                     */
/******************************************************************************/
static int
GA_Type(int NumberType)
{
  int gatype;
  switch (NumberType) {
  case NC_INT:
    gatype = C_INT;
    break;
  case NC_FLOAT:
    gatype = C_FLOAT;
    break;
  case NC_DOUBLE:
    gatype = C_DBL;
    break;
  case NC_BYTE:
  case NC_CHAR:
    gatype = C_CHAR;
    break;
  default:
    gatype = 0;
    ReportError("GAType", 40);
    break;
  }
  return gatype;
}

/******************************************************************************/
/*                          GA_Duplicate_type                                 */
/******************************************************************************/
/** 
 * 
 * 
 * @param ga_orig handle to existing GA to use as a model
 * @param ntype GA data type for new GA
 * 
 * @return handle to new ga with same dimensions and distribution as
 * @c ga_orig, but with a different type
 */
static int
GA_Duplicate_type(int oga, char *nname, int ntype)
{
  int nga;
  int otype;
  int ndim, dims[GA_MAX_DIM], chunk[GA_MAX_DIM];
  int nblock, map[GA_MAX_DIM];
  int i;
  
  NGA_Inquire(oga, &otype, &ndim, &dims[0]);

  /* if it's already the correct type, just duplicate */

  if (otype == ntype) {
    nga = GA_Duplicate(oga, nname);
  } else {
    /* FIXME: should actually copy the distribution here */
    for (i = 0; i < GA_MAX_DIM; ++i) chunk[i] = 1;
    nga = GA_Create_handle();
    GA_Set_data(nga, ndim, dims, ntype);
    GA_Set_array_name(nga, nname);
    GA_Set_chunk(nga, &chunk[0]);
    GA_Allocate(nga);
  }    
  return nga;
}
    
/******************************************************************************/
/*                             Distribute2DMatrix                             */
/******************************************************************************/
void
Distribute2DMatrix(void *MatrixZero, void *LocalMatrix, 
                   int NumberType, MAPSIZE *Map)
{
  int me;
  int ndim, dims[GA_MAX_DIM];
  int gNX, gNY;
  int i, j;
  int gatype, ga;
  int lo[2], hi[2], ld[2];

  me = ParallelRank();

  GlobalMapSize(Map, &gNX, &gNY);

  gatype = GA_Type(NumberType);
  ga = GA_Duplicate_type(Map->dist, "Distribute2DMatrix", GA_Type(NumberType));
  
  switch (gatype) {
  case (C_CHAR):
    break;
  default:
    GA_Print_distribution(ga);
    break;
  }
  
  
  if (me == 0) {

    lo[0] = 0;
    lo[1] = 0;
    hi[gaYdim] = gNY-1;
    hi[gaXdim] = gNX-1;
    ld[gaYdim] = gNY-1;
    ld[gaXdim] = gNX-1;
    NGA_Put(ga, &lo[0], &hi[0], MatrixZero, &ld[0]);
  }
  GA_Sync();

  lo[gaYdim] = Map->OffsetY;
  lo[gaXdim] = Map->OffsetX;
  hi[gaYdim] = lo[gaYdim] + Map->NY - 1;
  hi[gaXdim] = lo[gaXdim] + Map->NX - 1;
  ld[gaYdim] = Map->NY - 1;
  ld[gaXdim] = Map->NX - 1;
  NGA_Get(ga, &lo[0], &hi[0], LocalMatrix, &ld[0]);

  GA_Sync();
  GA_Destroy(ga);
}

/******************************************************************************/
/*                            Collect2DMatrix                                 */
/******************************************************************************/
void
Collect2DMatrix(void *MatrixZero, void *LocalMatrix, 
                int NumberType, MAPSIZE *Map)
{
  int me;
  int ndim, dims[GA_MAX_DIM];
  int gNX, gNY;
  int i, j;
  int gatype, ga;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  
  me = ParallelRank();
  
  GlobalMapSize(Map, &gNX, &gNY);
  
  gatype = GA_Type(NumberType);
  ga = GA_Duplicate_type(Map->dist, "Collect2DMatrix", GA_Type(NumberType));
  
  lo[gaYdim] = Map->OffsetY;
  lo[gaXdim] = Map->OffsetX;
  hi[gaYdim] = lo[gaYdim] + Map->NY;
  hi[gaXdim] = lo[gaXdim] + Map->NX;
  ld[gaYdim] = Map->NY;
  ld[gaXdim] = Map->NX;
  NGA_Put(ga, &lo[0], &hi[0], LocalMatrix, &ld[0]);
  GA_Sync();

  if (me == 0) {

    lo[gaYdim] = 0;
    lo[gaXdim] = 0;
    hi[gaYdim] = gNY;
    hi[gaXdim] = gNX;
    ld[gaYdim] = gNY;
    ld[gaXdim] = gNX;
    NGA_Get(ga, &lo[0], &hi[0], MatrixZero, &ld[0]);
  }
  GA_Sync();
  GA_Destroy(ga);
}


/******************************************************************************/
/*                              Read2DMatrix                                  */
/******************************************************************************/
/** 
 * 
 * 
 * @param FileName name of file to read
 * @param Matrix  @e local 2D array (NX, NY) to be filled
 * @param NumberType 
 * @param NY 
 * @param NX 
 * @param NDataSet 
 * @param VarName 
 * @param index 
 * 
 * @return 
 */
int 
Read2DMatrix(char *FileName, void *LocalMatrix, int NumberType, MAPSIZE *Map, 
             int NDataSet, char *VarName, int index)
{
  const char Routine[] = "Read2DMatrix";
  void *tmpArray;
  int gNX, gNY;
  int me;

  me = ParallelRank();

  GlobalMapSize(Map, &gNX, &gNY);

  if (me == 0) {
    if (!(tmpArray = (void *)calloc(gNY * gNX, SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    Read2DMatrixFmt(FileName, tmpArray, NumberType, gNY, gNX, NDataSet, VarName, index);
  }

  Distribute2DMatrix(tmpArray, LocalMatrix, NumberType, Map);

  if (me == 0) {
    free(tmpArray);
  }

  return 0;
}

/******************************************************************************/
/*                              Write2DMatrix                                  */
/******************************************************************************/
int
Write2DMatrix(char *FileName, void *LocalMatrix, int NumberType, 
              MAPSIZE *Map, MAPDUMP *DMap, int index)
{
  const char Routine[] = "Write2DMatrix";
  void *tmpArray;
  int gNX, gNY;
  int me;

  me = ParallelRank();

  GlobalMapSize(Map, &gNX, &gNY);

  if (me == 0) {
    if (!(tmpArray = (void *)calloc(gNY * gNX, SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
  } else {
    tmpArray = NULL;
  }
  Collect2DMatrix(tmpArray, LocalMatrix, NumberType, Map);

  if (me == 0) {
    Write2DMatrixFmt(FileName, tmpArray, NumberType, Map->NY, Map->NX, DMap, index);
    free(tmpArray);
  }
  return 0;
}
