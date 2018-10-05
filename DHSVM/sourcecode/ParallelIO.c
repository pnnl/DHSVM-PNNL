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
 * LAST CHANGE: 2018-10-05 13:21:03 d3g096
 * COMMENTS:
 */

/******************************************************************************/
/*                                INCLUDES                                    */
/******************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <ga.h>

#include "settings.h"
#include "MapSize.h"
#include "sizeofnt.h"
#include "DHSVMerror.h"
#include "ParallelDHSVM.h"
#include "fileio.h"

/* global function pointers */
extern void (*CreateMapFileFmt) (char *FileName, ...);
extern int (*Read2DMatrixFmt) (char *FileName, void *Matrix, int NumberType, int NY, int NX, int NDataSet, ...);
extern int (*Write2DMatrixFmt) (char *FileName, void *Matrix, int NumberType, int NY, int NX, ...);

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
/*                             Distribute2DMatrix                             */
/******************************************************************************/
static void
Distribute2DMatrix(void *MatrixZero, void *LocalMatrix, 
                   int NumberType, MAPSIZE *Map, int mirror)
{
  int me;
  int gNX, gNY;
  int gatype, ga;
  int lo[2], hi[2], ld[2];

  me = ParallelRank();

  gNX = Map->gNX;
  gNY = Map->gNY;

  gatype = GA_Type(NumberType);
  ga = GA_Duplicate_type(Map->dist, "Distribute2DMatrix", gatype);
  
  /* switch (gatype) { */
  /* case (C_CHAR): */
  /*   break; */
  /* default: */
  /*   GA_Print_distribution(ga); */
  /*   break; */
  /* } */
  
  
  if (me == 0) {
    lo[0] = 0;
    lo[1] = 0;
    hi[gaYdim] = gNY-1;
    hi[gaXdim] = gNX-1;
    ld[gaXdim] = gNY;
    ld[gaYdim] = gNX;
    NGA_Put(ga, &lo[0], &hi[0], MatrixZero, &ld[0]);
  }
  ParallelBarrier();

  if (mirror) {
    lo[0] = 0;
    lo[1] = 0;
    hi[gaYdim] = gNY-1;
    hi[gaXdim] = gNX-1;
    ld[gaXdim] = gNY;
    ld[gaYdim] = gNX;
  } else {
    lo[gaYdim] = Map->OffsetY;
    lo[gaXdim] = Map->OffsetX;
    hi[gaYdim] = lo[gaYdim] + Map->NY - 1;
    hi[gaXdim] = lo[gaXdim] + Map->NX - 1;
    ld[gaXdim] = Map->NY;
    ld[gaYdim] = Map->NX;
  }
  NGA_Get(ga, &lo[0], &hi[0], LocalMatrix, &ld[0]);

  ParallelBarrier();
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
  int ga;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  int gNX, gNY;

  me = ParallelRank();
  ga = Collect2DMatrixGA(LocalMatrix, NumberType, Map);

  if (me == 0) {

    gNX = Map->gNX;
    gNY = Map->gNY;

    lo[gaYdim] = 0;
    lo[gaXdim] = 0;
    hi[gaYdim] = gNY-1;
    hi[gaXdim] = gNX-1;
    ld[gaXdim] = gNY;
    ld[gaYdim] = gNX;
    NGA_Get(ga, &lo[0], &hi[0], MatrixZero, &ld[0]);
  }
  ParallelBarrier();
  GA_Destroy(ga);
}

/******************************************************************************/
/*                             intRead2DMatrix                                */
/******************************************************************************/
static int 
intRead2DMatrix(char *FileName, void *LocalMatrix, int NumberType, MAPSIZE *Map, 
                int NDataSet, char *VarName, int index, int mirror)
{
  const char Routine[] = "Read2DMatrix";
  void *tmpArray;
  int gNX, gNY;
  int flag;
  int me;

  me = ParallelRank();

  gNX = Map->gNX;
  gNY = Map->gNY;

  if (me == 0) {
    if (!(tmpArray = (void *)calloc(gNY * gNX, SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
    flag = Read2DMatrixFmt(FileName, tmpArray, NumberType, gNY, gNX, NDataSet, VarName, index);
  } else {
    flag = 0;
  }
  Distribute2DMatrix(tmpArray, LocalMatrix, NumberType, Map, mirror);

  if (me == 0) {
    free(tmpArray);
  }
  GA_Brdcst(&flag, sizeof(int), 0);

  return flag;
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

  return intRead2DMatrix(FileName, LocalMatrix, NumberType, Map, 
                         NDataSet, VarName, index, 0);
}

/******************************************************************************/
/*                         Read2DMatrixAll                                    */
/******************************************************************************/
int 
Read2DMatrixAll(char *FileName, void *LocalMatrix, int NumberType, MAPSIZE *Map, 
             int NDataSet, char *VarName, int index)
{
  const char Routine[] = "Read2DMatrixAll";
  return intRead2DMatrix(FileName, LocalMatrix, NumberType, Map, 
                         NDataSet, VarName, index, 1);
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
  int me, flag;

  me = ParallelRank();

  gNX = Map->gNX;
  gNY = Map->gNY;

  if (me == 0) {
    if (!(tmpArray = (void *)calloc(gNY * gNX, SizeOfNumberType(NumberType))))
      ReportError((char *)Routine, 1);
  } else {
    tmpArray = NULL;
  }
  Collect2DMatrix(tmpArray, LocalMatrix, NumberType, Map);

  if (me == 0) {
    flag = Write2DMatrixFmt(FileName, tmpArray, NumberType, gNY, gNX, DMap, index);
    free(tmpArray);
  }
  GA_Brdcst(&flag, sizeof(int), 0);
  return flag;
}
