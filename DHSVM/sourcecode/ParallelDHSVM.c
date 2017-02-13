/*
 * SUMMARY:      ParallelDHSVM.c
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
 * LAST CHANGE: 2017-02-13 10:22:34 d3g096
 * COMMENTS:
 */



/******************************************************************************/
/*                                INCLUDES                                    */
/******************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <ga.h>

#include "sizeofnt.h"
#include "DHSVMerror.h"
#include "ParallelDHSVM.h"


const gaXdim = 1;
const gaYdim = 0;


/******************************************************************************/
/*                          ParallelInitialize                                */
/******************************************************************************/
void
ParallelInitialize(int *argc, char ***argv)
{
  int ierr;
  ierr = 0;
  GA_Initialize_args(argc, argv);
  if (!MA_init(MT_C_DBL, 5000, 5000)) {
    ReportError("ParallelInitialize: MA_init: ", 70);
    ierr += 1;
  }
  if (ierr > 0) exit(3);
}

/******************************************************************************/
/*                            ParallelRank                                    */
/******************************************************************************/
int
ParallelRank(void)
{
  return GA_Nodeid();
}

/******************************************************************************/
/*                            ParallelSize                                    */
/******************************************************************************/
int
ParallelSize(void)
{
  return GA_Nnodes();
}


/******************************************************************************/
/*                                GA_Type                                     */
/******************************************************************************/
int
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
int
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
/*                                GA_Put_one                                  */
/******************************************************************************/
void
GA_Put_one(int ga, MAPSIZE *Map, int x, int y, void *value)
{
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  lo[gaXdim] = Map->OffsetX + x;
  hi[gaXdim] = lo[gaXdim];
  lo[gaYdim] = Map->OffsetY + y;
  hi[gaYdim] = lo[gaYdim];
  ld[gaXdim] = 1;
  ld[gaYdim] = 1;
  NGA_Put(ga, &lo[0], &hi[0], value, &ld[0]);
}
/******************************************************************************/
/*                                GA_Get_one                                  */
/******************************************************************************/
void
GA_Get_one(int ga, MAPSIZE *Map, int x, int y, void *value)
{
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  lo[gaXdim] = Map->OffsetX + x;
  hi[gaXdim] = lo[gaXdim];
  lo[gaYdim] = Map->OffsetY + y;
  hi[gaYdim] = lo[gaYdim];
  ld[gaXdim] = 1;
  ld[gaYdim] = 1;
  NGA_Get(ga, &lo[0], &hi[0], value, &ld[0]);
}


/******************************************************************************/
/*                             DomainSummary                                  */
/******************************************************************************/
/* 
--------------------------------------------------------------------------
Proc       NX      NY OffsetX OffsetY         Xorig         Yorig NumCells
--------------------------------------------------------------------------
##### ####### ####### ####### ####### ##########.## ##########.## ########
--------------------------------------------------------------------------
*/
void
DomainSummary(MAPSIZE *global, MAPSIZE *local)
{
  static const char bar[] = 
    "---------------------------------------------------------------------------\n";
  static const char fmt[] = 
    "%6d %7d %7d %7d %7d %13.2f %13.2f %8d\n";
  static const char sfmt[] = 
    "%6s %7d %7d %7d %7d %13.2f %13.2f %8d\n";
  static const char b[] = " ";

  int me, nproc, p;

  ParallelBarrier();
  fflush(stdout);

  me = ParallelRank();
  nproc = ParallelSize();

  if (me == 0) {
    printf(bar);
    printf("Proc       NX      NY OffsetX OffsetY         Xorig         Yorig NumCells\n");
    printf(bar);
  }
  ParallelBarrier();
  fflush(stdout);
  for (p = 0; p < nproc; ++p) {
    if (me == p) {
      printf(fmt, p, 
             local->NX, local->NY, 
             local->OffsetX, local->OffsetY,
             local->Xorig, local->Yorig,
             local->NumCells);
    }
    ParallelBarrier();
    fflush(stdout);
  }
  if (me == 0) {
    printf(bar);
    printf(sfmt, "global",
           global->NX, global->NY, 
           global->OffsetX, global->OffsetY,
           global->Xorig, global->Yorig,
           global->NumCells);
    printf(bar);
  }
  ParallelBarrier();
  GA_Print_distribution(global->dist);
  if (me == 0) {
    printf(bar);
  }
  ParallelBarrier();
  fflush(stdout);
}

/******************************************************************************/
/*                            DomainDecomposition                             */
/******************************************************************************/
void
DomainDecomposition(MAPSIZE *global, MAPSIZE *local)
{
  int gaid; 
  int nproc, me, p;
  int dims[GA_MAX_DIM];
  int chunk[GA_MAX_DIM], ghosts[GA_MAX_DIM];
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM];
  int gNX, gNY;

  me = ParallelRank();
  nproc = ParallelSize();
  

  /* initialize the local domain */
  memcpy(local, global, sizeof(MAPSIZE));

  /* these should not be set in global, but be sure local has safe
     values */
  local->OrderedCells = NULL;
  local->NumCells = 0;
  
  /* create an appropriate sized GA the default way and use it to
     determine local shares of the domain. */

  dims[gaYdim] = global->NY;
  dims[gaXdim] = global->NX;
  
  chunk[gaYdim] = 1;
  chunk[gaXdim] = 1;

  ghosts[gaYdim] = 1;
  ghosts[gaXdim] = 1;
  
  gaid = NGA_Create_ghosts(C_INT, 2, dims, ghosts, "Domain Decompsition", chunk);
  if (gaid == 0) {
    ReportError("DomainDecomposition", 70);
  }
  /* GA_Print_distribution(gaid); */
  NGA_Distribution(gaid, me, lo, hi);
  global->dist = gaid;
  local->dist = gaid;
  local->Xorig = global->Xorig + lo[gaXdim]*global->DX;
  local->Yorig = global->Yorig + lo[gaYdim]*global->DY;
  local->OffsetX = lo[gaXdim];
  local->OffsetY = lo[gaYdim];
  local->NX = hi[gaXdim] - lo[gaXdim] + 1;
  local->NY = hi[gaYdim] - lo[gaYdim] + 1;

}

/******************************************************************************/
/*                            ParallelBarrier                                 */
/******************************************************************************/
void
ParallelBarrier()
{
  GA_Sync();
}

/******************************************************************************/
/*                            ParallelFinalize                                */
/******************************************************************************/
void
ParallelFinalize(void)
{
  int ierr;
  ierr = 0;
  GA_Terminate();
  ierr = MPI_Finalize();
  if (ierr != 0) {
    ReportError("ParallelInitialize: MA_init: ", 70);
  }
}


/******************************************************************************/
/*                              Global2Local                                  */
/******************************************************************************/
/** 
 * Compute process local cell indexes given the global indexes.  
 * 
 * @param Map 
 * @param globalx 
 * @param globaly 
 * @param localx 
 * @param localy 
 * 
 * @return 1 if global coordinate is in local domain
 */
int
Global2Local(MAPSIZE *Map, int globalx, int globaly, int *localx, int *localy)
{
  int tmpx, tmpy;
  int result;

  tmpx = globalx - Map->OffsetX;
  tmpy = globaly - Map->OffsetY;
  if (tmpx < 0 || tmpy < 0 || tmpx >= Map->NX || tmpy >= Map->NY) {
    result = 0;
  } else {
    *localx = tmpx;
    *localy = tmpy;
    result = 1;
  }
  return result;
}



/******************************************************************************/
/*                              Local2Global                                  */
/******************************************************************************/
/** 
 * This computes the global column and row given @em valid local indexes.
 * 
 * @param Map local domain description
 * @param localx @em valid local column index
 * @param localy @em valid local row index
 * @param globalx global column index
 * @param globaly global row index
 */
void
Local2Global(MAPSIZE *Map, int localx, int localy, int *globalx, int *globaly)
{
  *globalx = localx + Map->OffsetX;
  *globaly = localy + Map->OffsetY;
}
