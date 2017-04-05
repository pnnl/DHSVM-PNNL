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
 * LAST CHANGE: 2017-04-05 14:05:21 d3g096
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
#include <macdecls.h>

#include "sizeofnt.h"
#include "DHSVMerror.h"
#include "ParallelDHSVM.h"


const int gaXdim = 1;
const int gaYdim = 0;


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

/* -------------------------------------------------------------
   compare_int
   ------------------------------------------------------------- */
int 
compare_int(const void *i1ptr, const void *i2ptr)
{
  int i1 = *((int *)i1ptr);
  int i2 = *((int *)i2ptr);
  return (int) (i1 - i2);
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
  int ndim, dims[GA_MAX_DIM];
  int otype;

  NGA_Inquire(oga, &otype, &ndim, &dims[0]);

  /* if it's already the correct type, just duplicate */

  if (otype == ntype) {
    nga = GA_Duplicate(oga, nname);
  } else {
    int np;
    int *idx;
    int lo[GA_MAX_DIM], hi[GA_MAX_DIM];
    int nblk[GA_MAX_DIM], *mapc, *mapcptr;
    int p; 
    int i, d;

    np = ParallelSize();

    if (!(idx = (int *)calloc(np, sizeof(int)))) {
      ReportError("GA_Duplicate_type", 70);
    }
    if (!(mapc = (int *)calloc((GA_MAX_DIM-1)*np, sizeof(int)))) {
      ReportError("GA_Duplicate_type", 70);
    }

    mapcptr = mapc;
    for (d = 0; d < ndim; ++d) {
      for (p = 0; p < np; ++p) {
        NGA_Distribution(oga, p, &lo[0], &hi[0]);
        idx[p] = lo[d];
      }
      qsort(&idx[0], np, sizeof(int), compare_int);
      
      *mapcptr = 0;
      for (p = 0, i = 0; p < np ; ++p) {
        if (idx[p] != *mapcptr) {
          mapcptr++;
          *mapcptr = idx[p];
          i++;
        }
      }
      nblk[d] = i+1;
      mapcptr++;
    }
    
    nga = GA_Create_handle();
    GA_Set_array_name(nga, nname);
    GA_Set_data(nga, ndim, dims, ntype);
    GA_Set_irreg_distr(nga, mapc, nblk);
    GA_Allocate(nga);

    if (GA_Compare_distr(oga, nga)) {
      ReportError("GA_Duplicate_type: distributions differ", 70);
    }

    free(idx);
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
/*                                GA_Acc_one                                  */
/******************************************************************************/
void
GA_Acc_one(int ga, MAPSIZE *Map, int x, int y, void *value, void *alpha)
{
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  lo[gaXdim] = Map->OffsetX + x;
  hi[gaXdim] = lo[gaXdim];
  lo[gaYdim] = Map->OffsetY + y;
  hi[gaYdim] = lo[gaYdim];
  ld[gaXdim] = 1;
  ld[gaYdim] = 1;
  NGA_Acc(ga, &lo[0], &hi[0], value, &ld[0], alpha);
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
  int me;
  int dims[GA_MAX_DIM];
  int chunk[GA_MAX_DIM];
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM];

  me = ParallelRank();

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

  gaid = NGA_Create(C_FLOAT, 2, dims, "Domain Decompsition", chunk);
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

/******************************************************************************/
/*                                                                    */
/******************************************************************************/
static
float **
alloc_float_2d(int NY, int NX)
{
  static char Routine[] = "alloc_float_2d";
  int j;
  float **result;
  result = (float **)calloc(NY, sizeof(float *));
  if (result == NULL) {
    ReportError(Routine, 1);
  }
  result[0] = (float *) calloc(NX*NY, sizeof(float));
  if (result[0] == NULL) {
    ReportError(Routine, 1);
  }
  for (j = 1; j < NY; ++j) {
    result[j] = result[0] + j*NX;
  }
  return result;
}

/******************************************************************************/
/*                              GA_Alloc_patch                                */
/******************************************************************************/
void
GA_Alloc_patch(int ga, MAPSIZE *Map, GA_Patch *p)
{
  p->ixoff = 0;
  p->iyoff = 0;
  p->NX = Map->NX;
  p->NY = Map->NY;
  p->patch = alloc_float_2d(Map->NY, Map->NX);
}

/******************************************************************************/
/*                              GA_Alloc_patch_ghost                          */
/******************************************************************************/
void
GA_Alloc_patch_ghost(int ga, MAPSIZE *Map, GA_Patch *p)
{
  p->ixoff = 0;
  p->iyoff = 0;
  p->NX = Map->NX;
  p->NY = Map->NY;

  if (Map->OffsetX > 0) {
    p->NX += 1;
    p->ixoff = 1;
  }

  if (Map->OffsetY > 0) {
    p->NY += 1;
    p->iyoff = 1;
  }

  if (Map->OffsetX + Map->NX < Map->gNX) {
    p->NX += 1;
  } 
  if (Map->OffsetY + Map->NY < Map->gNY) {
    p->NY += 1;
  } 
  p->patch = alloc_float_2d(p->NY, p->NX);
}

/******************************************************************************/
/*                             fill_ga_dims                                   */
/******************************************************************************/
static void
fill_ga_dims(MAPSIZE *Map, GA_Patch *p, int *lo, int *hi, int *ld)
{
  lo[gaXdim] = Map->OffsetX - p->ixoff;
  lo[gaYdim] = Map->OffsetY - p->iyoff;
  hi[gaXdim] = lo[gaXdim] + p->NX - 1;
  hi[gaYdim] = lo[gaYdim] + p->NY - 1;
  ld[gaXdim] = p->NY;
  ld[gaYdim] = p->NX;
}

/******************************************************************************/
/*                               GA_Acc_patch                                 */
/******************************************************************************/
void
GA_Acc_patch(int ga, MAPSIZE *Map, GA_Patch *p)
{
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float alpha;
  alpha = 1.0;
  fill_ga_dims(Map, p, &lo[0], &hi[0], &ld[0]);
  NGA_Acc(ga, &lo[0], &hi[0], &(p->patch[0][0]), &ld[0], &alpha);
}

/******************************************************************************/
/*                               GA_Get_patch                                 */
/******************************************************************************/
void
GA_Get_patch(int ga, MAPSIZE *Map, GA_Patch *p)
{
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  fill_ga_dims(Map, p, &lo[0], &hi[0], &ld[0]);
  NGA_Get(ga, &lo[0], &hi[0], &(p->patch[0][0]), &ld[0]);
}

/******************************************************************************/
/*                               GA_Put_patch                                 */
/******************************************************************************/
void
GA_Put_patch(int ga, MAPSIZE *Map, GA_Patch *p)
{
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  fill_ga_dims(Map, p, &lo[0], &hi[0], &ld[0]);
  NGA_Put(ga, &lo[0], &hi[0], &(p->patch[0][0]), &ld[0]);
}


/******************************************************************************/
/*                              GA_Free_patch                                 */
/******************************************************************************/
void
GA_Free_patch(GA_Patch *p)
{
  free(p->patch[0]);
  free(p->patch);
}
