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
 * LAST CHANGE: 2017-04-27 12:45:08 d3g096
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

#include "constants.h"
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


/******************************************************************************/
/*                                compare_int                                    */
/******************************************************************************/
static int 
compare_int(const void *i1ptr, const void *i2ptr)
{
  int i1 = *((int *)i1ptr);
  int i2 = *((int *)i2ptr);
  return (int) (i1 - i2);
}

/******************************************************************************/
/*                             GA_Inquire_irreg_distr                               */
/******************************************************************************/
/* the nblk and mapc arrays are filled with the information requried
   by GA_Irreg_distr.  The length of nblk needs to be greater or equal
   to ga's number of dimensions. The length of mapc should be the
   number of processors times the number of dimensions */
static void
GA_Inquire_irreg_distr(int ga, int *mapc, int *nblk)
{
  int np;
  int *idx;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM];
  int *mapcptr;
  int p; 
  int i, d;
  int ndim, dims[GA_MAX_DIM];
  int gatype;

  NGA_Inquire(ga, &gatype, &ndim, &dims[0]);

  np = GA_Nnodes();
  if (!(idx = (int *)calloc(np, sizeof(int)))) {
    ReportError("GA_Inquire_irreg_distr", 70);
  }

  mapcptr = mapc;
  for (d = 0; d < ndim; ++d) {
    for (p = 0; p < np; ++p) {
      NGA_Distribution(ga, p, &lo[0], &hi[0]);
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
  free(idx);
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
    int nblk[GA_MAX_DIM], *mapc;

    np = ParallelSize();

    if (!(mapc = (int *)calloc((GA_MAX_DIM-1)*np, sizeof(int)))) {
      ReportError("GA_Duplicate_type", 70);
    }

    GA_Inquire_irreg_distr(oga, &mapc[0], &nblk[0]);

    nga = GA_Create_handle();
    GA_Set_array_name(nga, nname);
    GA_Set_data(nga, ndim, dims, ntype);
    GA_Set_irreg_distr(nga, mapc, nblk);
    GA_Allocate(nga);

    if (GA_Compare_distr(oga, nga)) {
      ReportError("GA_Duplicate_type: distributions differ", 70);
    }

    free(mapc);
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
/*                                GA_Acc_one                                  */
/******************************************************************************/
void
GA_Acc_one_global(int ga, MAPSIZE *Map, int x, int y, void *value, void *alpha)
{
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  lo[gaXdim] = x;
  hi[gaXdim] = lo[gaXdim];
  lo[gaYdim] = y;
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

  fflush(stdout);
  fflush(stderr);
  ParallelBarrier();

  me = ParallelRank();
  nproc = ParallelSize();

  if (me == 0) {
    printf(bar);
    printf("Proc       NX      NY OffsetX OffsetY         Xorig         Yorig NumCells\n");
    printf(bar);
  }
  fflush(stdout);
  ParallelBarrier();
  for (p = 0; p < nproc; ++p) {
    if (me == p) {
      printf(fmt, p, 
             local->NX, local->NY, 
             local->OffsetX, local->OffsetY,
             local->Xorig, local->Yorig,
             local->NumCells);
    }
    fflush(stdout);
    ParallelBarrier();
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
  fflush(stdout);
  ParallelBarrier();
  GA_Print_distribution(global->dist);
  if (me == 0) {
    printf(bar);
  }
  fflush(stdout);
  ParallelBarrier();
}

/******************************************************************************/
/*                             GA_Mapsize                                     */
/******************************************************************************/
static void
GA_Mapsize(MAPSIZE *global, MAPSIZE *local, int gaid)
{
  int me;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM];

  me = GA_Nodeid();

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
/*                            DomainDecomposition                             */
/******************************************************************************/
void
SimpleDomainDecomposition(MAPSIZE *global, MAPSIZE *local)
{
  int gaid; 
  int dims[GA_MAX_DIM];
  int chunk[GA_MAX_DIM];

  /* initialize the local domain */
  memcpy(local, global, sizeof(MAPSIZE));

  /* these should not be set in global, but be sure local has safe
     values */
  local->OrderedCells = NULL;
  local->NumCells = 0;
  local->AllCells = 0;
  
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

  GA_Mapsize(global, local, gaid);

}

/******************************************************************************/
/*                             find_splits                                    */
/******************************************************************************/
/* it's assumed that ga is 1-D float */
static void
find_splits(int ga, int nsplit, int *isplit)
{
  int me, np;
  int ga_mask, ga_sum;
  int gatype, ndim, dim[GA_MAX_DIM];
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM];
  int i, f, mylo, myhi, idx0;
  float value, *ga_data, *frac;

  me = GA_Nodeid();
  np = GA_Nnodes();

  NGA_Inquire(ga, &gatype, &ndim, &dim[0]);

  ga_mask = GA_Duplicate(ga, "find_splits Mask");
  GA_Zero(ga_mask);

  if (me == 0) {
    lo[0] = 0;
    hi[0] = 0;
    value = 1;
    NGA_Put(ga_mask, &lo[0], &hi[0], &value, NULL);
  }
  GA_Sync();

  ga_sum = GA_Duplicate(ga, "find_splits Sum");
  GA_Zero(ga_sum);

  GA_Scan_add(ga, ga_sum, ga_mask, 0, dim[0], 0);
  NGA_Select_elem(ga_sum, "max", &value, &lo[0]);
  value = 1.0/value;

  GA_Scale(ga_sum, &value);

  if (!(frac = (float *)calloc(nsplit, sizeof(float)))) {
    ReportError("find_splits", 70);
  }
  for (f = 0; f < nsplit; ++f) {
    isplit[f] = 0;
    frac[f] = ((double)f)/((double)nsplit);
  }

  NGA_Distribution(ga_sum, me, &lo[0], &hi[0]);
  mylo = lo[0];
  myhi = hi[0];
  NGA_Access(ga_sum, &mylo, &myhi, &ga_data, NULL);

  idx0 = mylo;
  for (f = 1; f < nsplit; ++f) {
    if (ga_data[idx0 - mylo] <= frac[f] && frac[f] <= ga_data[myhi - mylo]) {
      for (i = idx0 + 1; i <= myhi; ++i) {
        if (ga_data[i - mylo] > frac[f]) {
          isplit[f] = i - 1;
          break;
        }
      }
    }
  }
  NGA_Release(ga_sum, &mylo, &myhi);
  /* GA_Print(ga_sum); */

  GA_Igop(&isplit[0], nsplit, "+");
  
  GA_Destroy(ga_mask);
  GA_Destroy(ga_sum);
}
  

/******************************************************************************/
/*                            MaskedDomainDecomposition                       */
/******************************************************************************/
void 
MaskedDomainDecomposition(MAPSIZE *gmap, MAPSIZE *lmap, MAPSIZE *nmap, 
                          unsigned char *mask)
{
  static float one = 1.0;
  int me, np;
  int ga_xsum, ga_ysum;
  int gatype;
  int ndim, dims[GA_MAX_DIM];
  int nblk[GA_MAX_DIM], *mapc, *mapcptr;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float sum;
  int x, y, gx, gy, i;
  int p, d, b;
  char *nname;
  int nga;

  me = ParallelRank();
  np = ParallelSize();

  if (!(mapc = (int *)calloc(GA_MAX_DIM*np, sizeof(int)))) {
    ReportError("MaskedDomainDecomposition", 70);
  }
  NGA_Inquire(gmap->dist, &gatype, &ndim, &dims[0]);
  GA_Inquire_irreg_distr(gmap->dist, &mapc[0], &nblk[0]);

  if (nblk[gaYdim] > 1) {
    gy = gmap->NY;
    ga_ysum = NGA_Create(C_FLOAT, 1, &gmap->NY, "Sum along Y", NULL);
    GA_Zero(ga_ysum);

    for (y = 0, i = 0; y < lmap->NY; y++) {
      sum = 0.0;
      for (x = 0; x < lmap->NX; x++, i++) {
        if (INBASIN(mask[i])) {
          sum += 1.0;
        }
      }
      Local2Global(lmap, 0, y, &gx, &gy);
      NGA_Acc(ga_ysum, &gy, &gy, &sum, &ld[0], &one);
    }

    find_splits(ga_ysum, nblk[gaYdim], &mapc[0]);
    /* GA_Print(ga_ysum); */
    GA_Destroy(ga_ysum);
  } else {
    mapc[0] = 0;
  }
  
  if (nblk[gaXdim] > 1) {
    gx = gmap->NX;
    ga_xsum = NGA_Create(C_FLOAT, 1, &gmap->NX, "Sum along X", NULL);
    GA_Zero(ga_xsum);

    for (x = 0, i = 0; x < lmap->NX; x++) {
      sum = 0.0;
      for (y = 0; y < lmap->NY; y++, i++) {
        if (INBASIN(mask[i])) {
          sum += 1.0;
        }
      }
      Local2Global(lmap, x, 0, &gx, &gy);
      NGA_Acc(ga_xsum, &gx, &gx, &sum, &ld[0], &one);
    }

    find_splits(ga_xsum, nblk[gaXdim], &mapc[nblk[0]]);

    /* GA_Print(ga_xsum); */
    GA_Destroy(ga_xsum);
  } else {
    mapc[nblk[0]] = 0;
  }

  /*
  for (p = 0; p < np; ++p) {
    if (me == p) {
      if (p == 0) {
        fprintf(stdout, "MaskedDomainDecomposition:\n");
      }
      mapcptr = mapc;
      for (d = 0; d < 2; ++d) {
        fprintf(stdout, "%d: %d(%d): ", p, d, nblk[d]);
        for (b = 0; b < nblk[d]; ++b, ++mapcptr) {
          fprintf(stdout, "%d, ", *mapcptr);
        }
        fprintf(stdout, "\n");
        fflush(stdout);
      }
    }
    GA_Sync();
  }
  */

  memcpy(nmap, lmap, sizeof(MAPSIZE));

  nname = GA_Inquire_name(gmap->dist);

  nga = GA_Create_handle();
  GA_Set_array_name(nga, nname);
  GA_Set_data(nga, ndim, dims, gatype);
  GA_Set_irreg_distr(nga, mapc, nblk);
  GA_Allocate(nga);

  GA_Destroy(gmap->dist);
  GA_Mapsize(gmap, nmap, nga);
  
  free(mapc);
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

/******************************************************************************/
/*                            Collect2DMatrixGA                               */
/******************************************************************************/
int
Collect2DMatrixGA(void *LocalMatrix, int NumberType, MAPSIZE *Map)
{
  int me;
  int ndim, dims[GA_MAX_DIM];
  int gNX, gNY;
  int i, j;
  int ga, gatype;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  
  me = ParallelRank();
  
  gNX = Map->gNX;
  gNY = Map->gNY;
  
  gatype = GA_Type(NumberType);
  ga = GA_Duplicate_type(Map->dist, "Collect2DMatrix", gatype);
  /* GA_Print_distribution(ga); */
  
  lo[gaYdim] = Map->OffsetY;
  lo[gaXdim] = Map->OffsetX;
  hi[gaYdim] = lo[gaYdim] + Map->NY - 1;
  hi[gaXdim] = lo[gaXdim] + Map->NX - 1;
  ld[gaXdim] = Map->NY;
  ld[gaYdim] = Map->NX;
  NGA_Put(ga, &lo[0], &hi[0], LocalMatrix, &ld[0]);
  GA_Sync();
  /* GA_Print(ga); */
  return ga;
}

