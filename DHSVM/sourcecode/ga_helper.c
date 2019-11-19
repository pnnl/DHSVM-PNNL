/*
 * SUMMARY:      ga_helper.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    November 2018
 * DESCRIPTION:  Some very basic things to aid in Global Arrays usage.
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2019-11-19 09:24:36 d3g096
 * COMMENTS:
 */

#include <stdlib.h>
#include "sizeofnt.h"
#include "ga_helper.h"

const int gaXdim = 1;
const int gaYdim = 0;

/******************************************************************************/
/*                          ParallelInitialize                                */
/******************************************************************************/
void
ParallelInitialize(int *argc, char ***argv)
{
  int ierr;
  ierr = MPI_Init(argc, argv);
  if (ierr != 0) {
    ReportError("ParallelInitialize: MPI_Init: ", 70);
  }
  GA_Initialize();
  if (!MA_init(MT_C_DBL, 500000, 500000)) {
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
void
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

  /* TIMING_TASK_START("GA Creation", 4); */

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
  /* TIMING_TASK_END("GA Creation", 4); */
  return nga;
}


