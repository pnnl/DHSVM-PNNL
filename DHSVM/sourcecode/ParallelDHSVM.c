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
 * LAST CHANGE: 2017-01-23 12:38:59 d3g096
 * COMMENTS:
 */



/******************************************************************************/
/*                                INCLUDES                                    */
/******************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <ga.h>

#include "DHSVMerror.h"
#include "data.h"

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
/*                        ParallelDomainDecomposition                         */
/******************************************************************************/
void
DomainDecomposition(MAPSIZE *global, MAPSIZE *local)
{
  int gaid; 
  int nproc, me, p;
  int dims[2];
  int chunk[2];
  int lo[2], hi[2];

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

  dims[0] = global->NX;
  dims[1] = global->NY;
  
  chunk[0] = 1;
  chunk[1] = 1;
  
  gaid = NGA_Create(C_INT, 2, dims, "Domain Decompsition", chunk);
  if (gaid == 0) {
    ReportError("DomainDecomposition", 70);
  }
  GA_Print_distribution(gaid);
  NGA_Distribution(gaid, me, lo, hi);
  local->Xorig = global->Xorig + lo[0]*global->DX;
  local->Yorig = global->Yorig + lo[1]*global->DY;
  local->OffsetX = lo[0];
  local->OffsetY = lo[1];
  local->NX = hi[0] - lo[0]+ 1;
  local->NY = hi[1] - lo[1]+ 1;

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
  GA_Terminate();
  
}


