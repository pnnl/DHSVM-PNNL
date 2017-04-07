/*
 * SUMMARY:      ParallelDHSVM.h
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
 * LAST CHANGE: 2017-04-06 08:00:22 d3g096
 * COMMENTS:
 */

#ifndef _ParallelDHSVM_h_
#define _ParallelDHSVM_h_

#include "sizeofnt.h"
#include "data.h"

int Global2Local(MAPSIZE *Map, int globalx, int globaly, int *localx, int *localy);
void Local2Global(MAPSIZE *Map, int localx, int localy, int *globalx, int *globaly);
void ParallelInitialize(int *argc, char ***argv);
void MaskedDomainDecomposition(MAPSIZE *gmap, MAPSIZE *lmap, MAPSIZE* nmap, 
                               unsigned char *mask);
void SimpleDomainDecomposition(MAPSIZE *global, MAPSIZE *local);
void DomainSummary(MAPSIZE *global, MAPSIZE *local);
int ParallelRank(void);
int ParallelSize(void);
void ParallelBarrier(void);
int GA_Type(int NumberType);
int GA_Duplicate_type(int oga, char *nname, int ntype);
void GA_Put_one(int ga, MAPSIZE *Map, int x, int y, void *value);
void GA_Acc_one(int ga, MAPSIZE *Map, int x, int y, void *value, void *alpha);
void GA_Get_one(int ga, MAPSIZE *Map, int x, int y, void *value);

struct ga_patch_ {
  int NY, NX;                   /**< dimensions of the allocated patch */
  int ixoff;                    /**< add this to column index, when accessing path */
  int iyoff;                    /**< add this to row index, when accessing patch */
  float **patch;                /**< Allocated 2D[y][x] filled with current GA w/ ghosts  */
};
typedef struct ga_patch_ GA_Patch;
void GA_Alloc_patch(int ga, MAPSIZE *Map, GA_Patch *p);
void GA_Alloc_patch_ghost(int ga, MAPSIZE *Map, GA_Patch *p);
void GA_Get_patch(int ga, MAPSIZE *Map, GA_Patch *p);
void GA_Acc_patch(int ga, MAPSIZE *Map, GA_Patch *p);
void GA_Put_patch(int ga, MAPSIZE *Map, GA_Patch *p);
void GA_Free_patch(GA_Patch *p);

int Collect2DMatrixGA(void *LocalMatrix, int NumberType, MAPSIZE *Map);

void ParallelFinalize(void);

extern const int gaXdim, gaYdim;

#endif

