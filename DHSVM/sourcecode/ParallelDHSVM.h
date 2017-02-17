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
 * LAST CHANGE: 2017-02-16 10:52:36 d3g096
 * COMMENTS:
 */

#ifndef _ParallelDHSVM_h_
#define _ParallelDHSVM_h_

#include "sizeofnt.h"
#include "data.h"

int Global2Local(MAPSIZE *Map, int globalx, int globaly, int *localx, int *localy);
void Local2Global(MAPSIZE *Map, int localx, int localy, int *globalx, int *globaly);
void ParallelInitialize(int *argc, char ***argv);
void DomainDecomposition(MAPSIZE *global, MAPSIZE *local);
void DomainSummary(MAPSIZE *global, MAPSIZE *local);
int ParallelRank(void);
int ParallelSize(void);
void ParallelBarrier(void);
int GA_Type(int NumberType);
int GA_Duplicate_type(int oga, char *nname, int ntype);
void GA_Put_one(int ga, MAPSIZE *Map, int x, int y, void *value);
void GA_Acc_one(int ga, MAPSIZE *Map, int x, int y, void *value, void *alpha);
void GA_Get_one(int ga, MAPSIZE *Map, int x, int y, void *value);

void ParallelFinalize(void);

extern const int gaXdim, gaYdim;

#endif

