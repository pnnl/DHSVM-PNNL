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
 * LAST CHANGE: 2017-02-06 11:32:53 d3g096
 * COMMENTS:
 */

#ifndef _ParallelDHSVM_h_
#define _ParallelDHSVM_h_

#include "sizeofnt.h"
#include "data.h"

void ParallelInitialize(int *argc, char ***argv);
void DomainDecomposition(MAPSIZE *global, MAPSIZE *local);
void DomainSummary(MAPSIZE *global, MAPSIZE *local);
int ParallelRank(void);
int ParallelSize(void);
void ParallelBarrier(void);
int GA_Type(int NumberType);
int GA_Duplicate_type(int oga, char *nname, int ntype);
void ParallelFinalize(void);

extern const int gaXdim, gaYdim;

#endif

