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
 * LAST CHANGE: 2017-01-24 08:51:07 d3g096
 * COMMENTS:
 */

#ifndef _ParallelDHSVM_h_
#define _ParallelDHSVM_h_

#include "data.h"

void ParallelInitialize(int *argc, char ***argv);
void DomainDecomposition(MAPSIZE *global, MAPSIZE *local);
void DomainSummary(MAPSIZE *global, MAPSIZE *local);
int ParallelRank(void);
int ParallelSize(void);
void ParallelBarrier(void);
void ParallelFinalize(void);

#endif

