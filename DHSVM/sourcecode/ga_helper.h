/*
 * SUMMARY:      ga_helper.h
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
 * LAST CHANGE: 2019-11-19 09:23:10 d3g096
 * COMMENTS:
 */

#ifndef _ga_helper_h_
#define _ga_helper_h_

#include <mpi.h>
#include <ga.h>

extern const int gaXdim, gaYdim;

#ifdef __cplusplus
extern "C" {
#endif
void ParallelInitialize(int *argc, char ***argv);
int ParallelRank(void);
int ParallelSize(void);
void ParallelBarrier(void);
void ParallelFinalize(void);

int GA_Type(int NumberType);
int GA_Duplicate_type(int oga, char *nname, int ntype);
void GA_Inquire_irreg_distr(int ga, int *mapc, int *nblk);
#ifdef __cplusplus
}
#endif

#endif
