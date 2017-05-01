/*
 * SUMMARY:      array_alloc.h
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    May 2017
 * DESCRIPTION: A collection of routines to allocate 2D arrays in a
 * consistent (and hopefully efficient) manner
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2017-05-01 11:19:13 d3g096
 * COMMENTS:
 */

#ifndef _array_alloc_h_
#define _array_alloc_h_

#include <sys/types.h>

extern float **calloc_2D_float(int NY, int NX);
extern void free_2D_float(float **p);

#endif

