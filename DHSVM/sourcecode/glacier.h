/*
 * SUMMARY:      glacier.h - header file for DHSVM snow routines
 * USAGE:        Part of Glacier Model
 *
 * AUTHOR:       Bibi 
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              bibinaz@u.washington.edu
 * ORIG-DATE:    16-Sept-2011 
 * DESCRIPTION:  header file for DHSVM glacier routines
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: glacier.h,v 1.1.1.1 2011/09/16 03:45:00 bibinaz Exp $     
 */

#ifndef GLACIER_H
#define GLACIER_H

#include <stdarg.h>
#include "data.h"



void InitGlacierMap(MAPSIZE *Map, GLPIX ***GlacierMap);

void SetupIndexArrays(void);

//float Run_glacier(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, double *s_half, double *s_tmp,  double *s_inp,   double *h,  double *b, double *a_net, double *A_tri, double *B_tri, double *C_tri, double *E_tri, double *A_tmp, double *B_tmp, double *C_tmp, double *E_tmp, double *Dx_m, double *Dx_p,  double *Dy_m, double *Dy_p, double *s_out, double dt_yr);          
//void Run_glacier(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow,GLACIERPIX **GlacierMap, double dt_yr, double yr);
int calc_ela(MAPSIZE *Map, TOPOPIX **TopoMap,SNOWPIX **Snow, GLPIX **GlacierMap, DATE *Current, DUMPSTRUCT *Dump);
int RunGlacier(double *b, double *s_init, double *s_out, double yr_min, double yr_max, double dt_yr, double *b_dot, OPTIONSTRUCT * Options);
void main_gl(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, GLPIX **GlacierMap, double dt_year, double year_min, double year_max,DATE *Current, DUMPSTRUCT *Dump,OPTIONSTRUCT * Options);
void main_spinup(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, GLPIX **GlacierMap, double dt_year, double year_min, double year_max,DATE *Current, DUMPSTRUCT *Dump, OPTIONSTRUCT * Options);
void gl_massbalance(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow, GLPIX **GlacierMap, double dt_year, double year_min, double year_max,DATE *Current, DUMPSTRUCT *Dump);
//int GetBalance(double yr, double *s, double *b_dot_ppt, double *b_dot_melt, int N);
//float update_glmassbalance(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **Snow);
#endif
