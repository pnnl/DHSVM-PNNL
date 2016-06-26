#ifdef HAVE_GLACIER

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cs.h"
#include <time.h>
#include "constants.h"
#include "settings.h"
#include "massenergy.h"
#include "functions.h"
#include "snow.h"

int           GetBalance(double yr, double *s, double *b_dot_ppt, double *b_dot_melt, int N);
int           BuildDiffusivity(double *s, double *h, double *Dy_m, double *Dy_p, double *Dx_m, double *Dx_p);
void          BuildSparseRowColIndices(int *row, int *col);
void          BuildSparseArrayElements(double *s, double *b_dot, double *Dx_m, double *Dx_p, double *Dy_m, double *Dy_p, int *row, int *col, double *A_val, double *C_vec, double dt);
//static void   resid (int ok, cs *A, double *x, double *b, double *r);
//static double tic (void);
//static double btoc (double t);

//cs            *cs_compress (const cs *T); 

int RunGlacier(double *b, double *s_init, double *s_out, double yr_min, double yr_max, double dt_yr, 
               double *b_dot, OPTIONSTRUCT * Options)
{ 
  double  *s_inp;           // surface topography at start of time step (m)
  double  *h;               // ice thickness (m)
  //double  *h_out;           // LINT for debugging
  // double  *b_dot;           // net ice-equivalent glacier mass balance (m/yr)
  double  *b_dot_melt;      // ice-equivalent melt rate (m/yr) 
  double  *b_dot_ppt;       // ice-equivalent solid precipitation rate (m/yr)
  double  *Dx_m;            // nonlinear ice "diffusion" coefficient
  double  *Dx_p;            // nonlinear ice "diffusion" coefficient
  double  *Dy_m;            // nonlinear ice "diffusion" coefficient
  double  *Dy_p;            // nonlinear ice "diffusion" coefficient

  double  yr;

  extern int *ip_jc;
  extern int *im_jc;
  extern int *ic_jp;
  extern int *ic_jm;
  extern int *ic_jc;
  extern int *ip_jp;
  extern int *im_jp;
  extern int *ip_jm;
  extern int *im_jm;
  //extern int nx;
  //extern int ny;
  int  *row;
  int  *col;

  double *A_val;
  double *C;
  double *X;

  cs   *A_t;      // sparse matrix in triplet form
  cs   *A_c;      // sparse matrix in compressed column matrix form

  int *i_mask;    // ice mask
 
  extern int  N; 
  int    k, ok, order;

  double h_max,alpha_I ;
  clock_t tic = clock();	
  s_inp  = malloc(N*sizeof(double));
  h      = malloc(N*sizeof(double));
  //b_dot  = malloc(N*sizeof(double));
  b_dot_melt = malloc(N*sizeof(double));
  b_dot_ppt  = malloc(N*sizeof(double));
  Dx_m   = malloc(N*sizeof(double));
  Dx_p   = malloc(N*sizeof(double));
  Dy_m   = malloc(N*sizeof(double));
  Dy_p   = malloc(N*sizeof(double));

  C      = malloc(N*sizeof(double));
  X      = malloc(N*sizeof(double));        // ?? LINT ??
  
  i_mask = malloc(N*sizeof(int));
  
  A_t    = cs_spalloc(N, N, 5*N, 1, 1);
  // A_c    = cs_spalloc(N, N, 5*N, 1, 0);
 
  row    = A_t->i;
  col    = A_t->p;
  A_val  = A_t->x;

  A_t->nz = 5*N;                          // Here I assume nz = nzmax

  for (k=0; k<N; k++){
    s_inp[k] = s_init[k];
  }
 
  yr      = yr_min;

  BuildSparseRowColIndices(row, col);

  while(1)
    {
      //GetBalance(yr, s_inp, b_dot_ppt, b_dot_melt, N);   // update glacier mass balance forcings
   
    for (k=0; k<N; k++)
      {   
	// b_dot[k] = b_dot_ppt[k] + b_dot_melt[k];
	
	//if( b_dot[k] > 1)
	//if(  b_dot[k]> 0)
	// printf(" bdot =  %f",   b_dot[k]);
	h[k]     = s_inp[k]-b[k];
      }
        
    BuildDiffusivity(s_inp, h, Dy_m, Dy_p, Dx_m, Dx_p);
    BuildSparseArrayElements(s_inp, b_dot, Dx_m, Dx_p, Dy_m, Dy_p, row, col, A_val, C, dt_yr);
    A_t->x = A_val;
    A_c = cs_compress(A_t);     // change from triplet to compressed column format. space is allocated for A_c in cs_compress.c which needs to be freed after cs_cholsol.
    cs_dupl(A_c);              // sum and remove duplicate entries

    order = 1;                 // input parameter for cs_cholsol

    for (k=0; k<N; k++){
      X[k] = C[k];
      //printf("X= %f\n",x[k]);
    }
    ok     = cs_cholsol (order,A_c, X);
    cs_spfree(A_c);
    
    if (ok!=1)
      {
	printf("run_glacier(): Failure of cs_cholsol() Cholesky solver at time = %.2f yr\n", yr);
	return 0;
      }

    h_max = 0;
    alpha_I = 0.0;
    for (k=0; k<N; k++)
      {
      s_out[k] = X[k];
      if(s_out[k]<=b[k])
        {
	  s_out[k]  = b[k];
	  i_mask[k] = 0;       // evaluate the instantaneous ice mask (not used in the glacier model but needed for hydrological model)
        }
      else
        {
	  if (s_out[k]- b[k]>h_max)
	    h_max = s_out[k]-b[k];
	  i_mask[k] = 1;
        }
      s_inp[k] = s_out[k];
      }

    yr    = yr + dt_yr;
    
    for (k=0; k<N; k++)
      {
	if (s_out[k]>b[k])
	  alpha_I += 1;
      }

    alpha_I = 100*alpha_I/(double)N;
if (Options->Glacier == GLSPINUP){
    printf(" Ice covered area      = %8.3f percent\n", alpha_I);
    printf("Time %.2f yr: max(h) = %.3f m\n", yr, h_max);
}
    if(yr >=yr_max)
      break;
    }
  clock_t toc = clock();
  printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  free(s_inp); 
  free(h);
  //free(b_dot);
  free(b_dot_melt);
  free(b_dot_ppt);
  free(Dx_m);
  free(Dx_p);
  free(Dy_m);
  free(Dy_p);
  //free(row);
  //free(col);
  //free(A_val);
  free(C);
  free(X);
  free(i_mask);

  cs_spfree(A_t);    // Fails for unknown reasons. Hence rely on default free on exit
  // cs_spfree(A_c);    // Fails for unknown reasons. Hence rely on default free on exit
  free(ip_jc);
  free(im_jc);
  free(ic_jp);
  free(ic_jm);
  free(ic_jc);
  free(ip_jp);
  free(im_jp);
  free(ip_jm);
  free(im_jm);
  
  return 1;
}
#endif
