#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This function does not exploit the fact that the sparse matrix is symmetric
// This should be exploited once code as been demonstrated to work properly

void BuildSparseArrayElements(double *s, double *b_dot, double *Dx_m, double *Dx_p, double *Dy_m, double *Dy_p, int *row, int *col, double *A_val, double *C_vec, double dt)
  {
  extern int N;
  //extern int *ic_jc;
  extern int *ic_jm;
  extern int *ip_jc;
  extern int *ic_jp;       // Delete when symmetry exploited
  extern int *im_jc;       // Delete when symmetry exploited
  extern double OMEGA;
  int k, cnt=0;

  double *D_sum;

  D_sum = malloc(N*sizeof(double));
 
  for (k=0; k<N; k++)
    {
    D_sum[k] = Dx_m[k] + Dx_p[k] + Dy_m[k] + Dy_p[k];
    C_vec[k] = (1-OMEGA)*(Dx_m[k]*s[im_jc[k]] + Dx_p[k]*s[ip_jc[k]] + Dy_m[k]*s[ic_jm[k]] + Dy_p[k]*s[ic_jp[k]]) + (1/dt - (1-OMEGA)*D_sum[k])*s[k] + b_dot[k];
    }

  for (k=0; k<N; k++) 
    {
//  A_val[cnt]   = 1/dt + OMEGA*D_sum[k];    // reinstate when symmetry is exploited
//  A_val[cnt+1] = -OMEGA*Dy_m[k];           // reinstate when symmetry is exploited
//  A_val[cnt+2] = -OMEGA*Dx_p[k];           // reinstate when symmetry is exploited

    A_val[cnt]   = -OMEGA*Dx_m[k];
    A_val[cnt+1] = -OMEGA*Dy_p[k];
    A_val[cnt+2] = 1/dt + OMEGA*D_sum[k];
    A_val[cnt+3] = -OMEGA*Dy_m[k];
    A_val[cnt+4] = -OMEGA*Dx_p[k];

    cnt += 5;                                // change to cnt+=3 when symmetry is exploited
    }
  free(D_sum);
  }
