#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int BuildDiffusivity(double *s, double *h, double *Dy_m, double *Dy_p, double *Dx_m, double *Dx_p)
{
  double  *h_ic_jc;
  double  *h_ic_jp;
  double  *h_ip_jc;
  double  *h_ic_jm;
  double  *h_im_jc;
  double  *s_ic_jc;
  double  *s_ip_jc;
  double  *s_im_jc;
  double  *s_im_jm;
  double  *s_im_jp;
  double  *s_ic_jp;
  double  *s_ic_jm;
  double  *s_ip_jm;

  double  *h_ic_JC;
  double  *h_ic_JC_up;
  double  *h_IC_jc;
  double  *h_IC_jc_up;

  double  *ds_dx_ic_JC;
  double  *ds_dx_IC_jc;
  double  *ds_dy_ic_JC;
  double  *ds_dy_IC_jc;

  double  *S2_ic_JC;
  double  *S2_IC_jc;

  int k;

  extern int N;
  extern int *ic_jc;
  extern int *im_jc; 
  extern int *ip_jc;
  extern int *ic_jm;
  extern int *ic_jp;
  extern int *im_jm;
  extern int *im_jp;
  extern int *ip_jm;
  //extern int *ip_jp;

  //extern double  RHO;
  //extern double  g;
  extern double  A_tilde;
  extern double  C_tilde;
  //extern double  nm_half;
  extern double  np1;
  //extern double  mm_half;
  extern double  m1;
  extern double  K0_eps;

  extern double  dx;

  h_ic_jc = malloc(N*sizeof(double));
  h_ic_jp = malloc(N*sizeof(double));
  h_ip_jc = malloc(N*sizeof(double));
  h_ic_jm = malloc(N*sizeof(double));
  h_im_jc = malloc(N*sizeof(double));
  s_ic_jc = malloc(N*sizeof(double));
  s_ip_jc = malloc(N*sizeof(double));
  s_im_jc = malloc(N*sizeof(double));
  s_im_jm = malloc(N*sizeof(double));
  s_im_jp = malloc(N*sizeof(double));
  s_ic_jp = malloc(N*sizeof(double));
  s_ic_jm = malloc(N*sizeof(double));
  s_ip_jm = malloc(N*sizeof(double));
  h_ic_JC = malloc(N*sizeof(double));
  h_ic_JC_up = malloc(N*sizeof(double));
  h_IC_jc = malloc(N*sizeof(double));
  h_IC_jc_up = malloc(N*sizeof(double));
  ds_dx_ic_JC = malloc(N*sizeof(double));
  ds_dx_IC_jc = malloc(N*sizeof(double));
  ds_dy_ic_JC = malloc(N*sizeof(double));
  ds_dy_IC_jc = malloc(N*sizeof(double));
  S2_ic_JC    = malloc(N*sizeof(double));
  S2_IC_jc    = malloc(N*sizeof(double));

  // Find all of the relevant surface and thickness nodes that will be used to construct the
  // diffusion vectors.

  for (k=0; k<N; k++)
    {
    h_ic_jc[k] = h[ic_jc[k]];
    h_ic_jp[k] = h[ic_jp[k]];
    h_ip_jc[k] = h[ip_jc[k]];
    h_ic_jm[k] = h[ic_jm[k]];
    h_im_jc[k] = h[im_jc[k]];

    h_IC_jc[k] = (h_im_jc[k]+h_ic_jc[k])/2;
    h_ic_JC[k] = (h_ic_jm[k]+h_ic_jc[k])/2;

    s_ic_jc[k] = s[ic_jc[k]];
    s_ip_jc[k] = s[ip_jc[k]];
    s_im_jc[k] = s[im_jc[k]];
    s_ic_jp[k] = s[ic_jp[k]];
    s_ic_jm[k] = s[ic_jm[k]];
    s_ip_jm[k] = s[ip_jm[k]];
    s_im_jp[k] = s[im_jp[k]];
    s_im_jm[k] = s[im_jm[k]];

    // Calculate some surface differences

    ds_dx_IC_jc[k] = (s_ic_jc[k] - s_im_jc[k])/dx;
    ds_dy_IC_jc[k] = (s_ic_jp[k] + s_im_jp[k] - s_ic_jm[k] - s_im_jm[k])/(4*dx);
    ds_dx_ic_JC[k] = (s_ip_jc[k] + s_ip_jm[k] - s_im_jc[k] - s_im_jm[k])/(4*dx);
    ds_dy_ic_JC[k] = (s_ic_jc[k] - s_ic_jm[k])/dx;

    S2_ic_JC[k] = pow(ds_dx_ic_JC[k],2) + pow(ds_dy_ic_JC[k],2) + pow(K0_eps,2);
    S2_IC_jc[k] = pow(ds_dx_IC_jc[k],2) + pow(ds_dy_IC_jc[k],2) + pow(K0_eps,2);

    // Switched JSA correction (replaces tanh-smoothed method)

    h_IC_jc_up[k]   = h_im_jc[k];
    h_ic_JC_up[k]   = h_ic_jm[k];

    if (s_ic_jc[k]>s_im_jc[k])
      h_IC_jc_up[k] = h_ic_jc[k];

    if (s_ic_jc[k]>s_ic_jm[k])
      h_ic_JC_up[k] = h_ic_jc[k];

    // calculate the D(i,j+1) diffusion coefficient (note hardwired flow law and sliding law exponents)

    Dy_m[k] = A_tilde*h_ic_JC_up[k]*pow(h_ic_JC[k],np1)*S2_ic_JC[k]
            + C_tilde*h_ic_JC_up[k]*pow(h_ic_JC[k],m1)*pow(S2_ic_JC[k],0.5);

    // Calculate the D(i+1,j) diffusion coefficient (note hardwired flow law and sliding law exponents)

    Dx_m[k] = A_tilde*h_IC_jc_up[k]*pow(h_IC_jc[k],np1)*S2_IC_jc[k]
            + C_tilde*h_IC_jc_up[k]*pow(h_IC_jc[k],m1)*pow(S2_IC_jc[k],0.5);

    }
  
  for (k=0; k<N; k++)
    {
    Dy_p[k] = Dy_m[ic_jp[k]];
    Dx_p[k] = Dx_m[ip_jc[k]];
    }

  free(h_ic_jc);
  free(h_ic_jp);
  free(h_ip_jc);
  free(h_ic_jm);
  free(h_im_jc);
  free(s_ic_jc);
  free(s_ip_jc);
  free(s_im_jm);
  free(s_im_jc);
  free(s_im_jp);
  free(s_ic_jp);
  free(s_ic_jm);
  free(s_ip_jm);
  free(h_ic_JC);
  free(h_ic_JC_up);
  free(h_IC_jc);
  free(h_IC_jc_up);
  free(ds_dx_ic_JC);
  free(ds_dx_IC_jc);
  free(ds_dy_ic_JC);
  free(ds_dy_IC_jc);
  free(S2_ic_JC);
  free(S2_IC_jc);
}
