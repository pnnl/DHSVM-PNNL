#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int diffusion(double *s, double *h, double *Dy_m, double *Dy_p, double *Dx_m, double *Dx_p, double A_GLEN, double C_SLIDE)
{

  // Notes

  // 2. Hardwiring of ice flow exponents n_GLEN = 3 and m_SLIDE = 2      REMOVE THIS QUIRK 

  double  *h_ic_jc;
  double  *h_ic_jp;
  double  *h_ip_jc;
  double  *h_ic_jm;
  double  *h_im_jc;
  double  *s_ic_jc;
  double  *s_ip_jc;
  double  *s_im_jc;
  double  *s_ip_jp;
  double  *s_im_jp;
  double  *s_ic_jp;
  double  *s_ic_jm;
  double  *s_ip_jm;
  double  *h_ic_JP;
  double  *h_ic_JP_up;
  double  *h_IP_jc;
  double  *h_IP_jc_up;
  double  *ds_dx_ic_JP;
  double  *ds_dx_IP_jc;
  double  *ds_dy_ic_JP;
  double  *ds_dy_IP_jc;
  double  *S2_ic_JP;
  double  *S2_IP_jc;

  double  A_tilde;
  double  C_tilde;

  unsigned long k;

  extern unsigned long N;
  extern unsigned long *ic_jc;
  extern unsigned long *im_jc; 
  extern unsigned long *ip_jc;
  extern unsigned long *ic_jm;
  extern unsigned long *ic_jp;
  extern unsigned long *im_jm;
  extern unsigned long *im_jp;
  extern unsigned long *ip_jm;
  extern unsigned long *ip_jp;

  extern double  RHO;
  extern double  g;
  // extern double  A_GLEN;    // not used but should be
  // extern double  C_SLIDE;   // not used but should be
  extern double  K0_eps;

  extern double  dx;
 
  A_tilde = 2*A_GLEN*pow(RHO*g,3)/(5*pow(dx,2));               // hardwired exponents here
  C_tilde = C_SLIDE*pow(RHO*g,2)/pow(dx,2);

  if((h_ic_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_ic_jc\n");
    exit(0);
    }
  if((h_ic_jp = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_ic_jp\n");
    exit(0);
    }
  if((h_ip_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_ip_jc\n");
    exit(0);
    }
  if((h_ic_jm = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_ic_jm\n");
    exit(0);
    }
  if((h_im_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_im_jc\n");
    exit(0);
    }
  if((s_ic_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_ic_jc\n");
    exit(0);
    }
  if((s_ip_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_ip_jc\n");
    exit(0);
    }
  if((s_im_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_im_jc\n");
    exit(0);
    }
  if((s_ip_jp = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_ip_jp\n");
    exit(0);
    }
  if((s_im_jp = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_im_jp\n");
    exit(0);
    }
  if((s_ic_jp = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_ic_jp\n");
    exit(0);
    }
  if((s_ic_jm = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_ic_jm\n");
    exit(0);
    }
  if((s_ip_jm = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for s_ip_jm\n");
    exit(0);
    }
  if((h_ic_JP = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_ic_JP\n");
    exit(0);
    }
  if((h_ic_JP_up = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_ic_JP_up\n");
    exit(0);
    }
  if((h_IP_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_IP_jc\n");
    exit(0);
    }
  if((h_IP_jc_up = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for h_IP_jc_up\n");
    exit(0);
    }
  if((ds_dx_ic_JP = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for ds_dx_ic_JP\n");
    exit(0);
    }
  if((ds_dx_IP_jc = malloc(N*sizeof(double)))==NULL)
    {
    printf("diffusion(): Memory allocation failed for ds_dx_IP_jc\n");
    exit(0);
    }
  if((ds_dy_ic_JP = malloc(N*sizeof(double)))==NULL)
    {
      printf("diffusion(): Memory allocation failed for ds_dy_ic_JP\n");
      exit(0);
    }
  if((ds_dy_IP_jc = malloc(N*sizeof(double)))==NULL)
    {
      printf("diffusion(): Memory allocation failed for ds_dy_IP_jc\n");
      exit(0);
    }
  if((S2_ic_JP    = malloc(N*sizeof(double)))==NULL)
    {
      printf("diffusion(): Memory allocation failed for S2_ic_JP\n");
      exit(0);
    }
  
  S2_IP_jc    = malloc(N*sizeof(double));
  //{
  //   printf("diffusion(): Memory allocation failed for S2_IP_jc\n");
  //   exit(0);
  // }

  // Find all of the relevant surface and thickness nodes that will be used to construct the
  // diffusion vectors.

  for (k=0; k<N; k++)
    {
    h_ic_jc[k] = h[ic_jc[k]];
    h_ic_jp[k] = h[ic_jp[k]];
    h_ip_jc[k] = h[ip_jc[k]];
    h_ic_jm[k] = h[ic_jm[k]];
    h_im_jc[k] = h[im_jc[k]];

    h_IP_jc[k] = (h_ic_jc[k]+h_ip_jc[k])/2;
    h_ic_JP[k] = (h_ic_jc[k]+h_ic_jp[k])/2;

    s_ic_jc[k] = s[ic_jc[k]];
    s_ip_jc[k] = s[ip_jc[k]];
    s_im_jc[k] = s[im_jc[k]];
    s_ic_jp[k] = s[ic_jp[k]];
    s_ic_jm[k] = s[ic_jm[k]];
    s_ip_jm[k] = s[ip_jm[k]];
    s_ip_jp[k] = s[ip_jp[k]];
    s_im_jp[k] = s[im_jp[k]];

    // Calculate some surface differences

    ds_dx_IP_jc[k] = (s_ip_jc[k] - s_ic_jc[k])/dx;
    ds_dy_ic_JP[k] = (s_ic_jp[k] - s_ic_jc[k])/dx;

    ds_dx_ic_JP[k] = (s_ip_jc[k] + s_ip_jp[k] - s_im_jc[k] - s_im_jp[k])/(4*dx);
    ds_dy_IP_jc[k] = (s_ic_jp[k] - s_ic_jm[k] + s_ip_jp[k] - s_ip_jm[k])/(4*dx);

    S2_ic_JP[k] = pow(ds_dx_ic_JP[k],2) + pow(ds_dy_ic_JP[k],2) + pow(K0_eps,2);
    S2_IP_jc[k] = pow(ds_dx_IP_jc[k],2) + pow(ds_dy_IP_jc[k],2) + pow(K0_eps,2);

    // Switched JSA correction (replaces tanh-smoothed method)

    h_IP_jc_up[k]   = h_ic_jc[k];
    h_ic_JP_up[k]   = h_ic_jc[k];

    if(ds_dx_IP_jc[k]>0)
      h_IP_jc_up[k] = h_ip_jc[k];
 
    if(ds_dy_ic_JP[k]>0)
      h_ic_JP_up[k] = h_ic_jp[k];
 
    // calculate the D(i,j+1) diffusion coefficient (note hardwired flow law and sliding law exponents)

    Dy_p[k] = A_tilde*h_ic_JP_up[k]*pow(h_ic_JP[k],4)*S2_ic_JP[k]
            + C_tilde*h_ic_JP_up[k]*pow(h_ic_JP[k],2)*pow(S2_ic_JP[k],0.5);

    

    // Calculate the D(i+1,j) diffusion coefficient (note hardwired flow law and sliding law exponents)

    Dx_p[k] = A_tilde*h_IP_jc_up[k]*pow(h_IP_jc[k],4)*S2_IP_jc[k]
            + C_tilde*h_IP_jc_up[k]*pow(h_IP_jc[k],2)*pow(S2_IP_jc[k],0.5);
    // if(Dx_p[k] > 999999.0)
    //  printf("Dx_p=%f\n", Dx_p[k]);

    }
  
  for (k=0; k<N; k++)
    {
      Dy_m[k] = Dy_p[ic_jm[k]];
      Dx_m[k] = Dx_p[im_jc[k]];
     
    }

  free(h_ic_jc);
  free(h_ic_jp);
  free(h_ip_jc);
  free(h_ic_jm);
  free(h_im_jc);
  free(s_ic_jc);
  free(s_ip_jc);
  free(s_im_jc);
  free(s_ip_jp);
  free(s_im_jp);
  free(s_ic_jp);
  free(s_ic_jm);
  free(s_ip_jm);
  free(h_ic_JP);
  free(h_ic_JP_up);
  free(h_IP_jc);
  free(h_IP_jc_up);
  free(ds_dx_ic_JP);
  free(ds_dx_IP_jc);
  free(ds_dy_ic_JP);
  free(ds_dy_IP_jc);
  free(S2_ic_JP);
  free(S2_IP_jc);
  }
