#include <stdlib.h>
#include <stdio.h>  
void SetupIndexArrays(void)
{
  int i, j, cnt;

  extern int *ip_jc;
  extern int *im_jc;
  extern int *ic_jp;
  extern int *ic_jm;
  extern int *ic_jc;
  extern int *ip_jp;
  extern int *im_jp;
  extern int *ip_jm;
  extern int *im_jm;
  extern int nx;
  extern int ny;
  extern int N;

  int *ic;
  int *im;
  int *ip;

  int *jc;
  int *jm;
  int *jp;

  int k;

  ic_jc = malloc(N*sizeof(int));
 
  ip_jc = malloc(N*sizeof(int));
  im_jc = malloc(N*sizeof(int));
  ic_jp = malloc(N*sizeof(int));
  ic_jm = malloc(N*sizeof(int));

  ip_jp = malloc(N*sizeof(int));
  im_jp = malloc(N*sizeof(int));
  ip_jm = malloc(N*sizeof(int));
  im_jm = malloc(N*sizeof(int));

  ic = malloc(nx*sizeof(int));
  im = malloc(nx*sizeof(int));
  ip = malloc(nx*sizeof(int));

  jc = malloc(ny*sizeof(int));
  jm = malloc(ny*sizeof(int));
  jp = malloc(ny*sizeof(int));

  cnt = 0;
  
  for (j=0; j<ny; j++)
    jc[j] = j;

  for (j=0; j<ny-1; j++)
    jm[j] = j+1;
  jm[ny-1] = ny-1;

  jp[0] = 0;
  for (j=1; j<ny; j++)
    jp[j] = j-1;

  for (i=0; i<nx; i++)
   ic[i] = i;

  im[0] = 0;
  for (i=1; i<nx; i++)
    im[i] = i-1;

  for (i=0; i<nx-1; i++)
    ip[i] = i+1;
  ip[nx-1] = nx-1;

  for (k=0; k<N; k++)
    ic_jc[k] = k;

  k = 0;
 
  for (i=0; i<nx; i++)
    {
      for (j=0; j<ny; j++)
	{
	  ic_jp[k] = ny*ic[i] + jp[j];
	  ic_jm[k] = ny*ic[i] + jm[j];
	  im_jc[k] = ny*im[i] + jc[j];
	  ip_jc[k] = ny*ip[i] + jc[j];
      
	  im_jm[k] = ny*im[i] + jm[j];
	  im_jp[k] = ny*im[i] + jp[j];
      
	  ip_jm[k] = ny*ip[i] + jm[j];
	  ip_jp[k] = ny*ip[i] + jp[j];

	  k++;

	  
	}
    }
 
  free(ic);
  free(im);
  free(ip);
  free(jc);
  free(jm);
  free(jp);
}

