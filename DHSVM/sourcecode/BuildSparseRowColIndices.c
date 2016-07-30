#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This function does not exploit the fact that the sparse matrix is symmetric
// This should be exploited once code as been demonstrated to work properly

void BuildSparseRowColIndices(int *row, int *col)
  {
  extern int N;
  extern int *ic_jc;
  extern int *ic_jm;
  extern int *ip_jc;

  extern int *im_jc;   // delete when symmetry exploited
  extern int *ic_jp;   // delete when symmetry exploited

  int k, cnt=0;
  
  for (k=0; k<N; k++) 
    {
    row[cnt]   = ic_jc[k];
    row[cnt+1] = ic_jc[k]; 
    row[cnt+2] = ic_jc[k];
    row[cnt+3] = ic_jc[k];     // delete when symmetry exploited
    row[cnt+4] = ic_jc[k];     // delete when symmetry exploited

//  col[cnt]   = ic_jc[k];     // reinstate when symmetry exploited
//  col[cnt+1] = ic_jm[k];     // reinstate when symmetry exploited
//  col[cnt+2] = ip_jc[k];     // reinstate when symmetry exploited

    col[cnt]   = im_jc[k];     // delete following lines when symmetry exploited
    col[cnt+1] = ic_jp[k];
    col[cnt+2] = ic_jc[k];
    col[cnt+3] = ic_jm[k];
    col[cnt+4] = ip_jc[k];

    cnt       += 5;            // change to c+=3 when symmetry exploited
    }
  }
