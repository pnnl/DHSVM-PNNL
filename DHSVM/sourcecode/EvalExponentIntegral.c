/*
 * SUMMARY:      EvalExponentIntegral.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    2001
 * DESCRIPTION:  
 * DESCRIP-END.
 * FUNCTIONS:   
 * COMMENTS:
 * $Id: EvalExponentIntegral.c,v 1.5 2006/10/03 22:50:22 nathalie Exp $     
 */

#include <math.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"
#define EULER 0.57721566
#define MAXIT 100
#define FPMIN 1.0e-30
#define EPS 1.0e-7

float evalexpint(int n, float x)
{

  int i, ii, nm1;
  float a, b, c, d, del, fact, h, psi, ans;

  nm1 = n - 1;
  if (n < 0 || x < 0.0 || (fequal(x, 0.0) && (n == 0 || n == 1))) {
    printf("bad args to EvalExpInt \n");
    exit(-1);
  }
  else {
    if (n == 0)
      ans = exp(-x) / x;
    else {
      if (x == 0)
	ans = 1.0 / nm1;
      else {
	if (x > 1.0) {
	  b = x + n;
	  c = 1.0 / FPMIN;
	  d = 1.0 / b;
	  h = d;
	  for (i = 1; i <= MAXIT; i++) {
	    a = -i * (nm1 + i);
	    b += 2.0;
	    d = 1.0 / (a * d + b);
	    c = b + a / c;
	    del = c * d;
	    h *= del;
	    if (fabs(del - 1.0) < EPS) {
	      ans = h * exp(-x);
	      return ans;
	    }
	  }
	  printf("error in expint\n");
	  exit(-1);
	}
	else {
	  ans = (nm1 != 0 ? 1.0 / nm1 : -log(x) - EULER);
	  fact = 1.0;
	  for (i = 1; i <= MAXIT; i++) {
	    fact *= -x / i;
	    if (i != nm1)
	      del = -fact / (i - nm1);
	    else {
	      psi = -EULER;
	      for (ii = 1; ii <= nm1; ii++)
		psi += 1.0 / ii;
	      del = fact * (-log(x) + psi);
	    }
	    ans += del;
	    if (fabs(del) < fabs(ans) * EPS)
	      return ans;
	  }
	  printf("error in routine \n");
	  exit(-1);
	}
      }
    }
  }
  return ans;
}
