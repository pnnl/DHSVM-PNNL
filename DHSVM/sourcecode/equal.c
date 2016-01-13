/*
 * SUMMARY:      equal.c - Determine whether two float or two doubles are equal
 *                         to within machine precision
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Wed Feb  3 08:39:13 1999
 * DESCRIPTION:  
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     Adapted from netcdf library
 * $Id: equal.c,v 1.4 2003/07/01 21:26:29 olivier Exp $     
 */

/******************************************************************************/
/*				    INCLUDES                                  */
/******************************************************************************/
#ifdef TEST_EQUAL
#include <stdio.h>
#include <stdlib.h>
#endif
/* if float.h is not available on your system, comment out the following line,
   and uncomment the line after that */
#include <float.h>
/* #define NO_FLOAT_H */
#include <math.h>
#include "functions.h"

/******************************************************************************/
/*				GLOBAL VARIABLES                              */
/******************************************************************************/
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef FLT_DIGITS
#define FLT_DIGITS 7		/* default sig. digits for float data */
#endif
#ifndef DBL_DIGITS
#define DBL_DIGITS 15		/* default sig. digits for double data */
#endif
static double double_eps;
static float float_eps;
static int initeps = 0;

/******************************************************************************/
/******************************************************************************/
/*			           FUNCTIONS                                  */
/******************************************************************************/
/******************************************************************************/

static double double_epsilon(void);
static float float_epsilon(void);
static void init_epsilons(void);

/******************************************************************************/
/*				     dequal                                   */
/******************************************************************************/
unsigned char dequal(double a, double b)
{
  if (!initeps) {		/* make sure epsilons get initialized */
    init_epsilons();
    initeps = 1;
  }

  /* Two double values only need to be equal to within machine precision */
  if ((a > 0) == (b > 0) &&	/* prevents potential overflow */
      (ABSVAL(a - b) <= ABSVAL(double_eps * b)))
    return TRUE;
  else
    return FALSE;
}

/******************************************************************************/
/*				     fequal                                   */
/******************************************************************************/
unsigned char fequal(float a, float b)
{
  if (!initeps) {		/* make sure epsilons get initialized */
    init_epsilons();
    initeps = 1;
  }

  /* Two float values anly need to be equal to within machine precision */
  if ((a > 0) == (b > 0) &&	/* prevents potential overflow */
      (ABSVAL(a - b) <= ABSVAL(float_eps * b)))
    return TRUE;
  else
    return FALSE;
}

/******************************************************************************/
/*				 double_epsilon                               */
/******************************************************************************/
static double double_epsilon(void)
{
  double double_eps;
#ifndef NO_FLOAT_H
  double_eps = DBL_EPSILON;
#else /* NO_FLOAT_H */
  {
    double etop, ebot, eps;
    double one = 1.0;
    double two = 2.0;
    etop = 1.0;
    ebot = 0.0;
    eps = ebot + (etop - ebot) / two;
    while (eps != ebot && eps != etop) {
      double epsp1;

      epsp1 = one + eps;
      if (epsp1 > one)
	etop = eps;
      else
	ebot = eps;
      eps = ebot + (etop - ebot) / two;
    }
    double_eps = two * etop;
  }
#endif /* NO_FLOAT_H */
  return double_eps;
}

/******************************************************************************/
/*				 float_epsilon                                */
/******************************************************************************/
static float float_epsilon(void)
{
  float float_eps;
#ifndef NO_FLOAT_H
  float_eps = FLT_EPSILON;
#else /* NO_FLOAT_H */
  {
    float etop, ebot, eps;
    float one = 1.0;
    float two = 2.0;
    etop = 1.0;
    ebot = 0.0;
    eps = ebot + (etop - ebot) / two;
    while (eps != ebot && eps != etop) {
      float epsp1;

      epsp1 = one + eps;
      if (epsp1 > one)
	etop = eps;
      else
	ebot = eps;
      eps = ebot + (etop - ebot) / two;
    }
    float_eps = two * etop;
  }
#endif /* NO_FLOAT_H */
  return float_eps;
}

/******************************************************************************/
/*				 init_epsilons                                */
/******************************************************************************/
static void init_epsilons(void)
{
  float_eps = float_epsilon();
  double_eps = double_epsilon();
}

/*******************************************************************************
  Test main. Compile by typing:
  gcc -Wall -g -o test_equal -DTEST_EQUAL equal.c 
  then run the program by typing test_equal
*******************************************************************************/
#ifdef TEST_EQUAL
int main(void)
{
  char *str[] = { "FALSE", "TRUE" };
  double a;
  double b;
  float c;
  float d;

  printf("Testing double ...\n");
  a = 0.;
  b = 0.;
  printf("a = %15.13g, b = %15.13g, equal? %s\n", a, b, str[dequal(a, b)]);
  a = 0.;
  b = DBL_EPSILON;;
  printf("a = %15.13g, b = %15.13g, equal? %s\n", a, b, str[dequal(a, b)]);
  a = 1.;
  b = 1. + .5 * DBL_EPSILON;
  printf("a = %15.13g, b = %15.13g, equal? %s\n", a, b, str[dequal(a, b)]);
  a = 1.;
  b = 1. + DBL_EPSILON;
  printf("a = %15.13g, b = %15.13g, equal? %s\n", a, b, str[dequal(a, b)]);
  a = 1.;
  b = 1. + 1.5 * DBL_EPSILON;
  printf("a = %15.13g, b = %15.13g, equal? %s\n", a, b, str[dequal(a, b)]);
  a = 2.e17;
  b = (1. + .5 * DBL_EPSILON) * 2.e17;
  printf("a = %15.13g, b = %15.13g, equal? %s\n", a, b, str[dequal(a, b)]);
  a = 2.e17;
  b = (1. + 1.5 * DBL_EPSILON) * 2.e17;
  printf("a = %15.13g, b = %15.13g, equal? %s\n", a, b, str[dequal(a, b)]);

  printf("Testing float ...\n");
  c = 0.;
  d = 0.;
  printf("c = %7.5g, d = %7.5g, equal? %s\n", c, d, str[fequal(c, d)]);
  c = 0.;
  d = FLT_EPSILON;;
  printf("c = %7.5g, d = %7.5g, equal? %s\n", c, d, str[fequal(c, d)]);
  c = 1.;
  d = 1. + .5 * FLT_EPSILON;
  printf("c = %7.5g, d = %7.5g, equal? %s\n", c, d, str[fequal(c, d)]);
  c = 1.;
  d = 1. + FLT_EPSILON;
  printf("c = %7.5g, d = %7.5g, equal? %s\n", c, d, str[fequal(c, d)]);
  c = 1.;
  d = 1. + 1.5 * FLT_EPSILON;
  printf("c = %7.5g, d = %7.5g, equal? %s\n", c, d, str[fequal(c, d)]);
  c = 2.e17;
  d = (1. + .5 * FLT_EPSILON) * 2.e17;
  printf("c = %7.5g, d = %7.5g, equal? %s\n", c, d, str[fequal(c, d)]);
  c = 2.e17;
  d = (1. + 1.5 * FLT_EPSILON) * 2.e17;
  printf("c = %7.5g, d = %7.5g, equal? %s\n", c, d, str[fequal(c, d)]);

  return EXIT_SUCCESS;
}
#endif
