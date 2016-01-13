/*
 * SUMMARY:      brent.h - header file for the brent method
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for the brent method
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: brent.h,v 1.4 2003/07/01 21:26:26 olivier Exp $     
 */

#ifndef BRENT_H
#define BRENT_H

float RootBrent(int y, int x, float LowerBound, float UpperBound,
		float (*Function) (float Estimate, va_list ap), ...);

#define MACHEPS      3e-8	/* machine floating point precision (float) */

#define T            1e-5	/* tolerance */

#define MAXITER      1000	/* maximum number of allowed iterations */

#define MAXTRIES     5		/* maximum number of tries to bracket the 
				   root */

#define TSTEP        10		/* step to take in both directions if
				   attempting to bracket thr root  */
#endif
