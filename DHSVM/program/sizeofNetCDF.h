/*
 * SUMMARY:      sizeofnt.h - header for number types defined in DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header for number types defined in DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: sizeofnt.h,v 3.1.1 2012/02/04 Ning Exp $     
 */

#ifndef SIZEOFNETCDF_H
#define SIZEOFNETCDF_H

/* type info codes (taken from netcdf.h, although they are enumerated in that
   case) */
#ifndef HAVE_NETCDF
#define NC_BYTE 1		/* signed 1 byte integer */
#define NC_CHAR	2		/* ISO/ASCII character */
#define NC_SHORT 3		/* signed 2 byte integer */
#define NC_INT 4		/* signed 4 byte integer */
/* #define NC_LONG  *//* 8 bit integer not yet implemented in NetCDF 3.4, but
   anticipated in future versions */
#define	NC_FLOAT 5		/* single precision floating point number */
#define NC_DOUBLE 6		/* double precision floating point number */
#else
#include <netcdf.h>
#endif

size_t SizeOfNumberType(int NumberType);

#endif
