/*
 * SUMMARY:      sizeofnt.h - header for number types defined in DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD:     23-Aug-1996 at 17:26:48 by DHSVM Project Account
 * DESCRIPTION:  header for number types defined in DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 */

/* 	$Id: sizeofnt.h,v 1.2 1996/08/29 02:22:14 dhsvm Exp $	 */

#ifndef SIZEOFNT_H
#define SIZEOFNT_H

/* type info codes (taken from HDF.  DO NOT CHANGE, otherwise 
   compatability with HDF will be lost) */

#define NT_NONE        0    /* indicates that number type not set */

#define NT_FLOAT32     5
#define NT_FLOAT64     6

#define NT_INT8       20
#define NT_UINT8      21

#define NT_INT16      22
#define NT_UINT16     23
#define NT_INT32      24
#define NT_UINT32     25
#define NT_INT64      26
#define NT_UINT64     27

#define NT_UCHAR8      3
#define NT_CHAR8       4
#define NT_CHAR16     42
#define NT_UCHAR16    43

size_t SizeOfNumberType(int NumberType);

#endif
