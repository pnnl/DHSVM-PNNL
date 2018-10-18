/*
 * SUMMARY:      byte_swap.h
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    October 2018
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-10-17 14:56:37 d3g096
 * COMMENTS:
 */
 
#ifndef _byte_swap_h_
#define _byte_swap_h_

extern void byte_swap_short(short *buffer, int number_of_swaps);
extern void byte_swap_long(long *buffer, int number_of_swaps);


#endif

