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
 * LAST CHANGE: 2018-11-06 09:48:12 d3g096
 * COMMENTS:
 */
 
#ifndef _byte_swap_h_
#define _byte_swap_h_

#ifdef __cplusplus
extern "C" {
#endif
void byte_swap_short(short *buffer, int number_of_swaps);
void byte_swap_long(long *buffer, int number_of_swaps);
#ifdef __cplusplus
}
#endif


#endif

