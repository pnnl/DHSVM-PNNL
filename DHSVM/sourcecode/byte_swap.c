/*
 * SUMMARY:      byte_swap.c
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
 * LAST CHANGE: 2018-10-17 14:57:57 d3g096
 * COMMENTS:
 */

/******************************************************************************/
/*                          byte_swap_short                                   */
/******************************************************************************/
void 
byte_swap_short(short *buffer, int number_of_swaps)
{
  short *temp;
  int swap_loop;

  for (swap_loop = 0, temp = buffer; swap_loop < number_of_swaps;
       swap_loop++, temp++) {
    *temp = ((*temp & 0x00ff) << 8) | ((*temp & 0xff00) >> 8);
  }
}

/******************************************************************************/
/*                             byte_swap_long                                 */
/******************************************************************************/
void 
byte_swap_long(long *buffer, int number_of_swaps)
{
  long *temp;
  int swap_loop;

  for (swap_loop = 0, temp = buffer; swap_loop < number_of_swaps;
       swap_loop++, temp++) {
    *temp = ((*temp & 0x000000ff) << 24) | ((*temp & 0x0000ff00) << 8) |
      ((*temp & 0x00ff0000) >> 8) | ((*temp & 0xff000000) >> 24);
  }
}

