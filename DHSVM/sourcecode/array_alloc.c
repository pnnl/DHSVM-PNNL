/*
 * SUMMARY:      array_alloc.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    May 2017
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2017-05-01 11:48:05 d3g096
 * COMMENTS:
 */

#include <stdlib.h>


/******************************************************************************/
/*                            ALLOC_2D_TYPE                                   */
/******************************************************************************/
/** 
 * This should work the same regardless of type, hence the macro. 
 * 
 * @param theroutine name of calling routine
 * @param nx size of first index
 * @param ny size of second index
 * @param thetype element type
 * @param theresult @c "thetype **" pointing to the result
 */

#define ALLOC_2D_TYPE(theroutine, nx, ny, thetype, theresult)   \
  { int j;  \
    theresult = (thetype **)calloc(NY, sizeof(thetype *));        \
    if (theresult == NULL) ReportError(theroutine, 1);                  \
    theresult[0] = (thetype *) calloc(NX*NY, sizeof(thetype));          \
    if (theresult[0] == NULL) ReportError(theroutine, 1);               \
    for (j = 1; j < NY; ++j) theresult[j] = theresult[0] + j*NX;    \
  }

#define FREE_2D_TYPE(thearray)                  \
  free(thearray[0]);                            \
  free(thearray);                               \
  thearray = NULL;


/******************************************************************************/
/*                              calloc_2D_float                               */
/******************************************************************************/
float ** 
calloc_2D_float(int NY, int NX)
{
  static char Routine[] = "calloc_2D_float";
  float **result;
  ALLOC_2D_TYPE(Routine, NY, NX, float, result);
  return result;
}

/******************************************************************************/
/*                             free_2D_float                                  */
/******************************************************************************/
void
free_2D_float(float **p)
{
  FREE_2D_TYPE(p);
}
