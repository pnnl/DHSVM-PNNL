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
 * LAST CHANGE: 2018-05-07 12:46:50 d3g096
 * COMMENTS:
 */

#ifdef TEST_MAIN
#undef GPTL_TIMING
#endif

#include <stdlib.h>
#include <timing.h>


/******************************************************************************/
/*                            ALLOC_2D_TYPE                                   */
/******************************************************************************/
/** 
 * This should work the same regardless of type, hence the macro. 
 * 
 * @param nx size of first index
 * @param ny size of second index
 * @param thetype element type
 * @param theresult @c "thetype **" pointing to the result
 */

#define ALLOC_2D_TYPE(n1, n2, thetype, theresult)                       \
  {                                                                     \
    int i;                                                              \
    thetype *p;                                                         \
    TIMING_TASK_START("memory allocation");                             \
    theresult = NULL;                                                   \
    p = (thetype *) calloc(n1*n2, sizeof(thetype));                     \
    if (p != NULL)                                                      \
      {                                                                 \
        theresult = (thetype **)calloc(n1, sizeof(thetype *));          \
        if (theresult != NULL)                                          \
          {                                                             \
            for (i = 0; i < n1; ++i)                                    \
              {                                                         \
                size_t idx = i*n2;                                      \
                theresult[i] = &p[idx];                                 \
              }                                                         \
          }                                                             \
        else                                                            \
          {                                                             \
            free(p);                                                    \
          }                                                             \
      }                                                                 \
    TIMING_TASK_END("memory allocation");                               \
  }

#define FREE_2D_TYPE(thearray)                  \
  TIMING_TASK_START("memory allocation");       \
  free(thearray[0]);                            \
  free(thearray);                               \
  TIMING_TASK_END("memory allocation");         \
  thearray = NULL;


/******************************************************************************/
/*                              calloc_2D_float                               */
/******************************************************************************/
float ** 
calloc_2D_float(int NY, int NX)
{
  float **result;
  ALLOC_2D_TYPE(NY, NX, float, result);
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

/******************************************************************************/
/*                              calloc_2D_uint                               */
/******************************************************************************/
unsigned int ** 
calloc_2D_uint(int NY, int NX)
{
  unsigned int **result;
  ALLOC_2D_TYPE(NY, NX, unsigned int, result);
  return result;
}

/******************************************************************************/
/*                             free_2D_uint                                  */
/******************************************************************************/
void
free_2D_uint(unsigned int **p)
{
  FREE_2D_TYPE(p);
}


/******************************************************************************/
/*                              calloc_3D_uint                               */
/******************************************************************************/
unsigned int *** 
calloc_3D_uint(int N1, int N2, int N3)
{
  unsigned int ***result, *p;
  int i, j;
  result = NULL;
  p = (unsigned int *) calloc(N1*N2*N3, sizeof(unsigned int));
  if (p != NULL) {
    ALLOC_2D_TYPE(N1, N2, unsigned int *, result);
    if (result != NULL) {
      for (i = 0; i < N1; ++i) {
        for (j = 0; j < N2; ++j) {
          int idx = i*N2*N3 + j*N3;
          result[i][j] = &p[idx];
        }
      }
    } else {
      free(p);
    }
  }
  return result;
}

/******************************************************************************/
/*                             free_2D_uint                                  */
/******************************************************************************/
void
free_3D_uint(unsigned int ***p)
{
  free(p[0][0]);
  FREE_2D_TYPE(p);
}

/******************************************************************************/
/*                              calloc_3D_uchar                               */
/******************************************************************************/
unsigned char *** 
calloc_3D_uchar(int N1, int N2, int N3)
{
  unsigned char ***result, *p;
  int i, j;
  result = NULL;
  p = (unsigned char *) calloc(N1*N2*N3, sizeof(unsigned char));
  if (p != NULL) {
    ALLOC_2D_TYPE(N1, N2, unsigned char *, result);
    if (result != NULL) {
      for (i = 0; i < N1; ++i) {
        for (j = 0; j < N2; ++j) {
          int idx = i*N2*N3 + j*N3;
          result[i][j] = &p[idx];
        }
      }
    } else {
      free(p);
    }
  }
  return result;
}

/******************************************************************************/
/*                             free_2D_uchar                                  */
/******************************************************************************/
void
free_3D_uchar(unsigned char ***p)
{
  free(p[0][0]);
  FREE_2D_TYPE(p);
}

#ifdef TEST_MAIN

/******************************************************************************/
/*                           main program                                     */
/******************************************************************************/

int
main(int argc, char **argv)
{
  unsigned int **uint2d;
  unsigned int ***uint3d;
  const nx = 10, ny = 5, nz = 3;
  int i, j, k, n;

  uint2d = calloc_2D_uint(nx, ny);
  for (i = 0, n = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j, ++n) {
      uint2d[i][j] = n;
    }
  }
  for (i = 0; i < nx; ++i) {
    printf("%04d: ", i);
    for (j = 0; j < ny; ++j) {
      if (j > 0) printf(", ");
      printf("%4d", uint2d[i][j]);
    }
    printf("\n");
  }

  free_2D_uint(uint2d);

  uint3d = calloc_3D_uint(nx, ny, nz);
  for (i = 0, n = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j) {
      for (k = 0; k < nz; ++k, ++n) {
        uint3d[i][j][k] = n;
      }
    }
  }
  for (k = 0; k < nz; ++k) {
    printf("z = %d:\n", k);
    for (i = 0; i < nx; ++i) {
      printf("%04d: ", i);
      for (j = 0; j < ny; ++j) {
        if (j > 0) printf(", ");
        printf("%4d", uint3d[i][j][k]);
      }
      printf("\n");
    }
  }
  free_3D_uint(uint3d);
}


#endif
