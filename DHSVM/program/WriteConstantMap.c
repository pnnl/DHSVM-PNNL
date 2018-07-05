/*
 * SUMMARY:      WriteConstantMap.c
 * USAGE:        Part of DHSVM

 *               This program is used to make a map file that contains
 *               a single value for the entire map. Multiple maps,
 *               with the same value, can be included in the output. 
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    July 2018
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-07-03 12:17:35 d3g096
 * COMMENTS:
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libgen.h>

#include "sizeofnt.h"
#include "fileio.h"
#include "fifobin.h"

static char program[1024];

/* -------------------------------------------------------------
   Main Program
   ------------------------------------------------------------- */
int
main(int argc, char **argv)
{
  const char usage[] = 
    "usage: %s rows cols value outfile [steps]\n";
  FILE *outFile;                /* Output file */
  
  char outFilename[1024];   /* Name of output file */
  int nRows;                    /* Number of rows */
  int nCols;                    /* Number of columns */
  int value;                    /* Value to put in output map */
  int nSteps;
  int i, j, ierr;
  const int NumberType = NT_FLOAT32;
  const char DataLabel[] = "Just a constant value";
  const char Units[] = "none";

  float *Array;

  strncpy(program, basename(argv[0]), 1024);

  if (argc < 5) {
    fprintf(stderr, usage, program);
    return(3);
  }

  nRows = atoi(argv[1]);
  nCols = atoi(argv[2]);

  if (nRows <= 0 || nCols <= 0) {
    fprintf(stderr, "%s: error: invalid rows and/or columns (%s x %s)\n",
            argv[1], argv[2]);
    return(3);
  }
  value = atof(argv[3]);
  strncpy(outFilename, argv[4], 1024);

  if (argc > 5) {
    nSteps = atoi(argv[5]);
    if (nSteps <= 0) {
      fprintf(stderr, "%s: error: invalid number of steps (%s)\n",
              program, argv[5]);
      return(3);
    }
  } else {
    nSteps = 1;
  }

  CreateFileBin(outFilename, NULL);

  if (!(Array = (float*) calloc(nRows*nCols, sizeof(float)))) {
    fprintf(stderr, "%s: error: cannot allocate memory\n", program);
    return(2);
  }

  for (j = 0; j < nCols; ++j) {
    for (i = 0; i < nRows; ++i) {
      Array[j*nRows + i] = value;
    }
  }

  for (i = 0; i < nSteps; ++i) {
    Write2DMatrixBin(nCols, nRows, NumberType, DataLabel, Units, Array, outFilename);
  }

  return (0);
}
