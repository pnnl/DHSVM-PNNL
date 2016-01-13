/*
 * SUMMARY:      InitArray.c - Initialize array
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize array
 * DESCRIP-END.
 * FUNCTIONS:    InitCharArray()
 * COMMENTS:
 * $Id: InitArray.c,v 1.4 2003/07/01 21:26:15 olivier Exp $     
 */

#include <stdlib.h>
#include <stdio.h>
#include "functions.h"

void InitCharArray(char *Array, int Size)
{
  int i;			/* counter */

  for (i = 0; i < Size; i++)
    Array[i] = '\0';
}
