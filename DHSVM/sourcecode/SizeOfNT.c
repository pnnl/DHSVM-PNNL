/*
 * SUMMARY:      SizeOfNT() - Determine size of number type in bytes
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Determine size of number type in bytes
 * DESCRIP-END.
 * FUNCTIONS:    SizeOfNumberType()
 * COMMENTS:
 * $Id: SizeOfNT.c,v 1.4 2003/07/01 21:26:24 olivier Exp $     
 */

#include <stdlib.h>
#include <sys/types.h>
#include "DHSVMerror.h"
#include "sizeofnt.h"

/*****************************************************************************
  Function name: SizeOfNumberType()

  Purpose      : Determines the size of each number type.  Number type 
                 descriptors are "borrowed" from NetCDF.

  Required     :
    NumberType - descriptor of number type (for more info see comments in 
                 InitFileIO.c)

  Returns      : size of number type

  Modifies     : NA

  Comments     :
   Information about the number type has to be passed to the functions that
   deal with 2DMatrix.  The number types are adopted from the NetCDF standard.
*****************************************************************************/
size_t SizeOfNumberType(int NumberType)
{
  const char *Routine = "SizeOfNumberType";
  size_t ElemSize;

  switch (NumberType) {
  case NC_BYTE:
  case NC_CHAR:
    ElemSize = sizeof(char);
    break;
  case NC_SHORT:
    ElemSize = sizeof(short);
    break;
  case NC_INT:
/*   case NC_LONG: *//* NC_LONG is the same as NC_INT in NetCDF 3.4 */
    ElemSize = sizeof(int);
    break;
  case NC_FLOAT:
    ElemSize = sizeof(float);
    break;
  case NC_DOUBLE:
    ElemSize = sizeof(double);
    break;
  default:
    ElemSize = 0;
    ReportError((char *) Routine, 40);
    break;
  }

  return ElemSize;
}
