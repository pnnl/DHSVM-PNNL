/*
 * SUMMARY:      SizeOfNT() - Determine size of number type in bytes
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD:      6-Sep-1996 at 20:59:48 by Bart Nijssen
 * DESCRIPTION:  Determine size of number type in bytes
 * DESCRIP-END.
 * FUNCTIONS:    SizeOfNumberType()
 * COMMENTS:     
 */

#ifndef lint
static char vcid[] = "nil";
#endif /* lint */

#include <stdlib.h>
#include <sys/types.h>
#include "DHSVMerror.h"
#include "sizeofnt.h"

/*****************************************************************************
  Function name: SizeOfNumberType()

  Purpose      : Determines the size of each number type.  Number type 
                 descriptors are "borrowed" from HDF.

  Required     :
    NumberType - descriptor of number type (for more info see comments in 
                 InitFileIO.c)

  Returns      : size of number type

  Modifies     : NA

  Comments     :
   Information about the number type has to be passed to the functions that
   deal with 2DMatrix and 2DSlab.  The following number types are recognized
   (they are adopted from the HDF standard):
   - NT_UINT8
     NT_UCHAR8:   8 bit unsigned int (i.e. unsigned char)
   - NT_UINT16:  16 bit unsigned int (i.e. unsigned short on most machines)
   - NT_UINT32:  32 bit unsigned int (i.e. unsigned int on most machines)
   - NT_INT8
     NT_CHAR8:    8 bit int (i.e. char)
   - NT_INT16:   16 bit int (i.e. short on most machines)
   - NT_INT32:   32 bit int (i.e. int on most machines)
   - NT_FLOAT32: 32 bit float (i.e. float on most machines)
   - NT_FLOAT64: 64 bit float (i.e. double on most machines)
*****************************************************************************/
size_t SizeOfNumberType(int NumberType)
{
  const char *Routine = "SizeOfNumberType";
  size_t ElemSize;

  switch (NumberType) {
  case NT_CHAR8:
  case NT_INT8:
    ElemSize = sizeof(char);
    break;
  case NT_UCHAR8:
  case NT_UINT8:
    ElemSize = sizeof(unsigned char);
    break;
  case NT_INT16:
    ElemSize = sizeof(short);
    break;
  case NT_UINT16:
    ElemSize = sizeof(unsigned short);
    break;
  case NT_INT32:
    ElemSize = sizeof(int);
    break;
  case NT_UINT32:
    ElemSize = sizeof(unsigned int);
    break;
  case NT_FLOAT32:
    ElemSize = sizeof(float);
    break;
  case NT_FLOAT64:
    ElemSize = sizeof(double);
    break;
  case NT_NONE:
  default:
    ElemSize = 0;
    ReportError((char *) Routine, 40);
    break;
  }
  
  return ElemSize;
}
