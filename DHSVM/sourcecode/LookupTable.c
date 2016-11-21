/*
 * SUMMARY:      LookupTable.c
 * USAGE:        Part of DHSVM - initializes and allows lookup for tables with
 *               regularly spaced indices
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    18-Feb-97 at 17:46:13
 * DESCRIPTION:  initializes and allows lookup for tables with
 *               regularly spaced indices
 * DESCRIP-END.
 * FUNCTIONS:    init_float_table()
 *               float float_lookup(float x, FLOATTABLE *table)
 * COMMENTS:
 * $Id: LookupTable.c,v 1.4 2003/07/01 21:26:19 olivier Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "lookuptable.h"
#include "DHSVMerror.h"

/*****************************************************************************
  Function name: InitFloatTable()

  Purpose      : Initialize a table structure
                 
  Required     :
    unsigned long Size        - Number of entries in the lookup table
    float Offset              - Value of key for first entry in the table
    float Delta               - Key interval
    float (*Function)(float)  - Function used to fill the table entries
    FLOATTABLE *Table         - pointer to structure that holds the table

  Returns      : void

  Modifies     : FLOATTABLE *table

  Comments     :
*****************************************************************************/
void InitFloatTable(unsigned long Size, float Offset, float Delta,
		    float (*Function) (float), FLOATTABLE * Table)
{
  int i;
  float x;

  Table->Size = Size;
  Table->Offset = Offset;
  Table->Delta = Delta;

  Table->Data = calloc(Table->Size, sizeof(float));
  if (Table->Data == NULL)
    ReportError("InitFloatTable", 1);

  for (i = 0, x = Table->Offset + .5 * Table->Delta; i < Table->Size;
       i++, x += Table->Delta)
    Table->Data[i] = Function(x);
}

/*****************************************************************************
  Function name: FloatLookup()

  Purpose      : Lookup a table entry corresponding to key x
                 
  Required     : 
    float x           - key to be looked up
    FLOATTABLE *Table - Table structure that contains the entries

  Returns      : float

  Modifies     : None

  Comments     :
*****************************************************************************/
float FloatLookup(float x, FLOATTABLE * Table)
{
  int i;

  i = (int) (((x - Table->Offset) / Table->Delta));
  if (i < 0 || i >= Table->Size) {
    sprintf(errorstr, "FloatLookup: attempting lookup of value %f \n", x);
    ReportError(errorstr, 47);
  }

  return Table->Data[i];
}
