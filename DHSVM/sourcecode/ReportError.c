/*
 * SUMMARY:      ReportError.c - Report error and exit
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Display a context-dependent error message and exit
 * DESCRIP-END.
 * FUNCTIONS:    ReportError()
 *               ReportWarning()
 * COMMENTS:
 * $Id: ReportError.c,v 1.6 2004/08/24 23:21:48 tbohn Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"

static char *ErrorMessage[] = {
  "Cannot allocate enough memory in function:",	/* 1 */
  "Error while reading file:",	/* 2 */
  "Cannot open file:",		/* 3 */
  "File already exists, and should not be overwritten:",	/* 4 */
  "Error reading input in file:",	/* 5 */
  "No met locations specified in file:",	/* 6 */
  "Hourly timestep does not divide 24 in function:",	/* 7 */
  "No types defined in file: ",	/* 8 */
  "Invalid starting time and timestep combination.\nFirst timestep of the day will not coincide with midnight: ",	/* 9 */
  "Met Station Outside Basin Bounding Box:",	/* 10 */
  "Soil Porosity, Field Capacity, Wilting Point mismatch: ",	/* 11 */
  "Canopy Radiation Transmission less than 5%: ",	/* 12 */
  "Below canopy windspeed less than 5% of above canopy value: ",	/* 13 */
  "Unexpected error returned by:",	/* 14 */
  "Error in switch():",		/* 15 */
  "Absolute value of number larger than MAX_INT: ",	/* 16 */
  "Pixel to be dumped does not lie in the model area:",	/* 17 */
  "Pixel to be dumped does not lie within the modeled basin:",	/* 18 */
  "Not a valid ID for variable to be dumped:",	/* 19 */
  "Wrong layer number to be dumped:",	/* 20 */
  "Invalid resolution specifier in:",	/* 21 */
  "Invalid number of single maps to dump specified in:",	/* 22 */
  "Invalid date specified in:",	/* 23 */
  "Invalid time interval for map dumps in:",	/* 24 */
  "Invalid period for map dumps in:",	/* 25 */
  "Invalid variable ID in:",	/* 26 */
  "Invalid heat flux option specifier: ",	/* 27 */
  "Date in met file does not correspond to current model time:",	/* 28 */
  "Too many vegetation layers: ",	/* 29 */
  "Incorrect flag (valid flags: T, F): ",	/* 30 */
  "Radar or MM5 does not cover entire model area:",	/* 31 */
  "Unknown soil type in file:",	/* 32 */
  "Maximum number of iterations exceeded in RootBrent():",	/* 33 */
  "Root not bracketed in RootBrent():",	/* 34 */
  "Soil moisture profile is supersaturated: ",	/* 35 */
  "Grid NOT square, resulting in problems with flow width calculation:",	/* 36 */
  "Radar precipitation file starts later than start of run:",	/* 37 */
  "Invalid file format specifier:",	/* 38 */
  "Error setting file pointer:",	/* 39 */
  "Unknown number type:",	/* 40 */
  "Error writing to file:",	/* 41 */
  "Bad met interpolation combination",	/* 42 */
  "Invalid radiation type specifier:",	/* 43 */
  "Invalid channel/road network specifier:",	/* 44 */
  "Invalid flow gradient specifier:",	/* 45 */
  "Incorrect travel times in hydrograph file:",	/* 46 */
  "Table index out of bounds: ",	/* 47 */
  "Wind profile error:",	/* 48 */
  "Current version does not support more than two vegetation layers:",	/* 49 */
  "Section not found:",		/* 50 */
  "Key not found or invalid entry:",	/* 51 */
  "Invalid key or required key (no default value):",	/* 52 */
  "Invalid station number:",	/* 53 */
  "Cannot use precipitation lapse rate map if there is more than one met station:",	/* 54 */
  "Invalid Channel ID:",	/* 55 */
  "HAVE_NETCDF undefined during build, cannot read NetCDF format:",	/* 56 */
  "NETCDF error:",		/* 57 */
  "NETCDF Warning:",		/* 58 */
  "Wrong Y dimension:",		/* 59 */
  "Wrong X dimension:",		/* 60 */
  "Can't byteswap that elements size:",	/*61 */
  "No GLACIER land use class defined:", /*62*/
  "Ran out of data in surface routing file:", /* 63 */
  "Row column mismatch in surface routing file:", /* 64 */
  "Current version does not support this setup:", /* 65 */
  "Invalid Map->Resolution value for dumping map or image of variable ID:", /* 66 */
  "The options set in the input file do not support plotting variable ID:", /* 67 */
  "Buildup coefficients are out of range:", /* 68 */
  "Washoff coefficients are out of range\n Washoff Coeff must be > 0 & Washoff Expon must be within [-10, 10]:", /* 69 */
  "Invalid buildup function type:", /* 70 */
  "Invalid washoff function type:", /* 71 */
  "Water quality mass balance error greater than 10%: ",	/* 72 */
  "Runoff mass balance error greater than 10%: ",	/* 73 */
  "Number of pollutants for each land use category must be equal to the number of pollutants specified in [POLLUTANTS] section:", /* 74 */
  NULL
};

void ReportError(char *ErrorString, int ErrorCode)
{
  printf("%s %s\n", ErrorMessage[ErrorCode - 1], ErrorString);

  exit(ErrorCode);
}

void ReportWarning(char *ErrorString, int ErrorCode)
{
  fprintf(stderr, "%s %s\n", ErrorMessage[ErrorCode - 1], ErrorString);
}

/*******************************************************************************
  Test main. Compile by typing:
  gcc -DTEST_REPORTERROR -o test_error ReportError.c
  then run the program by typing test_error
*******************************************************************************/

#ifdef TEST_REPORTERROR

int main(void)
{
  char str[BUFSIZ + 1];
  int ErrorCode;

  ErrorCode = 1;
  while (ErrorMessage[ErrorCode - 1] != NULL) {
    sprintf(str, "Test -- Code %d", ErrorCode);
    ReportWarning(str, ErrorCode);
    ++ErrorCode;
  }
  fprintf(stderr, "\n");
  ReportWarning("", 3);
  fprintf(stderr,
	  "\nThe test is SUCCESSFUL if the next line is the same as the line "
	  "above:\n\n");
  ReportError("", 3);

  return EXIT_SUCCESS;
}
#endif
