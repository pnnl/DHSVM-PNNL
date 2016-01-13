/*
 * SUMMARY:      fifoNetCDF.h - header file for netcdf IO functions
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * Created:      Wed Jan 27 10:03:51 1999
 * DESCRIPTION:  header file for netcdf IO functions
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: fifoNetCDF.h,v3.1.2 2013/07/01 ning Exp $     
 */

#ifndef FIFONETCDF_H
#define FIFONETCDF_H

#define ATT_HISTORY   "history"
#define ATT_COMMENT   "comment"
#define ATT_MISSINGVALUE "missing_value"
#define ATT_LONGNAME  "long_name"
#define ATT_NAME      "name"
#define ATT_UNITS     "units"
#define ATT_FORMAT    "C_format"
#define TIME_DIM      "time"
#define X_DIM         "x"
#define Y_DIM         "y"

void CreateMapFileNetCDF(char *FileName, ...);
int Read2DMatrixNetCDF(char *FileName, void *Matrix, int NumberType, int NY,
		       int NX, int NDataSet, ...);
int Write2DMatrixNetCDF(char *FileName, void *Matrix, int NumberType, int NY,
			int NX, ...);

#endif
