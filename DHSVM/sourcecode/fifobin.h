/*
 * SUMMARY:      fifobin.h - header file for binary IO functions
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for binary IO functions
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: fifobin.h,v 1.4 2003/07/01 21:26:30 olivier Exp $     
 */

#ifndef FIFOBIN_H
#define FIFOBIN_H

void CreateMapFileBin(char *FileName, ...);
int Read2DMatrixBin(char *FileName, void *Matrix, int NumberType, int NY,
		    int NX, int NDataSet, ...); 
int Read2DMatrixByteSwapBin(char *FileName, void *Matrix, int NumberType,
			    int NY, int NX, int NDataSet, ...);
int Write2DMatrixBin(char *FileName, void *Matrix, int NumberType, int NY,
		     int NX, ...); 
int Write2DMatrixByteSwapBin(char *FileName, void *Matrix, int NumberType,
			     int NY, int NX, ...); 
void byte_swap_long(long *buffer, int number_of_swaps);
void byte_swap_short(short *buffer, int number_of_swaps);

#endif
