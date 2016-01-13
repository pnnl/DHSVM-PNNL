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

void CreateFileBin(char *FileName, char *FileLabel);
void MakeFileNameBin(char *Path, char *Str1, char *Str2, char *FileName);
int Read2DMatrixBin(int NY, int NX, int NumberType, int NDataSet,
                    void *Matrix, char *FileName);
int Write2DImageBin(int NY, int NX, char *DataLabel, void *Matrix,
                    char *FileName);
int Write2DMatrixBin(int NY, int NX, int NumberType, char *DataLabel,
                     char *Units, void *Matrix, char *FileName);


#endif
