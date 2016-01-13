/*
 * SUMMARY:      fileio.h - header file for file I/O
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD:     23-Aug-1996 at 17:20:49 by DHSVM Project Account
 * DESCRIPTION:  header file for file I/O
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 */

/* 	$Id: fileio.h,v 1.3 1996/08/29 02:22:14 dhsvm Exp $	 */

#ifndef FILEIO_H
#define FILEIO_H

/* define identifiers for different file formats */

#define BIN 1                   /* binary IO */
#define HDF 2                   /* hierarchical  data format */

void InitFileIO(int FileFormat);

/* function pointers for 2D file IO */

extern void (*CreateFile)(char *FileName, char *FileLabel);
extern void (*MakeFileName)(char *Path, char *Str1, char *Str2, char *FileName);
extern int (*Read2DMatrix)(int NY, int NX, int NumberType, int NDataSet, 
                           void *Matrix, char *FileName);
extern int (*Read2DSlab)(int NY, int NX, int NumberType, int NDataSet, 
                         void *Matrix, char *FileName);
extern int (*Write2DImage)(int NY, int NX, char *DataLabel, void *Matrix,  
                           char *FileName);
extern int (*Write2DMatrix)(int NY, int NX, int NumberType, char *DataLabel, 
                            char *Units, void *Matrix, char *FileName);

/* generic file functions */

void OpenFile(FILE **FilePtr, char *FileName, char *Mode, 
	      unsigned char OverWrite);
void MakeFileNameGen(char *Path, char *Str1, char *Str2, char *Str3,
                     char *FileName);


#endif

