/*
 * SUMMARY:      getinit.h - header file for GetInit.c
 * USAGE:        Not a stand-alone program
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     6-May-97 at 10:43:13
 $Id: getinit.h,v 1.4 2003/07/01 21:26:30 olivier Exp $
 */

#ifndef GETINIT_H
#define GETINIT_H

#ifndef BUFSIZE
#define BUFSIZE 255
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define SEPARATOR    '='
#define OPENCOMMENT  '#'
#define OPENSECTION  '['
#define CLOSESECTION ']'

typedef struct _INTINIENTRY {
  char *SectionName;
  char *KeyName;
  long *Var;
  long Default;
} INTINIENTRY;

typedef struct _INPUTSTRUCT *LISTPTR;

typedef struct _INPUTSTRUCT {
  char Str[BUFSIZE + 1];
  LISTPTR Next;
} INPUTSTRUCT;

typedef struct _DBLINIENTRY {
  char *SectionName;
  char *KeyName;
  double *Var;
  double Default;
} DBLINIENTRY;

typedef struct _STRINIENTRY {
  char *SectionName;
  char *KeyName;
  char VarStr[BUFSIZE + 1];
  char *Default;
} STRINIENTRY;

int CopyDouble(double *Value, char *Str, const int NValues);
int CopyFloat(float *Value, char *Str, const int NValues);
int CopyInt(int *Value, char *Str, const int NValues);
int CopyLong(long *Value, char *Str, const int NValues);
int CopyShort(short *Value, char *Str, const int NValues);
int CopyUChar(unsigned char *Value, char *Str, const int NValues);

int CountLines(FILE * InFile);

LISTPTR CreateNode(void);

void DeleteList(LISTPTR StartNode);

double GetInitDouble(const char *Section, const char *Key, double Default,
		     LISTPTR Input);

long GetInitLong(const char *Section, const char *Key, long Default,
		 LISTPTR Input);

unsigned long GetInitString(const char *Section, const char *Key,
			    const char *Default, char *ReturnBuffer,
			    unsigned long BufferSize, LISTPTR Input);

int IsEmptyStr(char *Str);

unsigned char IsKeyEntryPair(char *Buffer);

unsigned char IsSection(char *Buffer);

unsigned char LocateKey(const char *Key, char *Entry, LISTPTR Input);

LISTPTR LocateSection(const char *Section, LISTPTR Input);

void MakeKeyString(char *Buffer);

void ReadInitFile(char *TemplateFileName, LISTPTR * Input);

void Strip(char *Buffer);

#endif
