/* -------------------------------------------------------------
   file: Calendar.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created September 12, 1996 by  William A Perkins
   Last Change: Fri Oct 11 14:51:38 1996 by  William A Perkins <perk@doofus.pnl.gov>
   ------------------------------------------------------------- */

/* $Id: Calendar.h,v 1.3 1996/11/14 19:12:26 battelle Exp $ */

#ifndef _Calendar_h_
#define _Calendar_h_

#include <stdio.h>
#include "typenames.h"

/* -------------------------------------------------------------
   module type definitions
   ------------------------------------------------------------- */
typedef struct {
  int Year;
  int Month;
  int Day;
  int JDay;			/* Day of year(jan 1 = 1) */
  int Hour;
} DATE;

typedef struct {
  int Dt;			   /* Timestep (in sec) */
  DATE Start;		   /* Starting date of run */
  DATE End;			   /* Ending date of run */
  DATE Current;		   /* Current date in run */
  DATE StartRadar;	   /* Start radar file */
  DATE StartMM5;	   /* Start of MM5 files */
  int Step;			   /* Timestep since start */
  int DayStep;		   /* Time step since midnight */
  int NDaySteps;	   /* Number of timesteps per day */
  int NTotalSteps;	   /* Total number of steps in run */
} TIMESTRUCT;

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
uchar After(DATE Day1, DATE Day2);
uchar Before(DATE Day1, DATE Day2);
void CopyDate(DATE *Copy, DATE Original);
int DayOfYear(int Year, int Month, int Day);
uchar IsEqualTime(DATE Day1, DATE Day2);
uchar IsLeapYear(int Year);
uchar IsNewDay(int DayStep);
uchar IsNewMonth(DATE Day);
DATE NextDate(DATE Current, int Interval);
int NumberOfSteps(DATE Start, DATE End, int Interval);
void PrintDate(DATE Day, FILE *OutFile);
void SPrintDate(DATE Day, char *buffer);
int ScanDate(FILE *InFile, DATE *Day);
int SScanDate(char *Str, DATE *Day);

#endif
