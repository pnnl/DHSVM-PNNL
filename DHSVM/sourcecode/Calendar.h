/* -------------------------------------------------------------
   file: Calendar.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created September 12, 1996 by  William A Perkins
   $Id: Calendar.h,v 1.7 2004/07/08 19:53:03 colleen Exp $
   ------------------------------------------------------------- */

#ifndef _Calendar_h_
#define _Calendar_h_

#include <stdio.h>
#include "settings.h"

#define SECPMIN 60
#define SECPHOUR 3600
#define SECPDAY 86400
#define MINPHOUR 60
#define MINPDAY 1440
#define HOURPDAY 24
#define DAYPWEEK 7
#define DAYPYEAR 365
#define MONTHPYEAR 12

/* -------------------------------------------------------------
   module type definitions
   ------------------------------------------------------------- */
typedef struct {
  int Year;
  int Month;
  int Day;
  int Hour;
  int Min;
  int Sec;
  int JDay;			/* Day of year(jan 1 = 1) */
  double Julian;	/* Julian day */
} DATE;

typedef struct {
  int Dt;			/* Timestep (in sec) */
  DATE Start;		/* Starting date of run */
  DATE End;			/* Ending date of run */
  DATE *StartSed;	/* Starting date of sediment run */
  DATE *EndSed;		/* Ending date of sediment run */
  DATE *MWM;        /* Dates to run the mass wasting model */
  DATE MWMnext;     /* Date for next mass wasting model run*/
  DATE Current;		/* Current date in run */
  DATE StartRadar;	/* Start radar file */
  DATE StartMM5;	/* Start of MM5 files */
  int Step;			/* Timestep since start */
  int DayStep;		/* Time step since midnight */
  int NDaySteps;	/* Number of timesteps per day */
  int NTotalSteps;	/* Total number of steps in run */
  int NMWMTotalSteps;  /* Total number of times the mass wasting nmodel is run */
  int NSETotalSteps;   /* Total number of times the surface erosio nmodel is run */
} TIMESTRUCT;

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
uchar After(DATE * Day1, DATE * Day2);
uchar Before(DATE * Day1, DATE * Day2);
void CopyDate(DATE * Copy, DATE * Original);
int DayOfYear(int Year, int Month, int Day);
void IncreaseTime(TIMESTRUCT * Time);
void IncreaseVariableTime(TIMESTRUCT *Time, float VariableDT, TIMESTRUCT *NextTime);
int InitTime(TIMESTRUCT * Time, DATE * Start, DATE * End, DATE * StartRadar,
	     DATE * StartMM5, int Dt);
uchar IsEqualTime(DATE * Day1, DATE * Day2);
uchar IsLeapYear(int Year);
uchar IsNewDay(int DayStep);
uchar IsNewMonth(DATE * Now, int Interval);
DATE NextDate(DATE * Current, int Interval);
int NumberOfSteps(DATE * Start, DATE * End, int Interval);
void PrintDate(DATE * Day, FILE * OutFile);
void PrintRBMStartDate(int Dt, DATE * Day, FILE * OutFile);
void SPrintDate(DATE * Day, char *buffer);
int ScanDate(FILE * InFile, DATE * Day);
int SScanDate(char *Str, DATE * Day);
int SScanMonthDay(char *Str, DATE * Day);
void JulianDayToGregorian(double jd, int *y, int *m, int *d, int *h, int *mi,
			  double *sec);
double GregorianToJulianDay(int year, int mon, int day, int h, int mi,
			    double se);
int DayOfWeek(double j);

#endif
