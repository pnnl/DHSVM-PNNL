/*
 * SUMMARY:      Calendar.c - Generic functions to manipulate times and dates
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD: Fri Oct 11 14:50:46 1996 by  William A Perkins <perk@doofus.pnl.gov>
 * DESCRIPTION:  Generic functions to manipulate times and dates
 * DESCRIP-END.
 * FUNCTIONS:    DayOfYear() 
 *               IsLeapYear() 
 *               IsEqualTime() 
 *               ScanDate()
 *               NumberOfSteps() 
 *               NextDate() 
 *               CopyDate() 
 *               PrintDate() 
 *               IsNewMonth()
 *               IsNewDay()
 *               Before() 
 *               After()
 * COMMENTS:     
 */

#ifndef lint
static char vcid[] = "$Id: Calendar.c,v 1.6 1996/12/03 22:09:49 nijssen Exp $";
#endif /* lint */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "functions.h"

#define FAIL -1

/*****************************************************************************
  DayOfYear()
*****************************************************************************/
int DayOfYear(int Year, int Month, int Day) 
{
  int i;
  int Jday;
  int DaysPerMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if (IsLeapYear(Year))
    DaysPerMonth[1] = 29;
  else
    DaysPerMonth[1] = 28;

  for (i = 0, Jday = 0; i < (Month - 1); i++)
    Jday += DaysPerMonth[i];

  Jday += Day;

  return Jday;
}

/*****************************************************************************
  IsLeapYear()
*****************************************************************************/
uchar IsLeapYear(int Year) 
{
  if ((Year % 4 == 0 && Year % 100 != 0) || Year % 400 == 0) 
    return TRUE;
  return FALSE;
}

/*****************************************************************************
  IsEqualTime()
*****************************************************************************/
uchar IsEqualTime(DATE Day1, DATE Day2)
{
  if ((Day1.Year == Day2.Year) && (Day1.Month == Day2.Month) && 
      (Day1.Day == Day2.Day) && (Day1.Hour == Day2.Hour))
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  ScanDate()
*****************************************************************************/
int ScanDate(FILE *InFile, DATE *Day)
{
  char Str[NAMESIZE+1];

  if (fscanf(InFile, "%s", Str) != 1)
    return FALSE;

  return SScanDate(Str, Day);
}

/*****************************************************************************
  SScanDate()
*****************************************************************************/
int SScanDate(char *Str, DATE *Day)
{
  int i;
  int j;
  int Length;
  int Number[4];
  int DaysPerMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if (Str == NULL)
    return FALSE;

  Length = strlen(Str);

  for (i = Length-1, j = 0; i > 0; i--) {
    if (!isdigit(Str[i])) {
      Number[j] = atoi(&Str[i+1]);
      Str[i] = '\0';
      j++;
    }
  }

  Number[j] = atoi(Str);

  Day->Month = Number[3];
  Day->Day   = Number[2];
  Day->Year  = Number[1];
  Day->Hour  = Number[0];

  if (IsLeapYear(Day->Year))
    DaysPerMonth[1] = 29;
  else
    DaysPerMonth[1] = 28;

  if (Day->Month < 1 || Day->Month > MONTHSPYR)
    return FALSE;
  if (Day->Day < 1 || Day->Day > DaysPerMonth[Day->Month-1])
    return FALSE;
  if (Day->Hour < 0 || Day->Hour > 23)
    return FALSE;
  
  Day->JDay = DayOfYear(Day->Year, Day->Month, Day->Day); 

  return TRUE;
}

/*****************************************************************************
  NumberOfSteps()
*****************************************************************************/
int NumberOfSteps(DATE Start, DATE End, int Interval)
{
  int DaysPerYear;
  int NDays = 0;		/* Total number of days */
  int NSteps = 0;		/* Number of steps */
  int Day1;			/* Starting date */
  int Year;			/* Current year */

  if (Start.Year >= End.Year) {
    if (Start.Year > End.Year)
      return FAIL;
    else {
      if (Start.JDay >= End.JDay) {
	if (Start.JDay > End.JDay)
	  return FAIL;
	else {
	  if (Start.Hour >= End.Hour) {
	    if (Start.Hour > End.Hour)
	      return FAIL;
	    else
	      return 1;
	  }
	}
      }
    }
  }
  
  Day1 = Start.JDay;

  if (Start.Year < End.Year) {
    if (IsLeapYear(Start.Year))
      DaysPerYear = 366;
    else
      DaysPerYear = 365;
    NDays = DaysPerYear - Day1;
    Day1 = 0;
  }

  for (Year = Start.Year+1; Year < End.Year; Year++) {

    if (IsLeapYear(Year))
      DaysPerYear = 366;
    else
      DaysPerYear = 365;
    
    NDays += DaysPerYear;
  }

  NDays += End.JDay - Day1;
  NSteps = NDays * (int) HOURSPDAY;

  NSteps += End.Hour - Start.Hour;
  
  if (NSteps % Interval != 0)
    return FAIL;

  NSteps /= Interval;
  NSteps += 1;

  return NSteps;
}

/*****************************************************************************
  NextDate()
*****************************************************************************/
DATE NextDate(DATE Current, int Interval)
{
  int DaysPerMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int SumHours;
  DATE Next;                    /* Next date */

  CopyDate(&Next, Current);


  if ((SumHours = Current.Hour + Interval) >= HOURSPDAY) {

    Next.Hour = SumHours % Interval;
    Next.Day = Current.Day + SumHours/(int)HOURSPDAY;

    if (IsLeapYear(Current.Year))
      DaysPerMonth[1] = 29;
    else 
      DaysPerMonth[1] = 28;

    if (Next.Day > DaysPerMonth[Current.Month-1]) {
      Next.Month +=1;
      Next.Day -= DaysPerMonth[Current.Month-1];
    }
    else
      Next.Month = Current.Month;

    if (Next.Month > MONTHSPYR) {
      Next.Year = Current.Year + 1;
      Next.Month = 1;
    }
    else 
      Next.Year = Current.Year;
  }
  else
    Next.Hour = Current.Hour + Interval;

  Next.JDay = DayOfYear(Current.Year, Current.Month, Current.Day); 

  return Next;
}

/*****************************************************************************
  CopyDate()
*****************************************************************************/
void CopyDate(DATE *Copy, DATE Original)
{
  Copy->Year  = Original.Year;
  Copy->Month = Original.Month;
  Copy->Day   = Original.Day;
  Copy->JDay  = Original.JDay;
  Copy->Hour  = Original.Hour;
}

/*****************************************************************************
  PrintDate()
*****************************************************************************/
void PrintDate(DATE Day, FILE *OutFile)
{
  fprintf(OutFile, "%02d/%02d/%4d-%02dhr", Day.Month, Day.Day, 
	  Day.Year, Day.Hour);
} 

/* -------------------------------------------------------------
   SPrintDate
   Formats a DATE to a string
   ------------------------------------------------------------- */
void SPrintDate(DATE Day, char *buffer)
{
  sprintf(buffer, "%02d/%02d/%4d-%02dhr", Day.Month, Day.Day, 
	  Day.Year, Day.Hour);
}


/*****************************************************************************
  IsNewMonth()
*****************************************************************************/
uchar IsNewMonth(DATE Day)
{
  if (Day.Day == 1 && Day.Hour == 0)
    return TRUE;
  else 
    return FALSE;
}

/*****************************************************************************
  IsNewDay()
*****************************************************************************/
uchar IsNewDay(int DayStep)
{
  if (DayStep == 0)
    return TRUE;
  else 
    return FALSE;
}

/*****************************************************************************
  Before()
*****************************************************************************/
uchar Before(DATE Day1, DATE Day2)
{
  int Time1;
  int Time2;

  Time1 = ((Day1.Year * 365) + Day1.JDay) * (int)HOURSPDAY + Day1.Hour;
  Time2 = ((Day2.Year * 365) + Day2.JDay) * (int)HOURSPDAY + Day2.Hour;

  if (Time1 < Time2)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  After()
*****************************************************************************/
uchar After(DATE Day1, DATE Day2)
{
  int Time1;
  int Time2;

  Time1 = ((Day1.Year * 365) + Day1.JDay) * (int)HOURSPDAY + Day1.Hour;
  Time2 = ((Day2.Year * 365) + Day2.JDay) * (int)HOURSPDAY + Day2.Hour;

  if (Time1 > Time2)
    return TRUE;
  else
    return FALSE;
}
