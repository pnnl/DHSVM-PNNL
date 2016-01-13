/*
 * SUMMARY:      Calendar.c - Generic functions to manipulate times and dates
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
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
 *               IncreaseTime()
 *               InitTime()
 *               GregorianToJulianDay()
 *               JulianDayToGregorian()
 *               DayOfWeek()
 * COMMENTS:
 * $Id: Calendar.c,v 1.5 2004/02/17 20:40:58 jlanini Exp $     
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "functions.h"
#include "Calendar.h"
#include "DHSVMerror.h"

#define FAIL -1

/*****************************************************************************
  DayOfYear()
*****************************************************************************/
int DayOfYear(int Year, int Month, int Day)
{
  int i;
  int Jday;
  int DaysPerMonth[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

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
uchar IsEqualTime(DATE * Day1, DATE * Day2)
{
  return dequal(Day1->Julian, Day2->Julian);
}

/*****************************************************************************
  ScanDate()
*****************************************************************************/
int ScanDate(FILE * InFile, DATE * Day)
{
  char Str[NAMESIZE + 1];

  if (fscanf(InFile, "%s", Str) != 1)
    return FALSE;

  return SScanDate(Str, Day);
}

/*****************************************************************************
  SScanDate()
*****************************************************************************/
int SScanDate(char *DateStr, DATE * Day)
{
  char Str[BUFSIZE + 1];
  int i;
  int j;
  int Length;
  int Number[6];
  int DaysPerMonth[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

  if (Str == NULL)
    return FALSE;

  strcpy(Str, DateStr);

  Length = strlen(Str);

  for (i = Length - 1, j = 0; i > 0; i--) {
    if (!isdigit((int) Str[i])) {
      Number[j] = atoi(&Str[i + 1]);
      Str[i] = '\0';
      j++;
    }
  }

  Number[j] = atoi(Str);

  if (j < 2 || j > 5)
    return FALSE;

  Day->Month = Number[j--];
  Day->Day = Number[j--];
  Day->Year = Number[j--];
  if (j >= 0)
    Day->Hour = Number[j--];
  else
    Day->Hour = 0;
  if (j >= 0)
    Day->Min = Number[j--];
  else
    Day->Min = 0;
  if (j >= 0)
    Day->Sec = Number[j--];
  else
    Day->Sec = 0;

  if (IsLeapYear(Day->Year))
    DaysPerMonth[1] = 29;
  else
    DaysPerMonth[1] = 28;

  if (Day->Month < 1 || Day->Month > MONTHPYEAR)
    return FALSE;
  if (Day->Day < 1 || Day->Day > DaysPerMonth[Day->Month - 1])
    return FALSE;
  if (Day->Hour < 0 || Day->Hour > 23)
    return FALSE;
  if (Day->Min < 0 || Day->Min > 59)
    return FALSE;
  if (Day->Sec < 0 || Day->Sec > 59)
    return FALSE;

  Day->JDay = DayOfYear(Day->Year, Day->Month, Day->Day);
  Day->Julian =
    GregorianToJulianDay(Day->Year, Day->Month, Day->Day, Day->Hour, Day->Min,
			 Day->Sec);

  return TRUE;
}

/*****************************************************************************
  NumberOfSteps()
*****************************************************************************/
int NumberOfSteps(DATE * Start, DATE * End, int Interval)
{
  int NSteps = 0;		/* Number of steps */

  NSteps = (End->Julian - Start->Julian) * (SECPDAY / Interval);
  NSteps++;
  if (NSteps <= 0)
    return FAIL;

  return NSteps;
}

/*****************************************************************************
  NextDate()
*****************************************************************************/
DATE NextDate(DATE * Current, int Interval)
{
  DATE Next;			/* Next date */
  double Sec;

  Next.Julian = Current->Julian + ((double) Interval) / SECPDAY;
  JulianDayToGregorian(Next.Julian, &(Next.Year), &(Next.Month), &(Next.Day),
		       &(Next.Hour), &(Next.Min), &Sec);
  Next.Sec = (int) Sec;
  Next.JDay = DayOfYear(Next.Year, Next.Month, Next.Day);

  return Next;
}

/*****************************************************************************
  CopyDate()
*****************************************************************************/
void CopyDate(DATE * Copy, DATE * Original)
{
  Copy->Year = Original->Year;
  Copy->Month = Original->Month;
  Copy->Day = Original->Day;
  Copy->Hour = Original->Hour;
  Copy->Min = Original->Min;
  Copy->Sec = Original->Sec;
  Copy->JDay = Original->JDay;
  Copy->Julian = Original->Julian;
}

/*****************************************************************************
  PrintDate()
*****************************************************************************/
void PrintDate(DATE * Day, FILE * OutFile)
{
  fprintf(OutFile, "%02d/%02d/%4d-%02d:%02d:%02d", Day->Month, Day->Day,
	  Day->Year, Day->Hour, Day->Min, Day->Sec);
}

/* -------------------------------------------------------------
   SPrintDate
   Formats a DATE to a string
   ------------------------------------------------------------- */
void SPrintDate(DATE * Day, char *buffer)
{
  sprintf(buffer, "%02d.%02d.%4d-%02d:%02d:%02d", Day->Month, Day->Day,
	  Day->Year, Day->Hour, Day->Min, Day->Sec);
}

/*****************************************************************************
  IsNewMonth()
*****************************************************************************/
uchar IsNewMonth(DATE * Now, int Interval)
{
  int Year;
  int Month;
  int Day;
  int Hour;
  int Min;
  double Sec;
  double Julian;

  Julian = Now->Julian - ((double) Interval) / SECPDAY;
  JulianDayToGregorian(Julian, &Year, &Month, &Day, &Hour, &Min, &Sec);
  if (Month != Now->Month)
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
uchar Before(DATE * Day1, DATE * Day2)
{
  if (Day1->Julian < Day2->Julian)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  After()
*****************************************************************************/
uchar After(DATE * Day1, DATE * Day2)
{
  if (Day1->Julian > Day2->Julian)
    return TRUE;
  else
    return FALSE;
}

/*****************************************************************************
  IncreaseTime()
*****************************************************************************/
void IncreaseTime(TIMESTRUCT * Time)
{
  double Sec;

  (Time->Step)++;
  Time->DayStep = (Time->DayStep + 1) % Time->NDaySteps;
  Time->Current.Julian = Time->Start.Julian +
    ((double) Time->Step) * (((double) Time->Dt) / SECPDAY);
  JulianDayToGregorian(Time->Current.Julian, &(Time->Current.Year),
		       &(Time->Current.Month), &(Time->Current.Day),
		       &(Time->Current.Hour), &(Time->Current.Min), &Sec);
  Time->Current.Sec = (int) Sec;
  Time->Current.JDay = DayOfYear(Time->Current.Year, Time->Current.Month,
				 Time->Current.Day);
}

/*****************************************************************************
  IncreaseVariableTime()
*****************************************************************************/
void IncreaseVariableTime(TIMESTRUCT *Time, float VariableDT, TIMESTRUCT *NextTime)
{
  double Sec;

  Time->Current.Julian = Time->Current.Julian + VariableDT/SECPDAY;

  if(After(&(Time->Current), &(NextTime->Current)) || 
     IsEqualTime(&(Time->Current), &(NextTime->Current))) {

    (Time->Step)++;
    Time->DayStep = (Time->DayStep + 1) % Time->NDaySteps;
  }

  JulianDayToGregorian(Time->Current.Julian, &(Time->Current.Year),
		       &(Time->Current.Month), &(Time->Current.Day),
		       &(Time->Current.Hour), &(Time->Current.Min), 
		       &(Sec));
  Time->Current.Sec = (int) Sec;
  Time->Current.JDay = DayOfYear(Time->Current.Year, Time->Current.Month,
				 Time->Current.Day);
 
}

/******************************************************************************/
/*				   InitTime()                                 */
/******************************************************************************/
int InitTime(TIMESTRUCT * Time, DATE * Start, DATE * End, DATE * StartRadar,
	     DATE * StartMM5, int Dt)
{
  const char *Routine = "InitTime";
  int tmpsecond;

  Time->Dt = Dt;		/* timestep in seconds */

  /* Just to be sure, recalculate the Julian and JDay fields of all the dates */
  if (Start != NULL) {
    CopyDate(&(Time->Start), Start);
    Time->Start.JDay = DayOfYear(Time->Start.Year, Time->Start.Month,
				 Time->Start.Day);
    Time->Start.Julian = GregorianToJulianDay(Time->Start.Year,
					      Time->Start.Month,
					      Time->Start.Day,
					      Time->Start.Hour,
					      Time->Start.Min, Time->Start.Sec);
    CopyDate(&(Time->Current), &(Time->Start));
  }

  if (End != NULL) {
    CopyDate(&(Time->End), End);
    Time->End.JDay = DayOfYear(Time->End.Year, Time->End.Month, Time->End.Day);
    Time->End.Julian = GregorianToJulianDay(Time->End.Year, Time->End.Month,
					    Time->End.Day, Time->End.Hour,
					    Time->End.Min, Time->End.Sec);
  }
  if (StartRadar != NULL) {
    CopyDate(&(Time->StartRadar), StartRadar);
    Time->StartRadar.JDay =
      DayOfYear(Time->StartRadar.Year, Time->StartRadar.Month,
		Time->StartRadar.Day);
    Time->StartRadar.Julian =
      GregorianToJulianDay(Time->StartRadar.Year, Time->StartRadar.Month,
			   Time->StartRadar.Day, Time->StartRadar.Hour,
			   Time->StartRadar.Min, Time->StartRadar.Sec);
  }
  if (StartMM5 != NULL) {
    CopyDate(&(Time->StartMM5), StartMM5);
    Time->StartMM5.JDay =
      DayOfYear(Time->StartMM5.Year, Time->StartMM5.Month, Time->StartMM5.Day);
    Time->StartMM5.Julian =
      GregorianToJulianDay(Time->StartMM5.Year, Time->StartMM5.Month,
			   Time->StartMM5.Day, Time->StartMM5.Hour,
			   Time->StartMM5.Min, Time->StartMM5.Sec);
  }

  if (Start != NULL && End != NULL) {
    if (After(&(Time->Start), &(Time->End)))
      return FALSE;
    Time->Step = 0;
    Time->NDaySteps = SECPDAY / Time->Dt;
    if ((SECPDAY - Time->NDaySteps * Time->Dt) != 0)
      ReportError((char *) Routine, 7);
    /* the 0.5 in the next line is because Julian day is measured from midday,
       and  we want since midnight */
    tmpsecond = Round(((Time->Start.Julian - 0.5) -
		       floor(Time->Start.Julian - 0.5)) * SECPDAY);
    Time->DayStep = tmpsecond / Time->Dt;
    /* The following line insures that the first timestep of the day would
       coincide with midnight */
    if (tmpsecond - Time->DayStep * Time->Dt != 0)
      ReportError((char *) Routine, 9);
    Time->NTotalSteps = NumberOfSteps(&(Time->Start), &(Time->End), Time->Dt);
    if (Time->NTotalSteps == FAIL)
      return FALSE;
  }

  return TRUE;
}

/******************************************************************************/
/*			     GregorianToJulianDay()                           */
/******************************************************************************/
/*
** Takes a date, and returns a Julian day. A Julian day is the number of
** days since some base date  (in the very distant past).
** Handy for getting date of x number of days after a given Julian date
** (use jdate to get that from the Gregorian date).
** Author: Robert G. Tantzen, translator: Nat Howard
** Translated from the algol original in Collected Algorithms of CACM
** (This and jdate are algorithm 199).
** Copied from the xmgr auxiliary functions Mon Feb  1 14:16:40 1999
*/
double GregorianToJulianDay(int year, int mon, int day, int h, int mi,
			    double se)
{
  long m = mon;
  long d = day;
  long y = year;
  long c;
  long ya;
  long j;
  double seconds = h * 3600.0 + mi * 60 + se;

  if (m > 2)
    m -= 3;
  else {
    m += 9;
    --y;
  }
  c = y / 100L;
  ya = y - (100L * c);
  j = (146097L * c) / 4L + (1461L * ya) / 4L +
    (153L * m + 2L) / 5L + d + 1721119L;
  if (seconds < 12 * 3600.0) {
    j--;
    seconds += 12.0 * 3600.0;
  }
  else {
    seconds = seconds - 12.0 * 3600.0;
  }
  return (j + (seconds / 3600.0) / 24.0);
}

/******************************************************************************/
/*			     JulianDayToGregorian()                           */
/******************************************************************************/
/* Julian date converter. Takes a julian date (the number of days since
** some distant epoch or other), and returns an int pointer to static space.
** ip[0] = month;
** ip[1] = day of month;
** ip[2] = year (actual year, like 1977, not 77 unless it was  77 a.d.);
** ip[3] = day of week (0->Sunday to 6->Saturday)
** These are Gregorian.
** Copied from Algorithm 199 in Collected algorithms of the CACM
** Author: Robert G. Tantzen, Translator: Nat Howard
** Copied from the xmgr auxiliary functions Mon Feb  1 14:16:40 1999
*/
void JulianDayToGregorian(double jd, int *y, int *m, int *d, int *h, int *mi,
			  double *sec)
{
  static int ret[4];

  long j = jd;
  double tmp, frac = jd - j;

  /* The following four lines are added so that we round to the whole seconds, 
     Bart Nijssen Wed Feb  3 08:31:50 1999 */
  if (rint(frac * SECPDAY) != floor(frac * SECPDAY)) {
    jd += 1. / SECPDAY;
    j = jd;
    frac = jd - j;
  }
  /* end of added code */

  if (frac >= 0.5) {
    frac = frac - 0.5;
    j++;
  }
  else {
    frac = frac + 0.5;
  }

  ret[3] = (j + 1L) % 7L;
  j -= 1721119L;
  *y = (4L * j - 1L) / 146097L;
  j = 4L * j - 1L - 146097L * *y;
  *d = j / 4L;
  j = (4L * *d + 3L) / 1461L;
  *d = 4L * *d + 3L - 1461L * j;
  *d = (*d + 4L) / 4L;
  *m = (5L * *d - 3L) / 153L;
  *d = 5L * *d - 3 - 153L * *m;
  *d = (*d + 5L) / 5L;
  *y = 100L * *y + j;
  if (*m < 10)
    *m += 3;
  else {
    *m -= 9;
    *y += 1;
  }
  tmp = 3600.0 * (frac * 24.0);
  *h = (int) (tmp / 3600.0);
  tmp = tmp - *h * 3600.0;
  *mi = (int) (tmp / 60.0);
  *sec = tmp - *mi * 60.0;
}

/******************************************************************************/
/*				  DayOfWeek()                                 */
/******************************************************************************/
int DayOfWeek(double j)
{
  j += 0.5;
  return (int) (j + 1) % 7;
}

/*******************************************************************************
  Test main. Compile by typing:
  gcc -Wall -g -o test_calendar -DTEST_CALENDAR Calendar.c equal.c ReportError.c
  -lm
  then run the program by typing test_calendar
*******************************************************************************/
#ifdef TEST_CALENDAR
int main(int argc, char **argv)
{
  int Dt;
  DATE Start;
  DATE End;
  TIMESTRUCT Time;

  if (argc == 1) {
    printf("\nGive a timestep in seconds on the commend-line\n\n");
    exit(EXIT_FAILURE);
  }
  /* test program so no fancy I/O checking */
  Dt = atoi(argv[1]);

  SScanDate("2/3/1999", &Start);
  PrintDate(&Start, stdout);
  printf("\n");
  SScanDate("2/3/1999-0", &Start);
  PrintDate(&Start, stdout);
  printf("\n");
  SScanDate("2/3/1999-0:0", &Start);
  PrintDate(&Start, stdout);
  printf("\n");
  SScanDate("2/3/1999-0:0:0", &Start);
  PrintDate(&Start, stdout);
  printf("\n");
  SScanDate("3/4/2000", &End);
  PrintDate(&End, stdout);
  printf("\n");
  SScanDate("3/4/2000-16", &End);
  PrintDate(&End, stdout);
  printf("\n");
  SScanDate("3/4/2000-16:30", &End);
  PrintDate(&End, stdout);
  printf("\n");
  SScanDate("3/4/2000-16:30:0", &End);
  PrintDate(&End, stdout);
  printf("\n");

  if (!After(&End, &Start))
    printf("Error in After() function\n");
  if (!Before(&Start, &End))
    printf("Error in Before() function\n");

  printf("\nInitializing time structure\nStart: ");
  InitTime(&Time, &Start, &End, NULL, NULL, Dt);
  PrintDate(&(Time.Start), stdout);
  printf("\nEnd: ");
  PrintDate(&(Time.End), stdout);
  printf("\n");
  printf("Timestep: %d seconds\n", Time.Dt);
  printf("Number of timesteps per day: %d\n", Time.NDaySteps);
  printf("Number of timesteps in model run: %d\n", Time.NTotalSteps);

  while (Time.Step < Time.NTotalSteps) {
    if (IsNewMonth(&(Time.Current), Time.Dt))
      printf("Start of new month\n");
    if (IsNewDay(Time.DayStep))
      printf("Start of new day\n");
    PrintDate(&(Time.Current), stdout);
    printf("\n");
    IncreaseTime(&Time);
  }

  printf("\nInitializing time structure\nStart: ");
  SScanDate("2-27-2000-18:30:00", &Start);
  InitTime(&Time, &Start, &End, NULL, NULL, Dt);
  PrintDate(&(Time.Start), stdout);
  printf("\nEnd: ");
  PrintDate(&(Time.End), stdout);
  printf("\n");
  printf("Timestep: %d seconds\n", Time.Dt);
  printf("Number of timesteps per day: %d\n", Time.NDaySteps);
  printf("Number of timesteps in model run: %d\n", Time.NTotalSteps);

  printf("End of tests\n");

  return EXIT_SUCCESS;
}

#endif
/*****************************************************************************
  SScanMonthDay()
*****************************************************************************/
int SScanMonthDay(char *DateStr, DATE * Day)
{
  char Str[BUFSIZE + 1];
  int i;
  int j;
  int Length;
  int Number[6];
  int DaysPerMonth[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

  if (Str == NULL)
    return FALSE;

  strcpy(Str, DateStr);

  Length = strlen(Str);

  for (i = Length - 1, j = 0; i > 0; i--) {
    if (!isdigit((int) Str[i])) {
      Number[j] = atoi(&Str[i + 1]);
      Str[i] = '\0';
      j++;
    }
  }

  Number[j] = atoi(Str);

  if (j < 1 || j > 5)
    return FALSE;

  Day->Month = Number[j--];
  Day->Day = Number[j--];
  Day->Year = Number[j--];
  if (j >= 0)
    Day->Hour = Number[j--];
  else
    Day->Hour = 0;
  if (j >= 0)
    Day->Min = Number[j--];
  else
    Day->Min = 0;
  if (j >= 0)
    Day->Sec = Number[j--];
  else
    Day->Sec = 0;

  if (IsLeapYear(Day->Year))
    DaysPerMonth[1] = 29;
  else
    DaysPerMonth[1] = 28;

  if (Day->Month < 1 || Day->Month > MONTHPYEAR)
    return FALSE;
  if (Day->Day < 1 || Day->Day > DaysPerMonth[Day->Month - 1])
    return FALSE;
  if (Day->Hour < 0 || Day->Hour > 23)
    return FALSE;
  if (Day->Min < 0 || Day->Min > 59)
    return FALSE;
  if (Day->Sec < 0 || Day->Sec > 59)
    return FALSE;

  Day->JDay = DayOfYear(Day->Year, Day->Month, Day->Day);
  Day->Julian =
    GregorianToJulianDay(Day->Year, Day->Month, Day->Day, Day->Hour, Day->Min,
			 Day->Sec);

  return TRUE;
}

/*---------------Funtion rint (round to nearest integer)----------------------*/
/* float rint(float x) */
/* { */
/* return floor(x + 0.5); */
/* } */
