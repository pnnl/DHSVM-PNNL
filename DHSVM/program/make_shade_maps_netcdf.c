/*
 * SUMMARY:     make the dhsvm shade maps for a given dem 
 * USAGE:        
 *
 * AUTHOR:       Pascal Storck
 * E-MAIL:       pstorck@lightmail.com
 * ORIG-DATE:    June-2000
 * Last Change:  Feb-2013
 *               
 * DESCRIP-END.cd
 * $Id: make_dhsvm_shade_maps.c,v 3.1 2013/02/4 Ning Exp $  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fifoNetCDF.h"
#include "sizeofNetCDF.h"
#include "data.h"

#define DEGPRAD       57.29578  /* degree per radian */
#define MINPDEG        4.       /* minutes per degree longitude */
#undef PI
#define PI             3.14159265358979323846
#define RADPHOUR       0.2617994 /* radians per hour: Earth's Rotation 
                                    (2 PI rad/day) * (1 day/24 h) */
#define RADPDEG     (PI/180.0)   /* radians per degree */
#define SOLARCON    1360.       /* Solar constant (W/m^2) */
#define SECPMIN 60
#define SECPHOUR 3600
#define SECPDAY 86400
#define MINPHOUR 60
#define MINPDAY 1440
#define HOURPDAY 24
#define DAYPWEEK 7
#define DAYPYEAR 365
#define MONTHPYEAR 12
#define FALSE     0
#define TRUE      1


int GetNumber(char *numberStr);

float GetFloat(char *numberStr);

int CopyDouble(double *Value, char *Str, const int NValues);

int DayOfYear(int Year, int Month, int Day);

unsigned char IsLeapYear(int Year);

void SolarDay(int DayOfYear, float Longitude, float Latitude, 
              float StandardMeridian, float *NoonHour, float *Declination, 
              float *HalfDayLength, float *Sunrise, float *Sunset, 
              float *TimeAdjustment, float *SunEarthDist);

void SolarHour(float Latitude, float LocalHour, float Dt, float NoonHour, float *solar_hour,
               float Declination, float Sunrise, float Sunset, 
               float TimeAdjustment, float SunEarthDist, 
               float *SineSolarAltitude, int *DayLight, float *SolarTimeStep,
               float *SunMax, float *SolarAzimuth);

void CalcSlopeAspect(int nRows, int nCols, float dx, float **elev, float ***slope, float ***aspect);

void CalcHillShadeWithTerrainBlocking(int nColss,int nCols,float dx,float max_elev,float **elev,
				     float sal,float saz,float **slope,float **aspect,
				     float ***hillshade);

int main(int argc, char **argv)
{
  FILE   *demfile,*outfile,*outfile2;
  char   VarName[255];	
  int    flag;
  char   demfilename[255],outfilename[255],outfilename2[255];
  int    nRows;                    /* Number of rows */
  int    nCols;                    /* Number of columns */
  float  *temp;
  float  **elev,**slope,**aspect,**hillshade;
  unsigned char **outputgrid;
  void  *Array;
  int    i;
  int    ny,nx;
  float  lx,ly;
  double max_angle,angle;
  float  saz,sal;
  double theta;
  float  dx;
  float  x,y,sx,sy,dz,dist;
  float  mx,my;
  float  start_elev;
  float  max_elev,min_elev;
  float  target_row,target_col;
  int    stop_flag;
  float  a,b,c,d,e,f,g,h,j;
  float  dzdx,dzdy,rr;
  float  safe_distance;
  int    month,day,year,jday;
  float  hour;
  int    ihour;
  float  dt;
  int    count;
  float  outstep;
  int    stepsperday;
  int    daylight;
  int    sunlight;
  float  standardmeridian,latitude,longitude;
  float  noon_hour,  declination, halfdaylength, solar_hour;
  float  sunrise, sunset, timeadjustment, sunearthdistance;
  float  sinesolaraltitude, solartimestep, sunmax, solarazimuth;
  float  beam,diffuse;
  MAPSIZE Map;
  MAPDUMP DMap;

  if(argc < 14) {
    printf("usage is: make_dhsvm_shade_maps:  \n");
    printf("demfilename  \n");
    printf("outfilename  \n");
    printf("nrows, ncols \n");
    printf("cellsize (in the same units as the dem elevation)\n");
    printf("longitude and latitude of the site (dd)\n");
    printf("longitude of location for met file time stamp\n");
    printf("year month day output_time_step (hours)\n");
	printf("Xorig Yorig");
    exit(-1);
  }
  /* note: this program will loop over all time, starting at 0 and advancing */
  /* every hour */
  /* output images ranging from 0 to 255 are made at every output_time_step */
  /* these are in the proper format for DHSVM */

  strcpy(demfilename, argv[1]);   /* name of the nc_float dem input file - no header */
  strcpy(outfilename, argv[2]);   /* name of the nc_float hillshade output file        */
  nRows = GetNumber(argv[3]);
  nCols = GetNumber(argv[4]);
  dx = GetFloat(argv[5]); /* the cellsize of the dem */                            
  longitude=GetFloat(argv[6]) * RADPDEG;
  latitude=GetFloat(argv[7]) * RADPDEG;
  standardmeridian=GetFloat(argv[8]) * RADPDEG;
  year = GetNumber(argv[9]);
  month = GetNumber(argv[10]);
  day = GetNumber(argv[11]);
  outstep = GetFloat(argv[12]);
  /* extreme west coordinatee */
  if (!(CopyDouble(&Map.Xorig, argv[13], 1)))
	  exit (-1);;
  /* exterme north coordinate */
  if (!(CopyDouble(&Map.Yorig, argv[14], 1)))
	  exit (-1);;
  printf("calculating shade map for %d / %d / %d \n",month,day,year);

  Map.X = 0;
  Map.Y = 0;
  Map.OffsetX = 0;
  Map.OffsetY = 0;
  Map.NX = nCols;
  Map.NY = nRows;
  Map.DX = dx;
  Map.DY = dx;
  Map.DXY = (float) sqrt(Map.DX * Map.DX + Map.DY * Map.DY);

  strcpy(DMap.FileName, outfilename);
  DMap.ID = 304;				
  DMap.Layer = 1;
  DMap.Resolution = MAP_OUTPUT;	/* Full resolution maps */
  strcpy(DMap.Name, "Shade.Factor");
  strcpy(DMap.LongName, "Shade Factor");
  strcpy(DMap.Format, "%d");
  strcpy(DMap.FileLabel, "Shade Factor");
  strcpy(DMap.Units, "");
  DMap.NumberType = NC_BYTE;         
  DMap.MaxVal = 0;
  DMap.MinVal = 0;

  /* allocate memory */
  if (!((elev) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((elev)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }
  if (!(Array = (unsigned char *) calloc(nCols*nRows, sizeof(unsigned char))))
     exit(-1);
  if (!(temp = (float *) calloc(nRows*nCols, sizeof(float))))
    exit(-1);
  if (!((slope) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((slope)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }
  if (!((aspect) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((aspect)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }
  if (!((hillshade) = (float**) calloc(nRows, sizeof(float*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((hillshade)[ny] = (float*) calloc(nCols, sizeof(float))))
      exit(-1);
  }
  if (!((outputgrid) = (unsigned char**) calloc(nRows, sizeof(unsigned char*))))
    exit(-1);
  for (ny = 0; ny < nRows; ny++) {
    if (!((outputgrid)[ny] = (unsigned char*) calloc(nCols, sizeof(unsigned char))))
      exit(-1);
  }

  strcpy(VarName, "Basin.DEM");
  flag = Read2DMatrixNetCDF(demfilename, temp, NC_FLOAT, Map.NY, Map.NX, 0,
	       VarName, 0);  
  if (flag == 0){
	  for (ny = 0, i = 0; ny < Map.NY; ny++) {
		  for (nx = 0; nx < Map.NX; nx++, i++) {
			  elev[ny][nx] = temp[i]; }
	  }
  }
  else if (flag == 1){
	  for (ny = Map.NY - 1, i = 0; ny >= 0; ny--) {
		  for (nx = 0; nx < Map.NX; nx++, i++) {
			  elev[ny][nx] = temp[i]; }
	  }
  }
  else exit (-1);
  free(temp);

  max_elev = 0.0;
  for (ny = 0; ny < nRows; ny++) {
    for (nx = 0; nx < nCols; nx++) {
      if(elev[ny][nx]>max_elev) max_elev = elev[ny][nx];
    }
  }

  CalcSlopeAspect(nRows,nCols,dx,elev,&slope,&aspect);

  dt = outstep; // in hours
  stepsperday = (int)(24 / outstep);

  /* netcdf map properties */
  DMap.N = stepsperday;					
  if (!(DMap.DumpDate = (DATE *) calloc(DMap.N, sizeof(DATE))))
      exit(-1);

  CreateMapFileNetCDF(DMap.FileName, DMap.FileLabel, &Map);

  for(i = 0; i < stepsperday; i++){
    hour = (float)i * outstep;
    printf("working on hour %f \n",hour);
    jday = DayOfYear(year,month,day);
  
    SolarDay(jday, longitude, latitude,
           standardmeridian, &noon_hour, 
           &declination, &halfdaylength,
           &sunrise, &sunset, &timeadjustment, &sunearthdistance);

    SolarHour(latitude, hour+dt, dt, noon_hour, &solar_hour,
	    declination,sunrise, sunset,
	    timeadjustment, sunearthdistance,
	    &sinesolaraltitude, &daylight, &solartimestep,
	    &sunmax, &solarazimuth); 

    sal=asin(sinesolaraltitude);
    saz=solarazimuth;

    /*  printf("for %2d/%2d/%4d at met-file-time %5.2f and solar hour %5.2f \n",
	  month,day,year,hour+0.5*dt,solar_hour);
    printf(" sunrise is at %5.2f with sunset at %5.2f and solar alt: %f with azimuth %f \n",
    sunrise,sunset,sal*DEGPRAD,saz*DEGPRAD);*/

    CalcHillShadeWithTerrainBlocking(nRows,nCols,dx,max_elev,elev,
				   sal,saz,slope,aspect,&hillshade);

	DMap.DumpDate[i].Year = year;
    DMap.DumpDate[i].Month = month;
    DMap.DumpDate[i].Day = day;
    DMap.DumpDate[i].JDay = jday;
    DMap.DumpDate[i].Hour = (int)hour;

    /* at this point hillshade is between 0 and 255 */
    /* which is the standard arc-info for the hillshade command */
    /* we need to translate this to the proper dhsvm format */
    /* and output it as an unsigned char */
    for (ny = 0; ny < nRows; ny++) {
	    for (nx = 0; nx < nCols; nx++) {
		    if(sinesolaraltitude>0)
			    if(hillshade[ny][nx]/255/sinesolaraltitude>11.47)
				    outputgrid[ny][nx]=255;
			    else
				    outputgrid[ny][nx]=(unsigned char)(hillshade[ny][nx]/sinesolaraltitude/11.47);
		    else outputgrid[ny][nx]=0;
			((unsigned char *) Array)[ny * Map.NX + nx] = outputgrid[ny][nx];
	    }
    }
    Write2DMatrixNetCDF(DMap.FileName, Array, DMap.NumberType, Map.NY, Map.NX, &DMap, i);
  }

  return EXIT_SUCCESS;
}

/*****************************************************************************
  GetNumber()
*****************************************************************************/
int GetNumber(char *numberStr) 
{
  char *endPtr;
  int number = 0;

  number = (int) strtol(numberStr, &endPtr, 0);
  if (*endPtr != '\0'){
    printf("end Ptr is %s \n",endPtr);
 printf("problem extracting integer from %s \n",numberStr);
    exit(-1);
  }
  return number;
}

/*****************************************************************************
  GetFloat()
*****************************************************************************/
float GetFloat(char *numberStr) 
{
  char *endPtr;
  float number = 0;

  number = (float) strtod(numberStr, &endPtr);
  if (*endPtr != '\0'){
    printf("problem extracting float from %s \n",numberStr);
    exit(-1);
  }

  return number;
}
/*****************************************************************************
  CopyDouble()
*****************************************************************************/
int CopyDouble(double *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = strtod(Str, &EndPtr);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}

/*****************************************************************************
  CalcSolar()
*****************************************************************************/
/*
 * SUMMARY:      CalcSolar.c - inline solar calculations
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Batelle Pacific Northwest Laboratories
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Jul-96
 * Last Change: Wed Feb  3 23:11:44 1999 by Bart Nijssen <nijssen@u.washington.edu>
 * DESCRIPTION:  These functions make inline solar radiation calculations
 *               that take into slope and aspect of the pixel, but do not
 *               account for shadowing of neighbouring pixels.  These
 *               routines can be used instead of IPW, but it is up to the
 *               user to choose one method (IPW) or the other (inline
 *               calculations).  
 * DESCRIP-END.
 * FUNCTIONS:    SolarDay()
 *               SolarHour()
 *               SolarAngle()
 *               SolarConst()
 * COMMENTS:     
 */
/*****************************************************************************
  Function name: SolarDay()

  Purpose:	 This subroutine calculates daily solar values 

  Required:
   int DayOfYear	  - day of year (January 1 is 1)
   float Longitude        - site longitude (rad)
   float Latitude         - site latitude (rad)
   float StandardMeridian - longitude of time zone of standard meridian (rad)

  Returns: void

  Modifies:
    float *NoonHour       - true solar noon (hr)
    float *Declination    - solar Declination (rad)
    float *HalfDayLength  - half-day length (hr)
    float *Sunrise        - time of Sunrise (hr)
    float *Sunset         - time of Sunset (hr)
    float *TimeAdjustment - required adjustment to local time (hr)
    float *SunEarthDist   - distance from sun to earth

  Comments     : EXECUTE AT START OF EACH DAY 
*****************************************************************************/
void SolarDay(int DayOfYear, float Longitude, float Latitude, 
              float StandardMeridian, float *NoonHour, float *Declination, 
              float *HalfDayLength, float *Sunrise, float *Sunset, 
              float *TimeAdjustment, float *SunEarthDist)
{
  float B;                      /* coefficient for equation of time */
  float EqnOfTime;              /* adjustment for equation of time (min) */
  float LongitudeAdjust;        /* adjustment for longitude (min) */
  float CosineHalfDayLength;    /* cosine of the half-day length */
  
  
  /* note need to check if day light savings time calculate adjustment for 
     true solar time longitude adjustment add 4 min per degree away from 
     StandardMeridian (4 min/degree * 180 degree/pi radian) */
  LongitudeAdjust = (MINPDEG * DEGPRAD)*(StandardMeridian - Longitude);
  
  /* equation of time */
  B = (2.0*PI*(DayOfYear - 81)) /  364.;
  EqnOfTime = 9.87*sin(2*B) - 7.53*cos(B) - 1.5*sin(B);
  
  /* adjustment factor to convert local time to solar time
     solar time = local time + TimeAdjustment  */
  /* for example from GMT to the west coast of the us (PST)*/
  /* is a -8 hour shift, i.e. PST = GMT-8 */
  *TimeAdjustment = -(LongitudeAdjust + EqnOfTime) / MINPHOUR;
  
  /* work in solar time  */
  *NoonHour = 12.0;
  
  /* solar Declinationation  */
  *Declination = .4098*sin(2*PI*(284 + DayOfYear) / DAYPYEAR);
  
  /* half-day length  */    
  CosineHalfDayLength = - tan(Latitude)*tan(*Declination);
  if (CosineHalfDayLength >= 1.0)
    *HalfDayLength = PI;
  else
    *HalfDayLength = acos(CosineHalfDayLength);
  
  /* convert HalfDayLength from radians to Hours 
     1 radian = (180 deg / PI) * (1 hr / 15 degrees rotation) */
  *HalfDayLength = *HalfDayLength / RADPHOUR;
  
  /* solar time of Sunrise and Sunset  */
  *Sunrise = *NoonHour - *HalfDayLength;
  *Sunset = *NoonHour + *HalfDayLength;
  
  /* calculate the sun-earth distance */
  *SunEarthDist = 1.0 + 0.033*cos( RADPDEG * (360. * DayOfYear  / 365 ) );
}

/*****************************************************************************
  Function name: SolarHour()  
  
  Purpose: This subroutine calculates position of the sun as a function of
           the time of day, the length of time the sun is above the horizon,
           and the maximum radiation.

  Required:
    float Latitude		- site laditude (rad)
    float LocalHour		- local time (hr)
    float Dt			- length of current timestep (hr)
    float NoonHour		- true solar noon (hr)
    float Declination		- solar Declination (rad)
    float Sunrise		- time of Sunrise (hr)
    float Sunset		- time of Sunset (hr)
    float TimeAdjustment	- required adjustment to convert local time
				  to solar time (hr)
    float SunEarthDist          - distance from Sun to Earth 

  Returns: void

  Modifies:
    float *SineSolarAltitude - sine of sun's SolarAltitude 
    int *DayLight	     - FALSE: measured solar radiation and the sun is
                                      below the horizon.  
			       TRUE: sun is above the horizon
    float *SolarTimeStep     - fraction of the timestep the sun is above the 
                               horizon 
    float *SunMax            - calculated solar radiation at the top of the 
                               atmosphere (W/m^2) 
                                                 
  Comments     : EXECUTE AT START OF EACH TIMESTEP 
*****************************************************************************/
void SolarHour(float Latitude, float LocalHour, float Dt, float NoonHour, float *solar_hour,
               float Declination, float Sunrise, float Sunset, 
               float TimeAdjustment, float SunEarthDist, 
               float *SineSolarAltitude, int *DayLight, float *SolarTimeStep,
               float *SunMax, float *SolarAzimuth)
{
  float SolarAltitude;          /* SolarAltitude of sun from horizon (rads) */
  float SolarZenith;            /* sun zenith angle (rads) */
  float StartHour = 0;		/* currect Hour in solar time (hr) */
  float EndHour = 0;		/* mid-point of current solar Hour (hr) */
  float Hour;                   /* angle of current "halfhr" from solar noon
                                   (rads) */

  /* NOTE THAT HERE Dt IS IN HOURS NOT IN SECONDS */
  *SunMax = 0.0; 
  *SolarTimeStep=1.0;
  
  /* all calculations based on hour, the solar corrected local time */
  Hour = LocalHour + TimeAdjustment;
  if(Hour<0) Hour+=24;
  if(Hour>24) Hour-=24;
  
  *solar_hour=Hour;
  *DayLight = FALSE;
  if ((Hour > Sunrise) && ((Hour-Dt) < Sunset))
    *DayLight = TRUE;
  
    /*  if (*DayLight == TRUE) {*/
    /* sun is above the horizon */   
    if( Dt > 0.0 ) {  
      /* compute average solar SolarAltitude over the timestep */
      StartHour = (Hour - Dt > Sunrise) ? Hour - Dt : Sunrise;
      EndHour = (Hour < Sunset) ? Hour : Sunset;
      
      /*  convert to radians  */
      StartHour = RADPHOUR * (StartHour - NoonHour);                   
      EndHour = RADPHOUR * (EndHour - NoonHour);   
      *SolarTimeStep = EndHour - StartHour;     
      
      /*  determine the average geometry of the sun angle  */
      *SineSolarAltitude = sin(Latitude)*sin(Declination)
        + (cos(Latitude)*cos(Declination)*(sin(EndHour) - sin(StartHour))
           / *SolarTimeStep);
    } 
    else {
      Hour = RADPHOUR * ( Hour - NoonHour );
      *SineSolarAltitude = sin(Latitude)*sin(Declination)
        + (cos(Latitude)*cos(Declination)*cos(Hour));
    }

    SolarAltitude = asin(*SineSolarAltitude);
    SolarZenith = PI/2-SolarAltitude;
    
    *SolarAzimuth = ((sin(Latitude)*(*SineSolarAltitude) - sin(Declination))
		     / (cos(Latitude)*sin(SolarZenith)));
    
    if (*SolarAzimuth > 1.) 
      *SolarAzimuth = 1.;
    
    *SolarAzimuth=acos(-(*SolarAzimuth));

    if (Dt > 0.0) {
      if (fabs(EndHour) > fabs(StartHour)) {
	*SolarAzimuth = 2*PI - (*SolarAzimuth);
      }
    }
    else {
      if (Hour > 0) {
	*SolarAzimuth = 2*PI - (*SolarAzimuth);
      }
    }
    /*     printf("Solar Azimuth = %g\n", *SolarAzimuth*180./PI); */
    
    *SunMax = SOLARCON * SunEarthDist * *SineSolarAltitude;	
    /*  }*/
}

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
unsigned char IsLeapYear(int Year) 
{
  if ((Year % 4 == 0 && Year % 100 != 0) || Year % 400 == 0) 
    return TRUE;
  return FALSE;
}

/*****************************************************************************
  CalcSlopeAspect()
*****************************************************************************/
void CalcSlopeAspect(int nRows, int nCols, float dx, float **elev, float ***slope, float ***aspect)
{
  int ny,nx;
  float x,y,sx,sy,dz,dist;
  float a,b,c,d,e,f,g,h,j;
  float dzdx,dzdy,rr;

  for (ny = 0; ny < nRows; ny++) {
      for (nx = 0; nx < nCols; nx++) {
	if(nx==0 || ny == 0 || ny == nRows-1 || nx == nCols-1) {
	  (*slope)[ny][nx]=0;
	  (*aspect)[ny][nx]=0;
	}
	else
	  {
	    /* get the slope in rads based on the arc-info method */
	    /* see the arc-info help pages for more info */	  
	    /* also find the aspect using the arc-info method */
      
	    a=elev[ny-1][nx-1];
	    b=elev[ny-1][nx];
	    c=elev[ny-1][nx+1];
	    d=elev[ny][nx-1];
	    e=elev[ny][nx];
	    f=elev[ny][nx+1];
	    g=elev[ny+1][nx-1];
	    h=elev[ny+1][nx];
	    j=elev[ny+1][nx+1];
     	  
	    dzdx=((a+2*d+g)-(c+2*f+j))/(8*dx);
	    dzdy=((a+2*b+c)-(g+2*h+j))/(8*dx);
	    rr=sqrt(dzdx*dzdx+dzdy*dzdy);
	    (*slope)[ny][nx]=atan(rr);

	    if (dzdx == 0.0 && dzdy == 0.0) 
	      (*aspect)[ny][nx] = 0.0;
	    else 
	      (*aspect)[ny][nx] = atan2(dzdx, -dzdy);

	    /*aspect is calculated assuming that x is positive eastward */
	    /* and that y is positive southward */
	    /* aspects are returned as follows*/
	    /*   north = 0, east = 90, south = 180, west = -90 */
	    /*  sw = -135 */

	    /* to convert to 0-360 clockwise from north */ 
	    if((*aspect)[ny][nx]<0.0) (*aspect)[ny][nx] = 57.29578+(*aspect)[ny][nx];
	 
	  }
      }
  }
}

/*****************************************************************************
  CalcHillShadeWithTerrainBlocking()
*****************************************************************************/
void CalcHillShadeWithTerrainBlocking(int nRows,int nCols,float dx,float max_elev,float **elev,
				     float sal,float saz,float **slope,float **aspect,
				     float ***hillshade)
{
  int ny,nx;
  float lx,ly;
  double max_angle,angle;
  double theta;
  float x,y,sx,sy,dz,dist;
  float mx,my;
  float start_elev;
  int stop_flag;
  float safe_distance;

  if(sal>0){
    for (ny = 0; ny < nRows; ny++) {
      for (nx = 0; nx < nCols; nx++) {  
        (*hillshade)[ny][nx] =   255*(cos(sal)*sin(slope[ny][nx])
			       *cos(aspect[ny][nx]-saz)+sin(sal)
			       *cos(slope[ny][nx]));

	/* at this point hillshade can range from 0 to 255 */

	if ((*hillshade)[ny][nx]<0.0) (*hillshade)[ny][nx]=0.0;
      }
    }
    ly=(float)(nRows*dx-dx);
    lx=(float)(nCols*dx-dx);

    /* utilize a series of checks to speed up program */
    /* if one pixel in the given direction blocks, no need to check the rest in that direction */
    /* if we are beyond the distance that the max elev can't block, then don't check the rest  */
    /*   in that direction */
    /* here x increases eastward and y increases southward */
   for (ny = 0; ny < nRows; ny++) {
      for (nx = 0; nx < nCols; nx++) {
	    start_elev=elev[ny][nx];
	    if (start_elev>0) {
	      safe_distance=tan(sal)*(max_elev-start_elev);
	      sx=(float)nx*dx+0.5*dx;
	      sy=(float)ny*dx+0.5*dx;
	      x=sx;
	      y=sy;

	      while(x>dx && x < lx && y>dx && y<ly) {
	        x=x+((float)sin(saz))*dx;
	        y=y-((float)cos(saz))*dx;
	        dz=elev[(int)(y/dx)][(int)(x/dx)]-start_elev;
	        dist=sqrt((x-sx)*(x-sx)+(y-sy)*(y-sy));
	        if(dist > safe_distance){
	          x=0;
	          y=0;
	        }
	        if(dz>0) {
	           angle=atan((double)(dz/dist));
	           if(angle > sal) {
	              x=0;
	              y=0;
	              (*hillshade)[ny][nx]=0.0;
	           }	    
	         }
	       }
	     }
        }
      }
  }
  else
    {
      for (ny = 0; ny < nRows; ny++) {
	for (nx = 0; nx < nCols; nx++) {
	  (*hillshade)[ny][nx]=0.0;
	}
      }
    }
}
/*****************************************************************************/
