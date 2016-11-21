/*
 * SUMMARY:      FileIONetCDF.c - Functions for NetCDF IO
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-UPDATE:  Fen-2013
 * DESCRIPTION:  Functions for netcdf IO
 * DESCRIP-END.
 * FUNCTIONS:    CreateMapFileNetCDF()
 *               Read2DMatrixNetCDF()
 *               Write2DMatrixNetCDF()
 *               SizeOfNumberType()
 *
 * Modification was made to Read2DMatrix by Ning (2013) 

 * Comment      :First, make sure that the input NETCDF file to Read2DMatrix
                is in 3 dimensions. In this case, the 1st dimention should be time, 
		  the 2nd is y (lat), and 3rd is x (lon). 
		  Using the command 'ncdump -in.nc', we can see that the 
		  (descending/ascending) order in which the NetCDF file uses to stores 
		  coordinates. The order is important b/c it decides in which order    
		  the NetCDF files is read out and stored in a 2D matrix.
		  Simply put, in the input netcdf file, if the first variable value 
		  is the smallest x and y, which is the lower left corner in the spatial 
                sense, then it shouldn't be assigned to Matrix[0][0] which is spatially 
                located at the upper left corner. Instead it should be stored in lower left 
                corner. So the Read2DMatrix function is changed so that it reads x and y 
		  in .nc input file, and output a flag that will be used as an indicator
		  of whethe or not the matrix should be reversed. Note that this program now 
                handles two cases: 1) .nc generated from gdal tool, in which x and y are both 
                in an ascending order; 2) .nc generated from arcmap, in which x is in an 
                ascending and y in an descedning order, and no reverse needs to be made in this 
                case. Just be aware that a time dimention has to be added to arc generated 
		  .nc file. Check the tutorial for instructions. (Ning, Feb 2013)  
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <netcdf.h>
#include <time.h>
#include "settings.h"
#include "data.h"
#include "fifoNetCDF.h"
#include "DHSVMerror.h"
#include "sizeofnt.h"
#include "fileio.h"

static void nc_check_err(const int ncstatus, const int line, const char *file);
static int GenerateHistory(int argc, char **argv, char *History);
static int ncUpdateGlobalHistory(int argc, char **argv, int ncid);
char commandline[] = "Testing the NetCDF file format";

/*****************************************************************************
  Function name: MakeFileNameNetCDF()

  Purpose      : Create a new file name ending in extension ".nc" to
                 indicate that the file is binary

  Required     : 
    Path     - directory for files
    Str1     - first part of filename
    Str2     - second part of filename
    FileName - filename

  Returns      : void

  Modifies     : FileName

  Comments     : 
*****************************************************************************/
void MakeFileNameNetCDF(char *Path, char *Str1, char *Str2, char *FileName)
{
  const char *Routine = "MakeFileNameNetCDF";
  
  MakeFileNameGen(Path, Str1, Str2, ".nc", FileName);
}
/*******************************************************************************
  Function name: CreateMapFileNetCDF()

  Purpose      : Open and close a new file.  If the file already exists it 
                 will be overwritten.

  Required     : 
    FileName  - Name of the new file
    FileLabel - String describing file contents
    Map       - structure with information about spatial extent of model area

  Returns      : void

  Modifies     : 

  Comments     : NetCDF defines all the dimensions in the file before the file 
                 can be written to.  By default it creates the entire file
		 during when nc_endef() is called and fills all positions with
		 _FillValue.  This behavior is turned of here to speed up the
		 initialization process by or'ing  the  NC_NOFILL  flag  into
		 the  mode parameter of nc_create() 
*******************************************************************************/
void CreateMapFileNetCDF(char *FileName, ...)
{
  const char *Routine = "CreateMapFileNetCDF";
  va_list ap;
  MAPSIZE *Map = NULL;		/* pointer to structure with map info */
  char *FileLabel;		/* File label */
  double *Array;
  double missing_value[1];
  char *ptrstr[1];
  int x;
  int y;
  int varideast;
  int varidnorth;
  int varidtime;
  int ncstatus;
  int ncid;
  int dimids[3];		/* time, north, east */

  /****************************************************************************/
  /*                     HANDLE VARIABLE ARGUMENT LIST                        */
  /****************************************************************************/

  va_start(ap, FileName);
  FileLabel = va_arg(ap, char *);
  Map = va_arg(ap, MAPSIZE *);

  /* Go ahead and clobber any existing file */
  ncstatus = nc_create(FileName, NC_CLOBBER | NC_NOFILL, &ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /****************************************************************************/
  /*                              DEFINE MODE                                 */
  /****************************************************************************/

  ncstatus = nc_def_dim(ncid, TIME_DIM, NC_UNLIMITED, &dimids[0]);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_def_dim(ncid, Y_DIM, Map->NY, &dimids[1]);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_def_dim(ncid, X_DIM, Map->NX, &dimids[2]);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /* Define dimension variables and their attributes */
  /* time */
  ncstatus = nc_def_var(ncid, TIME_DIM, NC_DOUBLE, 1, &dimids[0], &varidtime);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varidtime, ATT_NAME, strlen(TIME_DIM),
			     TIME_DIM);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varidtime, ATT_LONGNAME, strlen(TIME_DIM),
			     TIME_DIM);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varidtime, ATT_UNITS, strlen("index"),
			     "index");
  nc_check_err(ncstatus, __LINE__, __FILE__);
/*   ncstatus = nc_put_att_text(ncid, varidtime, ATT_FORMAT, strlen("%g"), "%g"); */
/*   nc_check_err(ncstatus, __LINE__, __FILE__); */

  /* northing */
  ncstatus = nc_def_var(ncid, Y_DIM, NC_DOUBLE, 1, &dimids[1], &varidnorth);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varidnorth, ATT_NAME, strlen(Y_DIM), Y_DIM);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus =
    nc_put_att_text(ncid, varidnorth, ATT_LONGNAME, strlen("Northing"),
		    "Northing");
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varidnorth, ATT_UNITS, strlen("m"), "m");
  nc_check_err(ncstatus, __LINE__, __FILE__);
/*   ncstatus = nc_put_att_text(ncid, varidnorth, ATT_FORMAT, strlen("%g"), "%g"); */
/*   nc_check_err(ncstatus, __LINE__, __FILE__); */

  /* easting */
  ncstatus = nc_def_var(ncid, X_DIM, NC_DOUBLE, 1, &dimids[2], &varideast);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varideast, ATT_NAME, strlen(X_DIM), X_DIM);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varideast, ATT_LONGNAME, strlen("Easting"),
			     "Easting");
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_put_att_text(ncid, varideast, ATT_UNITS, strlen("m"), "m");
  nc_check_err(ncstatus, __LINE__, __FILE__);
/*   ncstatus = nc_put_att_text(ncid, varideast, ATT_FORMAT, strlen("%g"), "%g"); */
/*   nc_check_err(ncstatus, __LINE__, __FILE__); */

  /* Update the history attribute */
  ptrstr[0] = commandline;
  ncstatus = ncUpdateGlobalHistory(1, ptrstr, ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /* Insert the file label and other global attributes */
  ncstatus = nc_put_att_text(ncid, NC_GLOBAL, ATT_COMMENT, strlen(FileLabel),
			     FileLabel);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  missing_value[0] = NA;
  ncstatus = nc_put_att_double(ncid, NC_GLOBAL, ATT_MISSINGVALUE, NC_DOUBLE, 1,
			       missing_value);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /* exit the define mode */
  ncstatus = nc_enddef(ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /****************************************************************************/
  /*              WRITE X AND Y DIMENSIONS TO THE OUTPUT FILE                 */
  /****************************************************************************/

  Array = (double *) calloc(Map->NX, sizeof(double));
  if (Array == NULL)
    ReportError((char *) Routine, 1);
  for (x = 0; x < Map->NX; x++)
    Array[x] = Map->Xorig + x * Map->DX;
  ncstatus = nc_put_var_double(ncid, varideast, Array);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  free(Array);

  Array = (double *) calloc(Map->NY, sizeof(double));
  if (Array == NULL)
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NY; y++)
    Array[y] = Map->Yorig - y * Map->DY;
  ncstatus = nc_put_var_double(ncid, varidnorth, Array);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  free(Array);

  ncstatus = nc_close(ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);
}

/*******************************************************************************
  Function name: Read2DMatrixNetCDF()

  Purpose      : Function to read a 2D array from a file.

  Required     :
    FileName   - name of input file
    Matrix     - address of array data into
    NumberType - code for number type (taken from HDF, see comments at the
                 beginning of InitFileIO.c for more detail)
    NY         - Number of rows
    NX         - Number of columns
    NDataSet   - number of the dataset to read, i.e. the first matrix in a 
                 file is number 0, etc. (this is not used for the NetCDF file,
		 since we can retrieve the variable by name).
    VarName    - Name of variable to retrieve

  Returns      : Number of elements read

  Modifies     : Matrix

  Comments     : NOTE that we cannot modify anything other than the returned
                 Matrix, because we have to stay compatible with Read2DMatrixBin 
*******************************************************************************/
int Read2DMatrixNetCDF(char *FileName, void *Matrix, int NumberType, int NY,
		       int NX, int NDataSet, ...)
{
  const char *Routine = "Read2DMatrixNetCDF";
  char Str[BUFSIZE + 1];
  char dimname[NC_MAX_NAME + 1];
  char *VarName;
  int dimids[3];
  int ndims;
  int ncid;
  int ncstatus;
  nc_type TempNumberType;
  int varid;
  size_t count[3];
  size_t start[3] = { 0, 0, 0 };
  size_t dimlen;
  va_list ap;
  size_t index;			/* index of the time slice being dumped */
  int timid;
  size_t timelen;
  double time;
  double *Ycoord;
  double *Xcoord;  /* lat, lon variables */
  int	LatisAsc, LonisAsc, flag;    /* flag */
  int lon_varid, lat_varid;
  count[0] = 1;
  count[1] = NY;
  count[2] = NX;

  /****************************************************************************/
  /*                   GO THROUGH VARIABLE ARGUMENT LIST                      */
  /****************************************************************************/
  va_start(ap, NDataSet);
  VarName = va_arg(ap, char *);
  index = va_arg(ap, int);

  /****************************************************************************/
  /*                           QUERY NETDCF FILE                              */
  /****************************************************************************/

  ncstatus = nc_open(FileName, NC_NOWRITE, &ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /* check whether the variable exists and get its parameters */
  ncstatus = nc_inq_varid(ncid, VarName, &varid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  ncstatus = nc_inq_var(ncid, varid, NULL, &TempNumberType, &ndims, dimids,
			NULL);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (TempNumberType != NumberType) {
    sprintf(Str, "%s: nc_type for %s is different than expected.\n",
	    FileName, VarName);
    ReportWarning(Str, 58);
  }

  /* make sure that the x and y dimensions have the correct sizes */
  ncstatus = nc_inq_dim(ncid, dimids[1], dimname, &dimlen);  

  nc_check_err(ncstatus, __LINE__, __FILE__);

  ncstatus = nc_inq_varid(ncid, dimname, &lat_varid);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (dimlen != NY)
	  ReportError(VarName, 59);
  Ycoord = (double *) calloc(dimlen, sizeof(double));
  if (Ycoord == NULL)
    ReportError((char *) Routine, 1);
    /* Read the latitude coordinate variable data. */
  ncstatus = nc_get_var_double(ncid, lat_varid, Ycoord);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  /* A quick check if the lat, long are in a ascending order. 
  If so, matrix must be flipped so the first value in the matrix will be 
  assigned to the lower left corner cell that has lowest X (lon) & Y (lat) value. 
  (see more comments in the header of this C file). */
  LatisAsc = 1;
  if( Ycoord[0] > Ycoord[NY - 1] ) 
	  LatisAsc = 0;

  ncstatus = nc_inq_dim(ncid, dimids[2], dimname, &dimlen);

  nc_check_err(ncstatus, __LINE__, __FILE__);

  ncstatus = nc_inq_varid(ncid, dimname, &lon_varid);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (dimlen != NX)
    ReportError(VarName, 60);
  Xcoord = (double *) calloc(NX, sizeof(double));
  if (Xcoord == NULL)
    ReportError((char *) Routine, 1);
  /* Read the latitude coordinate variable data. */
  ncstatus = nc_get_var_double(ncid, lon_varid, Xcoord);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  LonisAsc = 1;
  if( Xcoord[0] > Xcoord[NX - 1] ) 
	  LonisAsc = 0;

  if (LonisAsc == 0){
	  printf("The current program does not handle the cases when longitude or X \
values in the .nc input in an descending order. You can either change the input \
.nc file format outside of this program. or you can easily modify this program to \
fit your needs. \n");
	  ReportError("Improper NetCDF input files", 58);
  }
  if ((LatisAsc == 0) & (LonisAsc == 1))
	  flag = 0;
  if ((LatisAsc == 1) & (LonisAsc == 1))
	  flag = 1;
  
  /* see whether the time dimension needs to be updated (the assumption is that
     the same index value refers to the same moment in time.  Since currently we
     make separate files for separate variables this is OK) */
  ncstatus = nc_inq_dimlen(ncid, dimids[0], &timelen);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (timelen < index + 1) {	/* need to add one to time */
    ncstatus = nc_inq_varid(ncid, TIME_DIM, &timid);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    time = (double) index;
    ncstatus = nc_put_var1_double(ncid, timid, &index, &time);
    nc_check_err(ncstatus, __LINE__, __FILE__);
  }
  start[0] = index;
  /****************************************************************************/
  /*                             READ VARIABLE                                */
  /****************************************************************************/

  switch (NumberType) {
  case NC_BYTE:
    ncstatus = nc_get_vara_uchar(ncid, varid, start, count, Matrix);
    break;
  case NC_CHAR:
    ncstatus = nc_get_vara_text(ncid, varid, start, count, Matrix);
    break;
  case NC_SHORT:
    ncstatus = nc_get_vara_short(ncid, varid, start, count, Matrix);
    break;
  case NC_INT:
    ncstatus = nc_get_vara_int(ncid, varid, start, count, Matrix);
    break;
    /* 8 bit integer not yet implemented in NetCDF 3.4, but anticipated in
       future versions */
    /*   case NC_LONG: */
    /*     ncstatus = nc_put_vara_long(ncid, varid, start, count, Matrix); */
    /*     break; */
  case NC_FLOAT:
    ncstatus = nc_get_vara_float(ncid, varid, start, count, Matrix);
    break;
  case NC_DOUBLE:
    ncstatus = nc_get_vara_double(ncid, varid, start, count, Matrix);
    break;
  default:
    ReportError((char *) Routine, 40);
    break;
  }
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /****************************************************************************/
  /*                                CLEAN UP                                  */
  /****************************************************************************/

  ncstatus = nc_close(ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  return flag;
}

/*******************************************************************************
  Function name: Write2DMatrixNetCDF()

  Purpose      : Function to write a 2D array to a file.  Data is appended to
                 the end of the file.

  Required     :   
    FileName   - name of output file
    Matrix     - address of array containing matrix elements
    NumberType - code for number type (taken from HDF, see comments at the
                 beginning of InitFileIO.c for more detail)
    NY         - Number of rows
    NX         - Number of columns

  Returns      : Number of elements written 

  Modifies     :

  Comments     :
*******************************************************************************/
int Write2DMatrixNetCDF(char *FileName, void *Matrix, int NumberType, int NY,
			int NX, ...)
{
  const char *Routine = "Write2DMatrixNetCDF";
  double time;
  int dimids[3];		/* time, north, east */
  size_t index;			/* index of the time slice being dumped */
  int ncid;
  int ncstatus;
  int timid;
  int varid;
  size_t count[3];
  size_t start[3] = { 0, 0, 0 };
  size_t timelen;
  va_list ap;
  MAPDUMP *DMap;

  count[0] = 1;
  count[1] = NY;
  count[2] = NX;

  /****************************************************************************/
  /*                   GO THROUGH VARIABLE ARGUMENT LIST                      */
  /****************************************************************************/
  va_start(ap, NX);
  DMap = va_arg(ap, MAPDUMP *);
  index = va_arg(ap, int);

  /****************************************************************************/
  /*                           QUERY NETDCF FILE                              */
  /****************************************************************************/

  ncstatus = nc_open(FileName, NC_WRITE, &ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /* get dimension ID's */
  nc_inq_dimid(ncid, TIME_DIM, &(dimids[0]));
  nc_check_err(ncstatus, __LINE__, __FILE__);
  nc_inq_dimid(ncid, Y_DIM, &(dimids[1]));
  nc_check_err(ncstatus, __LINE__, __FILE__);
  nc_inq_dimid(ncid, X_DIM, &(dimids[2]));
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /* see whether variable has been defined; if not defined, define it now */
  ncstatus = nc_inq_varid(ncid, DMap->Name, &varid);
  if (ncstatus == NC_ENOTVAR) {	/* Variable not defined */

    ncstatus = nc_redef(ncid);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    ncstatus = nc_def_var(ncid, DMap->Name, DMap->NumberType, 3, dimids,
			  &varid);
    nc_check_err(ncstatus, __LINE__, __FILE__);

    /* write variable attributes */
    ncstatus = nc_put_att_text(ncid, varid, ATT_NAME, strlen(DMap->Name),
			       DMap->Name);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    ncstatus = nc_put_att_text(ncid, varid, ATT_LONGNAME,
			       strlen(DMap->LongName), DMap->LongName);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    ncstatus = nc_put_att_text(ncid, varid, ATT_UNITS, strlen(DMap->Units),
			       DMap->Units);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    ncstatus = nc_put_att_text(ncid, varid, ATT_FORMAT, strlen(DMap->Format),
			       DMap->Format);
    nc_check_err(ncstatus, __LINE__, __FILE__);

    ncstatus = nc_enddef(ncid);
    nc_check_err(ncstatus, __LINE__, __FILE__);
  }
  else				/* Variable defined */
    nc_check_err(ncstatus, __LINE__, __FILE__);

  /* see whether the time dimension needs to be updated (the assumption is that
     the same index value refers to the same moment in time.  Since currently we
     make separate files for separate variables this is OK) */
  ncstatus = nc_inq_dimlen(ncid, dimids[0], &timelen);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (timelen < index + 1) {	/* need to add one to time */
    ncstatus = nc_inq_varid(ncid, TIME_DIM, &timid);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    time = (double) index;
    ncstatus = nc_put_var1_double(ncid, timid, &index, &time);
    nc_check_err(ncstatus, __LINE__, __FILE__);
  }
  start[0] = index;

  /****************************************************************************/
  /*                             WRITE VARIABLE                               */
  /****************************************************************************/

  switch (NumberType) {
  case NC_BYTE:
    ncstatus = nc_put_vara_uchar(ncid, varid, start, count, Matrix);
    break;
  case NC_CHAR:
    ncstatus = nc_put_vara_text(ncid, varid, start, count, Matrix);
    break;
  case NC_SHORT:
    ncstatus = nc_put_vara_short(ncid, varid, start, count, Matrix);
    break;
  case NC_INT:
    ncstatus = nc_put_vara_int(ncid, varid, start, count, Matrix);
    break;
    /* 8 bit integer not yet implemented in NetCDF 3.4, but anticipated in
       future versions */
    /*   case NC_LONG: */
    /*     ncstatus = nc_put_vara_long(ncid, varid, start, count, Matrix); */
    /*     break; */
  case NC_FLOAT:
    ncstatus = nc_put_vara_float(ncid, varid, start, count, Matrix);
    break;
  case NC_DOUBLE:
    ncstatus = nc_put_vara_double(ncid, varid, start, count, Matrix);
    break;
  default:
    ReportError((char *) Routine, 40);
    break;
  }
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /****************************************************************************/
  /*                                CLEAN UP                                  */
  /****************************************************************************/

  ncstatus = nc_close(ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  return NY * NX;
}

/*******************************************************************************
  Function name: nc_check_err()

  Purpose      : Check status returned by NetCDF functions.

  Required     : int ncstatus - error status returned by NetCDF functions
  
  Returns      : void

  Modifies     : void

  Comments     :
*******************************************************************************/
static void nc_check_err(const int ncstatus, const int line, const char *file)
{
  char str[BUFSIZE + 1];

  if (ncstatus != NC_NOERR) {
    sprintf(str, "%s, line: %d -- %s", file, line, nc_strerror(ncstatus));
    ReportError((char *) str, 57);
  }
}

/*******************************************************************************
  Function name: GenerateHistory

  Purpose      : Generate a string documenting when and by whom a change was
                 made to the NetCDF file

  Required     : int argc     - number of strings
                 char **argv  - pointer to array of strings
		 char *History - history string
  
  Returns      : 0

  Modifies     : History string

  Comments     : Originally part of another series of programs for NetCDF.  In
                 that case the command-line arguments were recorded in the
		 history. 
*******************************************************************************/
static int GenerateHistory(int argc, char **argv, char *History)
{
  int i;
  char TimeStr[BUFSIZ];
  time_t timer;

  /* Get date and time in local time */
  time(&timer);
  strftime(TimeStr, (size_t) BUFSIZ, "%b %d, %Y %X %z", localtime(&timer));
  strcpy(History, TimeStr);

  /* add user name */
  sprintf(History, "%s by %s:", History, getenv("LOGNAME"));

  /* add command-line arguments */
  for (i = 0; i < argc; i++)
    sprintf(History, "%s %s", History, argv[i]);

  return 0;
}

/*******************************************************************************
  Function name: ncUpdateGlobalHistory

  Purpose      : Update the global history attribute or create one if it does
                 not exist.

  Required     : int argc     - number of strings
                 char **argv  - pointer to array of strings
		 int ncid     - NetCDF file id
  
  Returns      : status (NC_NOERR if successful)

  Modifies     : History attribute

  Comments     : File has to be in define mode before calling ths function.
                 Originally part of another series of programs for NetCDF.  In
                 that case the command-line arguments were recorded in the
		 history. 
*******************************************************************************/
static int ncUpdateGlobalHistory(int argc, char **argv, int ncid)
{
  char update[BUFSIZ + 1];	/* update string */
  char *strOldHistory = NULL;	/* pointer to old history string */
  char *strNewHistory = NULL;	/* pointer to new and improved history string */
  int status = 0;		/* return status */
  size_t length;		/* string length */

  /* Create an update string */
  GenerateHistory(argc, argv, update);

  /* Query the NetCDF file to see if a history attribute is already present 
     The assumption here is that there is only one global history attribute,
     and that the attribute name is all in lower case.  The file is supposed to
     be in define mode  */

  status = nc_inq_att(ncid, NC_GLOBAL, ATT_HISTORY, NULL, &length);
  if (status == NC_ENOTATT) {	/* no history attibrute present */

    status = nc_put_att_text(ncid, NC_GLOBAL, ATT_HISTORY, strlen(update),
			     update);
    nc_check_err(status, __LINE__, __FILE__);
  }
  else if (status == NC_NOERR) {	/* history attribute present */

    strOldHistory = calloc(length + 1, sizeof(char));
    if (strOldHistory == NULL)
      return 1;
    status = nc_get_att_text(ncid, NC_GLOBAL, ATT_HISTORY, strOldHistory);
    nc_check_err(status, __LINE__, __FILE__);
    /* the NetCDF text is normally stored without the C string 
       terminator '\0'.  So check whether it is present at the end
       of the string, if not, append one. */
    if (strOldHistory[length - 1] != '\0')
      strOldHistory[length] = '\0';

    /* place the newest string at the beginning */
    strNewHistory = calloc(length + strlen(update) + 2, sizeof(char));
    if (strNewHistory == NULL)
      return 1;
    sprintf(strNewHistory, "%s\n%s", update, strOldHistory);

    status = nc_put_att_text(ncid, NC_GLOBAL, ATT_HISTORY,
			     strlen(strNewHistory), strNewHistory);
    nc_check_err(status, __LINE__, __FILE__);

    if (strOldHistory != NULL)
      free(strOldHistory);
    if (strNewHistory != NULL)
      free(strNewHistory);
  }

  else
    nc_check_err(status, __LINE__, __FILE__);

  return status;
}
