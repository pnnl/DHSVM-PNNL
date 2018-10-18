// -------------------------------------------------------------
// file: NetCDFInputMap2D.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 18, 2018 by William A. Perkins
// Last Change: 2018-10-18 08:44:16 d3g096
// -------------------------------------------------------------

#include <cstdio>
#include <vector>
#include <netcdf.h>

extern "C" {
#include "DHSVMerror.h"
}

#include "NetCDFInputMap2D.hpp"

// -------------------------------------------------------------
// nc_check_err
// -------------------------------------------------------------
static void
nc_check_err(const int ncstatus, const int line, const char *file)
{
  char str[BUFSIZE + 1];

  if (ncstatus != NC_NOERR) {
    sprintf(str, "%s, line: %d -- %s", file, line, nc_strerror(ncstatus));
    ReportError((char *) str, 57);
  }
}

// -------------------------------------------------------------
//  class NetCDFInputMap2D
// -------------------------------------------------------------

// -------------------------------------------------------------
// NetCDFInputMap2D:: constructors / destructor
// -------------------------------------------------------------
NetCDFInputMap2D::NetCDFInputMap2D(const std::string fname, const std::string vname,
                                   const int NumberType, const MAPSIZE *Map,
                                   const bool mirror)
  : SerialInputMap2D(fname, vname, NumberType, Map, mirror),
    my_ncid(0), my_varid(0), my_ndims(0), my_flip(0)
{
  // empty
}

NetCDFInputMap2D::~NetCDFInputMap2D(void)
{
  // empty
}

// -------------------------------------------------------------
// NetCDFInputMap2D::my_open
// -------------------------------------------------------------
void
NetCDFInputMap2D::my_open(void)
{
  int ncstatus;
  int TempNumberType;
  char msg[BUFSIZE + 1];
  
  ncstatus = nc_open(my_Name.c_str(), NC_NOWRITE, &my_ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  /// check whether the variable exists and get its parameters

  ncstatus = nc_inq_varid(my_ncid, my_VarName.c_str(), &my_varid);
  nc_check_err(ncstatus, __LINE__, __FILE__);

  ncstatus = nc_inq_var(my_ncid, my_varid, 0, &TempNumberType, &my_ndims, my_dimids, NULL);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (TempNumberType != my_NumberType) {
    sprintf(msg, "%s: nc_type for %s is different than expected.\n",
	    my_Name.c_str(), my_VarName.c_str());
    ReportWarning(msg, 58);
  }
  my_flip = my_check();
}


// -------------------------------------------------------------
// NetCDFInputMap2D::my_check
// make sure that the x and y dimensions have the correct sizes 
// -------------------------------------------------------------
int
NetCDFInputMap2D::my_check(void)
{
  int ncstatus;
  int flag(0);
  size_t dimlen;
  char dimname[NC_MAX_NAME + 1];
  int lat_varid, lon_varid;
  bool LatisAsc, LonisAsc;
  
  ncstatus = nc_inq_dim(my_ncid, my_dimids[1], dimname, &dimlen);  
  nc_check_err(ncstatus, __LINE__, __FILE__);
  
  ncstatus = nc_inq_varid(my_ncid, dimname, &lat_varid);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (dimlen != my_Map->gNY)
    ReportError(my_VarName.c_str(), 59);

  std::vector<double> Ycoord(dimlen);

  // Read the latitude coordinate variable data. 
  ncstatus = nc_get_var_double(my_ncid, lat_varid, &Ycoord[0]);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  
  /* A quick check if the lat, long are in a ascending order. 
     If so, matrix must be flipped so the first value in the matrix will be 
     assigned to the lower left corner cell that has lowest X (lon) & Y (lat) value. 
     (see more comments in the header of this C file). */

  LatisAsc = 1;
  if( Ycoord[0] > Ycoord[my_Map->gNY - 1] ) 
    LatisAsc = 0;

  ncstatus = nc_inq_dim(my_ncid, my_dimids[2], dimname, &dimlen);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  ncstatus = nc_inq_varid(my_ncid, dimname, &lon_varid);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  if (dimlen != my_Map->gNX)
    ReportError(my_VarName.c_str(), 60);

  std::vector<double> Xcoord(dimlen);

  /* Read the latitude coordinate variable data. */
  ncstatus = nc_get_var_double(my_ncid, lon_varid, &Xcoord[0]);
  nc_check_err(ncstatus, __LINE__, __FILE__);
  LonisAsc = true;
  if( Xcoord[0] > Xcoord[my_Map->gNX - 1] ) 
    LonisAsc = false;

  if (!LonisAsc){
    printf("The current program does not handle the cases when longitude or X \
values in the .nc input in an descending order. You can either change the input \
.nc file format outside of this program. or you can easily modify this program to \
fit your needs. \n");
    ReportError("Improper NetCDF input files", 58);
  }
  if (!LatisAsc && LonisAsc) flag = 0;
  if (LatisAsc && LonisAsc) flag = 1;

  return (flag);
}

// -------------------------------------------------------------
// NetCDFInputMap2D::my_close
// -------------------------------------------------------------
void
NetCDFInputMap2D::my_close(void)
{
  int ncstatus = nc_close(my_ncid);
  nc_check_err(ncstatus, __LINE__, __FILE__);
}

// -------------------------------------------------------------
// NetCDFInputMap2D::my_read_fmt
// -------------------------------------------------------------
int
NetCDFInputMap2D::my_read_fmt(const int& index, unsigned char *buffer)
{
  int ncstatus;
  size_t count[3];
  size_t start[3] = { index, 0, 0 };
  size_t timelen;

  count[0] = 1;
  count[1] = my_Map->gNY;
  count[2] = my_Map->gNX;

  /* see whether the time dimension needs to be updated (the assumption is that
     the same index value refers to the same moment in time.  Since currently we
     make separate files for separate variables this is OK) */
  ncstatus = nc_inq_dimlen(my_ncid, my_dimids[0], &timelen);
  nc_check_err(ncstatus, __LINE__, __FILE__);


  /****************************************************************************/
  /*                             READ VARIABLE                                */
  /****************************************************************************/

  switch (my_NumberType) {
  case NC_BYTE:
    ncstatus = nc_get_vara_uchar(my_ncid, my_varid, start, count, buffer);
    break;
  case NC_CHAR:
    ncstatus = nc_get_vara_text(my_ncid, my_varid, start, count, (char *)buffer);
    break;
  case NC_SHORT:
    ncstatus = nc_get_vara_short(my_ncid, my_varid, start, count, (short *)buffer);
    break;
  case NC_INT:
    ncstatus = nc_get_vara_int(my_ncid, my_varid, start, count, (int *)buffer);
    break;
    /* 8 bit integer not yet implemented in NetCDF 3.4, but anticipated in
       future versions */
    /*   case NC_LONG: */
    /*     ncstatus = nc_put_vara_long(my_ncid, my_varid, start, count, (void *)buffer); */
    /*     break; */
  case NC_FLOAT:
    ncstatus = nc_get_vara_float(my_ncid, my_varid, start, count, (float *)buffer);
    break;
  case NC_DOUBLE:
    ncstatus = nc_get_vara_double(my_ncid, my_varid, start, count, (double *)buffer);
    break;
  default:
    ReportError("NetCDFInputMap2D::my_read_fmt", 40);
    break;
  }
  nc_check_err(ncstatus, __LINE__, __FILE__);

  return (my_flip);
}


