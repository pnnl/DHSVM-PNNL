// -------------------------------------------------------------
// file: PNetCDFInputMap2D.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 30, 2018 by William A. Perkins
// Last Change: 2018-12-03 14:15:27 d3g096
// -------------------------------------------------------------

#include <iostream>
#include <netcdf.h>
#include <netcdf_par.h>
#include "ga_helper.h"
#include <ga-mpi.h>
#include "PNetCDFInputMap2D.hpp"

// -------------------------------------------------------------
//  class PNetCDFInputMap2D
// -------------------------------------------------------------

// -------------------------------------------------------------
// PNetCDFInputMap2D:: constructors / destructor
// -------------------------------------------------------------
PNetCDFInputMap2D::PNetCDFInputMap2D(const std::string fname, const std::string vname,
                                     const int NumberType, const MAPSIZE *Map,
                                     const bool mirror)
  : NetCDFInputMap2D(fname, vname, NumberType, Map, mirror)
{
  
}

PNetCDFInputMap2D::~PNetCDFInputMap2D(void)
{
}

// -------------------------------------------------------------
// PNetCDFInputMap2D::my_open
// -------------------------------------------------------------
void
PNetCDFInputMap2D::my_open(void)
{
  int me(ParallelRank());
  int ncstatus;
  int TempNumberType;
  char msg[BUFSIZE + 1];
  MPI_Comm comm(GA_MPI_Comm());
  MPI_Info info = MPI_INFO_NULL;
  int ierr(0);

  try {
    int themode(NC_NOWRITE | NC_MPIIO);
    ncstatus = nc_open_par(my_Name.c_str(), themode, comm, info, &my_ncid);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    /// check whether the variable exists and get its parameters

    ncstatus = nc_inq_varid(my_ncid, my_VarName.c_str(), &my_varid);
    nc_check_err(ncstatus, __LINE__, __FILE__);

    ncstatus = nc_inq_var(my_ncid, my_varid, 0, &TempNumberType, &my_ndims, my_dimids, NULL);
    nc_check_err(ncstatus, __LINE__, __FILE__);
    if (TempNumberType != my_NumberType) {
      sprintf(msg, "%s: nc_type for %s is different than expected.\n",
              my_Name.c_str(), my_VarName.c_str());
      std::cerr << msg << std::endl;
      // throw InputMap2D::exception(msg, 58);
    }

    // ncstatus = nc_var_par_access(my_ncid, my_varid, NC_INDEPENDENT);
    ncstatus = nc_var_par_access(my_ncid, my_varid, NC_COLLECTIVE);
    nc_check_err(ncstatus, __LINE__, __FILE__);

    my_flip = my_check();

    std::cerr << my_Name << ": " << my_VarName
              << ": successfully opened in parallel" << std::endl;
    
  } catch (const InputMap2D::exception& e) {
    std::cerr << me << ": " << my_Name << ": error: " << e.what() << std::endl;
    ierr++;
  } catch (...) {
    std::cerr << me << ": " << my_Name << ": unknown error" << std::endl;
  }
  GA_Igop(&ierr, 1, "+");
  if (ierr) {
    std::string msg(my_Name);
    msg += ": error: cannot open";
    throw InputMap2D::exception(msg, 3);
  }
}

// -------------------------------------------------------------
// PNetCDFInputMap2D::my_close
// -------------------------------------------------------------
void
PNetCDFInputMap2D::my_close(void)
{
  int ncstatus = nc_close(my_ncid);
  // Let's not worry about this
  // nc_check_err(ncstatus, __LINE__, __FILE__);
}

// -------------------------------------------------------------
// PNetCDFInputMap2D::my_indexes
// -------------------------------------------------------------
void
PNetCDFInputMap2D::my_indexes(size_t start[], size_t count[], const int& index)
{
  if (my_mirror) {
    start[0] = index;
    start[1] = 0;
    start[2] = 0;
    
    count[0] = 1;
    count[1] = my_Map->gNY;
    count[2] = my_Map->gNX;
  } else {
    start[0] = index;
    start[1] = my_Map->OffsetY;
    start[2] = my_Map->OffsetX;
    
    count[0] = 1;
    count[1] = my_Map->NY;
    count[2] = my_Map->NX;
  }
}


// -------------------------------------------------------------
// PNetCDFInputMap2D::my_read
// -------------------------------------------------------------
int
PNetCDFInputMap2D::my_read(const int& NDataSet, const int& index, void *LocalMatrix)
{
  const char Routine[] = "PNetCDFInputMap2D::my_read";
  int flag = my_read_fmt(NDataSet, index, static_cast<unsigned char *>(LocalMatrix));
  return flag;
}

  


