// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: NetCDFInputMap2D.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 18, 2018 by William A. Perkins
// Last Change: 2018-11-30 12:02:43 d3g096
// -------------------------------------------------------------


#ifndef _NetCDFInputMap2D_hpp_
#define _NetCDFInputMap2D_hpp_

#include "InputMap2D.hpp"

// -------------------------------------------------------------
//  class NetCDFInputMap2D
// -------------------------------------------------------------
class NetCDFInputMap2D
  : public SerialInputMap2D
{
protected:

  /// An error routine
  static void nc_check_err(const int& ncstatus, const int& line, 
			   const char *sfile, const char *dfile);

  /// NetCDF file handle
  int my_ncid;

  /// NetCDF variable handle
  int my_varid;

  /// Number of dimensions for this variable (should be 3)
  int my_ndims;

  /// The dimension ids for this variable
  int my_dimids[3];

  /// Is the y-dimension flipped?
  int my_flip;

  /// Open the input map file (specialized)
  void my_open();

  /// Performe necessary checks on the file and variable
  int my_check();

  /// Close the input map file (specialized)
  void my_close();

  /// Fill array index ranges
  virtual void my_indexes(size_t start[], size_t count[], const int& index);

  /// format specific read
  int my_read_fmt(const int& unused_index, const int& index, unsigned char *buffer);

public:

  /// Default constructor.
  NetCDFInputMap2D(const std::string fname, const std::string vname,
                   const int NumberType, const MAPSIZE *Map, const bool mirror);

  /// Destructor
  ~NetCDFInputMap2D(void);
};

#endif
