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
// Last Change: 2018-10-18 08:37:39 d3g096
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

  /// format specific read
  int my_read_fmt(const int& index, unsigned char *buffer);

public:

  /// Default constructor.
  NetCDFInputMap2D(const std::string fname, const std::string vname,
                   const int NumberType, const MAPSIZE *Map, const bool mirror);

  /// Destructor
  ~NetCDFInputMap2D(void);
};

#endif
