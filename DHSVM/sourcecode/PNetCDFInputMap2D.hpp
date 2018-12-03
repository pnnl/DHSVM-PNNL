// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: PNetCDFInputMap2D.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 30, 2018 by William A. Perkins
// Last Change: 2018-11-30 15:40:45 d3g096
// -------------------------------------------------------------


#ifndef _PNetCDFInputMap2D_hpp_
#define _PNetCDFInputMap2D_hpp_

#include "NetCDFInputMap2D.hpp"

// -------------------------------------------------------------
//  class PNetCDFInputMap2D
// -------------------------------------------------------------
class PNetCDFInputMap2D
  : public NetCDFInputMap2D
{
protected:

  /// Open the input map file (specialized)
  void my_open();

  /// Close the input map file (specialized)
  void my_close();
  
  /// Fill array index ranges
  void my_indexes(size_t start[], size_t count[], const int& index);

  /// Read a map into the memory specified (specialized)
  int my_read(const int& NDataSet, const int& index, void *LocalMatrix);


public:

  /// Default constructor.
  PNetCDFInputMap2D(const std::string fname, const std::string vname,
                    const int NumberType, const MAPSIZE *Map, const bool mirror);

  /// Destructor
  ~PNetCDFInputMap2D(void);
};



#endif
