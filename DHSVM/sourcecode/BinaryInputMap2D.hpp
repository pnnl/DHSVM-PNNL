// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: BinaryInputMap2D.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 17, 2018 by William A. Perkins
// Last Change: 2018-10-18 11:40:31 d3g096
// -------------------------------------------------------------


#ifndef _BinaryInputMap2D_hpp_
#define _BinaryInputMap2D_hpp_

#include <cstdio>

#include "InputMap2D.hpp"

// -------------------------------------------------------------
//  class BinaryInputMap2D
// -------------------------------------------------------------
class BinaryInputMap2D
  : public SerialInputMap2D
{
protected:

  /// The file descriptor for this map
  FILE *my_fd;

  /// Open the input map file (specialized)
  void my_open();

  /// Close the input map file (specialized)
  void my_close();

  /// format specific read
  int my_read_fmt(const int& index, unsigned char *LocalMatrix);

public:

  /// Default constructor.
  BinaryInputMap2D(const std::string fname, const std::string vname,
                   const int NumberType, const MAPSIZE *Map, const bool mirror);

  /// Destructor
  ~BinaryInputMap2D(void);
};

// -------------------------------------------------------------
//  class ByteSwapInputMap2d
// -------------------------------------------------------------
class ByteSwapInputMap2d
  : public BinaryInputMap2D
{
protected:

  /// format specific read
  int my_read_fmt(const int& index, unsigned char *LocalMatrix);

public:

  /// Default constructor.
  ByteSwapInputMap2d(const std::string fname, const std::string vname,
                     const int NumberType, const MAPSIZE *Map, const bool mirror);

  /// Destructor
  ~ByteSwapInputMap2d(void);
};




#endif
