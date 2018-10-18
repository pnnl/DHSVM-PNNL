// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: InputMap2D.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 15, 2018 by William A. Perkins
// Last Change: 2018-10-18 06:23:32 d3g096
// -------------------------------------------------------------


#ifndef _InputMap2D_hpp_
#define _InputMap2D_hpp_

#include <string>

#include "MapSize.h"

// -------------------------------------------------------------
//  class InputMap2D
// -------------------------------------------------------------
class InputMap2D {
protected:

  /// The input file name
  const std::string my_Name;

  /// The name of the variable
  const std::string my_VarName;

  /// The number type
  const int my_NumberType;

  /// Describes the @em local map extent
  const MAPSIZE *my_Map;

  /// Is the input map set mirrored on all processes?
  const int my_mirror;

  /// The last index read
  int my_last_index;

  /// Open the input map file (specialized)
  virtual void my_open() = 0;

  /// Close the input map file (specialized)
  virtual void my_close()  = 0;

  /// Read a map into the memory specified (specialized)
  virtual int my_read(const int& index, void *LocalMatrix) = 0;

  /// Undefined protected copy constructor to avoid unwanted copies.
  InputMap2D(const InputMap2D& old);

public:

  /// Default constructor.
  InputMap2D(const std::string fname, const std::string vname, const int NumberType,
             const MAPSIZE *Map, const bool mirror)
    : my_Name(fname), my_VarName(vname), my_NumberType(NumberType),
      my_Map(Map), my_mirror(mirror), my_last_index(-1)
  {}

  /// Destructor
  ~InputMap2D(void)
  {
    try {
      this->close();
    } catch (...) {
      // destructors cannot throw
    }
  }

  /// Open the input map file
  void open()
  {
    this->my_open();
  }

  /// Close the input map file
  void close()
  {
    this->my_close();
  }

  /// Read a map into the memory specified (if 1 returned flip y)
  int read(const int& index, void *LocalMatrix)
  {
    int flag(this->my_read(index, LocalMatrix));
    my_last_index = index;
  }
};

// -------------------------------------------------------------
//  class SerialInputMap2D
// -------------------------------------------------------------
class SerialInputMap2D
  : public InputMap2D
{
protected:

  /// format specific read
  virtual int my_read_fmt(const int& index, unsigned char *buffer) = 0;

  /// Read a map into the memory specified (specialized)
  int my_read(const int& index, void *LocalMatrix);

  /// Routine to distribute the map once read
  void my_distribute(unsigned char *buf0, void *LocalMatrix);

public:

  /// Default constructor.
  SerialInputMap2D(const std::string fname, const std::string vname,
                   const int NumberType, const MAPSIZE *Map, const bool mirror);

  /// Destructor
  ~SerialInputMap2D(void);
};




#endif
