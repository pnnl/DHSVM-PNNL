// -------------------------------------------------------------
// file: BinaryInputMap2D.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 17, 2018 by William A. Perkins
// Last Change: 2018-11-21 10:54:23 d3g096
// -------------------------------------------------------------


#include "BinaryInputMap2D.hpp"
#include "ParallelDHSVM.h"
#include "sizeofnt.h"
#include "ga_helper.h"
#include "byte_swap.h"


// -------------------------------------------------------------
//  class BinaryInputMap2D
// -------------------------------------------------------------

// -------------------------------------------------------------
// BinaryInputMap2D:: constructors / destructor
// -------------------------------------------------------------
BinaryInputMap2D::BinaryInputMap2D(const std::string fname, const std::string vname,
                                   const int NumberType, const MAPSIZE *Map,
                                   const bool mirror)
  : SerialInputMap2D(fname, vname, NumberType, Map, mirror),
    my_fd(NULL)
{
  // empty
}

BinaryInputMap2D::~BinaryInputMap2D(void)
{
  this->close();
}

// -------------------------------------------------------------
// BinaryInputMap2D:my_open
// -------------------------------------------------------------
void
BinaryInputMap2D::my_open()
{
  int me(ParallelRank());
  if (me == 0) {
    my_fd = fopen(my_Name.c_str(), "r");
    if (!my_fd) {
      std::string msg(my_Name);
      msg += ": cannot open";
      throw InputMap2D::exception(msg, 3);
    }
  } else {
    my_fd = NULL;
  }
}

// -------------------------------------------------------------
// BinaryInputMap2D:my_close
// -------------------------------------------------------------
void
BinaryInputMap2D::my_close()
{
  int me(ParallelRank());
  if (me == 0) {
    if (my_fd != NULL) {
      fclose(my_fd);
    }
  }
  my_fd = NULL;
}

// -------------------------------------------------------------
// BinaryInputMap2D::my_read_fmt
// -------------------------------------------------------------
int
BinaryInputMap2D::my_read_fmt(const int& index, const int& unused_index, unsigned char *LocalMatrix)
{
  int NX(this->my_Map->gNX);
  int NY(this->my_Map->gNY);
  size_t ElemSize(SizeOfNumberType(this->my_NumberType));
  int origin, offset;

  if (my_last_index < 0 || my_last_index > index) {
    origin = SEEK_SET;
    offset = NX*NY*ElemSize*index;
  } else {
    origin = SEEK_CUR;
    offset = NX*NY*ElemSize*(index - my_last_index - 1);
  }
  
  if (fseek(my_fd, offset, origin)) {
    std::string msg(my_Name);
    msg += ": fseek error";
    throw InputMap2D::exception(msg, 39);
  }
  
  int NElements(fread(LocalMatrix, ElemSize, NY * NX, this->my_fd));

  if (NElements != NY * NX) {
    std::string msg(my_Name);
    msg += ": fread returned wrong size";
    throw InputMap2D::exception(msg, 2);
  }
                
  return NElements;
}


// -------------------------------------------------------------
//  class ByteSwapInputMap2d
// -------------------------------------------------------------

// -------------------------------------------------------------
// ByteSwapInputMap2d:: constructors / destructor
// -------------------------------------------------------------
ByteSwapInputMap2d::ByteSwapInputMap2d(const std::string fname, const std::string vname,
                     const int NumberType, const MAPSIZE *Map, const bool mirror)
  : BinaryInputMap2D(fname, vname, NumberType, Map, mirror)
{}


ByteSwapInputMap2d::~ByteSwapInputMap2d(void)
{
}

// -------------------------------------------------------------
// ByteSwapInputMap2d::my_read_fmt
// -------------------------------------------------------------
int
ByteSwapInputMap2d::my_read_fmt(const int& index, const int& unused_index, unsigned char *LocalMatrix)
{
  int NElements =
    BinaryInputMap2D::my_read_fmt(index, unused_index, LocalMatrix);
  
  size_t ElemSize(SizeOfNumberType(this->my_NumberType));

  if (ElemSize == 4) {
    byte_swap_long((long int *)LocalMatrix, NElements);
  }
  else if (ElemSize == 2) {
    byte_swap_short((short *)LocalMatrix, NElements);
  }
  else if (ElemSize != 1) {
    std::string msg(my_Name);
    msg += ": unknown element size";
    throw InputMap2D::exception(msg, 61);
  }
}
