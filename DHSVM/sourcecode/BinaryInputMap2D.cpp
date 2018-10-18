// -------------------------------------------------------------
// file: BinaryInputMap2D.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created October 17, 2018 by William A. Perkins
// Last Change: 2018-10-17 15:01:23 d3g096
// -------------------------------------------------------------


extern "C" {
#include "sizeofnt.h"
#include "ParallelDHSVM.h"
#include "DHSVMerror.h"
#include "byte_swap.h"
}

#include "BinaryInputMap2D.hpp"


// I don't want to include all of fileio.h
extern "C" 
void OpenFile(FILE **FilePtr, const char *FileName, const char *Mode,
	      unsigned char OverWrite);


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
  // empty
}

// -------------------------------------------------------------
// BinaryInputMap2D:my_open
// -------------------------------------------------------------
void
BinaryInputMap2D::my_open()
{
  OpenFile(&(my_fd), my_Name.c_str(), "r", 0);
}

// -------------------------------------------------------------
// BinaryInputMap2D:my_close
// -------------------------------------------------------------
void
BinaryInputMap2D::my_close()
{
  fclose(my_fd);
  my_fd = NULL;
}

// -------------------------------------------------------------
// BinaryInputMap2D::my_read_fmt
// -------------------------------------------------------------
int
BinaryInputMap2D::my_read_fmt(const int& index, void *LocalMatrix)
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
  if (fseek(my_fd, offset, origin)) ReportError(this->my_Name.data(), 39);

  int NElements(fread(LocalMatrix, ElemSize, NY * NX, this->my_fd));

  if (NElements != NY * NX) ReportError(this->my_Name.data(), 2);
                
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
ByteSwapInputMap2d::my_read_fmt(const int& index, void *LocalMatrix)
{
  int NElements =
    BinaryInputMap2D::my_read_fmt(index, LocalMatrix);
  
  size_t ElemSize(SizeOfNumberType(this->my_NumberType));

  if (ElemSize == 4) {
    byte_swap_long((long int *)LocalMatrix, NElements);
  }
  else if (ElemSize == 2) {
    byte_swap_short((short *)LocalMatrix, NElements);
  }
  else if (ElemSize != 1) {
    ReportError(this->my_Name.c_str(), 61);
  }
}
