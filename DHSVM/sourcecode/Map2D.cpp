/*
 * SUMMARY:      Map2D.cpp
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    October 2018
 * DESCRIPTION:  Routines for input of 2D map data
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-10-23 13:50:17 d3g096
 * COMMENTS:
 */

#include <memory>

extern "C" {
#include "Map2D.h"
#include "ParallelDHSVM.h"
#include "DHSVMerror.h"
}

#include "BinaryInputMap2D.hpp"

#ifdef HAVE_NETCDF
#include "NetCDFInputMap2D.hpp"
#endif

// -------------------------------------------------------------
// functor that acts like a InputMap2D factory
// -------------------------------------------------------------
// -------------------------------------------------------------
//  class InputMap2DFactory
// -------------------------------------------------------------
class InputMap2DFactory {
protected:

  /// The format we are using
  int format;

  /// Protected copy constructor to avoid unwanted copies.
  InputMap2DFactory(const InputMap2DFactory& old);

public:

  /// Default constructor.
  InputMap2DFactory(const int& theformat)
    : format(theformat)
  {}

  /// Destructor
  ~InputMap2DFactory(void)
  {}

  InputMap2D *operator() (const std::string& fname, const std::string& vname,
                          const int& NumberType, const MAPSIZE *Map,
                          const bool& mirror);
};


// -------------------------------------------------------------
// InputMap2DFactory::operator()
// -------------------------------------------------------------
InputMap2D *
InputMap2DFactory::operator() (const std::string& fname, const std::string& vname,
                               const int& NumberType, const MAPSIZE *Map,
                               const bool& mirror)
{
  InputMap2D *result(NULL);
  switch (format) {
  case (BIN):
    result = new BinaryInputMap2D(fname, vname, NumberType, Map, mirror);
    break;
  case (BYTESWAP):
    result = new ByteSwapInputMap2d(fname, vname, NumberType, Map, mirror);
    break;
#ifdef HAVE_NETCDF
  case (NETCDF):
    result = new NetCDFInputMap2D(fname, vname, NumberType, Map, mirror);
    break;
#endif
  default:
    
    break;
  }
  return result;
}

static std::auto_ptr<InputMap2DFactory> input_factory(new InputMap2DFactory(BIN));

// -------------------------------------------------------------
// Map2DInit
// -------------------------------------------------------------
extern "C"
void
Map2DInit(int FileFormat)
{
  const char Routine[] = "Map2DInit";
  switch (FileFormat) {
  case (BIN):
  case (BYTESWAP):
    break;
  case (NETCDF):
#ifndef HAVE_NETCDF
    ReportError((char *) Routine, 56);
#endif
    break;
  default:
    ReportError((char *) Routine, 38);
  }
  input_factory.reset(new InputMap2DFactory(FileFormat));
}

// -------------------------------------------------------------
// Read2DMatrix
// -------------------------------------------------------------
extern "C"
int
Read2DMatrix(const char *FileName, void *Matrix, int NumberType, 
             MAPSIZE *Map, int NDataSet, const char *VarName, int index)
{
  int flag;
  bool mirror(false);
  std::auto_ptr<InputMap2D>
    map((*input_factory)(FileName, VarName, NumberType, Map, mirror));
  map->open();
  flag = map->read(NDataSet, index, Matrix);
  map->close();
  return flag;
}

// -------------------------------------------------------------
// Read2DMatrixAll
// -------------------------------------------------------------
extern "C"
int Read2DMatrixAll(const char *FileName, void *Matrix, int NumberType, 
                    MAPSIZE *Map, int NDataSet, const char *VarName, int index)
{
  int flag;
  std::auto_ptr<InputMap2D>
    map((*input_factory)(FileName, VarName, NumberType, Map, true));
  map->open();
  flag = map->read(NDataSet, index, Matrix);
  map->close();
  return flag;
}


// -------------------------------------------------------------
// InputMap2DAlloc
// -------------------------------------------------------------
extern "C"
void *
InputMap2DAlloc(const char* fname, const char* vname,
                int NumberType, MAPSIZE *Map, int mirror)
{
  InputMap2D *result;
  result = (*input_factory)(fname, vname, NumberType, Map, mirror);
  return static_cast<void *>(result);
}

// -------------------------------------------------------------
// InputMap2DOpen
// -------------------------------------------------------------
int
InputMap2DOpen(void *map2d)
{
  static_cast<InputMap2D *>(map2d)->open();
  return 0;
}

// -------------------------------------------------------------
// InputMap2DRead
// -------------------------------------------------------------
int
InputMap2DRead(void *map2d, int NDataSet, int index, void *ldata)
{
  return static_cast<InputMap2D *>(map2d)->read(NDataSet, index, ldata);
}

// -------------------------------------------------------------
// InputMap2DClose
// -------------------------------------------------------------
int
InputMap2DClose(void *map2d)
{
  static_cast<InputMap2D *>(map2d)->close();
  return 0;
}

// -------------------------------------------------------------
// InputMap2DFree
// -------------------------------------------------------------
extern "C"
void
InputMap2DFree(void *map2d)
{
  delete static_cast<InputMap2D *>(map2d);
}

