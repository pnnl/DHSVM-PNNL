/*
 * SUMMARY:      InputMap2D.cpp
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    October 2018
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-11-01 06:39:43 d3g096
 * COMMENTS:
 */

#include <exception>

extern "C" {
#include "DHSVMerror.h"
#include "ParallelDHSVM.h"
}

#include "InputMap2D.hpp"

// -------------------------------------------------------------
//  class SerialInputMap2D
// -------------------------------------------------------------

// -------------------------------------------------------------
// SerialInputMap2D:: constructors / destructor
// -------------------------------------------------------------
SerialInputMap2D::SerialInputMap2D(const std::string fname, const std::string vname,
                                   const int NumberType, const MAPSIZE *Map, const bool mirror)
  : InputMap2D(fname, vname, NumberType, Map, mirror)
{
  int gatype(GA_Type(this->my_NumberType));
  my_ga = GA_Duplicate_type(this->my_Map->dist, "Distribute2DMatrix", gatype);
}

SerialInputMap2D::~SerialInputMap2D(void)
{
  GA_Destroy(my_ga);
}

// -------------------------------------------------------------
// SerialInputMap2D::my_distribute
// -------------------------------------------------------------
void
SerialInputMap2D::my_distribute(unsigned char *buf0, void *LocalMatrix)
{
  int me(ParallelRank());
  int gNX(this->my_Map->gNX), gNY(this->my_Map->gNY);
  int lo[2], hi[2], ld[2];

  if (me == 0) {
    lo[0] = 0;
    lo[1] = 0;
    hi[gaYdim] = gNY-1;
    hi[gaXdim] = gNX-1;
    ld[gaXdim] = gNY;
    ld[gaYdim] = gNX;
    NGA_Put(my_ga, &lo[0], &hi[0], (void *)buf0, &ld[0]);
  }
  ParallelBarrier();

  if (this->my_mirror) {
    lo[0] = 0;
    lo[1] = 0;
    hi[gaYdim] = gNY-1;
    hi[gaXdim] = gNX-1;
    ld[gaXdim] = gNY;
    ld[gaYdim] = gNX;
  } else {
    lo[gaYdim] = this->my_Map->OffsetY;
    lo[gaXdim] = this->my_Map->OffsetX;
    hi[gaYdim] = lo[gaYdim] + this->my_Map->NY - 1;
    hi[gaXdim] = lo[gaXdim] + this->my_Map->NX - 1;
    ld[gaXdim] = this->my_Map->NY;
    ld[gaYdim] = this->my_Map->NX;
  }
  NGA_Get(my_ga, &lo[0], &hi[0], LocalMatrix, &ld[0]);

  ParallelBarrier();
}

// -------------------------------------------------------------
// SerialInputMap2D::my_read
// -------------------------------------------------------------
int
SerialInputMap2D::my_read(const int& NDataSet, const int& index, void *LocalMatrix)
{
  const char Routine[] = "SerialInputMap2D::my_read";
  int me(ParallelRank());
  int gNX(this->my_Map->gNX);
  int gNY(this->my_Map->gNY);

  size_t bufsize(gNX*gNY*SizeOfNumberType(this->my_NumberType));
  unsigned char *tmpArray;

  int flag(0);

  if (me == 0) {
    tmpArray = new unsigned char[bufsize];
    flag = my_read_fmt(NDataSet, index, tmpArray);
  }
    
  my_distribute(tmpArray, LocalMatrix);
    
  if (me == 0) {
    delete[] tmpArray;
  }
    
  GA_Brdcst(&flag, sizeof(int), 0);
  return flag;
}
