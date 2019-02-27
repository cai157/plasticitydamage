//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DRANKFOURTENSOR_H
#define DRANKFOURTENSOR_H

// MOOSE includes
#include "DataIO.h"

#include "libmesh/tensor_value.h"
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h"
#include "RankFourTensor.h"

// Forward declarations
class MooseEnum;
class RankTwoTensor;
class RankFourTensor;
class DRankFourTensor;
class DRankTwoTensor;

class DRankFourTensor:public RankFourTensor
{
public:

  void Drotate(const RealTensorValue & R, const RealTensorValue &);
    
protected:
    
  template <class T>
  friend void dataStore(std::ostream &, T &, void *);
    
  template <class T>
  friend void dataLoad(std::istream &, T &, void *);
    
  friend class RankTwoTensor;
  friend class DRankTwoTensor;
  friend class RankThreeTensor;

};

#endif // RANKFOURTENSOR_H
