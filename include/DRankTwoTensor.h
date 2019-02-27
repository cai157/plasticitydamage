//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DRANKTWOTENSOR_H
#define DRANKTWOTENSOR_H

#include "Moose.h"

#include "RankFourTensor.h"

// Any requisite includes here
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "RankTwoTensor.h"

#include <vector>
#include "MooseRandom.h"

// Forward declarations
class RankTwoTensor;
class RankFourTensor;
class DRankFourTensor;

class DRankTwoTensor:public RankTwoTensor
{

private:
    
  template <class T>
  friend void dataStore(std::ostream &, T &, void *);

  template <class T>
  friend void dataLoad(std::istream &, T &, void *);
  friend class RankFourTensor;
  friend class RankThreeTensor;
  friend class DRankFourTensor;
};

#endif // DRANKTWOTENSOR_H
