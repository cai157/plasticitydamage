//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DRankFourTensor.h"

// MOOSE includes
#include "RankTwoTensor.h"
#include "MooseEnum.h"
#include "MooseException.h"
#include "MooseUtils.h"
#include "MatrixTools.h"
#include "MaterialProperty.h"
#include "PermutationTensor.h"

#include "libmesh/utility.h"

// C++ includes
#include <iomanip>
#include <ostream>

void
DRankFourTensor::Drotate(const RealTensorValue & R, const RealTensorValue & dRdphi)
{
  DRankFourTensor old = *this;
    
  unsigned int index = 0;
  for (unsigned int i = 0; i < N; ++i)
  {
    for (unsigned int j = 0; j < N; ++j)
    {
      for (unsigned int k = 0; k < N; ++k)
      {
        for (unsigned int l = 0; l < N; ++l)
        {
            
          Real sum = 0.0;
          unsigned int index2 = 0;
          
          for (unsigned int m = 0; m < N; ++m)
          {
            const Real a = dRdphi(i, m);
            for (unsigned int n = 0; n < N; ++n)
            {
              const Real ab = a * R(j, n);
              for (unsigned int o = 0; o < N; ++o)
              {
                const Real abc = ab * R(k, o);
                for (unsigned int p = 0; p < N; ++p)
                sum += abc * R(l, p) * old._vals[index2++];
              }
             }
            }
            
            index2 = 0;
            for (unsigned int m = 0; m < N; ++m)
            {
                const Real a = R(i, m);
                for (unsigned int n = 0; n < N; ++n)
                {
                    const Real ab = a * dRdphi(j, n);
                    for (unsigned int o = 0; o < N; ++o)
                    {
                        const Real abc = ab * R(k, o);
                        for (unsigned int p = 0; p < N; ++p)
                            sum += abc * R(l, p) * old._vals[index2++];
                    }
                }
            }

            index2 = 0;
            for (unsigned int m = 0; m < N; ++m)
            {
                const Real a = R(i, m);
                for (unsigned int n = 0; n < N; ++n)
                {
                    const Real ab = a * R(j, n);
                    for (unsigned int o = 0; o < N; ++o)
                    {
                        const Real abc = ab * dRdphi(k, o);
                        for (unsigned int p = 0; p < N; ++p)
                            sum += abc * R(l, p) * old._vals[index2++];
                    }
                }
            }
            

            index2 = 0;
            for (unsigned int m = 0; m < N; ++m)
            {
                const Real a = R(i, m);
                for (unsigned int n = 0; n < N; ++n)
                {
                    const Real ab = a * R(j, n);
                    for (unsigned int o = 0; o < N; ++o)
                    {
                        const Real abc = ab * R(k, o);
                        for (unsigned int p = 0; p < N; ++p)
                            sum += abc * dRdphi(l, p) * old._vals[index2++];
                    }
                }
            }
            

        
            _vals[index++] = sum; // 4 pieces
          }
        }
      }
    }
 }
    
