/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTETHERMALDEFORMGRAD_H
#define COMPUTETHERMALDEFORMGRAD_H

#include "RankTwoTensor.h"
#include "ComputeVolumetricDeformGrad.h"

class ComputeThermalDeformGrad;

template <>
InputParameters validParams<ComputeThermalDeformGrad>();


/**
 * ComputeThermalDeformGrad is the class to compute volumetric deformation gradient
 * Modification based on pre-multiplication to a deformation gradient
 * Can be used to form a chain of volumetric corections on deformation
 */
class ComputeThermalDeformGrad : public ComputeVolumetricDeformGrad
{
public:
  ComputeThermalDeformGrad(const InputParameters & parameters);

protected:
  virtual void createVolumetricDeformGrad();
  
  const bool _has_temp;
  const VariableValue & _temperature;
  const Real _alpha;
  bool _has_stress_free_temp;
  Real _stress_free_temp;
};

#endif //COMPUTETHERMALDEFORMGARD_H
