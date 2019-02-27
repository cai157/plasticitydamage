/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTETHERMALDEFORMGRADANISOTROPIC_H
#define COMPUTETHERMALDEFORMGRADANISOTROPIC_H

#include "RankTwoTensor.h"
#include "ComputeVolumetricDeformGradThermal.h"

#include "ElementPropertyReadFile.h"
#include "RotationTensor.h"

class ComputeThermalDeformGradAnisotropic;

template <>
InputParameters validParams<ComputeThermalDeformGradAnisotropic>();


/**
 * ComputeThermalDeformGrad is the class to compute volumetric deformation gradient
 * Modification based on pre-multiplication to a deformation gradient
 * Can be used to form a chain of volumetric corections on deformation
 */
class ComputeThermalDeformGradAnisotropic :
public ComputeVolumetricDeformGradThermal
{
public:
  ComputeThermalDeformGradAnisotropic(const InputParameters & parameters);

protected:
  virtual void createThermalDeformGrad();
  
  const bool _has_temp;
  const VariableValue & _temperature;
  Real _stress_free_temp;
//  const std::vector<Real> & _thermal_expansion_coeff;
  const Real _thermal_expansion_coeff1;
  const Real _thermal_expansion_coeff2;
  const Real _thermal_expansion_coeff3;

  bool _has_stress_free_temp;
  MaterialProperty<RankTwoTensor> & _thermal_deform_grad;
  const ElementPropertyReadFile * _read_prop_user_object_eigenstrain;
  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;
 };

#endif //COMPUTETHERMALDEFORMGARDANISOTROPIC_H
