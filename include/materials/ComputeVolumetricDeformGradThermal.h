//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEVOLUMETRICDEFORMGRADTHERMAL_H
#define COMPUTEVOLUMETRICDEFORMGRADTHERMAL_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"

class ComputeVolumetricDeformGradThermal;

template <>
InputParameters validParams<ComputeVolumetricDeformGradThermal>();

/**
 * ComputeVolumetricDeformGrad is the class to compute volumetric deformation gradient
 * Modification based on pre-multiplication to a deformation gradient
 * Can be used to form a chain of volumetric corections on deformation
 */
class ComputeVolumetricDeformGradThermal : public DerivativeMaterialInterface<Material>
{
public:
  ComputeVolumetricDeformGradThermal(const InputParameters & parameters);

protected:
//  virtual void initQpStatefulProperties();
 // virtual void computeQpProperties();
 // virtual void createVolumetricDeformGrad();

//  const MaterialProperty<RankTwoTensor> & _pre_deform_grad;
 // MaterialProperty<RankTwoTensor> & _volumetric_deform_grad;
//  MaterialProperty<RankTwoTensor> & _post_deform_grad;
};

#endif // COMPUTEVOLUMETRICDEFORMGARDTHERMAL_H
