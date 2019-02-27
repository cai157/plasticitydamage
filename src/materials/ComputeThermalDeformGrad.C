/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeThermalDeformGrad.h"

template<>
InputParameters validParams<ComputeThermalDeformGrad>()
{
  InputParameters params = validParams<ComputeVolumetricDeformGrad>();
  params.addClassDescription("Computes thermal deformation gradient and adjusts the total deformation gradient");
  params.addCoupledVar("temp", "Coupled Temperature");
  params.addParam<Real>("stress_free_temperature", "The stress-free temperature.  If not specified, the initial temperature is used.");
  params.addParam<Real>("thermal_expansion", "The thermal expansion coefficient.");
  return params;
}

ComputeThermalDeformGrad::ComputeThermalDeformGrad(const InputParameters & parameters) :
    ComputeVolumetricDeformGrad(parameters),
    _has_temp(isCoupled("temp")),
    _temperature(_has_temp ? coupledValue("temp") : _zero),
    _alpha(parameters.isParamValid("thermal_expansion") ? getParam<Real>("thermal_expansion") : 0.),
    _has_stress_free_temp(false),
    _stress_free_temp(0.0)
{
  if (parameters.isParamValid("stress_free_temperature"))
  {
    _has_stress_free_temp = true;
    _stress_free_temp = getParam<Real>("stress_free_temperature");
    if (!_has_temp)
      mooseError("Cannot specify stress_free_temperature without coupling to temperature");
  }
}

void
ComputeThermalDeformGrad::createVolumetricDeformGrad()
{
  _volumetric_deform_grad[_qp].zero();
  Real theta = std::exp(_alpha * (_temperature[_qp] - _stress_free_temp));
  _volumetric_deform_grad[_qp].addIa(1.0 * theta);
}
