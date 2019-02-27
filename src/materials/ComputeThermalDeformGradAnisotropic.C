/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeThermalDeformGradAnisotropic.h"

template<>
InputParameters validParams<ComputeThermalDeformGradAnisotropic>()
{
  InputParameters params = validParams<ComputeVolumetricDeformGradThermal>();

  params.addClassDescription("Computes thermal deformation gradient and adjusts the total deformation gradient");
  params.addCoupledVar("temp", "Coupled Temperature");
  params.addParam<Real>("stress_free_temp", 0.0 , "The stress-free temperature.  If not specified, the initial temperature is used.");
  params.addParam<UserObjectName>("read_prop_user_object_eigenstrain","The ElementReadPropertyFile GeneralUserObject to read element specific property values from file");
  //params.addParam<std::vector<Real>>("thermal_expansion_coeff", "thermal expansion coefficient");
  params.addParam<Real>("thermal_expansion_coeff1", "The thermal expansion coefficient1");
  params.addParam<Real>("thermal_expansion_coeff2", "The thermal expansion coefficient2");
  params.addParam<Real>("thermal_expansion_coeff3", "The thermal expansion coefficient3");

  return params;
}

ComputeThermalDeformGradAnisotropic::ComputeThermalDeformGradAnisotropic(const InputParameters & parameters)
   :ComputeVolumetricDeformGradThermal(parameters),
    _has_temp(isCoupled("temp")),
    _temperature(_has_temp ? coupledValue("temp") : _zero),
//  _alpha(parameters.isParamValid("thermal_expansion") ? getParam<std::vector<Real>>("thermal_expansion") : std::vector<Real>(3,0)),
//    _stress_free_temp(0.0),
    _stress_free_temp(getParam<Real>("stress_free_temp")),
//    _thermal_expansion_coeff(getParam<std::vector<Real>>("Thermal_expansion_coeff")),

    _thermal_expansion_coeff1(parameters.isParamValid("thermal_expansion_coeff1") ? getParam<Real>("thermal_expansion_coeff1") : 0.),
    _thermal_expansion_coeff2(parameters.isParamValid("thermal_expansion_coeff2") ? getParam<Real>("thermal_expansion_coeff2") : 0.),
    _thermal_expansion_coeff3(parameters.isParamValid("thermal_expansion_coeff3") ? getParam<Real>("thermal_expansion_coeff3") : 0.),

//    _thermal_expansion_coeff(getMaterialProperty<std::vector<Real>>("thermal_expansion_coeff")),
    _has_stress_free_temp(false),
    _thermal_deform_grad(declareProperty<RankTwoTensor>("thermal_deform_grad")),
    _read_prop_user_object_eigenstrain(isParamValid("read_prop_user_object_eigenstrain") ? & getUserObject<ElementPropertyReadFile>("read_prop_user_object_eigenstrain") : NULL),
    _Euler_angles_mat_prop(declareProperty<RealVectorValue>("Euler_angles_eigenstrain"))


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
ComputeThermalDeformGradAnisotropic::createThermalDeformGrad()
{
 //_thermal_deform_grad[_qp].zero();
  _thermal_deform_grad[_qp].addIa(1.0);
    
//  Real theta1 = std::exp(_thermal_expansion_coeff[0] * (_temperature[_qp] - _stress_free_temp));
//  Real theta2 = std::exp(_thermal_expansion_coeff[1] * (_temperature[_qp] - _stress_free_temp));
//  Real theta3 = std::exp(_thermal_expansion_coeff[2] * (_temperature[_qp] - _stress_free_temp));
    
  Real theta1 = std::exp(_thermal_expansion_coeff1 * (_temperature[_qp] - _stress_free_temp));
  Real theta2 = std::exp(_thermal_expansion_coeff2 * (_temperature[_qp] - _stress_free_temp));
  Real theta3 = std::exp(_thermal_expansion_coeff3 * (_temperature[_qp] - _stress_free_temp));
    
  RankTwoTensor I1(1,0,0, 0,0,0);
  RankTwoTensor I2(0,1,0, 0,0,0);
  RankTwoTensor I3(0,0,1, 0,0,0);
    
    std::cout << "I100:" << I1(0,0) << std::endl;
    std::cout << "I111:" << I1(1,1) << std::endl;
    std::cout << "I122:" << I1(2,2) << std::endl;
    std::cout << "I110:" << I1(1,0) << std::endl;
    std::cout << "I112:" << I1(1,2) << std::endl;
    std::cout << "I120:" << I1(2,0) << std::endl;

    
  RankTwoTensor ft = theta1 * I1 + theta2 * I2 + theta3 * I3;

  if (_read_prop_user_object_eigenstrain)
  {
     _Euler_angles_mat_prop[_qp](0) = _read_prop_user_object_eigenstrain->getData(_current_elem, 0);
     _Euler_angles_mat_prop[_qp](1) = _read_prop_user_object_eigenstrain->getData(_current_elem, 1);
     _Euler_angles_mat_prop[_qp](2) = _read_prop_user_object_eigenstrain->getData(_current_elem, 2);
  }
//  else
//    _Euler_angles_mat_prop[_qp] = _Euler_angles_eigenstrain;
 
  RotationTensor R(_Euler_angles_mat_prop[_qp]);
  R.update(_Euler_angles_mat_prop[_qp]);
  RankTwoTensor crysrot_ft = R.transpose();
  ft.rotate(crysrot_ft);
  
    std::cout << "theta1" << theta1 << std::endl;
    std::cout << "ft00:" << ft(0,0) << std::endl;
    std::cout << "ft11:" << ft(1,1) << std::endl;
    std::cout << "ft22:" << ft(2,2) << std::endl;
    std::cout << "ft10:" << ft(1,0) << std::endl;
    std::cout << "ft12:" << ft(1,2) << std::endl;
    std::cout << "ft20:" << ft(2,0) << std::endl;

    
  _thermal_deform_grad[_qp] = ft;

 }
