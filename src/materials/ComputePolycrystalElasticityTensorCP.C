/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputePolycrystalElasticityTensorCP.h"
#include "RotationTensor.h"
#include "EulerAngleProvider.h"
//#include "FiniteStrainCrystalPlasticityGG.h"

template<>
InputParameters validParams<ComputePolycrystalElasticityTensorCP>()
{
  InputParameters params = validParams<ComputeElasticityTensorBase>();
  params.addClassDescription("Compute an evolving elasticity tensor coupled to a grain growth phase field model.");
  params.addRequiredParam<UserObjectName>("grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addRequiredParam<UserObjectName>("grain_tracker_crysrot", "Name of GrainTracker user object that provides RankTwoTensors");
  params.addParam<Real>("length_scale", 1.0e-9, "Lengthscale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
  params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

ComputePolycrystalElasticityTensorCP::ComputePolycrystalElasticityTensorCP(const InputParameters & parameters) :
    ComputeElasticityTensorBase(parameters),
    _length_scale(getParam<Real>("length_scale")),
    _pressure_scale(getParam<Real>("pressure_scale")),
    _grain_tracker(getUserObject<GrainDataTracker<RankFourTensor>>("grain_tracker")),
    _grain_tracker_crysrot(getUserObject<GrainDataTracker<EulerAngles>>("grain_tracker_crysrot")),
    _op_num(coupledComponents("v")),
    _vals(_op_num),

    _D_elastic_tensor(_op_num),
    _D_plastic_energy(_op_num),

    _crysrot(declareProperty<RankTwoTensor>("crysrot")),
    _crysrot_old(declarePropertyOld<RankTwoTensor>("crysrot")),

    _elasticity_tensor_old(declarePropertyOld<RankFourTensor>(_elasticity_tensor_name)),

    _ave_plastic_energy_name("ave_plastic_energy"),
    _ave_plastic_energy(declareProperty<Real>(_ave_plastic_energy_name)),

    _plastic_energy(declareProperty<Real>("plastic_energy")),
    _plastic_energy_old(declarePropertyOld<Real>("plastic_energy")),

    _slip_strain_out_old(getMaterialPropertyOld<std::vector<Real>>("slip_strain_out")),

    _JtoeV(6.24150974e18)

{
  // Loop over variables (ops)
  for (auto op_index = decltype(_op_num)(0); op_index < _op_num; ++op_index)
  {
    // Initialize variables
    _vals[op_index] = &coupledValue("v", op_index);

    // declare elasticity tensor derivative properties
    _D_elastic_tensor[op_index] = &declarePropertyDerivative<RankFourTensor>(_elasticity_tensor_name, getVar("v", op_index)->name());
    _D_plastic_energy[op_index] = &declarePropertyDerivative<Real>(_ave_plastic_energy_name, getVar("v", op_index)->name());
      
  }
    
 }


void ComputePolycrystalElasticityTensorCP::initQpStatefulProperties()
{

    _crysrot[_qp].zero();
    _crysrot[_qp].addIa(1.0);
    
    _crysrot_old[_qp].zero();
    _crysrot_old[_qp].addIa(1.0);
    
    _elasticity_tensor[_qp].zero();
    _elasticity_tensor_old[_qp].zero();
    
    _plastic_energy[_qp] = 0.0;
    _plastic_energy_old[_qp] = 0.0;
}


void
ComputePolycrystalElasticityTensorCP::computeQpElasticityTensor()
{
    
  // Calculate elasticity tensor
  _elasticity_tensor[_qp].zero();
  _elasticity_tensor_old[_qp].zero();
    
  _crysrot[_qp].zero();
  _crysrot[_qp].addIa(1.0);
  
  _crysrot_old[_qp].zero();
  _crysrot_old[_qp].addIa(1.0);
    
  _plastic_energy[_qp] = 0.0;
    
  _ave_plastic_energy[_qp]= 0.0;
    
  EulerAngles angles;
  RealVectorValue angle2;
  angle2.zero();
   
  // Calculate plasticity energy
    
  Real aab;
  Real k = 7.2;
  int nss = _slip_strain_out_old[_qp].size() - 1;
    
  Real slip_strain_bar[_qp];
  slip_strain_bar[_qp] = 0.0;
    
  Real cap_a[_qp];
    
  for (unsigned int i = 0; i < nss; ++i){
      slip_strain_bar[_qp] += std::abs(_slip_strain_out_old[_qp][i]);
  }

  slip_strain_bar[_qp] = slip_strain_bar[_qp] / nss;

  if (slip_strain_bar[_qp] == 0)
      cap_a[_qp] = 0.0; //any value that won't give inf
  else
      cap_a[_qp] = 1.0 / (k * slip_strain_bar[_qp]) * 1.0E-6; //unit: MPa

    
  for (unsigned int i = 0; i < nss; ++i){
    for (unsigned int j = 0; j < nss; ++j){
        
        //mooseWarning("CP before calculation plastic_energy = " , plastic_energy[_qp]);
        
        if (i == j)
            aab = 1.0;
        else
            aab = 1.4;
            
        _plastic_energy[_qp] = cap_a[_qp] * std::pow(std::abs(_slip_strain_out_old[_qp][i] * _slip_strain_out_old[_qp][j]) * aab, 3/4) + _plastic_energy_old[_qp];

    }
  }

/*
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  Real sum_h = 0.0;
    
  for (auto op_index = beginIndex(op_to_grains); op_index < op_to_grains.size(); ++op_index)
  {

    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;
 
    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;
      
    // Sum all rotated elasticity tensors
    _elasticity_tensor[_qp] += _grain_tracker.getData(grain_id) * h;
      
    // Sum all plastic energy
    _ave_plastic_energy[_qp] += plastic_energy[_qp] * h;
    sum_h += h;
  }
    
  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _elasticity_tensor[_qp] /= sum_h;
  _ave_plastic_energy[_qp] /= sum_h;

 */
    
  const auto & op_to_grains_euler = _grain_tracker_crysrot.getVarToFeatureVector(_current_elem->id());
  Real sum_h = 0.0;

  for (auto op_index_euler = beginIndex(op_to_grains_euler); op_index_euler < op_to_grains_euler.size(); ++op_index_euler)
    {
        auto grain_id_euler = op_to_grains_euler[op_index_euler];
        if (grain_id_euler == FeatureFloodCount::invalid_id)
            continue;
        
        // Interpolation factor for elasticity tensors
        Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index_euler])[_qp] - 0.5))) / 2.0;
        
        // Sum all rotated elasticity tensors
        angles = _grain_tracker_crysrot.getData(grain_id_euler) ;
        angle2 += RealVectorValue(angles) * h;
        
        _elasticity_tensor[_qp] += _grain_tracker.getData(grain_id_euler) * h;
        
        // Sum all plastic energy
        _ave_plastic_energy[_qp] += _plastic_energy[_qp] * h;
       
        sum_h += h;
       // mooseWarning("CP in loop sum_h = " , sum_h);
    }
    
    const Real tol = 1.0e-10;
    sum_h = std::max(sum_h, tol);
    angle2 /= sum_h;
    
    RotationTensor R(angle2);
    R.update(angle2);
    _crysrot[_qp] = R.transpose();

    _elasticity_tensor[_qp] /= sum_h;
    _ave_plastic_energy[_qp] /= sum_h;
    
    //mooseWarning("CP sum_h = " , sum_h);

  // Calculate elasticity tensor derivative: Cderiv = dhdopi/sum_h * (Cop - _Cijkl)
  for (auto op_index = decltype(_op_num)(0); op_index < _op_num; ++op_index)
  {
    (*_D_elastic_tensor[op_index])[_qp].zero();
    (*_D_plastic_energy[op_index])[_qp] = 0.0;

  }

  for (auto op_index = beginIndex(op_to_grains_euler); op_index < op_to_grains_euler.size(); ++op_index)
  {
      
    auto grain_id = op_to_grains_euler[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    RankFourTensor & C_deriv = (*_D_elastic_tensor[op_index])[_qp];
  
    C_deriv = (_grain_tracker.getData(grain_id) - _elasticity_tensor[_qp]) * dhdopi / sum_h;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    C_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
      
    //mooseWarning("CP before plastic_energy2 = " , plastic_energy[_qp]);
      
    Real & P_deriv = (*_D_plastic_energy[op_index])[_qp];
    P_deriv = (_plastic_energy[_qp] - _ave_plastic_energy[_qp]) * dhdopi / sum_h;
      
   // mooseWarning("CP before plastic_energy3 = " , plastic_energy[_qp]);
   // mooseWarning("CP before P_deriv = " , P_deriv);
      
    P_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
   // mooseWarning("CP after P_deriv = " , P_deriv);

    }

}
