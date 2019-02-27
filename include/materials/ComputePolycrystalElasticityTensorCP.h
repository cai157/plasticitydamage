/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H
#define COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H

#include "ComputeElasticityTensorBase.h"
#include "GrainDataTracker.h"
#include "RankTwoTensor.h"
#include "EulerAngleProvider.h"
//#include "FiniteStrainCrystalPlasticityGG.h"


//Forward Declarations
class ComputePolycrystalElasticityTensorCP;
class EulerAngleProvider;

template<>
InputParameters validParams<ComputePolycrystalElasticityTensorCP>();

/**
 * Compute an evolving elasticity tensor coupled to a grain growth phase field model.
 */
class ComputePolycrystalElasticityTensorCP : public ComputeElasticityTensorBase
{
public:
  ComputePolycrystalElasticityTensorCP(const InputParameters & parameters);
    
protected:
  virtual void initQpStatefulProperties();

  virtual void computeQpElasticityTensor();
  Real _length_scale;
  Real _pressure_scale;

  /// Grain tracker object
  const GrainDataTracker<RankFourTensor> & _grain_tracker;
  const GrainDataTracker<EulerAngles> & _grain_tracker_crysrot;

  /// Number of order parameters
  unsigned int _op_num;

  /// Order parameters
  std::vector<const VariableValue *> _vals;

  /// vector of elasticity tensor material properties
  std::vector< MaterialProperty<RankFourTensor> *> _D_elastic_tensor;
  std::vector< MaterialProperty<Real> *> _D_plastic_energy;
    
  MaterialProperty<RankTwoTensor> & _crysrot;
  MaterialProperty<RankTwoTensor> & _crysrot_old;
  MaterialProperty<RankFourTensor> & _elasticity_tensor_old;
  std::string _ave_plastic_energy_name;
  MaterialProperty<Real> & _ave_plastic_energy;
    
  MaterialProperty<Real> & _plastic_energy;
  MaterialProperty<Real> & _plastic_energy_old;
    
  const MaterialProperty<std::vector<Real> > & _slip_strain_out_old;
   
  /// Conversion factor from J to eV
  const Real _JtoeV;
};

#endif //COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H
