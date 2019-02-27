//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POLYCRYSTALPLASTICDRIVINGFORCEACTION_H
#define POLYCRYSTALPLASTICDRIVINGFORCEACTION_H

#include "Action.h"
#include "DerivativeMaterialPropertyNameInterface.h"

// Forward Declarations
class PolycrystalPlasticDrivingForceAction;

template <>
InputParameters validParams<PolycrystalPlasticDrivingForceAction>();
/**
 * Action that adds the Plastic driving force for each order parameter
 */
class PolycrystalPlasticDrivingForceAction : public Action,
                                             public DerivativeMaterialPropertyNameInterface
{
public:
  PolycrystalPlasticDrivingForceAction(const InputParameters & params);

  virtual void act();

private:
  /// Number of order parameters used in the model
  const unsigned int _op_num;

  /// Base name for the order parameters
  std::string _var_name_base;
  std::string _base_name;
  std::string _elasticity_tensor_name;
  std::string _ave_plastic_energy_name;
};

#endif // POLYCRYSTALPLASTICDRIVINGFORCEACTION_H
