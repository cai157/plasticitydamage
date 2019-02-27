/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "GrainTrackerCrysrot.h"
#include "EulerAngleProvider.h"
#include "RotationTensor.h"


template<>
InputParameters validParams<GrainTrackerCrysrot>()
{
  InputParameters params = validParams<GrainTracker>();
  params.addParam<bool>("random_rotations", true, "Generate random rotations when the Euler Angle provider runs out of data (otherwise error out)");
  params.addRequiredParam<UserObjectName>("euler_angle_provider", "Name of Euler angle provider user object");
  return params;
}

GrainTrackerCrysrot::GrainTrackerCrysrot(const InputParameters & parameters) :
    GrainDataTracker<EulerAngles>(parameters),
    _random_rotations(getParam<bool>("random_rotations")),

    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
}

EulerAngles
GrainTrackerCrysrot::newGrain(unsigned int new_grain_id)
{
  EulerAngles angles;
      
  if (new_grain_id < _euler.getGrainNum())
    {
    angles = _euler.getEulerAngles(new_grain_id);
    }
  else
  {
    if (_random_rotations)
      angles.random();
    else{
      mooseError("GrainTrackerElasticity has run out of grain rotation data.");
      _console << "else in GrainTrackerCrysrot" << std::endl;
    }
  }
    
   //RealVectorValue angle2 = RealVectorValue(angles);
    
  return angles;
}
