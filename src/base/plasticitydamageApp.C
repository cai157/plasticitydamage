#include "plasticitydamageApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "ComputeThermalExpansionEigenstrainAnisotropic.h"
#include "FiniteStrainCrystalPlasticityGG.h"
#include "ComputePolycrystalElasticityTensorCP.h"
#include "GrainTrackerCrysrot.h"
#include "ComputeThermalDeformGrad.h"
#include "ComputeThermalDeformGradAnisotropic.h"
#include "FiniteStrainCrystalPlasticityThermal.h"
#include "ComputeVolumetricDeformGradThermal.h"
#include "ACGrGrPlasticDrivingForce.h"
#include "PolycrystalPlasticDrivingForceAction.h"


template <>
InputParameters
validParams<plasticitydamageApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

plasticitydamageApp::plasticitydamageApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  plasticitydamageApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  plasticitydamageApp::associateSyntax(_syntax, _action_factory);
}

plasticitydamageApp::~plasticitydamageApp() {}

// External entry point for dynamic application loading
extern "C" void
plasticitydamageApp__registerApps()
{
  plasticitydamageApp::registerApps();
}
void
plasticitydamageApp::registerApps()
{
  registerApp(plasticitydamageApp);
}

// External entry point for dynamic object registration
extern "C" void
plasticitydamageApp__registerObjects(Factory & factory)
{
  plasticitydamageApp::registerObjects(factory);
}
void
plasticitydamageApp::registerObjects(Factory & factory)
{

    registerMaterial(ComputeThermalExpansionEigenstrainAnisotropic);
    registerMaterial(FiniteStrainCrystalPlasticityGG);
    registerMaterial(ComputePolycrystalElasticityTensorCP);
    registerPostprocessor(GrainTrackerCrysrot);
    registerMaterial(ComputeThermalDeformGradAnisotropic);
    registerMaterial(ComputeThermalDeformGrad);
    registerMaterial(FiniteStrainCrystalPlasticityThermal);
    registerMaterial(ComputeVolumetricDeformGradThermal);
    registerKernel(ACGrGrPlasticDrivingForce);

}

// External entry point for dynamic syntax association
extern "C" void
plasticitydamageApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  plasticitydamageApp::associateSyntax(syntax, action_factory);
}
void
plasticitydamageApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
    Registry::registerActionsTo(action_factory, {"plasticitydamageApp"});

    registerSyntax("PolycrystalPlasticDrivingForceAction", "Kernels/PolycrystalPlasticDrivingForce");
}
