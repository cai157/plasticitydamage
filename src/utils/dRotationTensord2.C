//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "dRotationTensord2.h"
#include "libmesh/libmesh.h"

dRotationTensord2::dRotationTensord2(Axis axis, Real angle) { update2(axis, angle); }

dRotationTensord2::dRotationTensord2(const RealVectorValue & euler_angles) { update2(euler_angles); }

void
dRotationTensord2::update2(Axis axis, Real angle)
{
  zero();

  RealVectorValue a;
  a(axis) = 1.0;

  const Real s = std::sin(angle * libMesh::pi / 180.0);
  const Real c = std::cos(angle * libMesh::pi / 180.0);

  // assemble row wise
  _coords[0] = a * RealVectorValue(1.0, -c, -c);
  _coords[1] = a * RealVectorValue(0.0, 0.0, s);
  _coords[2] = a * RealVectorValue(0.0, -s, 0.0);

  _coords[3] = a * RealVectorValue(0.0, 0.0, -s);
  _coords[4] = a * RealVectorValue(-c, 1.0, -c);
  _coords[5] = a * RealVectorValue(s, 0.0, 0.0);

  _coords[6] = a * RealVectorValue(0.0, s, 0.0);
  _coords[7] = a * RealVectorValue(-s, 0.0, 0.0);
  _coords[8] = a * RealVectorValue(-c, -c, 1.0);
}

void
dRotationTensord2::update2(const RealVectorValue & euler_angles)
{
  const Real phi_1 = euler_angles(0) * (libMesh::pi / 180.0);
  const Real Phi = euler_angles(1) * (libMesh::pi / 180.0);
  const Real phi_2 = euler_angles(2) * (libMesh::pi / 180.0);

  const Real c1 = std::cos(phi_1);
  const Real c2 = std::cos(Phi);
  const Real c3 = std::cos(phi_2);

  const Real s1 = std::sin(phi_1);
  const Real s2 = std::sin(Phi);
  const Real s3 = std::sin(phi_2);
    
  const Real dc1 = -std::sin(phi_1);
  const Real dc2 = -std::sin(Phi);
  const Real dc3 = -std::sin(phi_2);
    
  const Real ds1 = std::cos(phi_1);
  const Real ds2 = std::cos(Phi);
  const Real ds3 = std::cos(phi_2);


  // doing a Z1, X2, Z3 rotation
  // RealTensorValue is formed row-wise

  _coords[0] = -dc2 * s1 * s3;  // R11d2
  _coords[3] = -dc2 * c3 * s1;  // R12d2
  _coords[6] = s1 * ds2;         // R13d2

  _coords[1] = c1 * dc2 * s3; // R21d2
  _coords[4] = c1 * dc2 * c3; // R22d2
  _coords[7] = -c1 * ds2;     // R23d2

  _coords[2] = ds2 * s3; // R31d2
  _coords[5] = c3 * ds2; // R32d2
  _coords[8] = dc2;      // R33d2
}
