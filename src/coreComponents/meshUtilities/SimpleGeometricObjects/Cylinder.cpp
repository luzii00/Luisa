/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * Cylinder.cpp
 *
 *  Created on: Aug 4, 2017
 *      Author: settgast
 */

#include "Cylinder.hpp"

namespace geosx
{
using namespace dataRepository;

Cylinder::Cylinder( const std::string& name, ManagedGroup * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_point1{0.0,0.0,0.0},
  m_point2{0.0,0.0,0.0},
  m_radius{0.0}
{
  RegisterViewWrapper( viewKeyStruct::point1String, &m_point1, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Center point of one (upper or lower) face of the cylinder");

  RegisterViewWrapper( viewKeyStruct::point2String, &m_point2, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Center point of the other face of the cylinder");

  RegisterViewWrapper( viewKeyStruct::radiusString, &m_radius, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Radius of the cylinder");

}

Cylinder::~Cylinder()
{}


bool Cylinder::IsCoordInObject( const R1Tensor& coord ) const
{
  bool rval = false;

  R1Tensor axisVector, coord_minus_point1, projection, distance;
  realT height;

  axisVector = m_point2;
  axisVector -= m_point1;
  height = axisVector.Normalize();

  coord_minus_point1 = coord;
  coord_minus_point1 -= m_point1;

  projection = axisVector;
  projection *= Dot( axisVector, coord_minus_point1 );

  distance  = coord_minus_point1;
  distance -= projection;

  if (distance.L2_Norm()<m_radius && projection.L2_Norm()<height)
  {
    rval = true;
  }

  return rval;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Cylinder, std::string const &, ManagedGroup * const )

} /* namespace geosx */
