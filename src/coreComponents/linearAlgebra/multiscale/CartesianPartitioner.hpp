/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CartesianPartitioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_CARTESIANPARTITIONER_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_CARTESIANPARTITIONER_HPP_

#include "linearAlgebra/multiscale/PartitionerBase.hpp"

namespace geosx
{
namespace multiscale
{

class CartesianPartitioner : public PartitionerBase
{
public:

  explicit CartesianPartitioner( LinearSolverParameters::Multiscale::Coarsening params );

  virtual localIndex generate( multiscale::MeshLevel const & mesh,
                               arrayView1d< localIndex > const & partition ) override;

  virtual void setCoarseData( multiscale::MeshLevel & coarseMesh ) const override;

private:

  LinearSolverParameters::Multiscale::Coarsening m_params;

  integer m_cartNumPart[3]{};

};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_CARTESIANPARTITIONER_HPP_
