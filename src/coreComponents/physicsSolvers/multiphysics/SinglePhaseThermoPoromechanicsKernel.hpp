/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseThermoPoromechanicsKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASETHERMOPOROMECHANICSKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASETHERMOPOROMECHANICSKERNEL_HPP_

#include "SinglePhaseThermoPoromechanicsSolver.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

namespace geosx
{

namespace ThermoPoromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase thermo-poromechanics problems.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### SinglePhasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static single-phase ThermoPoromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhase :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Number of nodes per element which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;


  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  SinglePhase( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               localIndex const targetRegionIndex,
               SUBREGION_TYPE const & elementSubRegion,
               FE_TYPE const & finiteElementSpace,
               CONSTITUTIVE_TYPE & inputConstitutiveType,
               arrayView1d< globalIndex const > const & inputDispDofNumber,
               string const & inputFlowDofKey,
               arrayView1d< globalIndex const > const & temperatureDofNumber,
               globalIndex const rankOffset,
               CRSMatrixView< real64, globalIndex const > const & inputMatrix,
               arrayView1d< real64 > const & inputRhs,
               real64 const (&inputGravityVector)[3],
               arrayView1d< string const > const fluidModelNames,
               real64 const dt ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDispDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement()),
    m_uhat( nodeManager.incrementalDisplacement()),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).density() ),
    m_fluidDensityOld( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::densityOldString() ) ),
    m_dFluidDensity_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).dDensity_dPressure() ),
    m_flowDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey )),
    m_fluidPressure( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString() ) ),
    m_deltaFluidPressure( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString() ) ),
    m_thermalDofNumber( temperatureDofNumber ),
    m_temperature( nodeManager.template getReference< array1d< real64 > >( SinglePhaseThermoPoromechanicsSolver::viewKeyStruct::temperatureString() ) ),
    m_deltaTemperature( nodeManager.template getReference< array1d< real64 > >( SinglePhaseThermoPoromechanicsSolver::viewKeyStruct::newDeltaTemperatureString() ) ),
    m_thermalDiffusion( 1e-1 ),
    m_dt( dt )
  {}


  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    static constexpr int numDispDofPerElem =  Base::StackVariables::numRows;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            u_local(),
            uhat_local(),
            localFlowResidual{ 0.0 },
      localDispFlowJacobian{ {0.0} },
      localFlowDispJacobian{ {0.0} },
      localFlowFlowJacobian{ {0.0} },
      localFlowDofIndex{ 0 },
      localThermalResidual{ 0.0 },
      localThermalJacobian{ {0.0} },
      localThermalDofIndex{ 0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    real64 localFlowResidual[1];
    real64 localDispFlowJacobian[numDispDofPerElem][1];
    real64 localFlowDispJacobian[1][numDispDofPerElem];
    real64 localFlowFlowJacobian[1][1];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localFlowDofIndex[1];

    real64 localThermalResidual[numNodesPerElem];
    real64 localThermalJacobian[numNodesPerElem][numNodesPerElem];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localThermalDofIndex[numNodesPerElem];


    /// C-array storage for the local nodal temperature.
    real64 localTemperature[numNodesPerElem];

    /// C-array storage for the local nodal temperature increment.
    real64 localDeltaTemperature[numNodesPerElem];

  };


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the SinglePhase implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
#endif
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }

      stack.localThermalDofIndex[a] = m_thermalDofNumber[localNodeIndex];
      stack.localTemperature[a] = m_temperature[localNodeIndex];
      stack.localDeltaTemperature[a] = m_deltaTemperature[localNodeIndex];
    }

    stack.localFlowDofIndex[0] = m_flowDofNumber[k];
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Get displacement: (i) basis functions (N), (ii) basis function
    // derivatives (dNdX), and (iii) determinant of the Jacobian transformation
    // matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Evaluate total stress tensor
    real64 strainIncrement[6] = {0.0};
    real64 totalStress[6];
    real64 dPorosity_dPressure;
    real64 dPorosity_dVolStrainIncrement;
    real64 dTotalStress_dPressure[6] = {0.0};

    // --- Update total stress tensor
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );

    m_constitutiveUpdate.smallStrainUpdate( k,
                                            q,
                                            m_fluidPressure[k],
                                            m_deltaFluidPressure[k],
                                            strainIncrement,
                                            totalStress,
                                            dPorosity_dPressure,
                                            dPorosity_dVolStrainIncrement,
                                            dTotalStress_dPressure,
                                            stiffness );

    real64 const porosityNew = m_constitutiveUpdate.getPorosity( k, q );
    real64 const porosityOld = m_constitutiveUpdate.getOldPorosity( k, q );

    // Evaluate body force vector
    real64 bodyForce[3] = { m_gravityVector[0],
                            m_gravityVector[1],
                            m_gravityVector[2]};

    if( m_gravityAcceleration > 0.0 )
    {
      real64 mixtureDensity = ( 1.0 - porosityNew ) * m_solidDensity( k, q ) + porosityNew * m_fluidDensity( k, q );
      mixtureDensity *= detJxW;
      bodyForce[0] *= mixtureDensity;
      bodyForce[1] *= mixtureDensity;
      bodyForce[2] *= mixtureDensity;
    }

    // Assemble local jacobian and residual

    // --- Momentum balance
    for( localIndex i=0; i<6; ++i )
    {
      totalStress[i] *= -detJxW;
    }

    FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                     totalStress,
                                     N,
                                     bodyForce,
                                     reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );

    stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian );

    for( integer a = 0; a < numNodesPerElem; ++a )
    {
      stack.localDispFlowJacobian[a*3+0][0] -= dNdX[a][0] * dTotalStress_dPressure[0] * detJxW;
      stack.localDispFlowJacobian[a*3+0][0] -= dNdX[a][2] * dTotalStress_dPressure[4] * detJxW;
      stack.localDispFlowJacobian[a*3+0][0] -= dNdX[a][1] * dTotalStress_dPressure[5] * detJxW;

      stack.localDispFlowJacobian[a*3+1][0] -= dNdX[a][1] * dTotalStress_dPressure[1] * detJxW;
      stack.localDispFlowJacobian[a*3+1][0] -= dNdX[a][2] * dTotalStress_dPressure[3] * detJxW;
      stack.localDispFlowJacobian[a*3+1][0] -= dNdX[a][0] * dTotalStress_dPressure[5] * detJxW;

      stack.localDispFlowJacobian[a*3+2][0] -= dNdX[a][2] * dTotalStress_dPressure[2] * detJxW;
      stack.localDispFlowJacobian[a*3+2][0] -= dNdX[a][1] * dTotalStress_dPressure[3] * detJxW;
      stack.localDispFlowJacobian[a*3+2][0] -= dNdX[a][0] * dTotalStress_dPressure[4] * detJxW;
    }

    if( m_gravityAcceleration > 0.0 )
    {
      // Assumptions: ( i) dMixtureDens_dVolStrain contribution is neglected
      //              (ii) grains are assumed incompressible

      real64 const dMixtureDens_dPressure = dPorosity_dPressure * ( -m_solidDensity( k, q ) + m_fluidDensity( k, q ) )
                                            + porosityNew * m_dFluidDensity_dPressure( k, q );
      for( integer a = 0; a < numNodesPerElem; ++a )
      {
        stack.localDispFlowJacobian[a*3+0][0] += N[a] * dMixtureDens_dPressure * m_gravityVector[0] * detJxW;
        stack.localDispFlowJacobian[a*3+1][0] += N[a] * dMixtureDens_dPressure * m_gravityVector[1] * detJxW;
        stack.localDispFlowJacobian[a*3+2][0] += N[a] * dMixtureDens_dPressure * m_gravityVector[2] * detJxW;
      }
    }

    // --- Mass balance accumulation
    for( integer a = 0; a < numNodesPerElem; ++a )
    {
      stack.localFlowDispJacobian[0][a*3+0] += dPorosity_dVolStrainIncrement * m_fluidDensity( k, q ) * dNdX[a][0] * detJxW;
      stack.localFlowDispJacobian[0][a*3+1] += dPorosity_dVolStrainIncrement * m_fluidDensity( k, q ) * dNdX[a][1] * detJxW;
      stack.localFlowDispJacobian[0][a*3+2] += dPorosity_dVolStrainIncrement * m_fluidDensity( k, q ) * dNdX[a][2] * detJxW;
    }

    stack.localFlowResidual[0] += ( porosityNew * m_fluidDensity( k, q ) - porosityOld * m_fluidDensityOld( k ) ) * detJxW;
    stack.localFlowFlowJacobian[0][0] += ( dPorosity_dPressure * m_fluidDensity( k, q ) + porosityNew * m_dFluidDensity_dPressure( k, q ) ) * detJxW;

    // Thermal conditributon

/**
    // Test an unit jacobian
    real64 tmpmatr[8][8] = { { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
      { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
      { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }
    };
*/
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      for( localIndex b=0; b<numNodesPerElem; ++b )
      {
        real64 tmpValue1 = N[a] * N[b] * detJxW / m_dt;
        real64 tmpValue2 = LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJxW * m_thermalDiffusion;

        // Update local Jacobian
        stack.localThermalJacobian[a][b] += tmpValue1 + tmpValue2; //= tmpmatr[a][b];

        // Update local Rhs
        stack.localThermalResidual[a] += tmpValue1 * stack.localDeltaTemperature[b];//= a+b+1.0;//
        stack.localThermalResidual[a] += tmpValue2 * stack.localTemperature[b];
      }
    }


  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

    constexpr int nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[numDofPerTestSupportPoint * localNode + dim],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidual[numDofPerTestSupportPoint * localNode + dim] );
        maxForce = fmax( maxForce, fabs( stack.localResidual[numDofPerTestSupportPoint * localNode + dim] ) );

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localFlowDofIndex,
                                                                                stack.localDispFlowJacobian[numDofPerTestSupportPoint * localNode + dim],
                                                                                1 );

      }
    }


    localIndex const dof = LvArray::integerConversion< localIndex >( stack.localFlowDofIndex[0] - m_dofRankOffset );
    if( 0 <= dof && dof < m_matrix.numRows() )
    {
      m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof,
                                                                      stack.localRowDofIndex,
                                                                      stack.localFlowDispJacobian[0],
                                                                      nUDof );
      m_matrix.template addToRow< serialAtomic >( dof,
                                                  stack.localFlowDofIndex,
                                                  stack.localFlowFlowJacobian[0],
                                                  1 );
      RAJA::atomicAdd< serialAtomic >( &m_rhs[dof], stack.localFlowResidual[0] );
    }


    // Assemblage thermal term
    for( int a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const thermalDof = LvArray::integerConversion< localIndex >( stack.localThermalDofIndex[a] - m_dofRankOffset );

      if( 0 <= thermalDof && thermalDof < m_matrix.numRows() )
      {
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( thermalDof,
                                                                                stack.localThermalDofIndex,
                                                                                stack.localThermalJacobian[a],
                                                                                numNodesPerElem );


        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[thermalDof], stack.localThermalResidual[a] );
      }
    }

    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

  /// The rank global densities
  arrayView2d< real64 const > const m_solidDensity;
  arrayView2d< real64 const > const m_fluidDensity;
  arrayView1d< real64 const > const m_fluidDensityOld;
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_fluidPressure;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > const m_deltaFluidPressure;

  /// The global degree of freedom number of the thermal field.
  arrayView1d< globalIndex const > const m_thermalDofNumber;

  /// The rank-global temperature array.
  arrayView1d< real64 const > const m_temperature;

  /// The rank-global incremental temperature array.
  arrayView1d< real64 const > const m_deltaTemperature;

  /// Thermal diffusion coefficient
  real64 const m_thermalDiffusion;

  /// Time step
  real64 const m_dt;
};

using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhase,
                                                               arrayView1d< globalIndex const > const &,
                                                               string const &,
                                                               arrayView1d< globalIndex const > const &,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const &,
                                                               arrayView1d< real64 > const &,
                                                               real64 const (&)[3],
                                                               arrayView1d< string const > const,
                                                               real64 const >;

} /* namespace ThermoPoromechanicsKernels */

} /* namespace geosx */

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASETHERMOPOROMECHANICSKERNEL_HPP_
