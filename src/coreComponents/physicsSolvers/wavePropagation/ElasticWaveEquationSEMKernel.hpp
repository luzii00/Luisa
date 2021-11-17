/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ElasticWaveEquationSEMKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"


namespace geosx
{

/// Namespace to contain the elastic wave kernels.
namespace ElasticWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Convert a mesh element point coordinate into a coordinate on the reference element
   * @param coords coordinate of the point
   * @param coordsOnRefElem to contain the coordinate computed in the reference element
   * @param elementIndex index of the element containing the coords
   * @param faceNodeIndices array of face of the element
   * @param elemsToNodes map to obtaint global nodes from element index
   * @param X array of mesh nodes coordinates
   * @return true if coords is inside the element num index
   */
  template< typename FE_TYPE >
  GEOSX_HOST_DEVICE
  static bool
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        real64 const (&elemCenter)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arraySlice1d< localIndex const > const elemsToFaces,
                                        ArrayOfArraysView< localIndex const > const & facesToNodes,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                                        real64 (& coordsOnRefElem)[3] )
  {
    bool const isInsidePolyhedron =
      computationalGeometry::isPointInsidePolyhedron( X,
                                                      elemsToFaces,
                                                      facesToNodes,
                                                      elemCenter,
                                                      coords );
    if( isInsidePolyhedron )
    {
      real64 xLocal[FE_TYPE::numNodes][3]{};
      for( localIndex a = 0; a < FE_TYPE::numNodes; ++a )
      {
        LvArray::tensorOps::copy< 3 >( xLocal[a], X[ elemsToNodes[a] ] );
      }

      // coordsOnRefElem = invJ*(coords-coordsNode_0)
      real64 invJ[3][3]{};
      FE_TYPE::invJacobianTransformation( 0, xLocal, invJ );

      for( localIndex i = 0; i < 3; ++i )
      {
        // init at (-1,-1,-1) as the origin of the referential elem
        coordsOnRefElem[i] = -1.0;
        for( localIndex j = 0; j < 3; ++j )
        {
          coordsOnRefElem[i] += invJ[i][j] * (coords[j] - xLocal[0][j]);
        }
      }
      return true;
    }
    return false;

  }

  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] numNodesPerElem number of nodes per element
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] facesToNodes map from faces to nodes
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsLocal flag indicating whether the source is local or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceNodeConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverNodeConstants constant part of the receiver term
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const numNodesPerElem,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const & facesToNodes,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsLocal,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants )
  {

    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      real64 const center[3] = { elemCenter[k][0],
                                 elemCenter[k][1],
                                 elemCenter[k][2] };

      // Step 1: locate the sources, and precompute the source term

      /// loop over all the source that haven't been found yet
      for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
      {
        if( sourceIsLocal[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          real64 coordsOnRefElem[3]{};
          bool const sourceFound =
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             center,
                                                             elemsToNodes[k],
                                                             elemsToFaces[k],
                                                             facesToNodes,
                                                             X,
                                                             coordsOnRefElem );
          if( sourceFound )
          {
            sourceIsLocal[isrc] = 1;
            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }
          }
        }
      } // end loop over all sources


      // Step 2: locate the receivers, and precompute the receiver term

      /// loop over all the receivers that haven't been found yet
      for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
      {
        if( receiverIsLocal[ircv] == 0 )
        {
          real64 const coords[3] = { receiverCoordinates[ircv][0],
                                     receiverCoordinates[ircv][1],
                                     receiverCoordinates[ircv][2] };

          real64 coordsOnRefElem[3]{};
          bool const receiverFound =
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             center,
                                                             elemsToNodes[k],
                                                             elemsToFaces[k],
                                                             facesToNodes,
                                                             X,
                                                             coordsOnRefElem );
          if( receiverFound )
          {
            receiverIsLocal[ircv] = 1;

            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][a] = elemsToNodes[k][a];
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // end loop over receivers
    } );

  }
};

template< typename FE_TYPE >
struct MassMatrixKernel
{

  MassMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the mass and damping matrices
   * @param[in] size the number of cells in the subRegion
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] density cell-wise density
   * @param[out] mass diagonal of the mass matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real64 const > const density,
          arrayView1d< real64 > const mass )
  {
    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 xLocal[ numNodesPerElem ][ 3 ];
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      real64 N[ numNodesPerElem ];
      real64 gradN[ numNodesPerElem ][ 3 ];

      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        FE_TYPE::calcN( q, N );
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        for( localIndex a = 0; a < numNodesPerElem; ++a )
        {
          real64 const localIncrement = density[k] * detJ * N[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes[k][a]], localIncrement );
        }
      }
    } ); // end loop over element
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

template< typename FE_TYPE >
struct DampingMatrixKernel
{
  DampingMatrixKernel( FE_TYPE const & finiteElement )
   : m_finiteElement( finiteElement )
   {}
  /**
   * @brief Launches the precomputation of the damping matrix
   * @param[in]  X coordinates of the nodes
   * @param[in]  elemsToNodes map from element to nodes
   * @param[in]  targetSet Set of faces targetted here
   * @param[in]  faceToElemIndex map from face to element index
   * @param[in]  faceToNodes map from faces to nodes
   * @param[in]  faceNormal Normal to the faces
   * @param[in]  density cell-wise density
   * @param[in]  velocityVp cell-wise P-wave velocity
   * @param[in]  velocityVs cell-wise S-wave velocity
   * @param[out] damping_x diagonal of the damping matrix in x direction
   * @param[out] damping_y diagonal of the damping matrix in y direction
   * @param[out] damping_z diagonal of the damping matrix in z direction
   */

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          SortedArrayView< localIndex const > const  targetSet,
          arrayView2d< localIndex const > const  faceToElemIndex,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView2d< real64 const > const faceNormal,
          arrayView1d< real64 const > const density,
          arrayView1d< real64 > const velocityVp,
          arrayView1d< real64 > const velocityVs,
          arrayView1d< real64 > const damping_x,
          arrayView1d< real64 > const damping_y,
          arrayView1d< real64 > const damping_z)
  {
    localIndex const numNodesPerFace = 4;

    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
    constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

            
    for( localIndex i = 0; i < targetSet.size(); ++i )
    {
      localIndex const kf = targetSet[ i ];
      localIndex const k = faceToElemIndex[kf][0];
    
      real64 xLocal[numNodesPerElem][3];

      for( localIndex a=0; a< numNodesPerElem; ++a )
      {
        for( localIndex b=0; b<3; ++b )
        {
          xLocal[a][b] = X( elemsToNodes( k, a ), b );
        }
      }

      real64 N[numNodesPerElem];
      real64 gradN[ numNodesPerElem ][ 3 ];

      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {
        FE_TYPE::calcN( q, N );
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        real64 invJ[3][3]={{0}};
        FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

        for( localIndex a=0; a < numNodesPerFace; ++a )
        {
          /// compute ds=||detJ*invJ*normalFace_{kfe}||
          real64 tmp[3]={0};
          real64 ds = 0.0;
          for( localIndex b=0; b<3; ++b )
          {
            for( localIndex j = 0; j < 3; ++j )
            {
              tmp[b] += invJ[j][b]*faceNormal[kf][j];
            }
            ds +=tmp[b]*tmp[b];
          }
          ds = std::sqrt( ds );
  
          localIndex numNodeGl = facesToNodes[kf][a];
      
                 
          // Damping in x=xpos direction
          if ( faceNormal[kf][0] > 0.0 )
          {
            real64 const alpha_x = density[k] * velocityVp[k];
            real64 const alpha_y = density[k] * velocityVs[k];
            real64 const alpha_z = density[k] * velocityVs[k];
            damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
            damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
            damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
          }

          // Damping in x=xneg direction
          if ( faceNormal[kf][0] < 0.0 )
          {
            real64 const alpha_x = density[k] * velocityVp[k];
            real64 const alpha_y = density[k] * velocityVs[k];
            real64 const alpha_z = density[k] * velocityVs[k];
            damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
            damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
            damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
          }
          //Damping in y=ypos direction
          if ( faceNormal[kf][1] > 0.0 )
          {
            real64 const alpha_x = density[k] * velocityVs[k];
            real64 const alpha_y = density[k] * velocityVp[k];
            real64 const alpha_z = density[k] * velocityVs[k];
            damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
            damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
            damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
          }

          //Damping in y=yneg direction
          if ( faceNormal[kf][1] < 0.0 )
          {
            real64 const alpha_x = density[k] * velocityVs[k];
            real64 const alpha_y = density[k] * velocityVp[k];
            real64 const alpha_z = density[k] * velocityVs[k];
            damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
            damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
            damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];    
          }

          //Damping in z=zpos direction
          if ( faceNormal[kf][2] > 0.0 )
          {
            real64 const alpha_x = density[k] * velocityVs[k];
            real64 const alpha_y = density[k] * velocityVs[k];
            real64 const alpha_z = density[k] * velocityVp[k];
            damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
            damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
            damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
          }   
          //Damping in z=zneg direction
          if ( faceNormal[kf][2] < 0.0 )
          {
            real64 const alpha_x = density[k] * velocityVs[k];
            real64 const alpha_y = density[k] * velocityVs[k];
            real64 const alpha_z = density[k] * velocityVp[k];
            damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
            damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
            damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
          }   
        }
      }
    }  

  }
  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};


/**
 * @brief Implements kernels for solving the acoustic wave equations
 *   explicit central FD method and SEM
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitElasticDisplacementSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                                         CONSTITUTIVE_TYPE,
                                                                         FE_TYPE,
                                                                         1,
                                                                         1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   *   elements to be processed during this kernel launch.
   */
   ExplicitElasticDisplacementSEM( NodeManager & nodeManager,
                                   EdgeManager const & edgeManager,
                                   FaceManager const & faceManager,
                                   localIndex const targetRegionIndex,
                                   SUBREGION_TYPE const & elementSubRegion,
                                   FE_TYPE const & finiteElementSpace,
                                   CONSTITUTIVE_TYPE & inputConstitutiveType,

                                   real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_ux_np1( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >() ),
    m_uy_np1( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >() ),
    m_uz_np1( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >() ),
    m_stressxx( elementSubRegion.template getExtrinsicData< extrinsicMeshData::Stresstensor_xx>() ),
    m_stressyy( elementSubRegion.template getExtrinsicData< extrinsicMeshData::Stresstensor_yy>() ),
    m_stresszz( elementSubRegion.template getExtrinsicData< extrinsicMeshData::Stresstensor_zz>() ),
    m_stressxy( elementSubRegion.template getExtrinsicData< extrinsicMeshData::Stresstensor_xy>() ),
    m_stressxz( elementSubRegion.template getExtrinsicData< extrinsicMeshData::Stresstensor_xz>() ),
    m_stressyz( elementSubRegion.template getExtrinsicData< extrinsicMeshData::Stresstensor_yz>() ),
    m_mass( nodeManager.getExtrinsicData< extrinsicMeshData::MassVector > ()),
    m_damping_x( nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_x > ()),
    m_damping_y( nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_y > ()),
    m_damping_z( nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_z > ()),
    m_dt( dt )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }



  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitAcousticSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      xLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
  };
  //***************************************************************************


  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    /// numDofPerTrialSupportPoint = 1
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticSEM Description
   * Calculates stiffness vector
   *
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    real64 N[numNodesPerElem];
    real64 gradN[ numNodesPerElem ][ 3 ];

    FE_TYPE::calcN( q, N );
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, gradN );

    real64 uelemx[numNodesPerElem] = {{0.0}};
    real64 uelemy[numNodesPerElem] = {{0.0}};
    real64 uelemz[numNodesPerElem] = {{0.0}};
    real64 flowx[numNodesPerElem] = {{0.0}};
    real64 flowy[numNodesPerElem] = {{0.0}};
    real64 flowz[numNodesPerElem] = {{0.0}};

    //Calcul de u pour le premier ordre
    // Pré-multiplication par la masse globale
    for (localIndex i = 0; i < numNodesPerElem; i++)
    {         
      uelemx[i] = (m_mass[m_elemsToNodes[k][i]]-0.5*m_dt*m_damping_x[m_elemsToNodes[k][i]])*m_ux_np1[m_elemsToNodes[k][i]];
      uelemy[i] = (m_mass[m_elemsToNodes[k][i]]-0.5*m_dt*m_damping_y[m_elemsToNodes[k][i]])*m_uy_np1[m_elemsToNodes[k][i]];
      uelemz[i] = (m_mass[m_elemsToNodes[k][i]]-0.5*m_dt*m_damping_z[m_elemsToNodes[k][i]])*m_uz_np1[m_elemsToNodes[k][i]];
    }
          
    // Intégrale de Volume
    for (localIndex j = 0; j < numNodesPerElem; j++)
    {
      for (localIndex i = 0; i < numNodesPerElem; i++)
      {
        real64 dfx = detJ*gradN[i][0]*N[j];
        real64 dfy = detJ*gradN[i][1]*N[j];
        real64 dfz = detJ*gradN[i][2]*N[j];
        flowx[i] -= m_stressxx[k][j]*dfx + m_stressxy[k][j]*dfy + m_stressxz[k][j]*dfz;
        flowy[i] -= m_stressxy[k][j]*dfx + m_stressyy[k][j]*dfy + m_stressyz[k][j]*dfz;
        flowz[i] -= m_stressxz[k][j]*dfx + m_stressyz[k][j]*dfy + m_stresszz[k][j]*dfz;
      }

    }

    // Time update
    for (localIndex i = 0; i < numNodesPerElem; i++)
    {
      uelemx[i]+=m_dt*flowx[i];
      uelemy[i]+=m_dt*flowy[i];
      uelemz[i]+=m_dt*flowz[i];
    }

    // Mult by inverse mass matrix
    for (localIndex i = 0; i < numNodesPerElem; i++)
    {
    m_ux_np1[m_elemsToNodes[k][i]] = uelemx[i]/(m_mass[m_elemsToNodes[k][i]]+0.5*m_dt*m_damping_x[m_elemsToNodes[k][i]]);
    m_uy_np1[m_elemsToNodes[k][i]] = uelemy[i]/(m_mass[m_elemsToNodes[k][i]]+0.5*m_dt*m_damping_y[m_elemsToNodes[k][i]]);
    m_uz_np1[m_elemsToNodes[k][i]] = uelemz[i]/(m_mass[m_elemsToNodes[k][i]]+0.5*m_dt*m_damping_z[m_elemsToNodes[k][i]]);
    }

  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array in x direction.
  arrayView1d< real64 > const m_ux_np1;

  /// The array containing the nodal displacement array in y direction.
  arrayView1d< real64 > const m_uy_np1;

  /// The array containing the nodal displacement array in z direction.
  arrayView1d< real64 > const m_uz_np1;

  /// The array containing the xx component of the stresstensor
  arrayView2d < real64 const > const m_stressxx;

  /// The array containing the yy component of the stresstensor
  arrayView2d < real64 const > const m_stressyy;

  /// The array containing the zz component of the stresstensor
  arrayView2d < real64 const > const m_stresszz;

  /// The array containing the xy component of the stresstensor
  arrayView2d < real64 const > const m_stressxy;

  /// The array containing the xz component of the stresstensor
  arrayView2d < real64 const > const m_stressxz;

  /// The array containing the yz component of the stresstensor
  arrayView2d < real64 const > const m_stressyz;

  /// The array containing the diagonal of the mass matrix
  arrayView1d< real64 const > const m_mass;

  /// The array containing the diagonal of the damping matrix (x-part)
  arrayView1d< real64 const> const m_damping_x;

  /// The array containing the diagonal of the damping matrix (y-part)
  arrayView1d< real64 const > const m_damping_y;

  /// The array containing the diagonal of the damping matrix (z-part)
  arrayView1d< real64 const> const m_damping_z;

  /// The time increment for this time integration step.
  real64 const m_dt;


};

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitElasticStressSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                                   CONSTITUTIVE_TYPE,
                                                                   FE_TYPE,
                                                                   1,
                                                                   1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   *   elements to be processed during this kernel launch.
   */
   ExplicitElasticStressSEM( NodeManager & nodeManager,
                             EdgeManager const & edgeManager,
                             FaceManager const & faceManager,
                             localIndex const targetRegionIndex,
                             SUBREGION_TYPE const & elementSubRegion,
                             FE_TYPE const & finiteElementSpace,
                             CONSTITUTIVE_TYPE & inputConstitutiveType,
                             arrayView1d < real64 > const & inputMu,
                             arrayView1d < real64 > const & inputLambda,
                             arrayView2d < real64 > const & inputStressxx,
                             arrayView2d < real64 > const & inputStressyy,
                             arrayView2d < real64 > const & inputStresszz,
                             arrayView2d < real64 > const & inputStressxy,
                             arrayView2d < real64 > const & inputStressxz,
                             arrayView2d < real64 > const & inputStressyz,
                             real64 const & dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_ux_np1( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >() ),
    m_uy_np1( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >() ),
    m_uz_np1( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >() ),
    m_stressxx(inputStressxx),
    m_stressyy(inputStressyy),
    m_stresszz(inputStresszz),
    m_stressxy(inputStressxy),
    m_stressxz(inputStressxz),
    m_stressyz(inputStressyz),
    m_density( elementSubRegion.template getExtrinsicData< extrinsicMeshData::MediumDensity >() ),
    m_velocityVp( elementSubRegion.template getExtrinsicData< extrinsicMeshData::MediumVelocityVp >() ),
    m_velocityVs( elementSubRegion.template getExtrinsicData< extrinsicMeshData::MediumVelocityVs >() ),
    m_lambda(inputLambda),
    m_mu(inputMu),
    m_rhs( nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS > ()),
    m_dt( dt )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }



  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitElasticSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      xLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
  };
  //***************************************************************************


  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    /// numDofPerTrialSupportPoint = 1
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticSEM Description
   * Calculates stiffness vector
   *
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    
    m_mu[k] = m_density[k] * m_velocityVs[k] * m_velocityVs[k];
    m_lambda[k] = m_density[k] * m_velocityVp[k] * m_velocityVp[k] - 2.0*m_mu[k];

    real64 N[numNodesPerElem];
    real64 gradN[ numNodesPerElem ][ 3 ];


    real64 uelemxx[numNodesPerElem] = {{0.0}};
    real64 uelemyy[numNodesPerElem] = {{0.0}};
    real64 uelemzz[numNodesPerElem] = {{0.0}};
    real64 uelemxy[numNodesPerElem] = {{0.0}};
    real64 uelemxz[numNodesPerElem] = {{0.0}};
    real64 uelemyz[numNodesPerElem] = {{0.0}};
    real64 auxx[numNodesPerElem] = {{0.0}};
    real64 auyy[numNodesPerElem] = {{0.0}};
    real64 auzz[numNodesPerElem] = {{0.0}};
    real64 auxy[numNodesPerElem] = {{0.0}};
    real64 auxz[numNodesPerElem] = {{0.0}};
    real64 auyz[numNodesPerElem] = {{0.0}};


    FE_TYPE::calcN( q, N );
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, gradN );

    // // Computation of the stress tensor
    // //Pre-mult by mass matrix
    for (localIndex i = 0; i < numNodesPerElem; i++)
    {
      uelemxx[i] = detJ*m_stressxx[k][i];
      uelemyy[i] = detJ*m_stressyy[k][i];
      uelemzz[i] = detJ*m_stresszz[k][i];
      uelemxy[i] = detJ*m_stressxy[k][i];
      uelemxz[i] = detJ*m_stressxz[k][i];
      uelemyz[i] = detJ*m_stressyz[k][i];
    }

    //Volume integral
    for (localIndex j = 0; j < numNodesPerElem; j++)
    {
      for (localIndex i = 0; i < numNodesPerElem; i++)
      {
        real64 dfx2 = detJ*gradN[j][0]*N[i];
        real64 dfy2 = detJ*gradN[j][1]*N[i];
        real64 dfz2 = detJ*gradN[j][2]*N[i];
        auxx[i]+= dfx2*m_ux_np1[m_elemsToNodes[k][j]];
        auyy[i]+= dfy2*m_uy_np1[m_elemsToNodes[k][j]];
        auzz[i]+= dfz2*m_uz_np1[m_elemsToNodes[k][j]];
        auxy[i]+= dfx2*m_uy_np1[m_elemsToNodes[k][j]]+dfy2*m_ux_np1[m_elemsToNodes[k][j]];
        auxz[i]+= dfx2*m_uz_np1[m_elemsToNodes[k][j]]+dfz2*m_ux_np1[m_elemsToNodes[k][j]];
        auyz[i]+= dfy2*m_uz_np1[m_elemsToNodes[k][j]]+dfz2*m_uy_np1[m_elemsToNodes[k][j]];
      }

    }

    //Time integration
    for (localIndex i = 0; i < numNodesPerElem; i++)
    {
      real64 diag = m_lambda[k]*(auxx[i]+auyy[i]+auzz[i]);
      uelemxx[i]+= m_dt*(diag+2*m_mu[k]*auxx[i]);
      uelemyy[i]+= m_dt*(diag+2*m_mu[k]*auyy[i]);
      uelemzz[i]+= m_dt*(diag+2*m_mu[k]*auzz[i]);
      uelemxy[i]+= m_dt*m_mu[k]*auxy[i];
      uelemxz[i]+= m_dt*m_mu[k]*auxz[i];
      uelemyz[i]+= m_dt*m_mu[k]*auyz[i];
    }

    // Multiplication by inverse mass matrix
    for (localIndex i = 0; i < numNodesPerElem; i++)
    {
      m_stressxx[k][i] = uelemxx[i]/(detJ);
      m_stressyy[k][i] = uelemyy[i]/(detJ);
      m_stresszz[k][i] = uelemzz[i]/(detJ);
      m_stressxy[k][i] = uelemxy[i]/(detJ);
      m_stressxz[k][i] = uelemxz[i]/(detJ);
      m_stressyz[k][i] = uelemyz[i]/(detJ);
    }


    //Source injection

    for (localIndex i = 0; i < numNodesPerElem; i++)
    {
      m_stressxx[k][i]+= m_dt*(m_rhs[m_elemsToNodes[k][i]]/(detJ));
      m_stressyy[k][i]+= m_dt*(m_rhs[m_elemsToNodes[k][i]]/(detJ));
      m_stresszz[k][i]+= m_dt*(m_rhs[m_elemsToNodes[k][i]]/(detJ));

    }
   

  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array in x direction.
  arrayView1d< real64 const > const m_ux_np1;

  /// The array containing the nodal displacement array in y direction.
  arrayView1d< real64 const > const m_uy_np1;

  /// The array containing the nodal displacement array in z direction.
  arrayView1d< real64 const > const m_uz_np1;

  /// The array containing the xx component of the stresstensor
  arrayView2d < real64 > const m_stressxx;

  /// The array containing the yy component of the stresstensor
  arrayView2d < real64 > const m_stressyy;

  /// The array containing the zz component of the stresstensor
  arrayView2d < real64 > const m_stresszz;

  /// The array containing the xy component of the stresstensor
  arrayView2d < real64 > const m_stressxy;

  /// The array containing the xz component of the stresstensor
  arrayView2d < real64 > const m_stressxz;

  /// The array containing the yz component of the stresstensor
  arrayView2d < real64 > const m_stressyz;

  /// The array containing the density of the medium 
  arrayView1d< real64 const > const m_density;

  /// The array containing the P-wavespeed
  arrayView1d < real64 const > const m_velocityVp;

  /// The array containing the S-wavespeed
  arrayView1d < real64 const > const m_velocityVs;

  /// The array containing one of the Lamé coefficient: Lambda
  arrayView1d < real64  > const m_lambda;

  /// The array containing one of the Lamé coefficient: Mu
  arrayView1d< real64 > const m_mu;

  /// The array containing the RHS
  arrayView1d< real64  > const m_rhs;

  /// The time increment for this time integration step.
  real64 const m_dt;


};


/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitElasticDisplacementSEMFactory = finiteElement::KernelFactory< ExplicitElasticDisplacementSEM, real64>;

/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitElasticStressSEMFactory = finiteElement::KernelFactory< ExplicitElasticStressSEM,
                                                                      arrayView1d< real64 > const,
                                                                      arrayView1d< real64 > const,
                                                                      arrayView2d < real64 > const,
                                                                      arrayView2d < real64 > const,
                                                                      arrayView2d < real64 > const,
                                                                      arrayView2d < real64 > const,
                                                                      arrayView2d < real64 > const,
                                                                      arrayView2d < real64 > const,
                                                                      real64 >;                                                                 

} // namespace ElasticWaveEquationSEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMKERNEL_HPP_