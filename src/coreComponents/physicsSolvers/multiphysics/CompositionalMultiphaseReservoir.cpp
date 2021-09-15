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
 * @file CompositionalMultiphaseReservoir.cpp
 *
 */


#include "CompositionalMultiphaseReservoir.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanicsSolver.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

CompositionalMultiphaseReservoir::CompositionalMultiphaseReservoir( const string & name,
                                                                    Group * const parent ):
  ReservoirSolverBase( name, parent )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirFVM;
}

CompositionalMultiphaseReservoir::~CompositionalMultiphaseReservoir()
{}

void CompositionalMultiphaseReservoir::postProcessInput()
{
  ReservoirSolverBase::postProcessInput();

  GEOSX_THROW_IF( !m_flowSolver,
                  catalogName() << " " << getName()
                                << ": In postProcessInput, the flow solver has not been set yet",
                  InputError );

  // check that the flow solver is compatible with CompositionalMultiphaseReservoir
  bool const isSupported =
    dynamicCast< CompositionalMultiphaseFVM const * >( m_flowSolver ) ||
    dynamicCast< CompositionalMultiphaseHybridFVM const * >( m_flowSolver ) ||
    dynamicCast< MultiphasePoromechanicsSolver const * >( m_flowSolver );
  GEOSX_THROW_IF( !isSupported,
                  catalogName() << " " << getName()
                                << ": the solver of type " << m_flowSolver->catalogName()
                                << ", named " << m_flowSolver->getName()
                                << " cannot be used with " << catalogName(),
                  InputError );

  if( dynamicCast< CompositionalMultiphaseFVM const * >( m_flowSolver ) ||
      dynamicCast< CompositionalMultiphaseHybridFVM const * >( m_flowSolver ) )
  {
    m_wellSolver->setFlowSolverName( m_flowSolverName );
  }
  else
  {
    MultiphasePoromechanicsSolver const * solver = dynamicCast< MultiphasePoromechanicsSolver const * >( m_flowSolver );
    m_wellSolver->setFlowSolverName( solver->getFlowSolver()->getName() );
  }

}

void CompositionalMultiphaseReservoir::initializePostInitialConditionsPreSubGroups()
{
  ReservoirSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  // loop over the wells
  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion & subRegion )
  {

    array1d< array1d< arrayView3d< real64 const > > > permeability;
    if( dynamicCast< CompositionalMultiphaseFVM const * >( m_flowSolver ) )
    {
      CompositionalMultiphaseFVM const * solver =
        dynamicCast< CompositionalMultiphaseFVM const * >( m_flowSolver );
      permeability =
        elemManager.constructMaterialArrayViewAccessor< real64, 3 >( PermeabilityBase::viewKeyStruct::permeabilityString(),
                                                                     solver->targetRegionNames(),
                                                                     solver->permeabilityModelNames() );
    }
    else if( dynamicCast< CompositionalMultiphaseHybridFVM const * >( m_flowSolver ) )
    {
      CompositionalMultiphaseHybridFVM const * solver =
        dynamicCast< CompositionalMultiphaseHybridFVM const * >( m_flowSolver );
      permeability =
        elemManager.constructMaterialArrayViewAccessor< real64, 3 >( PermeabilityBase::viewKeyStruct::permeabilityString(),
                                                                     solver->targetRegionNames(),
                                                                     solver->permeabilityModelNames() );
    }
    else
    {
      MultiphasePoromechanicsSolver const * solver = dynamicCast< MultiphasePoromechanicsSolver const * >( m_flowSolver );
      permeability =
        elemManager.constructMaterialArrayViewAccessor< real64, 3 >( PermeabilityBase::viewKeyStruct::permeabilityString(),
                                                                     solver->getFlowSolver()->targetRegionNames(),
                                                                     solver->getFlowSolver()->permeabilityModelNames() );
    }

    // compute the Peaceman index (if not read from XML)
    PerforationData * const perforationData = subRegion.getPerforationData();
    perforationData->computeWellTransmissibility( meshLevel, subRegion, permeability );
  } );

  // bind the stored reservoir views to the current domain
  resetViews( domain );

  // set the MGR recipe
  if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirHybridFVM;
  }
}

void CompositionalMultiphaseReservoir::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                   DofManager const & dofManager,
                                                                   SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  // TODO: remove this and just call SolverBase::setupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  localIndex const resNDOF = m_wellSolver->numDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->numDofPerWellElement();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );
  string const resDofKey  = dofManager.getKey( m_wellSolver->resElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resDofNumber =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

  globalIndex const rankOffset = dofManager.rankOffset();

  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion const & subRegion )
  {
    PerforationData const * const perforationData = subRegion.getPerforationData();

    // get the well degrees of freedom and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get the well element indices corresponding to each perforation
    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );


    // Insert the entries corresponding to reservoir-well perforations
    // This will fill J_WR, and J_RW
    forAll< serialPolicy >( perforationData->size(), [=] ( localIndex const iperf )
    {
      stackArray1d< globalIndex, maxNumDof > eqnRowIndicesRes( resNDOF );
      stackArray1d< globalIndex, maxNumDof > eqnRowIndicesWell( wellNDOF );
      stackArray1d< globalIndex, maxNumDof > dofColIndicesRes( resNDOF );
      stackArray1d< globalIndex, maxNumDof > dofColIndicesWell( wellNDOF );

      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];
      localIndex const iwelem = perfWellElemIndex[iperf];

      for( localIndex idof = 0; idof < resNDOF; ++idof )
      {
        eqnRowIndicesRes[idof] = resDofNumber[er][esr][ei] + idof - rankOffset;
        dofColIndicesRes[idof] = resDofNumber[er][esr][ei] + idof;
      }

      for( localIndex idof = 0; idof < wellNDOF; ++idof )
      {
        eqnRowIndicesWell[idof] = wellElemDofNumber[iwelem] + idof - rankOffset;
        dofColIndicesWell[idof] = wellElemDofNumber[iwelem] + idof;
      }

      for( localIndex i = 0; i < eqnRowIndicesRes.size(); ++i )
      {
        if( eqnRowIndicesRes[i] >= 0 && eqnRowIndicesRes[i] < pattern.numRows() )
        {
          for( localIndex j = 0; j < dofColIndicesWell.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesRes[i], dofColIndicesWell[j] );
          }
        }
      }

      for( localIndex i = 0; i < eqnRowIndicesWell.size(); ++i )
      {
        if( eqnRowIndicesWell[i] >= 0 && eqnRowIndicesWell[i] < pattern.numRows() )
        {
          for( localIndex j = 0; j < dofColIndicesRes.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesWell[i], dofColIndicesRes[j] );
          }
        }
      }
    } );
  } );
}

void CompositionalMultiphaseReservoir::assembleCouplingTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const dt,
                                                              DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  using TAG = WellSolverBase::SubRegionTag;
  using COFFSET = CompositionalMultiphaseWell::ColOffset;
  using ROFFSET = CompositionalMultiphaseWell::RowOffset;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  localIndex constexpr MAX_NUM_COMP = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr MAX_NUM_DOF  = MAX_NUM_COMP + 1;

  localIndex const NC      = m_wellSolver->numFluidComponents();
  localIndex const resNDOF = m_wellSolver->numDofPerResElement();

  string const resDofKey = dofManager.getKey( m_wellSolver->resElementDofName() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const resDofNumberAccessor =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );
  ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber =
    resDofNumberAccessor.toNestedViewConst();
  globalIndex const rankOffset = dofManager.rankOffset();

  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion const & subRegion )
  {

    PerforationData const * const perforationData = subRegion.getPerforationData();

    // get the degrees of freedom
    string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get well variables on perforations
    arrayView2d< real64 const > const & compPerfRate =
      perforationData->getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::compPerforationRateString() );
    arrayView3d< real64 const > const & dCompPerfRate_dPres =
      perforationData->getReference< array3d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dPresString() );
    arrayView4d< real64 const > const & dCompPerfRate_dComp =
      perforationData->getReference< array4d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dCompString() );

    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

    // loop over the perforations and add the rates to the residual and jacobian
    forAll< parallelDevicePolicy<> >( perforationData->size(), [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
    {
      // local working variables and arrays
      stackArray1d< localIndex, 2 * MAX_NUM_COMP > eqnRowIndices( 2 * NC );
      stackArray1d< globalIndex, 2 * MAX_NUM_DOF > dofColIndices( 2 * resNDOF );

      stackArray1d< real64, 2 * MAX_NUM_COMP > localPerf( 2 * NC );
      stackArray2d< real64, 2 * MAX_NUM_COMP * 2 * MAX_NUM_DOF > localPerfJacobian( 2 * NC, 2 * resNDOF );

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const resOffset = resDofNumber[er][esr][ei];
      globalIndex const wellElemOffset = wellElemDofNumber[iwelem];

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        eqnRowIndices[TAG::RES *NC + ic]  = LvArray::integerConversion< localIndex >( resOffset - rankOffset ) + ic;
        eqnRowIndices[TAG::WELL *NC + ic] = LvArray::integerConversion< localIndex >( wellElemOffset - rankOffset ) + ROFFSET::MASSBAL + ic;
      }
      for( localIndex jdof = 0; jdof < resNDOF; ++jdof )
      {
        dofColIndices[TAG::RES *resNDOF + jdof]  = resOffset + jdof;
        dofColIndices[TAG::WELL *resNDOF + jdof] = wellElemOffset + COFFSET::DPRES + jdof;
      }

      // populate local flux vector and derivatives
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        localPerf[TAG::RES *NC + ic]  =  dt * compPerfRate[iperf][ic];
        localPerf[TAG::WELL *NC + ic] = -dt * compPerfRate[iperf][ic];

        for( localIndex ke = 0; ke < 2; ++ke )
        {
          localIndex const localDofIndexPres = ke * resNDOF;
          localPerfJacobian[TAG::RES *NC + ic][localDofIndexPres]  =  dt * dCompPerfRate_dPres[iperf][ke][ic];
          localPerfJacobian[TAG::WELL *NC + ic][localDofIndexPres] = -dt * dCompPerfRate_dPres[iperf][ke][ic];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
            localPerfJacobian[TAG::RES *NC + ic][localDofIndexComp]  =  dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
            localPerfJacobian[TAG::WELL *NC + ic][localDofIndexComp] = -dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
          }
        }
      }

      for( localIndex i = 0; i < localPerf.size(); ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                            dofColIndices.data(),
                                                                            localPerfJacobian[i].dataIfContiguous(),
                                                                            2 * resNDOF );
          atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localPerf[i] );
        }
      }
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseReservoir, string const &, Group * const )

} /* namespace geosx */
