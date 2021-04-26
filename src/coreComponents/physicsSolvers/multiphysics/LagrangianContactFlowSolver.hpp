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
 * @file LagrangianContactFlowSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTFLOWSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTFLOWSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactSolver.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{

class LagrangianContactSolver;
class FlowSolverBase;

class LagrangianContactFlowSolver : public SolverBase
{
public:
  LagrangianContactFlowSolver( const std::string & name,
                               Group * const parent );

  ~LagrangianContactFlowSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "LagrangianContactWithFlow"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  setNextDt( real64 const & currentDt,
             real64 & nextDt ) override;


  virtual real64
  explicitStep( real64 const & time_n,
                real64 const & dt,
                integer const cycleNumber,
                DomainPartition & domain ) override;

  virtual real64
  nonlinearImplicitStep( real64 const & time_n,
                         real64 const & dt,
                         integer const cycleNumber,
                         DomainPartition & domain ) override;

  virtual bool
  lineSearch( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain,
              DofManager const & dofManager,
              CRSMatrixView< real64, globalIndex const > const & localMatrix,
              arrayView1d< real64 > const & localRhs,
              arrayView1d< real64 const > const & localSolution,
              real64 const scaleFactor,
              real64 & lastResidual ) override;

  void updateOpeningForFlow( DomainPartition & domain );

  void assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs );

  void assembleStabilization( DomainPartition const & domain,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs );

  real64 splitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * contactSolverNameString() { return "contactSolverName"; }
    constexpr static char const * flowSolverNameString() { return "flowSolverName"; }
    constexpr static char const * stabilizationNameString() { return "stabilizationName"; }

    constexpr static char const * defaultConductivityString() { return "defaultConductivity"; }
  } LagrangianContactFlowSolverViewKeys;

protected:
  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

  /**
   * @Brief add the nnz induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLenghts the nnz in each row
   */
  void
  addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                  DofManager const & dofManager,
                                  arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void
  addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      SparsityPatternView< globalIndex > const & pattern ) const;

private:

  string m_contactSolverName;
  string m_flowSolverName;

  LagrangianContactSolver * m_contactSolver;
  FlowSolverBase * m_flowSolver;

  string m_stabilizationName;

  real64 m_defaultConductivity;

  integer m_activeSetIter = 0;

  string const m_tractionKey = LagrangianContactSolver::viewKeyStruct::tractionString();
  string const m_fractureStateKey = LagrangianContactSolver::viewKeyStruct::fractureStateString();
  string const m_localJumpKey = LagrangianContactSolver::viewKeyStruct::localJumpString();
  string const m_pressureKey = FlowSolverBase::viewKeyStruct::pressureString();
  string const m_deltaPressureKey = FlowSolverBase::viewKeyStruct::deltaPressureString();

  real64 m_initialResidual[4] = {0.0, 0.0, 0.0, 0.0};

  real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();

  void createPreconditioner( DomainPartition const & domain );

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_LAGRANGIANCONTACTFLOWSOLVER_HPP_ */