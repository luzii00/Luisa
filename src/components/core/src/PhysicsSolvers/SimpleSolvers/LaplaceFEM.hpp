/*
 * NewtonianMechanics.hpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#ifndef SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_
#define SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_

#include "PhysicsSolvers/SolverBase.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"


struct stabledt
{
  double m_maxdt;
};

namespace ML_Epetra
{ class MultiLevelPreconditioner; }

namespace geosx
{
namespace dataRepository
{
class ManagedGroup;
}
class BoundaryConditionBase;
class FiniteElementBase;
class DomainPartition;

class LaplaceFEM : public SolverBase
{
public:
  LaplaceFEM( const std::string& name,
                                ManagedGroup * const parent );


  virtual ~LaplaceFEM();

  static string CatalogName() { return "LaplaceFEM"; }

  virtual void FillDocumentationNode() override final;

  virtual void FillOtherDocumentationNodes( dataRepository::ManagedGroup * const group ) override final;

  virtual void InitializePreSubGroups( dataRepository::ManagedGroup * const problemManager ) override final;

  virtual void ReadXML_PostProcess() override final;

  virtual void TimeStep( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         dataRepository::ManagedGroup * domain ) override;

  void TimeStepExplicit( real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         DomainPartition * domain );

  void TimeStepQuasiStatic( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition& domain );

  real64 TimeStepImplicit( real64 const & time_n,
                           real64 const & dt,
                           integer const cycleNumber,
                           DomainPartition * const domain );

  void TimeStepImplicitSetup( real64 const& time_n,
                              real64 const& dt,
                              DomainPartition * const domain );

  void TimeStepImplicitComplete( real64 const & time,
                                 real64 const & dt,
                                 DomainPartition * const domain );

  void SetupSystem ( DomainPartition * const domain,
                     systemSolverInterface::EpetraBlockSystem * const blockSystem );

  void SetSparsityPattern( DomainPartition const * const domain,
                           Epetra_FECrsGraph * const sparsity );

  void SetNumRowsAndTrilinosIndices( ManagedGroup * const domain,
                                     localIndex & numLocalRows,
                                     localIndex & numGlobalRows,
                                     localIndex_array& localIndices,
                                     localIndex offset );

  void SetupMLPreconditioner( DomainPartition const & domain,
                              ML_Epetra::MultiLevelPreconditioner* MLPrec );


  real64 Assemble ( DomainPartition * const domain,
                    systemSolverInterface::EpetraBlockSystem * const blockSystem,
                    real64 const time,
                    real64 const dt );

  realT CalculateElementResidualAndDerivative( real64 const density,
                                               FiniteElementBase const * const fe,
                                               const Array2dT<R1Tensor>& dNdX,
                                               const realT* const detJ,
                                               R2SymTensor const * const refStress,
                                               array<R1Tensor> const & u,
                                               array<R1Tensor> const & uhat,
                                               array<R1Tensor> const & uhattilde,
                                               array<R1Tensor> const & vtilde,
                                               realT const dt,
                                               Epetra_SerialDenseMatrix& dRdU,
                                               Epetra_SerialDenseVector& R,
                                               real64 c[6][6] );

  void ApplyDirichletBC_implicit( ManagedGroup * object,
                                  BoundaryConditionBase const * const bc,
                                  lSet const & set,
                                  real64 const time_n,
                                  systemSolverInterface::EpetraBlockSystem & blockSystem );

  void ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                            real64 const scalingFactor,
                            localIndex const dofOffset,
                            dataRepository::ManagedGroup * const nodeManager );




  enum class timeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };

  struct viewKeyStruct
  {
    dataRepository::ViewKey trilinosIndex = { "trilinosIndex_LaplaceFEM" };
    dataRepository::ViewKey ghostRank = { "ghostRank" };
    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };
  } viewKeys;

  struct groupKeyStruct
  {
    dataRepository::GroupKey systemSolverParameters = { "SystemSolverParameters" };
  } groupKeys;


  systemSolverInterface::LinearSolverWrapper m_linearSolverWrapper;
  systemSolverInterface::EpetraBlockSystem m_linearSystem;

  SystemSolverParameters * getSystemSolverParameters() {return this->GetGroup<SystemSolverParameters>(groupKeys.systemSolverParameters); }

private:

  stabledt m_stabledt;
  timeIntegrationOption m_timeIntegrationOption;
  LaplaceFEM();

};


} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
