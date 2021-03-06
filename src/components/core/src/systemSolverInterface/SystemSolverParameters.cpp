/*
 * SystemSolverParameters.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#include "SystemSolverParameters.hpp"

namespace geosx
{

SystemSolverParameters::SystemSolverParameters( std::string const & name,
                                                ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  RegisterViewWrapper<integer>( viewKeys.verbosity );
  RegisterViewWrapper<real64>( viewKeys.krylovTol );
  RegisterViewWrapper<integer>( viewKeys.numKrylovIter );
  RegisterViewWrapper<integer>( viewKeys.kspace );
  RegisterViewWrapper<real64>( viewKeys.ilut_fill );
  RegisterViewWrapper<real64>( viewKeys.ilut_drop )->data();
  RegisterViewWrapper<integer>( viewKeys.useMLPrecond );
  RegisterViewWrapper<integer>( viewKeys.useInnerSolver );
  RegisterViewWrapper<integer>( viewKeys.scalingOption );
  RegisterViewWrapper<integer>( viewKeys.useBicgstab );
  RegisterViewWrapper<integer>( viewKeys.useDirectSolver );
  RegisterViewWrapper<real64>( viewKeys.KrylovResidualInit );
  RegisterViewWrapper<real64>( viewKeys.KrylovResidualFinal );
  RegisterViewWrapper<integer>( viewKeys.useNewtonSolve );
  RegisterViewWrapper<real64>( viewKeys.newtonTol );
  RegisterViewWrapper<integer>( viewKeys.maxIterNewton );
}

SystemSolverParameters::~SystemSolverParameters()
{
  // TODO Auto-generated destructor stub
}



void SystemSolverParameters::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("SystemSolverParamters");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Parameters for linear/non-linear system solver");


  docNode->AllocateChildNode( viewKeys.verbosity.Key(),
                              viewKeys.verbosity.Key(),
                              -1,
                              "integer",
                              "integer",
                              "verbosity level",
                              "verbosity level",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.solverType.Key(),
                              viewKeys.solverType.Key(),
                              -1,
                              "string",
                              "string",
                              "verbosity level",
                              "verbosity level",
                              "Klu",
                              "",
                              0,
                              1,
                              0 );


  docNode->AllocateChildNode( viewKeys.krylovTol.Key(),
                              viewKeys.krylovTol.Key(),
                              -1,
                              "real64",
                              "real64",
                              "verbosity level",
                              "verbosity level",
                              "1.0e-6",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.numKrylovIter.Key(),
                              viewKeys.numKrylovIter.Key(),
                              -1,
                              "integer",
                              "integer",
                              "verbosity level",
                              "verbosity level",
                              "100",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.useNewtonSolve.Key(),
                              viewKeys.useNewtonSolve.Key(),
                              -1,
                              "integer",
                              "integer",
                              "verbosity level",
                              "verbosity level",
                              "0",
                              "",
                              0,
                              1,
                              0 );



  docNode->AllocateChildNode( viewKeys.ilut_drop.Key(),
                              viewKeys.ilut_drop.Key(),
                              -1,
                              "real64",
                              "real64",
                              "verbosity level",
                              "verbosity level",
                              "0.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.ilut_fill.Key(),
                              viewKeys.ilut_fill.Key(),
                              -1,
                              "real64",
                              "real64",
                              "verbosity level",
                              "verbosity level",
                              "3.0",
                              "",
                              0,
                              1,
                              0 );
//  real64 m_krylovTol;          // Solver convergence criteria
//  integer  m_numKrylovIter;
//  integer  m_kspace;             // Number of krylov vectors before GMRES
// restart
//  real64 m_ilut_fill;          // Fill factor for ILUT preconditioner
//  real64 m_ilut_drop;          // Drop tolerance for ILUT preconditioner
//  bool   m_useMLPrecond;       // Use ML preconditioner
//  bool   m_useInnerSolver;     // Use row scaling
//  integer  m_scalingOption;      // Use row scaling
//  bool   m_useBicgstab;        // Use bicgstab instead of gmres
//  bool   m_useDirectSolver;    // Use Direct solver
//  real64 m_KrylovResidualInit;
//  real64 m_KrylovResidualFinal;
//
//  bool   m_useNewtonSolve;    // Use Newton-Raphson iterations
//  real64 m_newtonTol;
//  integer  m_maxIterNewton;     // Maximum number of Newton-Raphson iterations
}

} /* namespace geosx */
