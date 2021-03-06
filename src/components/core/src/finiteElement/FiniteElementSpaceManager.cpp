/*
 * FiniteElementSpaceManager.cpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#include "FiniteElementSpaceManager.hpp"
#include "FiniteElementSpace.hpp"

namespace geosx
{
using namespace dataRepository;

FiniteElementSpaceManager::FiniteElementSpaceManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{}

FiniteElementSpaceManager::~FiniteElementSpaceManager()
{
  // TODO Auto-generated destructor stub
}

void FiniteElementSpaceManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("FiniteElementSpace");
  docNode->setSchemaType("UniqueNode");
}

void FiniteElementSpaceManager::CreateChild( string const & childKey, string const & childName )
{
  // These objects should probably not be registered on managed group...
  std::unique_ptr<ManagedGroup> fem = ManagedGroup::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup( childName, std::move(fem) );
}



} /* namespace geosx */
