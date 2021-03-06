/*
 * MeshGeneratorBase.h
 *
 *  Created on: Oct 18, 2017
 *      Author: sherman
 */

#ifndef MESHGENERATORBASE_H_
#define MESHGENERATORBASE_H_

#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

namespace dataRepository
{}

class NodeManager;
class DomainPartition;

class MeshGeneratorBase : public dataRepository::ManagedGroup
{
public:
  explicit MeshGeneratorBase( std::string const & name,
                              ManagedGroup * const parent );

  virtual ~MeshGeneratorBase();

  static string CatalogName() { return "MeshGeneratorBase"; }

  MeshGeneratorBase() = default;
  MeshGeneratorBase( MeshGeneratorBase const & ) = default;
  MeshGeneratorBase( MeshGeneratorBase &&) = default;
  MeshGeneratorBase& operator=( MeshGeneratorBase const & ) = default;
  MeshGeneratorBase& operator=( MeshGeneratorBase&& ) = default;

  virtual void GenerateElementRegions( DomainPartition& domain ) = 0;

  virtual void GenerateMesh( dataRepository::ManagedGroup * const domain ) = 0;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) = 0;

  virtual void GetElemToNodesRelationInBox ( const std::string& elementType,
                                             const int index[],
                                             const int& iEle,
                                             int nodeIDInBox[],
                                             const int size) = 0;

  virtual void RemapMesh ( dataRepository::ManagedGroup * const domain ) = 0;

  int m_delayMeshDeformation;

  using CatalogInterface = cxx_utilities::CatalogInterface< MeshGeneratorBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

};
}

#endif /* MESHGENERATORBASE_H_ */
