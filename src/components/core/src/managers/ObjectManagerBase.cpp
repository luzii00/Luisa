/*
 * ObjectManagerBase.cpp
 *
 *  Created on: Sep 15, 2016
 *      Author: settgast1
 */

#include "ObjectManagerBase.hpp"

namespace geosx
{
using namespace dataRepository;

ObjectManagerBase::ObjectManagerBase( std::string const & name,
                                      ManagedGroup * const parent ):
  ManagedGroup(name,parent),
  m_sets(keys::sets,this),
  m_localToGlobalMap(),
  m_globalToLocalMap()
//  m_localToGlobalMap( RegisterViewWrapper< globalIndex_array >("localToGlobal")->reference() ),
//  m_globalToLocalMap( RegisterViewWrapper< map<globalIndex,localIndex> >("globalToLocal")->reference() )
{

  RegisterViewWrapper(viewKeyStruct::localToGlobalMapString, &m_localToGlobalMap, false );
  RegisterViewWrapper(viewKeyStruct::globalToLocalMapString, &m_globalToLocalMap, false );

  RegisterGroup( keys::sets, &m_sets, false );
  RegisterViewWrapper(viewKeyStruct::isExternalString, &m_isExternal, false );
  RegisterViewWrapper(viewKeyStruct::ghostRankString, &m_ghostRank, false );
  this->RegisterGroup(groupKeys.neighborData);
}
//ObjectManagerBase::ObjectManagerBase( std::string const & name,
//                                      ManagedGroup * const parent,
//                                      cxx_utilities::DocumentationNode *
// docNode ):
//    ManagedGroup(name,parent,docNode),
//    m_localToGlobalMap( RegisterViewWrapper< globalIndex_array
// >("localToGlobal")->reference() )
//{
//
//
//  this->RegisterGroup<ManagedGroup>("Sets");
//  this->RegisterViewWrapper< array<integer> >("isExternal");
//}


ObjectManagerBase::~ObjectManagerBase()
{}



ObjectManagerBase::CatalogInterface::CatalogType& ObjectManagerBase::GetCatalog()
{
  static ObjectManagerBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void ObjectManagerBase::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( viewKeys.ghostRank.Key(),
                              viewKeys.ghostRank.Key(),
                              -1,
                              "integer_array",
                              "integer_array",
                              "Array that indicates whether or not an index is a ghost. ",
                              "Array that indicates whether or not an index is a ghost. "
                              "If it is not a ghost the value will be -1. If it "
                              "is a ghost, then the value will be the owning rank.",
                              "-1",
                              "",
                              1,
                              0,
                              0 );

  docNode->AllocateChildNode( viewKeys.domainBoundaryIndicator.Key(),
                              viewKeys.domainBoundaryIndicator.Key(),
                              -1,
                              "integer_array",
                              "integer_array",
                              "List containing the element regions of the faces",
                              "List containing the element regions of the faces",
                              "0",
                              "",
                              1,
                              0,
                              0 );

//  docNode->AllocateChildNode( viewKeys.globalToLocalMap.Key(),
//                              viewKeys.globalToLocalMap.Key(),
//                              -1,
//                              "localIndex_map",
//                              "localIndex_map",
//                              "Map from globalIndex to localIndex. ",
//                              "Map from globalIndex to localIndex. ",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );

//
//  docNode->AllocateChildNode( viewKeys.localToGlobalMap.Key(),
//                              viewKeys.localToGlobalMap.Key(),
//                              -1,
//                              "globalIndex_array",
//                              "globalIndex_array",
//                              "Array that maps from localIndex to globalIndex. ",
//                              "Array that maps from localIndex to globalIndex. ",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );

}

void ObjectManagerBase::InitializePostSubGroups( ManagedGroup * const )
{
  m_ghostRank=-1;
}


void ObjectManagerBase::ConstructSetFromSetAndMap( const lSet& inputSet,
                                                   const lArray2d& map,
                                                   const std::string& newSetName )
{

  ManagedGroup * sets = GetGroup(std::string("Sets"));
  lSet& newset = sets->RegisterViewWrapper<lSet>(newSetName)->reference();
  newset.clear();

  localIndex mapSize = map.size(1);
  for( localIndex ka=0 ; ka<size() ; ++ka )
  {
    arrayView1d<localIndex const> const sublist = map[ka];
    int addToSet = 0;
    for( int a=0 ; a<mapSize ; ++a )
    {
      if( inputSet.count( sublist[a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == mapSize )
    {
      newset.insert( ka );
    }
  }
}

void ObjectManagerBase::ConstructSetFromSetAndMap( const lSet& inputSet,
                                                   const array<localIndex_array>& map,
                                                   const std::string& newSetName )
{

  ManagedGroup * sets = GetGroup(std::string("Sets"));
  lSet& newset = sets->RegisterViewWrapper<lSet>(newSetName)->reference();
  newset.clear();

  for( localIndex ka=0 ; ka<size() ; ++ka )
  {
    localIndex addToSet = 0;
    localIndex mapSize = map[ka].size();
    for( int a=0 ; a<mapSize ; ++a )
    {
      if( inputSet.count( map[ka][a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == mapSize )
    {
      newset.insert( ka );
    }
  }
}

void ObjectManagerBase::ConstructLocalListOfBoundaryObjects( localIndex_array& objectList ) const
{
  const array<integer>& isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( k );
    }
  }
}

void ObjectManagerBase::ConstructGlobalListOfBoundaryObjects( globalIndex_array& objectList ) const
{
  const array<integer>& isDomainBoundary = this->getReference<integer_array>(viewKeys.domainBoundaryIndicator);
  for( localIndex k=0 ; k<size() ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( this->m_localToGlobalMap[k] );
    }
  }
  std::sort( objectList.begin(), objectList.end() );
}

void ObjectManagerBase::ConstructGlobalToLocalMap()
{
  m_globalToLocalMap.clear();
  for( localIndex k=0 ; k<size() ; ++k )
  {
    m_globalToLocalMap[m_localToGlobalMap[k]] = k;
  }

}






int ObjectManagerBase::PackSize( array<string> const & wrapperNames,
                            localIndex_array const & packList,
                            integer const recursive ) const
{
  int packedSize = 0;
  buffer_unit_type * junk;
  packedSize += this->PackPrivate<false>( junk,
                                          wrapperNames,
                                          packList,
                                          recursive );

  return packedSize;
}




int ObjectManagerBase::Pack( buffer_unit_type * & buffer,
                             array<string> const & wrapperNames,
                             localIndex_array const & packList,
                             integer const recursive ) const
{
  int packedSize = 0;

  packedSize += this->PackPrivate<true>( buffer,
                                          wrapperNames,
                                          packList,
                                          recursive );

  return packedSize;
}


int ObjectManagerBase::UnpackGlobalMaps( buffer_unit_type const *& buffer,
                                         localIndex_array & packList,
                                         integer const recursive )
{
  int unpackedSize = 0;
  string groupName;
  unpackedSize += CommBufferOps::Unpack( buffer, groupName );
  GEOS_ASSERT( groupName==this->getName(), "ObjectManagerBase::Unpack(): group names do not match")

  string localToGlobalString;
  unpackedSize += CommBufferOps::Unpack( buffer, localToGlobalString);
  GEOS_ASSERT( localToGlobalString==viewKeyStruct::localToGlobalMapString, "ObjectManagerBase::Unpack(): label incorrect")

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  int sendingRank;
  unpackedSize += CommBufferOps::Unpack( buffer, sendingRank );

  int numUnpackedIndices;
  unpackedSize += CommBufferOps::Unpack( buffer, numUnpackedIndices );

  localIndex_array unpackedLocalIndices;
  unpackedLocalIndices.resize(numUnpackedIndices);

  globalIndex_array globalIndices;
  unpackedSize += CommBufferOps::Unpack( buffer, globalIndices );
  localIndex numNewIndices = 0;
  globalIndex_array newGlobalIndices;
  localIndex const oldSize = this->size();
  for( localIndex a=0 ; a<numUnpackedIndices ; ++a )
  {
    // check to see if the object already exists by checking for the global
    // index in m_globalToLocalMap. If it doesn't, then add the object
    map<globalIndex,localIndex>::iterator iterG2L = m_globalToLocalMap.find(globalIndices[a]);
    if( iterG2L == m_globalToLocalMap.end() )
    {
      // object does not exist on this domain
      const localIndex newLocalIndex = oldSize + numNewIndices;

      // add the global index of the new object to the globalToLocal map
      m_globalToLocalMap[globalIndices[a]] = newLocalIndex;

      unpackedLocalIndices(a) = newLocalIndex;

      newGlobalIndices.push_back( globalIndices[a] );

      ++numNewIndices;

      GEOS_ASSERT( packList.size()==0,
                   "ObjectManagerBase::Unpack(): packList specified, "
                   "but a new globalIndex is unpacked")
    }
    else
    {
      // object already exists on this domain
      // get the local index of the node
      localIndex b = iterG2L->second;
      unpackedLocalIndices(a) = b;
      if( ( sendingRank < rank && m_ghostRank[b] == -1) || ( sendingRank < m_ghostRank[b] ) )
      {
        m_ghostRank[b] = sendingRank;
      }
    }
  }
  newGlobalIndices.resize(numNewIndices);

  // figure out new size of object container, and resize it
  const localIndex newSize = oldSize + numNewIndices;
  this->resize( newSize );

  // add the new indices to the maps.
  for( int a=0 ; a<numNewIndices ; ++a )
  {
    localIndex const b = oldSize + a;
    m_localToGlobalMap[b] = newGlobalIndices(a);
    m_ghostRank[b] = sendingRank;
  }



//    for( localIndex a=0 ; a<numUnpackedIndices ; ++a )
//    {
//      localIndex const b = m_globalToLocalMap.at(globalIndices[a]);
//      if( m_ghostRank[b]==-1 )
//      {
//        if( sendingRank < rank )
//        {
//          m_ghostRank[b] = sendingRank;
//        }
//      }
//      else
//      {
//        if( sendingRank < m_ghostRank[b] )
//        {
//          m_ghostRank[b] = sendingRank;
//        }
//
//      }
//    }

    packList = unpackedLocalIndices;


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += CommBufferOps::Unpack( buffer, subGroups );
    GEOS_ASSERT( subGroups=="SubGroups", "ManagedGroup::Unpack(): group names do not match")

    decltype( this->GetSubGroups().size()) numSubGroups;
    unpackedSize += CommBufferOps::Unpack( buffer, numSubGroups );
    GEOS_ASSERT( numSubGroups==this->GetSubGroups().size(), "ManagedGroup::Unpack(): incorrect number of subGroups")

    for( auto const & index : this->GetSubGroups() )
    {
      string subGroupName;
      unpackedSize += CommBufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->GetGroup<ObjectManagerBase>(subGroupName)->UnpackGlobalMaps(buffer,packList,recursive);
    }
  }

  return unpackedSize;
}

int ObjectManagerBase::Unpack( buffer_unit_type const *& buffer,
                               localIndex_array & packList,
                               integer const recursive )
{
  int unpackedSize = 0;
  string groupName;
  unpackedSize += CommBufferOps::Unpack( buffer, groupName );
  GEOS_ASSERT( groupName==this->getName(), "ObjectManagerBase::Unpack(): group names do not match")

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  int sendingRank;
  unpackedSize += CommBufferOps::Unpack( buffer, sendingRank );

  int numUnpackedIndices;
  unpackedSize += CommBufferOps::Unpack( buffer, numUnpackedIndices );

//  integer readIncludeGlobalIndices;
//  unpackedSize += CommBufferOps::Unpack( buffer, readIncludeGlobalIndices);

  string wrappersLabel;
  unpackedSize += CommBufferOps::Unpack( buffer, wrappersLabel);
  GEOS_ASSERT( wrappersLabel=="Wrappers", "ObjectManagerBase::Unpack(): wrapper label incorrect")


//  if( readIncludeGlobalIndices )
//  {
//    localIndex_array unpackedLocalIndices;
//    unpackedLocalIndices.resize(numUnpackedIndices);
//
//  }

  int numWrappers;
  unpackedSize += CommBufferOps::Unpack( buffer, numWrappers);
  for( localIndex a=0 ; a<numWrappers ; ++a )
  {
    string wrapperName;
    unpackedSize += CommBufferOps::Unpack( buffer, wrapperName );
    ViewWrapperBase * const wrapper = this->getWrapperBase(wrapperName);
    wrapper->Unpack(buffer,packList);
  }


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += CommBufferOps::Unpack( buffer, subGroups );
    GEOS_ASSERT( subGroups=="SubGroups", "ManagedGroup::Unpack(): group names do not match")

    decltype( this->GetSubGroups().size()) numSubGroups;
    unpackedSize += CommBufferOps::Unpack( buffer, numSubGroups );
    GEOS_ASSERT( numSubGroups==this->GetSubGroups().size(), "ManagedGroup::Unpack(): incorrect number of subGroups")

    for( auto const & index : this->GetSubGroups() )
    {
      string subGroupName;
      unpackedSize += CommBufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->GetGroup(subGroupName)->Unpack(buffer,packList,recursive);
    }
  }

  return unpackedSize;
}

template< bool DOPACK >
int ObjectManagerBase::PackPrivate( buffer_unit_type * & buffer,
                                    array<string> const & wrapperNames,
                                    localIndex_array const & packList,
                                    integer const recursive ) const
{
  int packedSize = 0;
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, this->getName() );

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, rank );


  int const numPackedIndices = packList.size()==0 ? this->size() : packList.size();
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, numPackedIndices );

  packedSize += CommBufferOps::Pack<DOPACK>( buffer, string("Wrappers") );


  array<string> wrapperNamesForPacking;
  if( wrapperNames.size()==0 )
  {
    set<localIndex> exclusionList;
    ViewPackingExclusionList(exclusionList);
    wrapperNamesForPacking.resize( this->wrappers().size() );
    localIndex count = 0;
    for( localIndex k=0 ; k<this->wrappers().size() ; ++k )
    {
      if( exclusionList.count(k) == 0)
      {
        wrapperNamesForPacking[count++] = wrappers().values()[k].first;
      }
    }
    wrapperNamesForPacking.resize(count);
  }
  else
  {
    wrapperNamesForPacking = wrapperNames;
  }

  packedSize += CommBufferOps::Pack<DOPACK,int>( buffer, wrapperNamesForPacking.size() );
  for( auto const & wrapperName : wrapperNamesForPacking )
  {
    dataRepository::ViewWrapperBase const * const wrapper = this->getWrapperBase(wrapperName);
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, wrapperName );
    if( packList.empty() )
    {
      if(DOPACK)
      {
        packedSize += wrapper->Pack( buffer );
      }
      else
      {
        packedSize += wrapper->PackSize();
      }
    }
    else
    {
      if(DOPACK)
      {
        packedSize += wrapper->Pack( buffer, packList );
      }
      else
      {
        packedSize += wrapper->PackSize( packList );
      }
    }
  }


  if( recursive > 0 )
  {
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, string("SubGroups") );
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, this->GetSubGroups().size() );
    for( auto const & keyGroupPair : this->GetSubGroups() )
    {
      packedSize += CommBufferOps::Pack<DOPACK>( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->Pack( buffer, wrapperNames, packList, recursive );
    }
  }

  return packedSize;
}
template int ObjectManagerBase::PackPrivate<true>( buffer_unit_type * & ,array<string> const & , localIndex_array const &, integer const ) const;
template int ObjectManagerBase::PackPrivate<false>( buffer_unit_type * & ,array<string> const & , localIndex_array const &, integer const ) const;

template< bool DOPACK >
int ObjectManagerBase::PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                                              localIndex_array const & packList,
                                              integer const recursive ) const
{
  int packedSize = 0;


  packedSize += CommBufferOps::Pack<DOPACK>( buffer, this->getName() );

  // this doesn't link without the string()...no idea why.
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::localToGlobalMapString) );

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, rank );

  int const numPackedIndices = packList.size()==0 ? this->size() : packList.size();
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, numPackedIndices );


  globalIndex_array globalIndices;
  globalIndices.resize(numPackedIndices);
  if( packList.size()==0 )
  {
    for( localIndex a=0 ; a<numPackedIndices ; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[a];
    }
  }
  else
  {
    for( localIndex a=0 ; a<numPackedIndices ; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[packList[a]];
    }
  }
  packedSize += CommBufferOps::Pack<DOPACK>( buffer, globalIndices );


//  packedSize += CommBufferOps::Pack<DOPACK>( buffer, recursive );
  if( recursive > 0 )
  {
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, string("SubGroups") );
    packedSize += CommBufferOps::Pack<DOPACK>( buffer, this->GetSubGroups().size() );
    for( auto const & keyGroupPair : this->GetSubGroups() )
    {
      packedSize += CommBufferOps::Pack<DOPACK>( buffer, keyGroupPair.first );
      ObjectManagerBase const * const subObjectManager = ManagedGroup::group_cast<ObjectManagerBase const *>(keyGroupPair.second);
      if( subObjectManager )
      {
        packedSize += subObjectManager->PackGlobalMapsPrivate<DOPACK>( buffer, packList, recursive );
      }
    }
  }

  return packedSize;
}
template int ObjectManagerBase::PackGlobalMapsPrivate<true>( buffer_unit_type * &, localIndex_array const &, integer const ) const;
template int ObjectManagerBase::PackGlobalMapsPrivate<false>( buffer_unit_type * &, localIndex_array const &, integer const ) const;

void ObjectManagerBase::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  exclusionList.insert(this->getWrapperIndex(this->viewKeys.localToGlobalMapString));
  exclusionList.insert(this->getWrapperIndex(this->viewKeys.globalToLocalMapString));
  exclusionList.insert(this->getWrapperIndex(this->viewKeys.ghostRankString));
}


localIndex ObjectManagerBase::GetNumberOfGhosts() const
{
  return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i>-1;} );
}

localIndex ObjectManagerBase::GetNumberOfLocalIndices() const
{
  return std::count_if( m_ghostRank.begin(), m_ghostRank.end(), [](integer i)->localIndex {return i==-1;} );
}

void ObjectManagerBase::SetReceiveLists(  )
{

  map<int,localIndex_array>  receiveIndices;
  for( localIndex a=0 ; a<size() ; ++a )
  {
    if( m_ghostRank[a] > -1 )
    {
      receiveIndices[m_ghostRank[a]].push_back(a);
    }
  }

  for( map<int,localIndex_array>::const_iterator iter=receiveIndices.begin() ; iter!=receiveIndices.end() ; ++iter )
  {
    ManagedGroup * const neighborData = GetGroup(groupKeys.neighborData)->GetGroup( std::to_string( iter->first ) );

    localIndex_array & nodeAdjacencyList = neighborData->getReference<localIndex_array>( viewKeys.ghostsToReceive );
    nodeAdjacencyList = iter->second;
  }

}


} /* namespace geosx */
