/*
 * CommunicationTools.hpp
 *
 *  Created on: Jan 6, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_
#define SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_

#include "common/DataTypes.hpp"
#include "mpi.h"
#include<set>
namespace geosx
{


class ObjectManagerBase;
class NeighborCommunicator;
class MeshLevel;

class CommunicationTools
{
public:
  CommunicationTools();
  ~CommunicationTools();

  static void AssignGlobalIndices( ObjectManagerBase & object,
                                   ObjectManagerBase const & compositionObject,
                                   array<NeighborCommunicator> & neighbors );

  static void FindGhosts( MeshLevel * const meshLevel,
                          array<NeighborCommunicator> & neighbors );

  static int MPI_Size( MPI_Comm const & comm );
  static int MPI_Rank( MPI_Comm const & comm );

  static std::set<int> & getFreeCommIDs();
  static int reserveCommID();
  static void releaseCommID( int & ID );

  static void FindMatchedPartitionBoundaryObjects( ObjectManagerBase * const group,
                                            array<NeighborCommunicator> & allNeighbors );

  static void SynchronizeFields( const std::map<string, array<string> >& fieldNames,
                                 MeshLevel * const mesh,
                                 array<NeighborCommunicator> & allNeighbors );

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMUNICATIONTOOLS_HPP_ */
