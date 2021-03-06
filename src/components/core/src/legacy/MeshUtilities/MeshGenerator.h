/*
 * MeshGenerator.h
 *
 *  Created on: Nov 19, 2012
 *      Author: settgast1
 */

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include "IO/ticpp/TinyXMLParser.h"
#include "MPI_Communications/SpatialPartition.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "../IO/ticpp/HierarchicalDataNode.h.old"



class MeshGenerator
{
public:
  MeshGenerator();
  virtual ~MeshGenerator();

  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  void GenerateElementRegions( PhysicalDomainT * domain );

  void GenerateMesh( SpatialPartition& partition,
                     PhysicalDomainT * domain );

  void GenerateNodesets( TICPP::HierarchicalDataNode& hdn,
                         NodeManager * nodeManager );

  void GetElemToNodesRelationInBox ( const std::string& elementType,
                                     const int index[],
                                     const int& iEle,
                                     int nodeIDInBox[],
                                     const int size);

  void RemapMesh ( PhysicalDomainT * domain );


  int m_delayMeshDeformation;

private:


  int m_dim;
//  array<real64> m_vertices[3];
//  array<integer> m_nElems[3];
//  array<real64> m_nElemScaling[3];
//  array<real64> m_nElemBias[3];
//
//  array<string> m_regionNames;

  realT m_min[3]; // Minimum extent of mesh dimensions
  realT m_max[3]; // Maximum extent of mesh dimensions

  //int m_numElems[3];
  array<integer> m_firstElemIndexForBlock[3];
  array<integer> m_lastElemIndexForBlock[3];



//  realT m_wExtensionMin[3];
//  realT m_wExtensionMax[3];
//  int m_nExtensionLayersMin[3];
//  int m_nExtensionLayersMax[3];
//  realT m_extendedMin[3];
//  realT m_extendedMax[3]; // This is the domain size after we apply n layers
// of elements which are of the same size as the core elements.  We will move
// these nodes to where they should be later when we finish the meshing.
  int m_numElemsTotal[3];
//  realT m_commonRatioMin[3];
//  realT m_commonRatioMax[3];



  array<string> m_elementType;

  array<integer> m_numElePerBox;

  int m_trianglePattern;   // In pattern 0, half nodes have 4 edges and the
                           // other half have 8; for Pattern 1, every node has
                           // 6.

  realT m_fPerturb;
  int m_randSeed;

  int m_mapToRadial = 0;
  int meshAxis;
  float meshTheta;
  float meshPhi;
  float meshRout;
  float meshRact;

  realT m_skewAngle;
  R1Tensor m_skewCenter;

  std::string m_meshDx, m_meshDy, m_meshDz;

  inline globalIndex NodeGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) + index[1]*(m_numElemsTotal[2]+1) + index[2];
    return rval;
  }

  inline globalIndex ElemGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*m_numElemsTotal[1]*m_numElemsTotal[2] + index[1]*m_numElemsTotal[2] + index[2];
    return rval;
  }

  inline R1Tensor NodePosition( const int a[3], int trianglePattern )
  {
    R1Tensor X;
    realT xInterval(0);

    int xPosIndex = 0;
    if (trianglePattern == 1)
    {
      int startingIndex = 0;
      int endingIndex = 0;
      unsigned int block = 0;
      for( block=0 ; block<m_nElems[0].size() ; ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[0][block];
      }
      xPosIndex = endingIndex;
    }

    for( int i=0 ; i<3 ; ++i )
    {

      int startingIndex = 0;
      int endingIndex = 0;
      unsigned int block = 0;
      for( block=0 ; block<m_nElems[i].size() ; ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[i][block];
        if( a[i]>=startingIndex && a[i]<=endingIndex )
        {
          break;
        }
      }
      realT min = m_vertices[i][block];
      realT max = m_vertices[i][block+1];


      X[i] = min + (max-min) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

      if ((m_nElemBias[i][block] != 0) & (m_nElems[i][block]>1))
      {
        if (fabs(m_nElemBias[i][block]) >= 1)
        {
          throw GPException("Mesh bias must between -1 and 1!");
        }

        realT len = max -  min;
        realT xmean = len / m_nElems[i][block];
        realT x0 = xmean * double( a[i] - startingIndex );
        realT chi = m_nElemBias[i][block]/(xmean/len - 1.0);
        realT dx = -x0*chi + x0*x0*chi/len;
        X[i] += dx;
      }

      // This is for creating regular triangle pattern
      if (i==0) xInterval = (max-min) / m_nElems[i][block];
      if (trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 && a[0] != xPosIndex)
        X[0] -= xInterval * 0.5;
    }

    return X;
  }

  inline R1Tensor ElemCenterPosition( const int k[3] )
  {
    R1Tensor X;

    for( int i=0 ; i<3 ; ++i )
    {
      X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }

    return X;
  }

public:
  inline bool isRadial()
  {
    bool rval = (m_mapToRadial > 0);
    return rval;
  }

};

#endif /* MESHGENERATOR_H_ */
