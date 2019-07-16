

#ifndef FINITE_ELEMENT_SHAPE_KERNEL
#define FINITE_ELEMENT_SHAPE_KERNEL

#include "common/DataTypes.hpp"

namespace geosx
{

class FiniteElementShapeKernel
{
  constexpr localIndex numNodes = 8;
  constexpr real64 weight = 8.0 / numNodes;
  constexpr real64 quadratureFactor = 1.0 / sqrt(3.0);
  constexpr real64 parentCoords[3][numNodes] = { { -1,  1, -1,  1, -1,  1, -1,  1 },
                                                 { -1, -1,  1,  1, -1, -1,  1,  1 },
                                                 { -1, -1, -1, -1,  1,  1,  1,  1 } };


  inline real64 parentShapeFunctionValue( localIndex const a, real64 const xi0, real64 const xi1, real64 const x2 )
  {
    return 0.125 * ( 1 + xi0*parentCoords[0][a] ) * ( 1 + xi1*parentCoords[1][a] ) * ( 1 + xi2*parentCoords[2][a] ) ;
  }

  inline void parentShapeFunctionDerivatives( localIndex const a,
                                              real64 const xi0,
                                              real64 const xi1,
                                              real64 const x2,
                                              real64 dNdXi[3] )
  {
    dNdXi[0] = 0.125 * parentCoords[0][a] * ( 1 + xi1*parentCoords[0][a] ) * ( 1 + xi2*parentCoords[0][a] ) ;
    dNdXi[1] = 0.125 * ( 1 + xi0*parentCoords[0][a] ) * parentCoords[1][a] * ( 1 + xi2*parentCoords[2][a] ) ;
    dNdXi[2] = 0.125 * ( 1 + xi0*parentCoords[0][a] ) * ( 1 + xi1*parentCoords[1][a] ) * parentCoords[2][a] ;
  }



  inline real64 shapeFunctionDerivatives( real64 const xi0,
                                          real64 const xi1,
                                          real64 const xi2,
                                          real64 const X[3][numNodes],
                                          real64 dNdXi[3][numNodes] )
  {
    real64 J[3][3] = {0};

    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      real64 dNdXi[3];
      parentShapeFunctionDerivatives( a, xi0, xi1, xi2, dNdXi );
      J[0][0] += X[0][a] * dNdXi[0];
      J[0][1] += X[0][a] * dNdXi[1];
      J[0][2] += X[0][a] * dNdXi[2];

      J[1][0] += X[1][a] * dNdXi[0];
      J[1][1] += X[1][a] * dNdXi[1];
      J[1][2] += X[1][a] * dNdXi[2];

      J[2][0] += X[2][a] * dNdXi[0];
      J[2][1] += X[2][a] * dNdXi[1];
      J[2][2] += X[2][a] * dNdXi[2];
    }

    real64 detJ = inverse( J );

    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      real64 dNdXi[3];
      parentShapeFunctionDerivatives( a, xi0, xi1, xi2, dNdXi );
      dNdXi[0][a] = J[0][0] * dNdXi[0] + J[1][0] * dNdXi[1] + J[2][0] * dNdXi[2];
      dNdXi[1][a] = J[0][1] * dNdXi[0] + J[1][1] * dNdXi[1] + J[2][1] * dNdXi[2];
      dNdXi[2][a] = J[0][2] * dNdXi[0] + J[1][2] * dNdXi[1] + J[2][2] * dNdXi[2];
    }

    return detJ;
  }

  inline real64 shapeFunctionDerivatives( localIndex q,
                                          real64 const X[3][numNodes],
                                          real64 dNdXi[3][numNodes] )
  {
    return shapeFunctionDerivatives( quadratureFactor*parentCoords[0][q],
                                     quadratureFactor*parentCoords[1][q],
                                     quadratureFactor*parentCoords[2][q],
                                     X,
                                     dNdXi );
  }

  inline real64 inverse( real64 J[3][3] )
  {
    real64 const o1 = J[1][1]*J[2][2] - J[1][2]*J[2][1];
    real64 const o2 = J[0][2]*J[2][1] - J[0][1]*J[2][2];
    real64 const o3 = J[0][1]*J[1][2] - J[0][2]*J[1][1];
    real64 const o4 = J[1][2]*J[2][0] - J[1][0]*J[2][2];
    real64 const o5 = J[0][0]*J[2][2] - J[0][2]*J[2][0];
    real64 const o6 = J[0][2]*J[1][0] - J[0][0]*J[1][2];
    real64 const o7 = J[1][0]*J[2][1] - J[1][1]*J[2][0];
    real64 const o8 = J[0][1]*J[2][0] - J[0][0]*J[2][1];
    real64 const o9 = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    real64 const detJ = J[0][0] * o1 + J[1][0] * o2 + J[2][0] * o3 ;
    real64 const o10 = 1.0 / detJ;

    J[0][0] = o1 * o10;
    J[0][1] = o2 * o10;
    J[0][2] = o3 * o10;
    J[1][0] = o4 * o10;
    J[1][1] = o5 * o10;
    J[1][2] = o6 * o10;
    J[2][0] = o7 * o10;
    J[2][1] = o8 * o10;
    J[2][2] = o9 * o10;
  }

};

}

#endif //FINITE_ELEMENT_SHAPE_KERNEL