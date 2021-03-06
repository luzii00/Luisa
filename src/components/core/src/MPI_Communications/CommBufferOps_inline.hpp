/*
 * packer_inline.hpp
 *
 *  Created on: Jan 5, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMBUFFEROPS_INLINE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMBUFFEROPS_INLINE_HPP_

#include "common/integer_conversion.hpp"
namespace geosx
{
namespace CommBufferOps
{

template< bool DO_PACKING,
          typename T,
          typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer,
                     T const * const var,
                     INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  sizeOfPackedChars += length*sizeof(T);
  static_if( DO_PACKING )
  {
    memcpy( buffer, var, length*sizeof(T) );
    buffer += length*sizeof(T);
  });

  return sizeOfPackedChars;
}


template< bool DO_PACKING,
          typename T,
          typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer,
                     T const * const var,
                     INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ a ] );
  }

  return sizeOfPackedChars;
}

template< typename T,
          typename INDEX_TYPE >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        T * const var,
        INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += Unpack( buffer, length );

//  char * const ptr_var = reinterpret_cast<char *>(var);
  memcpy( var, buffer, length * sizeof(T) );
  buffer += length * sizeof(T);

  return sizeOfUnpackedChars;
}

template< typename T,
          typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer,
        T * const var,
        INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += Unpack( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a] );
  }

  return sizeOfUnpackedChars;
}



template< bool DO_PACKING,
          typename T,
          typename INDEX_TYPE >
localIndex
Pack( char*&  buffer,
      T const * const var,
      INDEX_TYPE const * const indices,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}

template< typename T,
          typename INDEX_TYPE >
localIndex
Unpack( char const *& buffer,
        T * const var,
        INDEX_TYPE const * const indices,
        INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += Unpack( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[ indices[a] ] );
  }

  return sizeOfUnpackedChars;
}





/**
 *
 * @param buffer
 * @param var
 * @return
 */
template< bool DO_PACKING, typename T >
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Pack( char*&  buffer, T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  static_if( DO_PACKING )
  {
    memcpy( buffer, &var, sizeOfPackedChars );
    buffer += sizeOfPackedChars;
  });
  return sizeOfPackedChars;
}

/**
 * @param buffer
 * @param var
 * @return size (in bytes) of unpacked data
 */
template< typename T>
typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
Unpack( char const *& buffer, T & var )
{
  localIndex const sizeOfUnpackedChars = sizeof(T);
  memcpy( &var,
          buffer,
          sizeOfUnpackedChars );
  buffer += sizeOfUnpackedChars;
  return sizeOfUnpackedChars;
}



template< bool DO_PACKING >
localIndex Pack( char*& buffer,  const std::string& var )
{
  string::size_type sizeOfPackedChars = var.size();

  Pack<DO_PACKING>( buffer, sizeOfPackedChars );

  static_if( DO_PACKING )
  {
    for( string::size_type i=0 ; i<sizeOfPackedChars ; ++i )
    {
      *buffer = var[i];
      buffer++;
    }
  });

  sizeOfPackedChars += sizeof( localIndex );
  return integer_conversion<localIndex>(sizeOfPackedChars);
}


template< bool DO_PACKING, typename T >
typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
Pack( char*& buffer, T const & var )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.Data(), var.Length());
  return sizeOfPackedChars;
}

template< typename T >
typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
Unpack( char const *& buffer, T & var )
{
  localIndex sizeOfUnpackedChars = 0;
  real64 * const pVar = var.Data();
  int length = var.Length();
  sizeOfUnpackedChars += Unpack( buffer, pVar, length );
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING,
          typename T >
localIndex Pack( char *& buffer, set<T> const & var )
{

  localIndex sizeOfPackedChars = 0;

//  const localIndex length = integer_conversion<localIndex>(var.size());
//
//  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );
//
//  for( typename std::set<T>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
//  {
//    sizeOfPackedChars += Pack<DO_PACKING>( buffer, *i);
//  }

  return sizeOfPackedChars;
}


template< typename T>
localIndex Unpack( char const *& buffer, set<T> & setToRead )
{
  setToRead.clear();

  localIndex sizeOfUnpackedChars = 0;
//
//  localIndex set_length;
//  sizeOfUnpackedChars += Unpack( buffer, set_length );
//
//
//  for( localIndex a=0 ; a<set_length ; ++a )
//  {
//    T temp;
//    sizeOfUnpackedChars += Unpack( buffer, temp );
//    setToRead.insert( temp );
//  }

  return sizeOfUnpackedChars;
}



//********************************************************************************************************************
template< bool DO_PACKING , typename T_KEY, typename T_VAL >
typename std::enable_if<
is_packable_map< map<T_KEY,T_VAL> >::value,
localIndex >::type
Pack( char*& buffer,
      std::map<T_KEY,T_VAL> const & var )
{

  localIndex sizeOfPackedChars = 0;

  const typename std::map<T_KEY,T_VAL>::size_type length = var.size();

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->first );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->second );
  }

  return sizeOfPackedChars;
}

template< typename T_KEY, typename T_VAL >
typename std::enable_if<
is_packable_map< map<T_KEY,T_VAL> >::value,
localIndex >::type
Unpack( char const *& buffer, std::map<T_KEY,T_VAL>& map )
{
  map.clear();

  localIndex sizeOfUnpackedChars = 0;

  typename std::map<T_KEY,T_VAL>::size_type map_length;
  sizeOfUnpackedChars += Unpack( buffer, map_length );


  for( typename std::map<T_KEY,T_VAL>::size_type a=0 ; a<map_length ; ++a )
  {
    T_KEY key;
    T_VAL value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value );

    map[key] = value;
    map.insert(std::make_pair(key, std::move(value)));
  }

  return sizeOfUnpackedChars;
}

template< bool DO_PACKING ,
          typename T_KEY,
          typename T_VAL,
          typename T_INDICES >
typename std::enable_if<
is_packable_map< map<T_KEY,T_VAL> >::value && is_packable_by_index<T_VAL>::value,
localIndex >::type
 Pack( char*& buffer,
       std::map<T_KEY,T_VAL> const & var,
       T_INDICES const & packIndices )
{

  localIndex sizeOfPackedChars = 0;

  const typename std::map<T_KEY,T_VAL>::size_type length = var.size();

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->first );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->second, packIndices );
  }

  return sizeOfPackedChars;
}

template< typename T_KEY,
          typename T_VAL,
          typename T_INDICES >
typename std::enable_if<
is_packable_map< map<T_KEY,T_VAL> >::value && is_packable_by_index<T_VAL>::value,
localIndex >::type
Unpack( char const *& buffer,
        std::map<T_KEY,T_VAL>& map,
        T_INDICES const & unpackIndices )
{
  map.clear();

  localIndex sizeOfUnpackedChars = 0;

  typename std::map<T_KEY,T_VAL>::size_type map_length;
  sizeOfUnpackedChars += Unpack( buffer, map_length );


  for( typename std::map<T_KEY,T_VAL>::size_type a=0 ; a<map_length ; ++a )
  {
    T_KEY key;
    T_VAL value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value, unpackIndices );

    map[key] = value;
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, int NDIM, typename INDEX_TYPE >
typename std::enable_if<
is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >::value,
localIndex >::type
Pack( char*& buffer,
      multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.dims(), NDIM );

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.strides(), NDIM );

  const localIndex length = var.size();
  localIndex sizeOfPackedArrayChars = length*sizeof(T);

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.data(), length );

  return sizeOfPackedChars;
}


template< typename T, int NDIM, typename INDEX_TYPE >
typename std::enable_if<
is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >::value,
localIndex >::type
Unpack( char const *& buffer,
        multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var )
{
  localIndex sizeOfUnpackedChars = 0;

  int numDimsRead = 0;
  INDEX_TYPE dims[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, dims, numDimsRead );

  if( numDimsRead != NDIM )
  {
    GEOS_ERROR( "error reading dims");
  }
  else
  {
    var.resize( NDIM, dims );
  }

  int numStrideRead = 0;
  INDEX_TYPE strides[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, strides, numStrideRead );

  if( numStrideRead != NDIM )
  {
    GEOS_ERROR( "error reading strides");
  }
  else
  {
    INDEX_TYPE const * const existingStrides = var.strides();
    for( int i=0 ; i<NDIM ; ++i )
    {
      GEOS_ASSERT( strides[i]==existingStrides[i], "CommBufferOps::Unpack(): strides are inconsistent." )
    }
  }

  localIndex numValuesRead;
  sizeOfUnpackedChars += Unpack( buffer, var.data(), numValuesRead );
  if( numValuesRead != var.size() )
  {
    GEOS_ERROR( "error reading data");
  }


  return sizeOfUnpackedChars;

}


template< bool DO_PACKING,
          typename T,
          int NDIM,
          typename T_indices,
          typename INDEX_TYPE >
typename std::enable_if<
is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >::value,
localIndex >::type
Pack( char*& buffer,
      multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var,
      const T_indices& indices )
{


  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.strides(), NDIM );


  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.data(indices[a]), var.strides()[0] );
  }
//  const localIndex length = indices.size();
//  localIndex sizeOfPackedArrayChars = length*sizeof(T);
//
//  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.data(), indices.data(), length );

  return sizeOfPackedChars;
}



template< typename T,
          int NDIM,
          typename T_indices,
          typename INDEX_TYPE >
localIndex Unpack( char const *& buffer,
                   multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var,
                   const T_indices& indices )
{
  localIndex sizeOfUnpackedChars = 0;


  int numStrideRead = 0;
  INDEX_TYPE strides[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, strides, numStrideRead );

  if( numStrideRead != NDIM )
  {
    GEOS_ERROR( "error reading strides");
  }

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex temp = strides[var.getSingleParameterResizeIndex()];
    sizeOfUnpackedChars += Unpack( buffer, var.data( indices[a] ), temp );
  }

  return sizeOfUnpackedChars;

}


template< bool DO_PACKING >
int Pack( char*& buffer,
          localIndex const * const var,
          localIndex const length,
          globalIndex_array const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  sizeOfPackedChars += length*sizeof(globalIndex);
  static_if( DO_PACKING )
  {
    for( localIndex a=0 ; a<length ; ++a )
    {
      memcpy( buffer, &(localToGlobalMap[var[a]]), sizeof(globalIndex) );
      buffer += sizeof(globalIndex);
    }
  });

  return sizeOfPackedChars;
}

inline
int Unpack( char const *& buffer,
            localIndex_array & var,
            map<globalIndex,localIndex> const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex length;
  sizeOfUnpackedChars += Unpack( buffer, length );
  var.resize(length);

  for( localIndex a=0 ; a<length ; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex );
    var[a] = globalToLocalMap.at(unpackedGlobalIndex);
  }

  return sizeOfUnpackedChars;
}



inline
int Unpack( char const *& buffer,
            localIndex * const var,
            localIndex const length,
            map<globalIndex,localIndex> const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex lengthUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, lengthUnpacked );

  GEOS_ASSERT( length==lengthUnpacked, "CommBufferOps::Unpack(): length incorrect." )

  for( localIndex a=0 ; a<length ; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex );
    var[a] = globalToLocalMap.at(unpackedGlobalIndex);
  }

  return sizeOfUnpackedChars;
}




//
//template< bool DO_PACKING >
//int Pack( char*& buffer,
//          multidimensionalArray::ManagedArray<localIndex,2,localIndex> const & var,
//          localIndex_array const & indices,
//          globalIndex_array const & localToGlobalMap )
//{
//  int sizeOfPackedChars;
//
//  return sizeOfPackedChars;
//}
//
//inline
//int Unpack( char const *& buffer,
//            multidimensionalArray::ManagedArray<localIndex,2,localIndex> & var,
//            localIndex_array const & indices,
//            globalIndex_array const & globalToLocalMap )
//{
//  int sizeOfUnpackedChars;
//
//  return sizeOfUnpackedChars;
//}


template< bool DO_PACKING >
int Pack( char*& buffer,
          array< localIndex_array > const & var,
          localIndex_array const & indices,
          globalIndex_array const & localToGlobalMap,
          globalIndex_array const & relatedObjectLocalToGlobalMap )
{
  int sizeOfPackedChars=0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobalMap[li] );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[li].data(), var[li].size(), relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

inline
int Unpack( char const *& buffer,
            array< localIndex_array > & var,
            localIndex_array const & indices,
            map<globalIndex,localIndex> const & globalToLocalMap,
            map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap )
{
  int sizeOfUnpackedChars=0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];

    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );
    // do a check here on the global Index??

    sizeOfUnpackedChars += Unpack( buffer, var[li], relatedObjectGlobalToLocalMap );
  }


  return sizeOfUnpackedChars;
}

template< bool DO_PACKING >
int Pack( char*& buffer,
          Array2dT<localIndex> const & var,
          localIndex_array const & indices,
          globalIndex_array const & localToGlobalMap,
          globalIndex_array const & relatedObjectLocalToGlobalMap )
{
  int sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, indices.size() );
  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, localToGlobalMap[li] );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[li].data(), var[li].size(), relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


inline
int Unpack( char const *& buffer,
            Array2dT<localIndex> & var,
            localIndex_array const & indices,
            map<globalIndex,localIndex> const & globalToLocalMap,
            map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap )
{
  int sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked );
  GEOS_ASSERT( numIndicesUnpacked==indices.size(), "CommBufferOps::Unpack(): Incorrect number of indices unpacked." )

  for( localIndex a=0 ; a<indices.size() ; ++a )
  {
    localIndex li = indices[a];

    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi );
    // do a check here on the global Index??

    sizeOfUnpackedChars += Unpack( buffer, var[li].data(), var[li].size(), relatedObjectGlobalToLocalMap );
  }


  return sizeOfUnpackedChars;
}


}
}



#endif /* SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMBUFFEROPS_INLINE_HPP_ */
