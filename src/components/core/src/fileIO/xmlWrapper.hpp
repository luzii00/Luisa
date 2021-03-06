/*
 * xmlWrapper.hpp
 *
 *  Created on: Jun 3, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_

#include "common/DataTypes.hpp"
#include "pugixml.hpp"
//#include <string>
#include <sstream>
#include "math/TensorT/TensorT.h"

namespace cxx_utilities
{
class DocumentationNode;
}
namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}


class xmlWrapper
{
public:
  using xmlDocument = pugi::xml_document;
  using xmlResult = pugi::xml_parse_result;
  using xmlNode = pugi::xml_node;

  xmlWrapper();
  virtual ~xmlWrapper();

  template< typename T >
  static void as_type( std::vector<T> & target, std::string value, std::string defValue );

//  template< typename T >
//  static T as_type( xmlNode const & node, std::string const name, T defValue
// );

  static R1Tensor as_type( xmlNode const & node, std::string const name, R1Tensor defValue );

  static void ReadAttributeAsType( dataRepository::ManagedGroup & group,
                                   cxx_utilities::DocumentationNode const & subDocNode,
                                   xmlNode const & targetNode );
};


template< typename T >
void xmlWrapper::as_type( std::vector<T> & target, std::string inputValue, std::string defValue )
{
  std::string csvstr = ( inputValue!="") ? inputValue : defValue;
  std::istringstream ss( csvstr );

  T value;

  while(ss.peek() == ',' || ss.peek() == ' ')
  {
    ss.ignore();
  }
  while( !((ss>>value).fail()) )
  {
    target.push_back( value );
    while(ss.peek() == ',' || ss.peek() == ' ')
    {
      ss.ignore();
    }
  }
}


//
//template< typename T >
//T xmlWrapper::as_type( xmlNode const & node, std::string const name, T
// defValue )
//{
//  T rval = defValue;
//  pugi::xml_attribute att = node.attribute( name.c_str() );
//
//  if( !att.empty() )
//  {
//
//  }
//
//  return rval;
//}



} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_ */
