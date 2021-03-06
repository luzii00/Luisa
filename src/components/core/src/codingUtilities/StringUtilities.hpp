/*
 * StringUtilities.hpp
 *
 *  Created on: Nov 12, 2014
 *      Author: rrsettgast
 */

#ifndef STRINGUTILITIES_HPP_
#define STRINGUTILITIES_HPP_

#include <cxxabi.h>
#include <cstring>
#include <memory>
//#include <vector>
#include <sstream>
#include <algorithm>
#include <map>

//#include "Common/Common.h"
//#include <string>

#include "../common/DataTypes.hpp"

namespace geosx
{
namespace stringutilities
{
template <class Type>
std::string toString(Type theVar);

template <class Type>
Type fromString(std::string theString);

//////////////////////////////////////////////////

//
// toString
//
/// convert a variable to a string
template <class Type>
std::string toString(Type theVar){
  std::ostringstream oss;
  oss << theVar;
  return oss.str();
}
// override the template for string->string
template <>
inline std::string toString<std::string>(std::string theVar)
{
  return theVar;
}

//
// fromString
//
/// returns a variable from a string
template <class Type>
Type fromString(std::string theString){
  std::istringstream iss(theString);
  Type theVar;
  iss >> theVar;
  return theVar;
}
// override the template for string->string
template <>
inline std::string fromString<std::string>(std::string theVar)
{
  return theVar;
}

// override the template for string->real64
// Allows unit manager to convert units
template <>
real64 fromString<real64>(std::string theVar);

// override the template for FieldType
//template <>
//inline FieldType fromString<FieldType>(std::string theString){
//  using namespace FieldInfo;
//  if(theString==FieldInfo::IntegerStr){
//    return integerField;
//  } else if(theString==FieldInfo::GlobalIndexStr){
//  return globalIndexField;
//  } else if(theString==FieldInfo::LocalIndexStr){
//  return localIndexField;
//  } else if(theString==FieldInfo::RealStr){
//    return realField;
//  } else if(theString==FieldInfo::R1TensorStr){
//    return R1TensorField;
//  } else if(theString==FieldInfo::R2TensorStr){
//    return R2TensorField;
//  } else if(theString==FieldInfo::R2SymTensorStr){
//    return R2SymTensorField;
//  }else {
//    throw GPException("Error fromString: unrecognized field type: " +
// theString +".");
//  }
//}

/// Convert a string to lowercase
void toLower(std::string& theString);
std::string lowercase(std::string theString);

/// Convert a string to uppercase
void toUpper(std::string& theString);
std::string uppercase(std::string theString);

/// Check for case insensitive equality between strings
bool ieq(std::string strA,std::string strB);

/// Overloaded function to check equality between strings and char arrays
/// Mainly used to avoid char*==char* mistakes
inline bool streq(const std::string& strA, const std::string& strB){
  return strA == strB;
}
inline bool streq(const std::string& strA, const char * strB){
  return strA == strB;
}
inline bool streq(const char * strA, const std::string& strB){
  return strA == strB;
}
inline bool streq(const char * strA, const char * strB){
  return !strcmp(strA, strB);
}

/// string is integer
inline bool strIsInt(std::string theString){
  std::istringstream iss(theString);
  int dummy;
  iss >> dummy;
  return !iss.fail();
}

/// Subdivide string by delimiters
string_array Tokenize(const std::string& str, const std::string& delimiters);

/// Subdivide string delimited by sequence of characters
string_array TokenizeSeq(const std::string& str, const std::string& seq);

/// Split string at first token
string_array Split(const std::string& str, const std::string& delimiters);

/// Remove comments from end of string
void RemoveComments(std::string& str, char d='%');

/// Remove all spaces ' ' from a string
inline void RemoveSpaces(std::string& aString){
  aString.erase(std::remove(aString.begin(), aString.end(), ' '), aString.end());
}

/// Expand string vector based on multiple tokens eg [a, b**3, c] => [a,b,b,b,c]
inline void ExpandMultipleTokens(string_array& sVector, const std::string& multipleToken="**"){
  localIndex n= sVector.size();
  string_array newVec;
  for( int i =0 ; i < n ; ++i)
  {
    string_array keyMult = TokenizeSeq(sVector[i], multipleToken);
    if( (keyMult.size() == 2) && strIsInt(keyMult[1]) )
    {
      int numMult = fromString<int>(keyMult[1]);
      for(int j=0 ; j < numMult ; ++j )
      {
        newVec.push_back(keyMult[0]);
      }
    }
    else
    {
      newVec.push_back(sVector[i]);
    }

  }
  sVector = newVec;
}

/// Trim whitespace from string
void TrimLeft(std::string& str, const std::string& d=" \t\n\r");
void TrimRight(std::string& str, const std::string& d=" \t\n\r");
void Trim(std::string& str, const std::string& d=" \t\n\r");

inline void Trim(string_array& strVect, const std::string& d=" \t\n\r"){
  for(int i =0 ; i < int(strVect.size()) ; ++i)
    Trim(strVect[i],d);
}

/// Replace parameters of form "$:NAME" in a string, returns true if a parameter
// is detected
/// @param lineStr - the string containing the parameters
/// @param parameterMap - map of parameter names and replacement strings,
/// @param prefix - Character sequence used to signal start of parameter ("$:"
// by default),
bool ReplaceParameters(std::string& lineStr, const std::map<std::string,std::string>& parameterMap,const std::string& prefix= "$:");



template< typename T >
inline void StringToType( std::vector<T>& destination, std::string const & source )
{
  std::istringstream ss( source );

  T value;

  while(ss.peek() == ',' || ss.peek() == ' ')
  {
    ss.ignore();
  }
  while( ss>>value )
  {
    destination.push_back( value );
    while(ss.peek() == ',' || ss.peek() == ' ')
    {
      ss.ignore();
    }
  }

}

template <>
inline void StringToType<std::string>( std::vector<std::string>& destination, std::string const & source )
{
  destination.push_back( source );
}


}
}

#endif /* STRINGUTILITIES_HPP_ */
