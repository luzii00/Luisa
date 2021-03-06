/*
 * Logger.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: settgast1
 */

#ifndef SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_

#include <string>
#include "common/GeosxConfig.hpp"
#include <sstream>
#ifdef USE_ATK
#include "slic/slic.hpp"
#include "slic/GenericOutputStream.hpp"
#endif


namespace geosx
{
void geos_abort( std::string message );

#if 0
#define GEOS_ERROR(msg) SLIC_ERROR(msg)
#else
#define GEOS_ERROR(msg) \
    std::cerr<<"***** GEOS_ERROR "<<std::endl;\
    std::cerr<<"***** FILE: "<<__FILE__<<std::endl;\
    std::cerr<<"***** LINE: "<<__LINE__<<std::endl;\
    std::ostringstream oss;\
    oss << msg;\
    geosx::geos_abort(oss.str());
#endif

#define GEOS_ASSERT( CONDITION, msg) \
  if( !(CONDITION) )\
  {\
    std::cerr<<"***** GEOS_ASSERT "<<std::endl;\
    std::cerr<<"***** FILE: "<<__FILE__<<std::endl;\
    std::cerr<<"***** LINE: "<<__LINE__<<std::endl;\
    std::ostringstream oss;\
    oss << msg;\
    geosx::geos_abort(oss.str());\
  }

}

#endif /* SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_ */
