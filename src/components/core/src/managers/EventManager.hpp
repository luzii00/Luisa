/*
 * EventManager.hpp
 *
 *  Created on: Oct 5, 2016
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace pugi
{
class xml_node;
}

namespace geosx
{

class EventManager : public dataRepository::ManagedGroup
{
public:
  EventManager( std::string const & name,
                ManagedGroup * const parent );

  virtual ~EventManager();

  void CreateSolverApplication(pugi::xml_node const & applicationNode);

  void ReadXML( pugi::xml_node const & problemNode );

  void CheckEventTiming();

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_EVENTMANAGER_HPP_ */
