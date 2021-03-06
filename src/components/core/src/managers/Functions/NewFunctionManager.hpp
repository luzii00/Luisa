/*
 * NewFunctionManager.hpp
 *
 *  Created on: July 6, 2017
 *      Author: sherman
 */

#ifndef NEWFUNCTIONMANAGER_HPP_
#define NEWFUNCTIONMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "FunctionBase.hpp"

namespace geosx
{



class NewFunctionManager : public dataRepository::ManagedGroup
{
public:
  NewFunctionManager( const std::string& name,
                      dataRepository::ManagedGroup * const parent );
  virtual ~NewFunctionManager();

  static NewFunctionManager * Instance()
  {
    static NewFunctionManager theFunctionManager("LastFunctionManagerOnEarth", nullptr);

    return &theFunctionManager;
  }

  static string CatalogName() { return "NewFunctionManager"; }
  virtual void FillDocumentationNode() override;
  virtual void CreateChild( string const & functionCatalogKey, string const & functionName ) override;
};


} /* namespace geosx */

#endif /* NEWFUNCTIONMANAGER_HPP_ */
