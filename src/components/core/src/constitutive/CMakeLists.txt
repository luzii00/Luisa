#
# Specify all headers
#
set(geosx_core_sources
  ${geosx_core_sources}
    constitutive/ConstitutiveManager.hpp
    constitutive/ConstitutiveBase.hpp
    constitutive/LinearElasticIsotropic.hpp
    PARENT_SCOPE )
#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    constitutive/ConstitutiveManager.cpp
    constitutive/ConstitutiveBase.cpp
    constitutive/LinearElasticIsotropic.cpp
    PARENT_SCOPE )
