
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    schema/SchemaUtilities.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    schema/SchemaUtilities.cpp
    PARENT_SCOPE )
