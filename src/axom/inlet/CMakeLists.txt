# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       Inlet
                        COMPONENTS SLIC Sidre PRIMAL )

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(inlet_headers
    Field.hpp
    Function.hpp
    Inlet.hpp
    InletVector.hpp
    Reader.hpp
    Container.hpp 
    Proxy.hpp
    VariantKey.hpp
    Verifiable.hpp
    VerifiableScalar.hpp
    Writer.hpp
    SphinxWriter.hpp
    JSONSchemaWriter.hpp
    ConduitReader.hpp
    YAMLReader.hpp
    JSONReader.hpp
    inlet_utils.hpp )

set(inlet_sources
    Field.cpp
    Function.cpp
    Inlet.cpp
    Container.cpp
    Proxy.cpp
    VariantKey.cpp
    VerifiableScalar.cpp
    SphinxWriter.cpp
    JSONSchemaWriter.cpp
    ConduitReader.cpp
    inlet_utils.cpp )

blt_list_append(TO inlet_headers ELEMENTS LuaReader.hpp IF SOL_FOUND)
blt_list_append(TO inlet_sources ELEMENTS LuaReader.cpp IF SOL_FOUND)

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
set(inlet_depends sidre fmt)

blt_list_append(TO inlet_depends ELEMENTS sol lua IF SOL_FOUND)

axom_add_library(NAME        inlet
                 SOURCES     ${inlet_sources}
                 HEADERS     ${inlet_headers}
                 DEPENDS_ON  ${inlet_depends}
                 FOLDER      axom/inlet )

axom_write_unified_header(NAME    inlet
                          HEADERS ${inlet_headers} )

axom_install_component(NAME    inlet
                       HEADERS ${inlet_headers} )


#------------------------------------------------------------------------------
# Add examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()
