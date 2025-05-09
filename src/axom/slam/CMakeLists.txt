# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# SLAM -- Set-theoretic Lightweight API for Meshes (low level mesh API)
#------------------------------------------------------------------------------
# Some additional flags that might be useful
#   To avoid warnings: -Wno-missing-braces
#   Other: -fno-strict-aliasing -Wextra
#   Vectorizing: -ftree-vectorize -fdump-tree-vect -ftree-vectorizer-verbose=7
#                -fdump-tree-optimized-details -march=native 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       SLAM
                        COMPONENTS SLIC)


#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(slam_headers
    # Utility files
    Utilities.hpp
    FieldRegistry.hpp
    ModularInt.hpp

    #SRM policies
    policies/CardinalityPolicies.hpp
    policies/SizePolicies.hpp
    policies/OffsetPolicies.hpp
    policies/StridePolicies.hpp
    policies/IndirectionPolicies.hpp
    policies/SubsettingPolicies.hpp
    policies/PolicyTraits.hpp
    policies/InterfacePolicies.hpp
    policies/SetInterfacePolicies.hpp
    policies/MapInterfacePolicies.hpp
    policies/BivariateSetInterfacePolicies.hpp

    # SRM Set headers
    IndirectionSet.hpp
    NullSet.hpp
    OrderedSet.hpp
    RangeSet.hpp
    ProductSet.hpp
    DynamicSet.hpp
    Set.hpp
    BivariateSet.hpp
    BitSet.hpp
    RelationSet.hpp

    # SRM Relation headers
    Relation.hpp
    StaticRelation.hpp
    DynamicVariableRelation.hpp
    DynamicConstantRelation.hpp

    # SRM Map headers
    Map.hpp
    MapBase.hpp
    BivariateMap.hpp
    SubMap.hpp
    DynamicMap.hpp

    # Topological mesh headers
    mesh_struct/IA.hpp
    mesh_struct/IA_impl.hpp
    )

set(slam_sources
    # SRM Set sources
    BitSet.cpp
    OrderedSet.cpp

    # SRM Relation sources

    # SRM Map sources

    # Other sources
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
axom_add_library(
    NAME        slam
    SOURCES     ${slam_sources}
    HEADERS     ${slam_headers}
    DEPENDS_ON  slic
    FOLDER      axom/slam )

axom_write_unified_header(NAME    slam
                          HEADERS ${slam_headers} )

axom_install_component(NAME    slam
                       HEADERS ${slam_headers} )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
  if (ENABLE_BENCHMARKS)
    add_subdirectory(benchmarks)
  endif()
endif()

if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
