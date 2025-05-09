# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Spin component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       SPIN
                        COMPONENTS PRIMAL SLIC SLAM)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set( spin_headers

     BVH.hpp
     Brood.hpp
     DenseOctreeLevel.hpp
     ImplicitGrid.hpp
     MortonIndex.hpp
     OctreeBase.hpp
     OctreeLevel.hpp
     RectangularLattice.hpp
     SparseOctreeLevel.hpp
     SpatialOctree.hpp
     UniformGrid.hpp

     ## internal
     internal/linear_bvh/RadixTree.hpp
     internal/linear_bvh/build_radix_tree.hpp
     internal/linear_bvh/bvh_traverse.hpp
     internal/linear_bvh/bvh_vtkio.hpp

     ## policy
     policy/LinearBVH.hpp
     policy/UniformGridStorage.hpp
   )

#------------------------------------------------------------------------------
# Specify spin dependencies
#------------------------------------------------------------------------------
set( spin_depends_on
     primal 
     slam)

blt_list_append( TO spin_depends_on ELEMENTS sparsehash IF SPARSEHASH_FOUND )
blt_list_append( TO spin_depends_on ELEMENTS RAJA IF RAJA_FOUND )
blt_list_append( TO spin_depends_on ELEMENTS umpire IF UMPIRE_FOUND )

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
axom_add_library(
    NAME                  spin
    HEADERS               ${spin_headers}
    FOLDER                axom/spin
    DEPENDS_ON            ${spin_depends_on})

axom_write_unified_header(NAME    spin
                          HEADERS ${spin_headers})

axom_install_component(NAME    spin
                       HEADERS ${spin_headers}
                       )

#------------------------------------------------------------------------------
# add examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

#------------------------------------------------------------------------------
# add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()
