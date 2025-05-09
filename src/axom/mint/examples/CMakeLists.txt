# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Mint component
#------------------------------------------------------------------------------

set( mint_examples
     mint_curvilinear_mesh.cpp
     mint_heat_equation.cpp
     mint_nbody_solver.cpp
     mint_particle_mesh.cpp
     mint_rectilinear_mesh.cpp
     mint_su2_mesh.cpp
     mint_unstructured_mixed_topology_mesh.cpp
     mint_unstructured_single_topology_mesh.cpp

     user_guide/mint_tutorial.cpp
     user_guide/mint_getting_started.cpp
   )

foreach( example ${mint_examples} )

   get_filename_component( example_name ${example} NAME_WE )

   axom_add_executable(
        NAME       ${example_name}_ex
        SOURCES    ${example}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON mint
        FOLDER     axom/mint/examples
        )

endforeach()
