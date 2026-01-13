# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Mint component
#------------------------------------------------------------------------------

set( mint_tests
    ## execution model tests
    mint_execution_cell_traversals.cpp
    mint_execution_node_traversals.cpp
    mint_execution_face_traversals.cpp
    mint_execution_xargs.cpp

    ## mesh tests
    mint_mesh.cpp
    mint_mesh_cell_types.cpp
    mint_mesh_blueprint.cpp
    mint_mesh_connectivity_array.cpp
    mint_mesh_coordinates.cpp
    mint_mesh_curvilinear_mesh.cpp
    mint_mesh_field.cpp
    mint_mesh_field_data.cpp
    mint_mesh_field_types.cpp
    mint_mesh_field_variable.cpp
    mint_mesh_particle_mesh.cpp
    mint_mesh_rectilinear_mesh.cpp
    mint_mesh_uniform_mesh.cpp
    mint_mesh_unstructured_mesh.cpp
    mint_mesh_face_relation.cpp

    ## fem tests
    mint_fem_shape_functions.cpp
    mint_fem_single_fe.cpp

    ## util tests
    mint_su2_io.cpp
    mint_util_write_vtk.cpp

    ## container tests
    mint_deprecated_mcarray.cpp
   )

foreach( test ${mint_tests} )

    get_filename_component( test_name ${test} NAME_WE )

    axom_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON mint gtest
        FOLDER     axom/mint/tests
        )

    # sets the number of threads for the omp policy; leave empty for non-omp policies
    set(_num_threads)
    if(AXOM_ENABLE_OPENMP)
        set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
    endif()

    axom_add_test( NAME    ${test_name}
                   COMMAND ${test_name}_test
                   NUM_OMP_THREADS ${_num_threads}
                 )

endforeach()
