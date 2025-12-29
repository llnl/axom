# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Unit tests for Spin component
#------------------------------------------------------------------------------

set(spin_tests
    spin_implicit_grid.cpp
    spin_morton.cpp
    spin_octree.cpp
    spin_rectangular_lattice.cpp
    spin_spatial_octree.cpp
    spin_uniform_grid.cpp
   )

set(spin_tests_depends
    spin
    gtest
    fmt
    )

foreach ( test ${spin_tests} )

    get_filename_component( test_name ${test} NAME_WE )

    axom_add_executable(
      NAME       ${test_name}_test
      SOURCES    ${test}
      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
      DEPENDS_ON ${spin_tests_depends}
      FOLDER     axom/spin/tests
      )

    axom_add_test(
      NAME    ${test_name}
      COMMAND ${test_name}_test
      )

endforeach()

if (AXOM_ENABLE_MINT)

    axom_add_executable(
      NAME       spin_bvh_test
      SOURCES    spin_bvh.cpp
      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
      DEPENDS_ON ${spin_tests_depends} mint
      FOLDER     axom/spin/tests
      )

    # sets the number of threads for the omp policy; leave empty for non-omp policies
    set(_num_threads)
    if(AXOM_ENABLE_OPENMP)
        set(_num_threads ${AXOM_TEST_NUM_OMP_THREADS})
    endif()

    axom_add_test(
      NAME            spin_bvh
      COMMAND         spin_bvh_test
      NUM_OMP_THREADS ${_num_threads}
      )

endif()
