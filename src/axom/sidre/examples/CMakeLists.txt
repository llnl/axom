# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sidre and Spio examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(example_sources
    sidre_data_vs_metadata.cpp
    sidre_stressgroups.cpp
    sidre_shocktube.cpp
)

if(AXOM_ENABLE_MPI)
    set (example_parallel_sources
         sidre_array.cpp
         sidre_createdatastore.cpp
         sidre_external_array.cpp
         sidre_generateindex.cpp)
else()
    list(APPEND example_sources sidre_createdatastore.cpp)
endif()

set(F_example_sources
    sidre_shocktube_F.f
)

if(APPLE)
    set(EXTRA_LIBS "")
elseif(UNIX)
    set(EXTRA_LIBS rt)
endif()

#------------------------------------------------------------------------------
# Add targets and tests for serial examples
#------------------------------------------------------------------------------
foreach(example_source ${example_sources})
    get_filename_component(exe_name ${example_source} NAME_WE)
    axom_add_executable(
        NAME       ${exe_name}_ex
        SOURCES    ${example_source}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  sidre ${EXTRA_LIBS}
        FOLDER      axom/sidre/examples )

    if(AXOM_ENABLE_TESTS)
        axom_add_test(
            NAME ${exe_name}
            COMMAND ${exe_name}_ex )
    endif()
endforeach()

#------------------------------------------------------------------------------
# Add targets and tests for parallel examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_MPI)
    foreach(example_source ${example_parallel_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        axom_add_executable(
            NAME       ${exe_name}_ex
            SOURCES    ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON sidre ${EXTRA_LIBS}
            FOLDER axom/sidre/examples)

        if(AXOM_ENABLE_TESTS)
            axom_add_test( NAME ${exe_name}
                           COMMAND ${exe_name}_ex
                           NUM_MPI_TASKS 1
                           )
        endif()
    endforeach()

    blt_add_test( NAME sidre_createdatastore_parallel
                  COMMAND sidre_createdatastore_ex
                  NUM_MPI_TASKS 4
                  )

    blt_add_test( NAME sidre_generateindex_parallel
                  COMMAND sidre_generateindex_ex
                  NUM_MPI_TASKS 4
                  )

endif()

#------------------------------------------------------------------------------
# Add MFEM Sidre DataCollection tests
#------------------------------------------------------------------------------
if(MFEM_FOUND AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)
    set(_mfem_dc_tests
        sidre_mfem_datacollection_materials.cpp
        sidre_mfem_datacollection_restart.cpp
        sidre_mfem_datacollection_vis.cpp)

    foreach(example_source ${_mfem_dc_tests})
        get_filename_component(exe_name ${example_source} NAME_WE)
        axom_add_executable(
            NAME        ${exe_name}_ex
            SOURCES     ${example_source}
            OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON  sidre ${EXTRA_LIBS}
            FOLDER      axom/sidre/examples)

        if(AXOM_ENABLE_TESTS)
            if(AXOM_ENABLE_MPI)
                if(MFEM_USE_MPI)
                    set(_nranks 2)
                else()
                    set(_nranks 1)
                endif()

                axom_add_test( NAME ${exe_name}
                               COMMAND ${exe_name}_ex
                               NUM_MPI_TASKS ${_nranks})
            else()
                axom_add_test( NAME ${exe_name}
                               COMMAND ${exe_name}_ex)
            endif()
        endif()

    endforeach()
endif()

#------------------------------------------------------------------------------
# Add targets and tests for Fortran examples
#------------------------------------------------------------------------------
if(ENABLE_FORTRAN)
    foreach(example_source ${F_example_sources})
        get_filename_component(exe_name ${example_source} NAME_WE)
        axom_add_executable(
            NAME ${exe_name}_ex
            SOURCES ${example_source}
            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
            DEPENDS_ON sidre ${EXTRA_LIBS}
            FOLDER axom/sidre/examples )

        if(AXOM_ENABLE_TESTS)
            axom_add_test(
                NAME ${exe_name}
                COMMAND ${exe_name}_ex )
        endif()
    endforeach()
endif()

#------------------------------------------------------------------------------
# Add Spio examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_MPI)
    add_subdirectory(spio)
endif()

#------------------------------------------------------------------------------
# Sidre version of lulesh 2.0
# Note: Disabled on windows builds due to Lulesh's usage of unix APIs.
#------------------------------------------------------------------------------
if(NOT WIN32)
  add_subdirectory(lulesh2)
endif()
