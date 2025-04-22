# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Spio examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(spio_example_sources
    IORead.cpp
    IOWrite.cpp
    )

set(spio_example_depends sidre ${EXTRA_LIBS})

#------------------------------------------------------------------------------
# Add targets and tests for spio examples
#------------------------------------------------------------------------------
foreach(src ${spio_example_sources})
    get_filename_component(exe_name ${src} NAME_WE)
    axom_add_executable(
        NAME       spio_${exe_name}_ex
        SOURCES    ${src}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${spio_example_depends}
        FOLDER     axom/sidre/examples
        )
endforeach()

if(SCR_FOUND AND AXOM_ENABLE_TESTS AND HDF5_FOUND)
    # Note: This example is a combination of a test and an example
    # since we're unable to run two SCR executables at the same time
    #
    # The test was written against the sidre_hdf5 protocol, so disable it
    # in configurations that do not include hdf5
    
    axom_add_executable(
        NAME       spio_IO_SCR_Checkpoint_ex
        SOURCES    IO_SCR_Checkpoint.cpp
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${spio_example_depends} gtest
        FOLDER     axom/sidre/examples
        )

    axom_add_test(
        NAME          spio_IO_SCR_Checkpoint
        COMMAND       spio_IO_SCR_Checkpoint_ex
        NUM_MPI_TASKS 4 )
endif()
