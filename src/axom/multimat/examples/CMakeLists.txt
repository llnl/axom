# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# MultiMat examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Single source standalone examples
#------------------------------------------------------------------------------
set(multimat_example_headers
    helper.hpp )

set(multimat_example_sources
    basic.cpp
    calculate.cpp
    calculate_gpu.cpp
    traversal.cpp )

if(CONDUIT_FOUND)
    set(multimat_basic_depends conduit::conduit)
endif()

foreach(example_source ${multimat_example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    set(deps "multimat;${multimat_${example_name}_depends}")
    axom_add_executable(
        NAME        multimat_${example_name}_ex
        HEADERS     ${multimat_example_headers}
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${deps}
        FOLDER      axom/multimat/examples
        )
    unset(deps)

    axom_add_test(
        NAME        multimat_${example_name}
        COMMAND     multimat_${example_name}_ex
        )
endforeach()
