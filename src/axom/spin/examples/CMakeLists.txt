# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Spin component
#------------------------------------------------------------------------------

set(spin_example_depends
        spin
        fmt
        )

axom_add_executable(
    NAME        spin_introduction_ex
    SOURCES     spin_introduction.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${spin_example_depends}
    FOLDER      axom/spin/examples
    )

if(CONDUIT_FOUND AND CLI11_FOUND)
    axom_add_executable(
        NAME        spin_space_filling_curve_ex
        SOURCES     spin_space_filling_curve.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${spin_example_depends} cli11 conduit::conduit
        FOLDER      axom/spin/examples)

    if(AXOM_ENABLE_TESTS)
        foreach(_dimension 2 3)
            foreach(_order morton hilbert)
                set(_test "spin_space_filling_curve_${_dimension}d_${_order}")
                axom_add_test(
                    NAME    ${_test}
                    COMMAND spin_space_filling_curve_ex
                            --dimension ${_dimension}
                            --order ${_order}
                            --level 2
                            --num-points 128
                            --seed 42
                            --protocol json
                            --output ${_test})
            endforeach()
        endforeach()
    endif()
endif()
