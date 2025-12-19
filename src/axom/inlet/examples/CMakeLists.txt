# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of examples
#------------------------------------------------------------------------------

set(inlet_examples
    nested_structs.cpp )

blt_list_append(
  TO
    inlet_examples
  ELEMENTS
    arrays.cpp
    documentation_generation.cpp
    fields.cpp
    lua_library.cpp
    containers.cpp
    user_defined_type.cpp
    verification.cpp
  IF
    SOL_FOUND )


#------------------------------------------------------------------------------
# Build each executable example
#------------------------------------------------------------------------------

foreach(example ${inlet_examples})
    get_filename_component( ex_name ${example} NAME_WE )
    axom_add_executable(NAME       inlet_${ex_name}_ex
                        SOURCES    ${example}
                        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet
                        FOLDER     axom/inlet/examples )
endforeach()

if (SOL_FOUND AND MFEM_FOUND)
    axom_add_executable(NAME       inlet_mfem_coefficient_ex
                        SOURCES    mfem_coefficient.cpp
                        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet
                        FOLDER     axom/inlet/examples )
endif()


#------------------------------------------------------------------------------
# Add examples as tests
#------------------------------------------------------------------------------

axom_add_test( NAME    inlet_nested_structs_ex
               COMMAND inlet_nested_structs_ex )

if (SOL_FOUND)
    axom_add_test( NAME    inlet_arrays_ex
                   COMMAND inlet_arrays_ex )

    axom_add_test( NAME    inlet_documentation_generation_ex
                   COMMAND inlet_documentation_generation_ex --enableDocs --file ${CMAKE_CURRENT_LIST_DIR}/example1.lua )

    axom_add_test( NAME    inlet_fields_ex
                   COMMAND inlet_fields_ex)

    axom_add_test( NAME    inlet_lua_library_ex
                   COMMAND inlet_lua_library_ex )

    axom_add_test( NAME    inlet_containers_ex
                   COMMAND inlet_containers_ex)

    axom_add_test( NAME    inlet_user_defined_type_ex
                   COMMAND inlet_user_defined_type_ex --file ${CMAKE_CURRENT_LIST_DIR}/example1.lua)

    axom_add_test( NAME    inlet_verification_ex
                   COMMAND inlet_verification_ex )

    if (MFEM_FOUND)
        axom_add_test( NAME    inlet_mfem_coefficient_ex
                       COMMAND inlet_mfem_coefficient_ex --file ${CMAKE_CURRENT_LIST_DIR}/mfem_coef.lua)
    endif()

endif()
