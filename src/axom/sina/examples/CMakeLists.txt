# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina examples
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# List of single source file examples
#------------------------------------------------------------------------------
set(sina_example_sources
    sina_basic.cpp
    sina_check_datum_type.cpp
    sina_create_datum.cpp
    sina_curve_set.cpp
    sina_document_assembly.cpp
    sina_file_object_creation.cpp
    sina_file_object_removal.cpp
    sina_local_id_relationship.cpp
    sina_query_record_for_files.cpp
    sina_query_records_relationships.cpp
    sina_relationship_assembly.cpp
    sina_set_datum_units_tags.cpp
    sina_tutorial.cpp
    sina_view_datum_types.cpp
    sina_view_datum_values.cpp
)

set(sina_example_depends sina conduit::conduit slic)

if (ENABLE_FORTRAN)
    blt_list_append( TO sina_example_sources ELEMENTS sina_fortran.f)
endif()

#------------------------------------------------------------------------------
# Add targets for Sina examples
#------------------------------------------------------------------------------
foreach(src ${sina_example_sources})
    get_filename_component(exe_name ${src} NAME_WE)

    axom_add_executable(
        NAME       ${exe_name}_ex
        SOURCES    ${src}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${sina_example_depends}
        FOLDER     axom/sina/examples
    )

    if (${exe_name} STREQUAL "sina_fortran" AND CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
        target_compile_options(${exe_name}_ex PRIVATE -qextname)
    endif()
endforeach()
