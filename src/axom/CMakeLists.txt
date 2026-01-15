# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Add Axom components
#------------------------------------------------------------------------------
# Note: Default state of components is controlled by the
# AXOM_ENABLE_ALL_COMPONENTS option.
#------------------------------------------------------------------------------

# Core is an essential part of Axom and cannot be turned off
axom_add_component(COMPONENT_NAME core DEFAULT_STATE ON)
if(DEFINED AXOM_ENABLE_CORE AND NOT AXOM_ENABLE_CORE)
    message(FATAL_ERROR "Core is an essential part of Axom and cannot be turned off")
endif()

# Lumberjack is a parallel message filtering and reduction library. It can be used
# by itself or as a SLIC Log Stream. It is not meant for serial (non-mpi) executables.
if(AXOM_ENABLE_MPI)
    axom_add_component(COMPONENT_NAME lumberjack
                       DEFAULT_STATE  ${AXOM_ENABLE_ALL_COMPONENTS})
else()
    message(STATUS "Axom Component Lumberjack turned off due to AXOM_ENABLE_MPI=OFF")
    set(AXOM_ENABLE_LUMBERJACK OFF CACHE BOOL "")
endif()

# Add components so that later ones depend on earlier ones
# (i.e. quest depends on mint, so quest follows mint)
axom_add_component(COMPONENT_NAME slic   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME sina   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME slam   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME primal DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME sidre  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME mint   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME spin   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME bump   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME mir    DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME inlet  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME klee   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME quest  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME multimat DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})

install(TARGETS              ${AXOM_COMPONENTS_ENABLED}
        EXPORT               axom-targets
        DESTINATION          lib)
install(EXPORT               axom-targets 
        NAMESPACE            axom::
        DESTINATION          lib/cmake)

# Create alias targets for projects that include us as a subdirectory
add_library(axom INTERFACE)
target_link_libraries(axom INTERFACE ${AXOM_COMPONENTS_ENABLED})

foreach(_comp ${AXOM_COMPONENTS_ENABLED})
    add_library(axom::${_comp} ALIAS ${_comp})
endforeach()


#------------------------------------------------------------------------------
# Generate export symbols
#------------------------------------------------------------------------------
include(GenerateExportHeader)

foreach(_component ${AXOM_COMPONENTS_ENABLED})
    # Skip header only libraries
    get_property(_targetType TARGET ${_component} PROPERTY TYPE)
    if(${_targetType} STREQUAL "INTERFACE_LIBRARY")
        continue()
    endif()

    set(_export_header ${PROJECT_BINARY_DIR}/include/axom/export/${_component}.h)
    generate_export_header( ${_component} 
                            PREFIX_NAME AXOM_
                            EXPORT_FILE_NAME ${_export_header})
    install(FILES ${_export_header} 
            DESTINATION include/axom/export )
endforeach()

#------------------------------------------------------------------------------
# Output some information about the configuration
#------------------------------------------------------------------------------
foreach(_component ${AXOM_COMPONENTS_FULL})
    if(${_component} IN_LIST AXOM_COMPONENTS_ENABLED)
        message(STATUS "Axom component ${_component} is ON")
    else()
        message(STATUS "Axom component ${_component} is OFF")
    endif()
endforeach()

if(AXOM_DATA_DIR)
    message(STATUS "AXOM_DATA_DIR: ${AXOM_DATA_DIR}")
else()
    message(STATUS "AXOM_DATA_DIR: <undefined>")
endif()

#------------------------------------------------------------------------------
# Fix FOLDER property for some targets
#------------------------------------------------------------------------------
# If Axom is the main project, set FOLDER property for top-level check targets
if ("${PROJECT_SOURCE_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
    set(_code_check_targets "docs" "doxygen_docs"
                            "check" "style" 
                            "clangformat_check" "clangformat_style" )
    foreach(_tgt ${_code_check_targets})
        if(TARGET ${_tgt})
            set_target_properties(${_tgt} PROPERTIES FOLDER "axom/code_checks")
        endif()
    endforeach()
endif()
