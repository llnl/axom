# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Doxygen documentation targets
#------------------------------------------------------------------------------

blt_add_doxygen_target( axom_doxygen )

# Required directory for doxygen to generate output to 
add_custom_target( build_doxygen_output_dir 
                   ALL 
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/html 
                 ) 

add_dependencies( axom_doxygen build_doxygen_output_dir )

