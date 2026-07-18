# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Lumberjack component
#------------------------------------------------------------------------------

axom_add_executable(NAME lumberjack_basic_ex
                    SOURCES basicExample.cpp
                    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                    DEPENDS_ON lumberjack 
                    FOLDER axom/lumberjack/examples )

axom_add_executable(NAME lumberjack_noncollective_messages_ex
                    SOURCES nonCollectiveMessagesExample.cpp
                    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                    DEPENDS_ON lumberjack 
                    FOLDER axom/lumberjack/examples )
