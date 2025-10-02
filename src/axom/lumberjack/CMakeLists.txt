# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Lumberjack component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(lumberjack_headers
    BinaryTreeCommunicator.hpp
    Combiner.hpp
    Communicator.hpp
    Lumberjack.hpp
    Message.hpp
    MPIUtility.hpp
    RootCommunicator.hpp
    NonCollectiveRootCommunicator.hpp
    TextEqualityCombiner.hpp
    TextTagCombiner.hpp
    LineFileTagCombiner.hpp
    )

set(lumberjack_sources
    BinaryTreeCommunicator.cpp
    Lumberjack.cpp
    Message.cpp
    MPIUtility.cpp
    RootCommunicator.cpp
    NonCollectiveRootCommunicator.cpp
    )


#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
axom_add_library(NAME        lumberjack
                 SOURCES     ${lumberjack_sources}
                 HEADERS     ${lumberjack_headers}
                 DEPENDS_ON  core mpi
                 FOLDER      axom/lumberjack
                 )

axom_write_unified_header(NAME    lumberjack
                          HEADERS ${lumberjack_headers})

axom_install_component(NAME      lumberjack
                       HEADERS   ${lumberjack_headers}
                       )


#------------------------------------------------------------------------------
# Add examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()
