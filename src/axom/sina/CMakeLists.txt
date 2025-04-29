# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina -- API to write data to Sina's common file format
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify necessary dependencies
#
# Note: Sina also optionally depends on Adiak when AXOM_USE_ADIAK=ON
#------------------------------------------------------------------------------
axom_component_requires(NAME       Sina
                        TPLS       Conduit )

#------------------------------------------------------------------------------
# Specify the sina headers/sources
#------------------------------------------------------------------------------
set(sina_headers
    core/ConduitUtil.hpp
    core/Curve.hpp
    core/CurveSet.hpp
    core/DataHolder.hpp
    core/Datum.hpp
    core/Document.hpp
    core/File.hpp
    core/ID.hpp
    core/Record.hpp
    core/Relationship.hpp
    core/Run.hpp
    )

set(sina_sources
    core/ConduitUtil.cpp
    core/Curve.cpp
    core/CurveSet.cpp
    core/DataHolder.cpp
    core/Datum.cpp
    core/Document.cpp
    core/File.cpp
    core/ID.cpp
    core/Record.cpp
    core/Relationship.cpp
    core/Run.cpp
    )


# Add Adiak header and source
blt_list_append( TO sina_headers ELEMENTS core/AdiakWriter.hpp IF AXOM_USE_ADIAK )
blt_list_append( TO sina_sources ELEMENTS core/AdiakWriter.cpp IF AXOM_USE_ADIAK )

# Add fortran interface for Sina
if(AXOM_USE_HDF5 AND ENABLE_FORTRAN)
  set(AXOM_USE_HDF5_FORTRAN ".true.")
else()
  set(AXOM_USE_HDF5_FORTRAN ".false.")
endif()

if (ENABLE_FORTRAN)
  configure_file(interface/sina_fortran_interface.f.in 
                 interface/sina_fortran_interface.f 
                 @ONLY)

  blt_list_append(TO sina_headers 
                  ELEMENTS interface/sina_fortran_interface.h)
  blt_list_append(TO sina_sources 
                  ELEMENTS interface/sina_fortran_interface.cpp
                           ${CMAKE_CURRENT_BINARY_DIR}/interface/sina_fortran_interface.f)
  
  install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/sina_functions.mod
                ${CMAKE_Fortran_MODULE_DIRECTORY}/hdf5_config.mod
          DESTINATION interface)
endif()

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(sina_depends
    core
    conduit::conduit
    )

blt_list_append( TO sina_depends ELEMENTS adiak::adiak IF AXOM_USE_ADIAK )

axom_add_library(NAME       sina
                 SOURCES    ${sina_sources}
                 HEADERS    ${sina_headers}
                 DEPENDS_ON ${sina_depends}
                 FOLDER     axom/sina)

axom_write_unified_header(NAME    sina
                          HEADERS ${sina_headers})

axom_install_component(NAME    sina
                       HEADERS ${sina_headers})

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if(AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
