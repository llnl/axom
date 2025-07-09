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
blt_list_append( TO sina_headers ELEMENTS core/AdiakWriter.hpp IF ADIAK_FOUND )
blt_list_append( TO sina_sources ELEMENTS core/AdiakWriter.cpp IF ADIAK_FOUND )

# Add fortran interface for Sina
if (ENABLE_FORTRAN)
  if(HDF5_FOUND)
    set(AXOM_USE_HDF5_FORTRAN ".true.")
  else()
    set(AXOM_USE_HDF5_FORTRAN ".false.")
  endif()

  axom_configure_file(
      interface/sina_fortran_interface.f.in
      ${CMAKE_CURRENT_BINARY_DIR}/interface/sina_fortran_interface.f)

  blt_list_append(TO sina_headers 
                  ELEMENTS interface/sina_fortran_interface.h)
  blt_list_append(TO sina_sources 
                  ELEMENTS interface/sina_fortran_interface.cpp
                           ${CMAKE_CURRENT_BINARY_DIR}/interface/sina_fortran_interface.f)
  endif()

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(sina_depends
    core
    conduit::conduit
    )

blt_list_append( TO sina_depends ELEMENTS adiak::adiak IF ADIAK_FOUND )

axom_add_library(NAME       sina
                 SOURCES    ${sina_sources}
                 HEADERS    ${sina_headers}
                 DEPENDS_ON ${sina_depends}
                 FOLDER     axom/sina)

axom_write_unified_header(NAME    sina
                          HEADERS ${sina_headers})

axom_install_component(NAME    sina
                       HEADERS ${sina_headers})
if (ENABLE_FORTRAN)
  install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/sina_functions.mod
                ${CMAKE_Fortran_MODULE_DIRECTORY}/sina_hdf5_config.mod
          DESTINATION lib/fortran)
endif()

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if(AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
