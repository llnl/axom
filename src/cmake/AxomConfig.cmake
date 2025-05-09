# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# AxomConfig.cmake - Create header of configuration options
#------------------------------------------------------------------------------

## Get Axom version information
message(STATUS "Configuring Axom version ${AXOM_VERSION_FULL}")


## Add a definition to the generated config file for each library dependency
## (optional and built-in) that we might need to know about in the code. We
## check for vars of the form <DEP>_FOUND or ENABLE_<DEP>
set(TPL_DEPS ADIAK C2C CALIPER CAMP CLI11 CONDUIT CUDA FMT HIP HDF5 LUA MFEM MPI OPENMP OPENCASCADE RAJA SCR SOL SPARSEHASH UMPIRE ZLIB)
foreach(dep ${TPL_DEPS})
    if( ${dep}_FOUND OR ENABLE_${dep} )
        set(AXOM_USE_${dep} TRUE  )
    endif()
endforeach()

# Handle MPI Fortran headers
if(AXOM_ENABLE_MPI AND ENABLE_FORTRAN)
  if(MPI_Fortran_USE_MPIF)
    set(AXOM_USE_MPIF_HEADER TRUE)
  endif()
endif()


## Add a configuration define for each enabled axom component
foreach(comp ${AXOM_COMPONENTS_ENABLED})
    string(TOUPPER ${comp} comp_uc)
    set(AXOM_USE_${comp_uc} TRUE)
endforeach()

## Add compile-time options to the config file
## Check for options of the form AXOM_ENABLE_<OPTION> and sets a corresponding 
## AXOM_USE_<OPTION> accordingly to use when generating the axom/config.hpp
foreach(option MPI3)
    if( AXOM_ENABLE_${option} )
      set(AXOM_USE_${option} TRUE)
    endif()
endforeach()

# Handle paths
convert_to_native_escaped_file_path(${PROJECT_SOURCE_DIR} AXOM_SRC_DIR_NATIVE)
convert_to_native_escaped_file_path(${PROJECT_BINARY_DIR} AXOM_BIN_DIR_NATIVE)
if(AXOM_DATA_DIR)
  convert_to_native_escaped_file_path(${AXOM_DATA_DIR} AXOM_DATA_DIR_NATIVE)
endif()

#------------------------------------------------------------------------------
# Compiler and language related configuration variables
#------------------------------------------------------------------------------

if(ENABLE_FORTRAN)

    file(WRITE ${PROJECT_BINARY_DIR}/c_loc_with_assumed_shape.f "
! This is expected to fail with gcc 4.7.1
! Error: Assumed-shape array 'arg' at (1) cannot be an argument to the
!      procedure 'c_loc' because it is not C interoperable
! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53945
      program main
      end program main
      subroutine test(arg, addr)
        use iso_c_binding
        integer(C_INT), target, intent(IN) :: arg(:)
        type(C_ptr), intent(OUT) :: addr
        addr = C_LOC(arg)
      end subroutine test
    ")
    try_compile(
        USE_C_LOC_WITH_ASSUMED_SHAPE
        ${PROJECT_BINARY_DIR}
        ${PROJECT_BINARY_DIR}/c_loc_with_assumed_shape.f
    )

endif(ENABLE_FORTRAN)

# The following MSVC pragmas caused problems for the Cray Fortran compiler,
# so only add them when using MSVC
if(COMPILER_FAMILY_IS_MSVC)
  set(AXOM_MSVC_PRAGMAS [=[
#if defined(_MSC_VER)
  /* Turn off warning about lack of DLL interface */
  #pragma warning(disable:4251)
  /* Turn off warning non-dll class is base for dll-interface class */
  #pragma warning(disable:4275)
  /* Turn off warning about identifier truncation */
  #pragma warning(disable:4786)  
#endif  /* defined(_MSC_VER) */
]=])
endif()


# Generate the configuration header file
axom_configure_file(
    ${PROJECT_SOURCE_DIR}/axom/config.hpp.in
    ${PROJECT_BINARY_DIR}/include/axom/config.hpp
)

install(FILES ${PROJECT_BINARY_DIR}/include/axom/config.hpp DESTINATION include/axom)

#------------------------------------------------------------------------------
# Generate axom-config.cmake for importing Axom into other CMake packages
#------------------------------------------------------------------------------

# Set up some paths, preserve existing cache values (if present)
set(AXOM_INSTALL_INCLUDE_DIR "include" CACHE STRING "")
set(AXOM_INSTALL_CONFIG_DIR "lib" CACHE STRING "")
set(AXOM_INSTALL_LIB_DIR "lib" CACHE STRING "")
set(AXOM_INSTALL_BIN_DIR "bin" CACHE STRING "")
set(AXOM_INSTALL_CMAKE_MODULE_DIR "${AXOM_INSTALL_CONFIG_DIR}/cmake" CACHE STRING "")

convert_to_native_escaped_file_path(${CMAKE_INSTALL_PREFIX} AXOM_INSTALL_PREFIX)
set(AXOM_INSTALL_PREFIX ${AXOM_INSTALL_PREFIX} CACHE STRING "" FORCE)


include(CMakePackageConfigHelpers)

# Add version helper
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/axom-config-version.cmake
    VERSION ${AXOM_VERSION_FULL}
    COMPATIBILITY AnyNewerVersion
)

# Set up cmake package config file
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/axom-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/axom-config.cmake
  INSTALL_DESTINATION 
    ${AXOM_INSTALL_CONFIG_DIR}
  PATH_VARS
    AXOM_INSTALL_INCLUDE_DIR
    AXOM_INSTALL_LIB_DIR
    AXOM_INSTALL_BIN_DIR
    AXOM_INSTALL_CMAKE_MODULE_DIR
  )

# Install config files
install(
  FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/axom-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/axom-config-version.cmake
  DESTINATION 
    ${AXOM_INSTALL_CMAKE_MODULE_DIR}
)

# Install BLT files that recreate BLT targets in downstream projects
blt_install_tpl_setups(DESTINATION ${AXOM_INSTALL_CMAKE_MODULE_DIR})

#------------------------------------------------------------------------------
# Create a list of exported targets so that other projects that include Axom
# can add it to their CMake export list
#------------------------------------------------------------------------------

# Add it to a temporary list before creating the cache variable to use list(APPEND)
set(_axom_exported_targets ${AXOM_COMPONENTS_ENABLED})

set(_optional_targets cli11 fmt hdf5 lua sol sparsehash)
foreach(_tar ${_optional_targets})
    string(TOUPPER ${_tar} _upper_tar)
    if(ENABLE_${_upper_tar} OR ${_upper_tar}_FOUND)
        list(APPEND _axom_exported_targets ${_tar})
    endif()
endforeach()

set(axom_exported_targets ${_axom_exported_targets} CACHE INTERNAL "")
