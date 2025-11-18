# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# SIDRE -- SI(mulation) D(ata) RE(pository)
#
# Also includes functionality for SPIO -- S(idre) P(arallel) I/O
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       Sidre
                        COMPONENTS SLIC
                        TPLS       Conduit )

#------------------------------------------------------------------------------
# Add the sidre sources
#------------------------------------------------------------------------------
set(sidre_headers
    core/Array.hpp
    core/Buffer.hpp
    core/ConduitMemory.hpp
    core/Group.hpp
    core/DataStore.hpp
    core/View.hpp
    core/Attribute.hpp
    core/AttrValues.hpp
    core/SidreTypes.hpp
    core/SidreDataTypeIds.h )

set(sidre_sources
    core/Buffer.cpp
    core/ConduitMemory.cpp
    core/Group.cpp
    core/DataStore.cpp
    core/View.cpp
    core/Attribute.cpp
    core/AttrValues.cpp )

# Add spio headers and sources when MPI is available
if(AXOM_ENABLE_MPI)
   list(APPEND sidre_headers
        spio/IOManager.hpp
        spio/IOBaton.hpp )
   list(APPEND sidre_sources
        spio/IOManager.cpp
        spio/IOBaton.cpp )
endif()

if(SHROUD_FOUND)
    add_subdirectory(interface)
    add_subdirectory(spio/interface)
endif()

if (ENABLE_FORTRAN)
    # Sidre's C API is only built to provide a Fortran interface.
    list(APPEND sidre_headers
        interface/sidre.h
        interface/SidreTypes.h
        # generated headers
        interface/c_fortran/wrapSidre.h
        interface/c_fortran/wrapDataStore.h
        interface/c_fortran/wrapGroup.h
        interface/c_fortran/wrapBuffer.h
        interface/c_fortran/wrapView.h
        interface/c_fortran/typesSidre.h
    )
    list(APPEND sidre_sources
        # generated source
        interface/c_fortran/wrapSidre.cpp
        interface/c_fortran/wrapDataStore.cpp
        interface/c_fortran/wrapGroup.cpp
        interface/c_fortran/wrapBuffer.cpp
        interface/c_fortran/wrapView.cpp
        interface/c_fortran/utilSidre.cpp
        interface/c_fortran/wrapfsidre.F
    )

    # Add spio's C API when MPI is available
    if(AXOM_ENABLE_MPI)
        list(APPEND sidre_headers
             spio/interface/c_fortran/wrapIOManager.h
             spio/interface/c_fortran/typesSPIO.h
        )
        list(APPEND sidre_sources
             spio/interface/c_fortran/wrapfspio.f
             spio/interface/c_fortran/wrapIOManager.cpp
             spio/interface/c_fortran/utilSPIO.cpp
        )
    endif()
endif()

if (MFEM_FOUND AND AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)
    list(APPEND sidre_headers core/MFEMSidreDataCollection.hpp)
    list(APPEND sidre_sources core/MFEMSidreDataCollection.cpp)
endif()

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(sidre_depends
    slic
    conduit::conduit
    )

blt_list_append(TO sidre_depends ELEMENTS hdf5 IF HDF5_FOUND)
blt_list_append(TO sidre_depends ELEMENTS sparsehash IF SPARSEHASH_FOUND)
if(AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION)
    blt_list_append(TO sidre_depends ELEMENTS mfem IF MFEM_FOUND)
    blt_list_append(TO sidre_depends ELEMENTS fmt)
endif()

# Include additional dependencies for spio when MPI is available
if(AXOM_ENABLE_MPI)
    list(APPEND sidre_depends conduit::conduit_mpi fmt)
    blt_list_append(TO sidre_depends ELEMENTS scr IF SCR_FOUND)
endif()


axom_add_library(NAME       sidre
                 SOURCES    ${sidre_sources}
                 HEADERS    ${sidre_headers}
                 DEPENDS_ON ${sidre_depends}
                 FOLDER     axom/sidre)

axom_write_unified_header(NAME    sidre
                          HEADERS ${sidre_headers})

axom_install_component(NAME    sidre
                       HEADERS ${sidre_headers}
                       )

if(ENABLE_FORTRAN AND AXOM_ENABLE_MPI)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_spio.mod DESTINATION lib/fortran)
endif()


if(NANOBIND_FOUND)
    nanobind_add_module(pysidre nanobind_sidre.cpp)
    target_link_libraries(pysidre PRIVATE sidre)

    # Use HIP executable linker flags for the python module
    # (CMake treats modules separately from executables,
    # executable flags not automatically applied)
    if(AXOM_ENABLE_HIP)
        # Make CMake compile the pysidre file with the HIP compiler.
        set_source_files_properties(nanobind_sidre.cpp PROPERTIES LANGUAGE HIP)

        string (REPLACE " " ";" MODULE_LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        target_link_options(pysidre PRIVATE ${MODULE_LINK_FLAGS})
    endif()

    install(TARGETS pysidre LIBRARY DESTINATION lib)
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
