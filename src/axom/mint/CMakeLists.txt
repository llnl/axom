# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Mint component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#
# Note: Mint also optionally depends on Sidre when AXOM_MINT_USE_SIDRE=ON
#------------------------------------------------------------------------------
axom_component_requires(NAME       Mint
                        COMPONENTS SLIC SLAM)

#------------------------------------------------------------------------------
# Add mint options
#------------------------------------------------------------------------------
cmake_dependent_option( AXOM_MINT_USE_SIDRE "Enables Sidre support in Mint" ON
                       "AXOM_ENABLE_SIDRE" OFF)


#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
## generate mint config
axom_configure_file ( core/config.hpp.in
                      ${PROJECT_BINARY_DIR}/include/axom/mint/config.hpp )

set( mint_headers

    ## Deprecated containers
    deprecated/MCArray.hpp

    ## exec
    execution/xargs.hpp
    execution/interface.hpp
    execution/internal/for_all_cells.hpp
    execution/internal/for_all_nodes.hpp
    execution/internal/for_all_faces.hpp
    execution/internal/helpers.hpp

    ## FEM classes
    fem/FEBasis.hpp
    fem/FEBasisTypes.hpp
    fem/FiniteElement.hpp
    fem/shape_functions/Lagrange.hpp
    fem/shape_functions/ShapeFunction.hpp

    ## Lagrange Basis
    fem/shape_functions/Lagrange.hpp
    fem/shape_functions/lagrange/lagrange_hexa_27.hpp
    fem/shape_functions/lagrange/lagrange_hexa_8.hpp
    fem/shape_functions/lagrange/lagrange_prism_6.hpp
    fem/shape_functions/lagrange/lagrange_pyra_5.hpp
    fem/shape_functions/lagrange/lagrange_quad_4.hpp
    fem/shape_functions/lagrange/lagrange_quad_9.hpp
    fem/shape_functions/lagrange/lagrange_tetra_4.hpp
    fem/shape_functions/lagrange/lagrange_tri_3.hpp

    ## Internal implementation
    mesh/internal/ConnectivityArray_typed_indirection.hpp
    mesh/internal/ConnectivityArrayHelpers.hpp
    mesh/internal/MeshHelpers.hpp

    ## Mesh data-model
    mesh/blueprint.hpp
    mesh/CellTypes.hpp
    mesh/ConnectivityArray.hpp
    mesh/CurvilinearMesh.hpp
    mesh/Field.hpp
    mesh/FieldAssociation.hpp
    mesh/FieldData.hpp
    mesh/FieldTypes.hpp
    mesh/FieldVariable.hpp
    mesh/Mesh.hpp
    mesh/MeshCoordinates.hpp
    mesh/MeshTypes.hpp
    mesh/ParticleMesh.hpp
    mesh/RectilinearMesh.hpp
    mesh/StructuredMesh.hpp
    mesh/UniformMesh.hpp
    mesh/UnstructuredMesh.hpp

    ## utils
    utils/vtk_utils.hpp
    utils/su2_utils.hpp
    utils/ExternalArray.hpp
    utils/ArrayWrapper.hpp
   )

blt_list_append(TO mint_headers ELEMENTS deprecated/SidreMCArray.hpp IF AXOM_MINT_USE_SIDRE)

set( mint_sources

    ## FEM classes
    fem/FiniteElement.cpp

    ## Internal implementation
    mesh/internal/MeshHelpers.cpp

    ## Mesh data-model
    mesh/blueprint.cpp
    mesh/CurvilinearMesh.cpp
    mesh/FieldData.cpp
    mesh/Mesh.cpp
    mesh/MeshCoordinates.cpp
    mesh/ParticleMesh.cpp
    mesh/RectilinearMesh.cpp
    mesh/StructuredMesh.cpp
    mesh/UniformMesh.cpp

    ## utils
    utils/vtk_utils.cpp
    utils/su2_utils.cpp
   )

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
set( mint_dependencies
     slic
     )

if (AXOM_MINT_USE_SIDRE)
   list( APPEND mint_dependencies sidre)
   blt_list_append( TO mint_dependencies ELEMENTS hdf5 IF HDF5_FOUND )
endif()

blt_list_append( TO mint_dependencies ELEMENTS RAJA IF RAJA_FOUND )

axom_add_library(
    NAME       mint
    SOURCES    ${mint_sources}
    HEADERS    ${mint_headers}
    DEPENDS_ON ${mint_dependencies}
    FOLDER     axom/mint
   )

axom_write_unified_header(NAME    mint
                          HEADERS config.hpp ${mint_headers})

axom_install_component(NAME    mint
                       HEADERS ${mint_headers}
                       )
install(FILES       ${PROJECT_BINARY_DIR}/include/axom/mint/config.hpp
        DESTINATION include/axom/mint
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
