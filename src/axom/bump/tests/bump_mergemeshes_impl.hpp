// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*! 
 * \file bump_mergemeshes_impl.hpp
 * \brief Implementation shared by the MergeMeshes execution-policy tests.
 */

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump.hpp"

#include "axom/bump/tests/blueprint_testing_helpers.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_relay_io_blueprint.hpp>

#include <iostream>
#include <algorithm>
#include <iostream>

// NOTE: Conduit 0.9.6 and later has macros but Axom is not quite on that version yet.
#include "conduit/conduit_config.h"
#define AXOM_CONDUIT_MAKE_VERSION_VALUE(MAJOR, MINOR, PATCH) \
  (((MAJOR) * 10000) + ((MINOR) * 100) + (PATCH))
#define AXOM_CONDUIT_VERSION_VALUE \
  AXOM_CONDUIT_MAKE_VERSION_VALUE(CONDUIT_VERSION_MAJOR, CONDUIT_VERSION_MINOR, CONDUIT_VERSION_PATCH)

namespace bump = axom::bump;
namespace utils = axom::bump::utilities;

namespace
{
void conduit_debug_err_handler(const std::string& s1, const std::string& s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1);
}
}  // namespace

//#define AXOM_DEBUG_MERGE_MESHES_TEST
#ifdef AXOM_DEBUG_MERGE_MESHES_TEST
inline void saveMesh(const conduit::Node& n_mesh, const std::string& fileRoot)
{
  #ifdef CONDUIT_RELAY_IO_HDF5_ENABLED
  const std::string protocol("hdf5");
  conduit::relay::io::save(n_mesh, fileRoot + ".yaml", "yaml");
  #else
  const std::string protocol("yaml");
  #endif
  // These save with a ".root" extension
  conduit::relay::io::blueprint::save_mesh(n_mesh, fileRoot, protocol);
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_mergemeshes
{
  static void test()
  {
    std::vector<std::string> matsetTypes {"unibuffer", "element_dominant", "material_dominant"};
    for(const auto& matsetType : matsetTypes)
    {
      for(int matflags = 3; matflags >= 1; matflags--)
      {
        SLIC_INFO(axom::fmt::format("test({}, {})", matsetType, matflags));
        test(matsetType, matflags);
      }
    }
  }

  static void test(const std::string& matsetType, int matflags)
  {
    conduit::Node hostMesh;
    create(hostMesh, matsetType, matflags);
#ifdef AXOM_DEBUG_MERGE_MESHES_TEST
    const auto preMergeFilename = axom::fmt::format("preMerge_{}_{}", matflags, matsetType);
    saveMesh(hostMesh, preMergeFilename);
#endif

    // host->device
    conduit::Node deviceMesh;
    utils::copy<ExecSpace>(deviceMesh, hostMesh);

    // The node names for input 1 in the final merged mesh.
    const axom::IndexType nodeMap[] = {1, 2, 5, 6, 9, 10, 13, 14, 16, 17};
    // The 2 nodes in input 1 that do not appear in input 0
    const axom::IndexType nodeSlice[] = {8, 9};
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<axom::IndexType> deviceNodeMap(10, 10, allocatorID);
    axom::Array<axom::IndexType> deviceNodeSlice(2, 2, allocatorID);
    axom::copy(deviceNodeMap.data(), nodeMap, 10 * sizeof(axom::IndexType));
    axom::copy(deviceNodeSlice.data(), nodeSlice, 2 * sizeof(axom::IndexType));

    // Set up inputs.
    // _bump_utilities_mergemeshes_begin
    std::vector<bump::MeshInput> inputs(2);
    inputs[0].m_input = deviceMesh.fetch_ptr("domain0000");

    inputs[1].m_input = deviceMesh.fetch_ptr("domain0001");
    inputs[1].m_nodeMapView = deviceNodeMap.view();
    inputs[1].m_nodeSliceView = deviceNodeSlice.view();
    // _bump_utilities_mergemeshes_end

    // Execute
    conduit::Node opts, deviceResult;
    opts["topology"] = "mesh";
    bump::MergeMeshesAndMatsets<ExecSpace> mm;
    mm.execute(inputs, opts, deviceResult);

    // device->host
    conduit::Node hostResult;
    utils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);
#ifdef AXOM_DEBUG_MERGE_MESHES_TEST
    const auto postMergeFilename = axom::fmt::format("postMerge_{}_{}", matflags, matsetType);
    saveMesh(hostResult, postMergeFilename);
#endif
    constexpr double tolerance = 1.e-7;
    conduit::Node expectedResult, info;
    result(expectedResult, matflags);
    bool success = false;
    try
    {
      success = compareConduit(expectedResult, hostResult, tolerance, info);
    }
    catch(const conduit::Error& e)
    {
      e.print();
    }
    if(!success)
    {
      info.print();
      printNode(hostResult);
    }
    EXPECT_TRUE(success);
  }

  static void create(conduit::Node& mesh, const std::string& matsetType, int matflags)
  {
    const char* yaml = R"xx(
domain0000:
  coordsets:
    coords:
      type: explicit
      values:
        x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
        y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 3., 3., 3., 3.]
  topologies:
    mesh:
      type: unstructured
      coordset: coords
      elements:
        shape: quad
        connectivity: [0,1,5,4, 4,5,9,8, 8,9,13,12, 2,3,7,6, 6,7,11,10, 10,11,15,14]
        sizes: [4,4,4, 4,4,4]
        offsets: [0,4,8,12,16,20]
  fields:
    nodal:
      topology: mesh
      association: vertex
      values: [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0]
    zonal:
      topology: mesh
      association: element
      values: [0,1,2, 3,4,5]
    zonal_mixed:
      topology: mesh
      association: element
      matset: mat
      values: [100.1, 101.1, 102.1, 103.25, 104.25, 105.25]
      # matset_values encodes 100+zone.mat
      matset_values: [100.1, 101.1, 102.1, 103.2,103.3, 104.2,104.3, 105.2,105.3]
  matsets:
    mat:
      material_map:
        A: 1
        B: 2
        C: 3
      topology: mesh
      material_ids: [1, 1, 1, 2,3, 2,3, 2,3]
      volume_fractions: [1., 1., 1., 0.5,0.5, 0.5,0.5, 0.5,0.5]
      indices: [0, 1, 2, 3,4, 5,6, 7,8]
      sizes: [1, 1, 1, 2, 2, 2]
      offsets: [0, 1, 2, 3, 5, 7]
domain0001:
  coordsets:
    coords:
      type: explicit
      values:
        x: [1., 2., 1., 2., 1., 2., 1., 2., 1.5, 1.5]
        y: [0., 0., 1., 1., 2., 2., 3., 3., 0.5, 1.5]
  topologies:
    mesh:
      type: unstructured
      coordset: coords
      elements:
        shape: mixed
        shape_map:
          quad: 3
          tri: 2
        connectivity: [0,8,2, 0,1,8, 1,3,8, 8,3,2, 2,9,4, 2,3,9, 3,5,9, 5,4,9, 4,5,7,6]
        sizes: [3,3,3,3, 3,3,3,3, 4]
        offsets: [0,3,6,9, 12,15,18,21, 24]
        shapes: [2,2,2,2, 2,2,2,2, 3]
  fields:
    nodal:
      topology: mesh
      association: vertex
      values: [1,1,1,1,1,1,1,1, 2,2]
    zonal:
      topology: mesh
      association: element
      values: [0,1,2,3, 4,5,6,7, 8]
    zonal_mixed:
      topology: mesh
      association: element
      matset: mat
      values: [-200.1, -201.15, -202.3, -203.15, -204.1, -205.15, -206.2, -207.15, -208.15]
      # matset_values encodes -200+zone.mat
      matset_values: [-200.1, -201.1,-201.2, -202.3, -203.1,-203.2, -204.1, -205.1,-205.2, -206.2, -207.1,-207.2, -208.1,-208.2]
  matsets:
    mat:
      material_map:
        A: 1
        B: 2
        #C: 3
      topology: mesh
      material_ids: [1, 1,2, 2, 1,2, 1, 1,2, 2, 1,2, 1,2]
      volume_fractions: [1., 0.5,0.5, 1., 0.5,0.5, 1., 0.5,0.5, 1., 0.5,0.5, 0.5,0.5]
      indices: [0, 1,2, 3, 4,5, 6, 7,8, 9, 10,11, 12,13]
      sizes: [1, 2, 1, 2, 1, 2, 1, 2, 2]
      offsets: [0, 1, 3, 4, 6, 7, 9, 10, 12]
)xx";
    mesh.parse(yaml);

    for(int dom = 0; dom < 2; dom++)
    {
      changeMatsetType(mesh[dom], matsetType);
    }
    applyMatFlags(mesh, matflags);
  }

  /// Remove mixed field and matset on domains according to matflags. This tests merging domains that are missing materials/fields.
  static void applyMatFlags(conduit::Node& mesh, int matflags)
  {
    for(int dom = 0; dom < 2; dom++)
    {
      if(!axom::utilities::bitIsSet(matflags, dom))
      {
        conduit::Node& domain = mesh[dom];
        domain["fields"].remove("zonal_mixed");
        domain.remove("matsets");
      }
    }
  }

  static void changeMatsetType(conduit::Node& domain, const std::string& matsetType)
  {
    // Change the material and field representations
    if(matsetType == "element_dominant")
    {
      conduit::Node domainCopy(domain);
      conduit::Node& srcMatset = domainCopy["matsets/mat"];
      conduit::Node& srcField = domainCopy["fields/zonal_mixed"];

      domain.remove("matsets/mat");
      domain.remove("fields/zonal_mixed");

      // These functions changed in Conduit 0.9.6
#if AXOM_CONDUIT_VERSION_VALUE < AXOM_CONDUIT_MAKE_VERSION_VALUE(0, 9, 6)
      conduit::blueprint::mesh::matset::to_multi_buffer_full(srcMatset, domain["matsets/mat"]);
      conduit::blueprint::mesh::field::to_multi_buffer_full(srcMatset,
                                                            srcField,
                                                            "mat",
                                                            domain["fields/zonal_mixed"]);
#else
      conduit::blueprint::mesh::matset::to_multi_buffer_by_element(srcMatset, domain["matsets/mat"]);
      conduit::blueprint::mesh::field::to_multi_buffer_by_element(srcMatset,
                                                                  srcField,
                                                                  "mat",
                                                                  domain["fields/zonal_mixed"]);
#endif
      // Make sure we preserve the material_map.
      if(srcMatset.has_child("material_map"))
      {
        domain["matsets/mat/material_map"].set(srcMatset["material_map"]);
      }
    }
    else if(matsetType == "material_dominant")
    {
      conduit::Node domainCopy(domain);
      conduit::Node& srcMatset = domainCopy["matsets/mat"];
      conduit::Node& srcField = domainCopy["fields/zonal_mixed"];

      domain.remove("matsets/mat");
      domain.remove("fields/zonal_mixed");

      conduit::blueprint::mesh::matset::to_multi_buffer_by_material(srcMatset, domain["matsets/mat"]);
      conduit::blueprint::mesh::field::to_multi_buffer_by_material(srcMatset,
                                                                   srcField,
                                                                   "mat",
                                                                   domain["fields/zonal_mixed"]);
      // Make sure we preserve the material_map.
      if(srcMatset.has_child("material_map"))
      {
        domain["matsets/mat/material_map"].set(srcMatset["material_map"]);
      }
    }
  }

  static void result(conduit::Node& mesh, int matflags)
  {
    // NOTE: We pass back different baselines for different matflags. The fields and matset change.
    //       It is simpler to just have a totally separate baseline to parse.

    // Result for matflags=3 - both input domains had the material and the mixed field.
    const char* yaml3 = R"xx(
coordsets: 
  coords: 
    type: "explicit"
    values: 
      x: [0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 1.5, 1.5]
      y: [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 0.5, 1.5]
topologies: 
  mesh: 
    type: "unstructured"
    coordset: "coords"
    elements: 
      connectivity: [0, 1, 5, 4, 4, 5, 9, 8, 8, 9, 13, 12, 2, 3, 7, 6, 6, 7, 11, 10, 10, 11, 15, 14, 1, 16, 5, 1, 2, 16, 2, 6, 16, 16, 6, 5, 5, 17, 9, 5, 6, 17, 6, 10, 17, 10, 9, 17, 9, 10, 14, 13]
      sizes: [4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4]
      offsets: [0, 4, 8, 12, 16, 20, 24, 27, 30, 33, 36, 39, 42, 45, 48]
      shape: "mixed"
      shape_map: 
        quad: 3
        tri: 2
      shapes: [3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3]
matsets: 
  mat: 
    topology: "mesh"
    volume_fractions: [1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5]
    material_ids: [0, 0, 0, 1, 2, 1, 2, 1, 2, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1]
    sizes: [1, 1, 1, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2]
    offsets: [0, 1, 2, 3, 5, 7, 9, 10, 12, 13, 15, 16, 18, 19, 21]
    indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    material_map: 
      A: 0
      B: 1
      C: 2
fields: 
  nodal: 
    association: "vertex"
    topology: "mesh"
    values: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2]
  zonal: 
    association: "element"
    topology: "mesh"
    values: [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 7, 8]
  zonal_mixed: 
    association: "element"
    topology: "mesh"
    matset: "mat"
    values: [100.1, 101.1, 102.1, 103.25, 104.25, 105.25, -200.1, -201.15, -202.3, -203.15, -204.1, -205.15, -206.2, -207.15, -208.15]
    matset_values: [100.1, 101.1, 102.1, 103.2, 103.3, 104.2, 104.3, 105.2, 105.3, -200.1, -201.1, -201.2, -202.3, -203.1, -203.2, -204.1, -205.1, -205.2, -206.2, -207.1, -207.2, -208.1, -208.2]
)xx";

    // Result for matflags=2 - domain 0 lacked the material and domain so we get default values where domain 0's data would be.
    const char* yaml2 = R"xx(
coordsets: 
  coords: 
    type: "explicit"
    values: 
      x: [0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 1.5, 1.5]
      y: [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 0.5, 1.5]
topologies: 
  mesh: 
    type: "unstructured"
    coordset: "coords"
    elements: 
      connectivity: [0, 1, 5, 4, 4, 5, 9, 8, 8, 9, 13, 12, 2, 3, 7, 6, 6, 7, 11, 10, 10, 11, 15, 14, 1, 16, 5, 1, 2, 16, 2, 6, 16, 16, 6, 5, 5, 17, 9, 5, 6, 17, 6, 10, 17, 10, 9, 17, 9, 10, 14, 13]
      sizes: [4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4]
      offsets: [0, 4, 8, 12, 16, 20, 24, 27, 30, 33, 36, 39, 42, 45, 48]
      shape: "mixed"
      shape_map: 
        quad: 3
        tri: 2
      shapes: [3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3]
matsets: 
  mat: 
    topology: "mesh"
    volume_fractions: [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 0.5, 0.5]
    material_ids: [2, 2, 2, 2, 2, 2, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1]
    sizes: [1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2]
    offsets: [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 13, 15, 16, 18]
    indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    material_map: 
      A: 0
      B: 1
      default: 2
fields: 
  nodal: 
    association: "vertex"
    topology: "mesh"
    values: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2]
  zonal: 
    association: "element"
    topology: "mesh"
    values: [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 7, 8]
  zonal_mixed: 
    association: "element"
    topology: "mesh"
    matset: "mat"
    values: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -200.1, -201.15, -202.3, -203.15, -204.1, -205.15, -206.2, -207.15, -208.15]
    matset_values: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -200.1, -201.1, -201.2, -202.3, -203.1, -203.2, -204.1, -205.1, -205.2, -206.2, -207.1, -207.2, -208.1, -208.2]
)xx";

    // Result for matflags=1 - domain 1 lacked the material and domain so we get default values where domain 1's data would be.
    const char* yaml1 = R"xx(
coordsets: 
  coords: 
    type: "explicit"
    values: 
      x: [0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 1.5, 1.5]
      y: [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 0.5, 1.5]
topologies: 
  mesh: 
    type: "unstructured"
    coordset: "coords"
    elements: 
      connectivity: [0, 1, 5, 4, 4, 5, 9, 8, 8, 9, 13, 12, 2, 3, 7, 6, 6, 7, 11, 10, 10, 11, 15, 14, 1, 16, 5, 1, 2, 16, 2, 6, 16, 16, 6, 5, 5, 17, 9, 5, 6, 17, 6, 10, 17, 10, 9, 17, 9, 10, 14, 13]
      sizes: [4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4]
      offsets: [0, 4, 8, 12, 16, 20, 24, 27, 30, 33, 36, 39, 42, 45, 48]
      shape: "mixed"
      shape_map: 
        quad: 3
        tri: 2
      shapes: [3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3]
matsets: 
  mat: 
    topology: "mesh"
    volume_fractions: [1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    material_ids: [0, 0, 0, 1, 2, 1, 2, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    sizes: [1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    offsets: [0, 1, 2, 3, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    material_map: 
      A: 0
      B: 1
      C: 2
      default: 3
fields: 
  nodal: 
    association: "vertex"
    topology: "mesh"
    values: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2]
  zonal: 
    association: "element"
    topology: "mesh"
    values: [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 7, 8]
  zonal_mixed: 
    association: "element"
    topology: "mesh"
    matset: "mat"
    values: [100.1, 101.1, 102.1, 103.25, 104.25, 105.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    matset_values: [100.1, 101.1, 102.1, 103.2, 103.3, 104.2, 104.3, 105.2, 105.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)xx";

    switch(matflags)
    {
    case 3:
      mesh.parse(yaml3);
      break;
    case 2:
      mesh.parse(yaml2);
      break;
    case 1:
      mesh.parse(yaml1);
      break;
    default:
      SLIC_ERROR("Unsupported matflags value.");
      break;
    }
  }
};
