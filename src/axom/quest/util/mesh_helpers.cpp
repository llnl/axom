// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "mesh_helpers.hpp"
#if defined(AXOM_USE_SIDRE)
  #include <axom/sidre.hpp>
#endif
#include <axom/slic.hpp>
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#if defined(AXOM_USE_CONDUIT)
  #include <conduit/conduit_blueprint_mesh.hpp>
#endif
#include <iostream>

namespace axom
{
namespace quest
{
namespace util
{
#ifdef AXOM_USE_MFEM
mfem::Mesh* make_cartesian_mfem_mesh_2D(const primal::BoundingBox<double, 2>& bbox,
                                        const NumericArray<int, 2>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling)
{
  constexpr int DIM = 2;
  const auto range = bbox.range();

  auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian2D(res[0],
                                                          res[1],
                                                          mfem::Element::QUADRILATERAL,
                                                          true,
                                                          range[0],
                                                          range[1],
                                                          reorder_space_filling));

  // Offset the mesh to lie w/in the bounding box
  const int NV = mesh->GetNV();
  for(int i = 0; i < NV; ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < DIM; ++d)
    {
      v[d] += bbox.getMin()[d];
    }
  }

  // Ensure that mesh has high order nodes
  mesh->SetCurvature(polynomial_order);

  return mesh;
}

mfem::Mesh* make_cartesian_mfem_mesh_3D(const primal::BoundingBox<double, 3>& bbox,
                                        const NumericArray<int, 3>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling)
{
  constexpr int DIM = 3;
  const auto range = bbox.range();

  auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian3D(res[0],
                                                          res[1],
                                                          res[2],
                                                          mfem::Element::HEXAHEDRON,
                                                          range[0],
                                                          range[1],
                                                          range[2],
                                                          reorder_space_filling));

  // Offset the mesh to lie w/in the bounding box
  const int NV = mesh->GetNV();
  for(int i = 0; i < NV; ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < DIM; ++d)
    {
      v[d] += bbox.getMin()[d];
    }
  }

  // Ensure that mesh has high order nodes
  mesh->SetCurvature(polynomial_order);

  return mesh;
}

#endif  // AXOM_USE_MFEM

#if defined(AXOM_USE_SIDRE)

axom::sidre::Group* make_structured_blueprint_box_mesh_3d(axom::sidre::Group* meshGrp,
                                                          const primal::BoundingBox<double, 3>& bbox,
                                                          const NumericArray<int, 3>& res,
                                                          const std::string& topologyName,
                                                          const std::string& coordsetName,
                                                          axom::runtime_policy::Policy runtimePolicy)
{
  auto* topoGrp = meshGrp->createGroup("topologies")->createGroup(topologyName);
  SLIC_ERROR_IF(topoGrp == nullptr,
                "Cannot allocate topology '" + topologyName + "' in blueprint mesh '" +
                  meshGrp->getName() + "'.  It already exists.");

  auto* coordsetGrp = meshGrp->createGroup("coordsets")->createGroup(coordsetName);
  SLIC_ERROR_IF(coordsetGrp == nullptr,
                "Cannot allocate coordset '" + coordsetName + "' in blueprint mesh '" +
                  meshGrp->getName() + "'.  It already exists.");

  topoGrp->createView("type")->setString("structured");
  topoGrp->createView("coordset")->setString(coordsetName);
  auto* dimsGrp = topoGrp->createGroup("elements/dims");

  constexpr int DIM = 3;

  auto ni = res[0];
  auto nj = res[1];
  auto nk = res[2];

  const axom::StackArray<axom::IndexType, DIM> vertsShape {res[0] + 1, res[1] + 1, res[2] + 1};
  const auto numVerts = vertsShape[0] * vertsShape[1] * vertsShape[2];

  dimsGrp->createViewScalar("i", ni);
  dimsGrp->createViewScalar("j", nj);
  dimsGrp->createViewScalar("k", nk);

  coordsetGrp->createView("type")->setString("explicit");
  auto* valuesGrp = coordsetGrp->createGroup("values");
  auto* xVu = valuesGrp->createViewAndAllocate("x", axom::sidre::DataTypeId::FLOAT64_ID, numVerts);
  auto* yVu = valuesGrp->createViewAndAllocate("y", axom::sidre::DataTypeId::FLOAT64_ID, numVerts);
  auto* zVu = valuesGrp->createViewAndAllocate("z", axom::sidre::DataTypeId::FLOAT64_ID, numVerts);
  #ifdef AXOM_USE_UMPIRE
  SLIC_ASSERT(axom::getAllocatorIDFromPointer(xVu->getVoidPtr()) == meshGrp->getDefaultAllocatorID());
  SLIC_ASSERT(axom::getAllocatorIDFromPointer(yVu->getVoidPtr()) == meshGrp->getDefaultAllocatorID());
  SLIC_ASSERT(axom::getAllocatorIDFromPointer(zVu->getVoidPtr()) == meshGrp->getDefaultAllocatorID());
  #endif

  const axom::MDMapping<DIM> vertMapping(vertsShape, axom::ArrayStrideOrder::COLUMN);
  axom::ArrayView<double, DIM> xView(xVu->getData(), vertsShape, vertMapping.strides());
  axom::ArrayView<double, DIM> yView(yVu->getData(), vertsShape, vertMapping.strides());
  axom::ArrayView<double, DIM> zView(zVu->getData(), vertsShape, vertMapping.strides());

  fill_cartesian_coords_3d(runtimePolicy, bbox, xView, yView, zView);

  #if defined(AXOM_DEBUG) && defined(AXOM_USE_CONDUIT)
  {
    conduit::Node info;
    bool isValid = verifyBlueprintMesh(meshGrp, info);
    SLIC_ASSERT_MSG(isValid, "Internal error: Generated mesh is invalid.");
  }
  #endif

  return meshGrp;
}

axom::sidre::Group* make_structured_blueprint_box_mesh_2d(axom::sidre::Group* meshGrp,
                                                          const primal::BoundingBox<double, 2>& bbox,
                                                          const NumericArray<int, 2>& res,
                                                          const std::string& topologyName,
                                                          const std::string& coordsetName,
                                                          axom::runtime_policy::Policy runtimePolicy)
{
  auto* topoGrp = meshGrp->createGroup("topologies")->createGroup(topologyName);
  SLIC_ERROR_IF(topoGrp == nullptr,
                "Cannot allocate topology '" + topologyName + "' in blueprint mesh '" +
                  meshGrp->getName() + "'.  It already exists.");

  auto* coordsetGrp = meshGrp->createGroup("coordsets")->createGroup(coordsetName);
  SLIC_ERROR_IF(coordsetGrp == nullptr,
                "Cannot allocate coordset '" + coordsetName + "' in blueprint mesh '" +
                  meshGrp->getName() + "'.  It already exists.");

  topoGrp->createView("type")->setString("structured");
  topoGrp->createView("coordset")->setString(coordsetName);
  auto* dimsGrp = topoGrp->createGroup("elements/dims");

  constexpr int DIM = 2;

  auto ni = res[0];
  auto nj = res[1];

  const axom::StackArray<axom::IndexType, DIM> vertsShape {res[0] + 1, res[1] + 1};
  const auto numVerts = vertsShape[0] * vertsShape[1];

  dimsGrp->createViewScalar("i", ni);
  dimsGrp->createViewScalar("j", nj);

  coordsetGrp->createView("type")->setString("explicit");
  auto* valuesGrp = coordsetGrp->createGroup("values");
  auto* xVu = valuesGrp->createViewAndAllocate("x", axom::sidre::DataTypeId::FLOAT64_ID, numVerts);
  auto* yVu = valuesGrp->createViewAndAllocate("y", axom::sidre::DataTypeId::FLOAT64_ID, numVerts);
  #ifdef AXOM_USE_UMPIRE
  SLIC_ASSERT(axom::getAllocatorIDFromPointer(xVu->getVoidPtr()) == meshGrp->getDefaultAllocatorID());
  SLIC_ASSERT(axom::getAllocatorIDFromPointer(yVu->getVoidPtr()) == meshGrp->getDefaultAllocatorID());
  #endif

  const axom::MDMapping<DIM> vertMapping(vertsShape, axom::ArrayStrideOrder::COLUMN);
  axom::ArrayView<double, DIM> xView(xVu->getData(), vertsShape, vertMapping.strides());
  axom::ArrayView<double, DIM> yView(yVu->getData(), vertsShape, vertMapping.strides());

  fill_cartesian_coords_2d(runtimePolicy, bbox, xView, yView);

  #if defined(AXOM_DEBUG) && defined(AXOM_USE_CONDUIT)
  {
    conduit::Node info;
    bool isValid = verifyBlueprintMesh(meshGrp, info);
    SLIC_ASSERT_MSG(isValid, "Internal error: Generated mesh is invalid.");
  }
  #endif

  return meshGrp;
}

axom::sidre::Group* make_unstructured_blueprint_box_mesh_3d(axom::sidre::Group* meshGrp,
                                                            const primal::BoundingBox<double, 3>& bbox,
                                                            const NumericArray<int, 3>& res,
                                                            const std::string& topologyName,
                                                            const std::string& coordsetName,
                                                            axom::runtime_policy::Policy runtimePolicy)
{
  make_structured_blueprint_box_mesh_3d(meshGrp, bbox, res, topologyName, coordsetName, runtimePolicy);
  convert_blueprint_structured_explicit_to_unstructured_3d(meshGrp, topologyName, runtimePolicy);
  return meshGrp;
}

axom::sidre::Group* make_unstructured_blueprint_box_mesh_2d(axom::sidre::Group* meshGrp,
                                                            const primal::BoundingBox<double, 2>& bbox,
                                                            const NumericArray<int, 2>& res,
                                                            const std::string& topologyName,
                                                            const std::string& coordsetName,
                                                            axom::runtime_policy::Policy runtimePolicy)
{
  make_structured_blueprint_box_mesh_2d(meshGrp, bbox, res, topologyName, coordsetName, runtimePolicy);
  convert_blueprint_structured_explicit_to_unstructured_2d(meshGrp, topologyName, runtimePolicy);
  return meshGrp;
}

void convert_blueprint_structured_explicit_to_unstructured_3d(axom::sidre::Group* meshGrp,
                                                              const std::string& topoName,
                                                              axom::runtime_policy::Policy runtimePolicy)
{
  if(runtimePolicy == axom::runtime_policy::Policy::seq)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_3d<axom::SEQ_EXEC>(meshGrp, topoName);
  }
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  if(runtimePolicy == axom::runtime_policy::Policy::omp)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_3d<axom::OMP_EXEC>(meshGrp, topoName);
  }
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(runtimePolicy == axom::runtime_policy::Policy::cuda)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_3d<axom::CUDA_EXEC<256>>(meshGrp,
                                                                                        topoName);
  }
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(runtimePolicy == axom::runtime_policy::Policy::hip)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_3d<axom::HIP_EXEC<256>>(meshGrp,
                                                                                       topoName);
  }
  #endif
}

void convert_blueprint_structured_explicit_to_unstructured_2d(axom::sidre::Group* meshGrp,
                                                              const std::string& topoName,
                                                              axom::runtime_policy::Policy runtimePolicy)
{
  if(runtimePolicy == axom::runtime_policy::Policy::seq)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_2d<axom::SEQ_EXEC>(meshGrp, topoName);
  }
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  if(runtimePolicy == axom::runtime_policy::Policy::omp)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_2d<axom::OMP_EXEC>(meshGrp, topoName);
  }
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(runtimePolicy == axom::runtime_policy::Policy::cuda)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_2d<axom::CUDA_EXEC<256>>(meshGrp,
                                                                                        topoName);
  }
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(runtimePolicy == axom::runtime_policy::Policy::hip)
  {
    convert_blueprint_structured_explicit_to_unstructured_impl_2d<axom::HIP_EXEC<256>>(meshGrp,
                                                                                       topoName);
  }
  #endif
}

template <typename ExecSpace>
void convert_blueprint_structured_explicit_to_unstructured_impl_3d(axom::sidre::Group* meshGrp,
                                                                   const std::string& topoName)
{
  AXOM_ANNOTATE_SCOPE("convert_to_unstructured");
  // Note: MSVC required `static` to pass DIM to the axom::for_all w/ C++14
  // this restriction might be lifted w/ C++17
  static constexpr int DIM = 3;

  const std::string& coordsetName =
    meshGrp->getView("topologies/" + topoName + "/coordset")->getString();

  sidre::Group* coordsetGrp = meshGrp->getGroup("coordsets")->getGroup(coordsetName);
  SLIC_ASSERT(std::string(coordsetGrp->getView("type")->getString()) == "explicit");

  axom::sidre::Group* topoGrp = meshGrp->getGroup("topologies")->getGroup(topoName);
  axom::sidre::View* topoTypeView = topoGrp->getView("type");
  SLIC_ASSERT(std::string(topoTypeView->getString()) == "structured");
  topoTypeView->setString("unstructured");
  topoGrp->createView("elements/shape")->setString("hex");

  axom::sidre::Group* topoElemGrp = topoGrp->getGroup("elements");
  axom::sidre::Group* topoDimsGrp = topoElemGrp->getGroup("dims");

  // Assuming no ghost, but we should eventually support ghosts.
  SLIC_ASSERT(!topoGrp->hasGroup("elements/offsets"));
  SLIC_ASSERT(!topoGrp->hasGroup("elements/strides"));

  const axom::StackArray<axom::IndexType, DIM> cShape {
    axom::IndexType(topoDimsGrp->getView("i")->getNode().to_value()),
    axom::IndexType(topoDimsGrp->getView("j")->getNode().to_value()),
    axom::IndexType(topoDimsGrp->getView("k")->getNode().to_value())};
  const axom::StackArray<axom::IndexType, DIM> vShape {cShape[0] + 1, cShape[1] + 1, cShape[2] + 1};

  const axom::IndexType cCount = cShape[0] * cShape[1] * cShape[2];

  const axom::StackArray<axom::IndexType, 2> connShape {cCount, 8};
  axom::sidre::View* connView =
    topoGrp->createViewWithShapeAndAllocate("elements/connectivity",
                                            axom::sidre::detail::SidreTT<axom::IndexType>::id,
                                            2,
                                            connShape.begin());
  axom::ArrayView<axom::IndexType, 2> connArrayView(
    static_cast<axom::IndexType*>(connView->getVoidPtr()),
    connShape);

  axom::MDMapping<DIM> cIdMapping(cShape, axom::ArrayStrideOrder::COLUMN);
  axom::MDMapping<DIM> vIdMapping(vShape, axom::ArrayStrideOrder::COLUMN);

  const axom::StackArray<const axom::StackArray<axom::IndexType, DIM>, 8> cornerOffsets {
    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
  axom::for_all<ExecSpace>(
    cCount,
    AXOM_LAMBDA(axom::IndexType iCell) {
      axom::StackArray<axom::IndexType, DIM> cIdx = cIdMapping.toMultiIndex(iCell);
      for(int n = 0; n < 8; ++n)
      {
        const auto& cornerOffset = cornerOffsets[n];
        axom::StackArray<axom::IndexType, DIM> vIdx;
        for(int d = 0; d < DIM; ++d)
        {
          vIdx[d] = cIdx[d] + cornerOffset[d];
        }
        axom::IndexType iVert = vIdMapping.toFlatIndex(vIdx);
        connArrayView(iCell, n) = iVert;
      }
    });

  const bool addExtraDataForMint = true;
  if(addExtraDataForMint)
  {
    AXOM_ANNOTATE_BEGIN("add_extra");
    /*
      Constructing a mint mesh from meshGrp fails unless we add some
      extra data.  Blueprint doesn't require this extra data.  (The mesh
      passes conduit's Blueprint verification.)  This should be fixed,
      or we should write better blueprint support utilities.
    */
    /*
      Make the coordinate arrays 2D to use mint::Mesh.
      For some reason, mint::Mesh requires the arrays to be
      2D, even though the second dimension is always 1.
    */
    axom::IndexType curShape[2];
    int curDim;
    auto* valuesGrp = coordsetGrp->getGroup("values");
    curDim = valuesGrp->getView("x")->getShape(2, curShape);
    assert(curDim == 1);
    const axom::IndexType vertsShape[2] = {curShape[0], 1};
    valuesGrp->getView("x")->reshapeArray(2, vertsShape);
    valuesGrp->getView("y")->reshapeArray(2, vertsShape);
    valuesGrp->getView("z")->reshapeArray(2, vertsShape);

    // Make connectivity array 2D for the same reason.
    auto* elementsGrp = topoGrp->getGroup("elements");
    auto* connView = elementsGrp->getView("connectivity");
    curDim = connView->getShape(2, curShape);
    constexpr axom::IndexType NUM_VERTS_PER_HEX = 8;
    SLIC_ASSERT(curDim == 1 || curDim == 2);
    if(curDim == 1)
    {
      SLIC_ASSERT(curShape[0] % NUM_VERTS_PER_HEX == 0);
      axom::IndexType connShape[2] = {curShape[0] / NUM_VERTS_PER_HEX, NUM_VERTS_PER_HEX};
      connView->reshapeArray(2, connShape);
    }

    // mint::Mesh requires connectivity strides, even though Blueprint doesn't.
    elementsGrp->createViewScalar("stride", NUM_VERTS_PER_HEX);

    // mint::Mesh requires field group, even though Blueprint doesn't.
    meshGrp->createGroup("fields");
    AXOM_ANNOTATE_END("add_extra");
  }

  #if defined(AXOM_DEBUG) && defined(AXOM_USE_CONDUIT)
  AXOM_ANNOTATE_BEGIN("validate_post");
  conduit::Node info;
  bool isValid = verifyBlueprintMesh(meshGrp, info);
  SLIC_ASSERT_MSG(isValid, "Internal error: Generated mesh is invalid.");
  AXOM_ANNOTATE_END("validate_post");
  #endif

  return;
}

template <typename ExecSpace>
void convert_blueprint_structured_explicit_to_unstructured_impl_2d(axom::sidre::Group* meshGrp,
                                                                   const std::string& topoName)
{
  AXOM_ANNOTATE_SCOPE("convert_to_unstructured");
  static constexpr int DIM = 2;
  static constexpr int NUM_VERTS_PER_QUAD = 4;

  const std::string& coordsetName =
    meshGrp->getView("topologies/" + topoName + "/coordset")->getString();

  sidre::Group* coordsetGrp = meshGrp->getGroup("coordsets")->getGroup(coordsetName);
  SLIC_ASSERT(std::string(coordsetGrp->getView("type")->getString()) == "explicit");

  axom::sidre::Group* topoGrp = meshGrp->getGroup("topologies")->getGroup(topoName);
  axom::sidre::View* topoTypeView = topoGrp->getView("type");
  SLIC_ASSERT(std::string(topoTypeView->getString()) == "structured");
  topoTypeView->setString("unstructured");
  topoGrp->createView("elements/shape")->setString("quad");

  axom::sidre::Group* topoElemGrp = topoGrp->getGroup("elements");
  axom::sidre::Group* topoDimsGrp = topoElemGrp->getGroup("dims");

  // Assuming no ghost, but we should eventually support ghosts.
  SLIC_ASSERT(!topoGrp->hasGroup("elements/offsets"));
  SLIC_ASSERT(!topoGrp->hasGroup("elements/strides"));

  const axom::StackArray<axom::IndexType, DIM> cShape {
    axom::IndexType(topoDimsGrp->getView("i")->getNode().to_value()),
    axom::IndexType(topoDimsGrp->getView("j")->getNode().to_value())};
  const axom::StackArray<axom::IndexType, DIM> vShape {cShape[0] + 1, cShape[1] + 1};

  const axom::IndexType cCount = cShape[0] * cShape[1];

  // 4 vertices per quad cell
  const axom::StackArray<axom::IndexType, 2> connShape {cCount, NUM_VERTS_PER_QUAD};
  axom::sidre::View* connView =
    topoGrp->createViewWithShapeAndAllocate("elements/connectivity",
                                            axom::sidre::detail::SidreTT<axom::IndexType>::id,
                                            2,
                                            connShape.begin());
  axom::ArrayView<axom::IndexType, 2> connArrayView(
    static_cast<axom::IndexType*>(connView->getVoidPtr()),
    connShape);

  axom::MDMapping<DIM> cIdMapping(cShape, axom::ArrayStrideOrder::COLUMN);
  axom::MDMapping<DIM> vIdMapping(vShape, axom::ArrayStrideOrder::COLUMN);

  const axom::StackArray<const axom::StackArray<axom::IndexType, DIM>, NUM_VERTS_PER_QUAD> cornerOffsets {
    {{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
  axom::for_all<ExecSpace>(
    cCount,
    AXOM_LAMBDA(axom::IndexType iCell) {
      axom::StackArray<axom::IndexType, DIM> cIdx = cIdMapping.toMultiIndex(iCell);
      for(int n = 0; n < NUM_VERTS_PER_QUAD; ++n)
      {
        const auto& cornerOffset = cornerOffsets[n];
        axom::StackArray<axom::IndexType, DIM> vIdx;
        for(int d = 0; d < DIM; ++d)
        {
          vIdx[d] = cIdx[d] + cornerOffset[d];
        }
        axom::IndexType iVert = vIdMapping.toFlatIndex(vIdx);
        connArrayView(iCell, n) = iVert;
      }
    });

  const bool addExtraDataForMint = true;
  if(addExtraDataForMint)
  {
    AXOM_ANNOTATE_SCOPE("add_extra");
    /*
      Constructing a mint mesh from meshGrp fails unless we add some
      extra data.  Blueprint doesn't require this extra data.  (The mesh
      passes conduit's Blueprint verification.)  This should be fixed,
      or we should write better blueprint support utilities.
    */
    /*
      Make the coordinate arrays 2D to use mint::Mesh.
      For some reason, mint::Mesh requires the arrays to be
      2D, even though the second dimension is always 1.
    */
    axom::IndexType curShape[2] = {};
    auto* valuesGrp = coordsetGrp->getGroup("values");
    SLIC_ASSERT(valuesGrp != nullptr);

    int curDim = valuesGrp->getView("x")->getShape(2, curShape);
    SLIC_ASSERT(curDim == 1);
    curShape[1] = 1;

    valuesGrp->getView("x")->reshapeArray(2, curShape);
    valuesGrp->getView("y")->reshapeArray(2, curShape);

    // Make connectivity array 2D for the same reason.
    auto* elementsGrp = topoGrp->getGroup("elements");
    auto* connView = elementsGrp->getView("connectivity");
    curDim = connView->getShape(2, curShape);
    SLIC_ASSERT(curDim == 2);
    AXOM_UNUSED_VAR(curDim);

    // mint::Mesh requires connectivity strides, even though Blueprint doesn't.
    constexpr axom::IndexType BIT_SPECIFIC_NUM_VERTS_PER_QUAD = 4;
    elementsGrp->createViewScalar("stride", BIT_SPECIFIC_NUM_VERTS_PER_QUAD);

    // mint::Mesh requires field group, even though Blueprint doesn't.
    meshGrp->createGroup("fields");
  }

  return;
}

  #if defined(AXOM_USE_CONDUIT)
bool verifyBlueprintMesh(const axom::sidre::Group* meshGrp, conduit::Node info)
{
  conduit::Node meshNode;
  meshGrp->createNativeLayout(meshNode);
  bool isValid = conduit::blueprint::mesh::verify(meshNode, info);
  if(!isValid)
  {
    info.print();
  }
  return isValid;
}
  #endif

#endif  // AXOM_USE_SIDRE

void fill_cartesian_coords_3d(axom::runtime_policy::Policy runtimePolicy,
                              const primal::BoundingBox<double, 3>& domainBox,
                              axom::ArrayView<double, 3>& xView,
                              axom::ArrayView<double, 3>& yView,
                              axom::ArrayView<double, 3>& zView)
{
  if(runtimePolicy == axom::runtime_policy::Policy::seq)
  {
    fill_cartesian_coords_3d_impl<axom::SEQ_EXEC>(domainBox, xView, yView, zView);
  }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  if(runtimePolicy == axom::runtime_policy::Policy::omp)
  {
    fill_cartesian_coords_3d_impl<axom::OMP_EXEC>(domainBox, xView, yView, zView);
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(runtimePolicy == axom::runtime_policy::Policy::cuda)
  {
    fill_cartesian_coords_3d_impl<axom::CUDA_EXEC<256>>(domainBox, xView, yView, zView);
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(runtimePolicy == axom::runtime_policy::Policy::hip)
  {
    fill_cartesian_coords_3d_impl<axom::HIP_EXEC<256>>(domainBox, xView, yView, zView);
  }
#endif
}

void fill_cartesian_coords_2d(axom::runtime_policy::Policy runtimePolicy,
                              const primal::BoundingBox<double, 2>& domainBox,
                              axom::ArrayView<double, 2>& xView,
                              axom::ArrayView<double, 2>& yView)
{
  if(runtimePolicy == axom::runtime_policy::Policy::seq)
  {
    fill_cartesian_coords_2d_impl<axom::SEQ_EXEC>(domainBox, xView, yView);
  }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  if(runtimePolicy == axom::runtime_policy::Policy::omp)
  {
    fill_cartesian_coords_2d_impl<axom::OMP_EXEC>(domainBox, xView, yView);
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(runtimePolicy == axom::runtime_policy::Policy::cuda)
  {
    fill_cartesian_coords_2d_impl<axom::CUDA_EXEC<256>>(domainBox, xView, yView);
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(runtimePolicy == axom::runtime_policy::Policy::hip)
  {
    fill_cartesian_coords_2d_impl<axom::HIP_EXEC<256>>(domainBox, xView, yView);
  }
#endif
}

template <typename ExecSpace>
void fill_cartesian_coords_3d_impl(const primal::BoundingBox<double, 3>& domainBox,
                                   axom::ArrayView<double, 3>& xView,
                                   axom::ArrayView<double, 3>& yView,
                                   axom::ArrayView<double, 3>& zView)
{
#if defined(AXOM_DEBUG)
  using XS = axom::execution_space<ExecSpace>;
  SLIC_ASSERT_MSG(
    XS::usesAllocId(xView.getAllocatorID()) && XS::usesAllocId(yView.getAllocatorID()) &&
      XS::usesAllocId(zView.getAllocatorID()),
    std::string("fill_cartesian_coords_3d_impl: alloc ids ") +
      std::to_string(xView.getAllocatorID()) + ", " + std::to_string(yView.getAllocatorID()) +
      " and " + std::to_string(zView.getAllocatorID()) +
      " are not all compatible with execution space " + XS::name());
#endif

  const auto& shape = xView.shape();
  SLIC_ASSERT(shape == yView.shape());
  SLIC_ASSERT(shape == zView.shape());

  // Mesh resolution
  const axom::NumericArray<axom::IndexType, 3> res {shape[0] - 1, shape[1] - 1, shape[2] - 1};

  // Mesh spacings.
  double dx = (domainBox.getMax()[0] - domainBox.getMin()[0]) / res[0];
  double dy = (domainBox.getMax()[1] - domainBox.getMin()[1]) / res[1];
  double dz = (domainBox.getMax()[2] - domainBox.getMin()[2]) / res[2];

  const auto& mapping = xView.mapping();
  SLIC_ASSERT(mapping == yView.mapping());
  SLIC_ASSERT(mapping == zView.mapping());
  auto order = mapping.getStrideOrder();
  if(int(order) & int(axom::ArrayStrideOrder::COLUMN))
  {
    axom::for_all<ExecSpace>(
      shape,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        xView(i, j, k) = domainBox.getMin()[0] + i * dx;
        yView(i, j, k) = domainBox.getMin()[1] + j * dy;
        zView(i, j, k) = domainBox.getMin()[2] + k * dz;
      });
  }
  else
  {
    axom::StackArray<axom::IndexType, 3> shapeKJI;
    shapeKJI[0] = shape[2];
    shapeKJI[1] = shape[1];
    shapeKJI[2] = shape[0];

    axom::for_all<ExecSpace>(
      shapeKJI,
      AXOM_LAMBDA(axom::IndexType k, axom::IndexType j, axom::IndexType i) {
        xView(i, j, k) = domainBox.getMin()[0] + i * dx;
        yView(i, j, k) = domainBox.getMin()[1] + j * dy;
        zView(i, j, k) = domainBox.getMin()[2] + k * dz;
      });
  }
}

template <typename ExecSpace>
void fill_cartesian_coords_2d_impl(const primal::BoundingBox<double, 2>& domainBox,
                                   axom::ArrayView<double, 2>& xView,
                                   axom::ArrayView<double, 2>& yView)
{
#if defined(AXOM_DEBUG)
  using XS = axom::execution_space<ExecSpace>;
  SLIC_ASSERT_MSG(
    XS::usesAllocId(xView.getAllocatorID()) && XS::usesAllocId(yView.getAllocatorID()),
    std::string("fill_cartesian_coords_2d_impl: alloc ids ") +
      std::to_string(xView.getAllocatorID()) + " and " + std::to_string(yView.getAllocatorID()) +
      " are not all compatible with execution space " + XS::name());
#endif

  const auto& shape = xView.shape();
  SLIC_ASSERT(shape == yView.shape());

  // Mesh resolution
  const axom::NumericArray<axom::IndexType, 2> res {shape[0] - 1, shape[1] - 1};

  // Mesh spacings.
  double dx = (domainBox.getMax()[0] - domainBox.getMin()[0]) / res[0];
  double dy = (domainBox.getMax()[1] - domainBox.getMin()[1]) / res[1];

  const auto& mapping = xView.mapping();
  SLIC_ASSERT(mapping == yView.mapping());
  auto order = mapping.getStrideOrder();
  if(int(order) & int(axom::ArrayStrideOrder::COLUMN))
  {
    axom::for_all<ExecSpace>(
      shape,
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) {
        xView(i, j) = domainBox.getMin()[0] + i * dx;
        yView(i, j) = domainBox.getMin()[1] + j * dy;
      });
  }
  else
  {
    axom::StackArray<axom::IndexType, 2> shapeJI;
    shapeJI[0] = shape[1];
    shapeJI[1] = shape[0];
    axom::for_all<ExecSpace>(
      shapeJI,
      AXOM_LAMBDA(axom::IndexType j, axom::IndexType i) {
        xView(i, j) = domainBox.getMin()[0] + i * dx;
        yView(i, j) = domainBox.getMin()[1] + j * dy;
      });
  }
}

}  // namespace util
}  // namespace quest
}  // namespace axom
