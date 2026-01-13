// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/detail/clipping/Plane3DClipper.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

Plane3DClipper::Plane3DClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Plane3D") : name)
{
  extractClipperInfo();
}

bool Plane3DClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                     axom::Array<LabelType>& labels)
{
  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();
  if(labels.size() < cellCount || labels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelCellsInOutImpl<axom::SEQ_EXEC>(shapeMesh, labels.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelCellsInOutImpl<axom::OMP_EXEC>(shapeMesh, labels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelCellsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, labels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelCellsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, labels.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

bool Plane3DClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::ArrayView<const axom::IndexType> cellIds,
                                    axom::Array<LabelType>& tetLabels)
{
  int allocId = shapeMesh.getAllocatorID();
  const auto cellCount = cellIds.size();
  const auto tetCount = cellCount * NUM_TETS_PER_HEX;
  if(tetLabels.size() < tetCount || tetLabels.getAllocatorID() != allocId)
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), tetCount, tetCount, allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelTetsInOutImpl<axom::SEQ_EXEC>(shapeMesh, cellIds, tetLabels.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelTetsInOutImpl<axom::OMP_EXEC>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelTetsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelTetsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

bool Plane3DClipper::specializedClipCells(quest::experimental::ShapeMesh& shapeMesh,
                                          axom::ArrayView<double> ovlap,
                                          conduit::Node& statistics)
{
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    specializedClipCellsImpl<axom::SEQ_EXEC>(shapeMesh, ovlap, statistics);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    specializedClipCellsImpl<axom::OMP_EXEC>(shapeMesh, ovlap, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    specializedClipCellsImpl<axom::CUDA_EXEC<256>>(shapeMesh, ovlap, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    specializedClipCellsImpl<axom::HIP_EXEC<256>>(shapeMesh, ovlap, statistics);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

bool Plane3DClipper::specializedClipCells(quest::experimental::ShapeMesh& shapeMesh,
                                          axom::ArrayView<double> ovlap,
                                          const axom::ArrayView<IndexType>& cellIds,
                                          conduit::Node& statistics)
{
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    specializedClipCellsImpl<axom::SEQ_EXEC>(shapeMesh, ovlap, cellIds, statistics);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    specializedClipCellsImpl<axom::OMP_EXEC>(shapeMesh, ovlap, cellIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    specializedClipCellsImpl<axom::CUDA_EXEC<256>>(shapeMesh, ovlap, cellIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    specializedClipCellsImpl<axom::HIP_EXEC<256>>(shapeMesh, ovlap, cellIds, statistics);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

bool Plane3DClipper::specializedClipTets(quest::experimental::ShapeMesh& shapeMesh,
                                         axom::ArrayView<double> ovlap,
                                         const axom::ArrayView<IndexType>& tetIds,
                                         conduit::Node& statistics)
{
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    specializedClipTetsImpl<axom::SEQ_EXEC>(shapeMesh, ovlap, tetIds, statistics);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    specializedClipTetsImpl<axom::OMP_EXEC>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    specializedClipTetsImpl<axom::CUDA_EXEC<256>>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    specializedClipTetsImpl<axom::HIP_EXEC<256>>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void Plane3DClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                         axom::ArrayView<LabelType> labels)
{
  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();
  auto vertCount = shapeMesh.getVertexCount();
  auto cellVolumes = shapeMesh.getCellVolumes();
  constexpr double EPS = 1e-10;

  const auto& vertCoords = shapeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];

  /*
    Compute whether vertices are inside shape.
  */
  axom::Array<bool> vertIsInside {ArrayOptions::Uninitialized(), vertCount, vertCount, allocId};
  auto vertIsInsideView = vertIsInside.view();
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vertIsInsideView.getAllocatorID()));

  auto plane = m_plane;
  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};
      double signedDist = plane.signedDistance(vert);
      vertIsInsideView[vertId] = signedDist > 0;
    });

  /*
   * Label cell by whether it has vertices inside, outside or both.
   */
  axom::ArrayView<const axom::IndexType, 2> connView = shapeMesh.getCellNodeConnectivity();
  SLIC_ASSERT(connView.shape()[1] == NUM_VERTS_PER_CELL_3D);

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      if(axom::utilities::isNearlyEqual(cellVolumes[cellId], 0.0, EPS))
      {
        labels[cellId] = LabelType::LABEL_OUT;
        return;
      }
      auto cellVertIds = connView[cellId];
      bool hasIn = vertIsInsideView[cellVertIds[0]];
      bool hasOut = !hasIn;
      for(int vi = 0; vi < NUM_VERTS_PER_CELL_3D; ++vi)
      {
        int vertId = cellVertIds[vi];
        bool isIn = vertIsInsideView[vertId];
        hasIn |= isIn;
        hasOut |= !isIn;
      }
      labels[cellId] = !hasOut ? LabelType::LABEL_IN
        : !hasIn               ? LabelType::LABEL_OUT
                               : LabelType::LABEL_ON;
    });

  return;
}

template <typename ExecSpace>
void Plane3DClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                        axom::ArrayView<const axom::IndexType> cellIds,
                                        axom::ArrayView<LabelType> tetLabels)
{
  auto cellCount = cellIds.size();
  auto meshTets = shapeMesh.getCellsAsTets();
  auto meshTetVolumes = shapeMesh.getTetVolumes();

  auto plane = m_plane;

  constexpr double EPS = 1e-10;

  /*
   * Label tet by whether it has vertices inside, outside or both.
   * Degenerate tets as outside, because they contribute no volume.
   */
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
      axom::IndexType cellId = cellIds[ci];

      const TetrahedronType* tetsForCell = &meshTets[cellId * NUM_TETS_PER_HEX];
      const double* tetVolumesForCell = &meshTetVolumes[cellId * NUM_TETS_PER_HEX];

      for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        const auto& tet = tetsForCell[ti];
        LabelType& tetLabel = tetLabels[ci * NUM_TETS_PER_HEX + ti];

        if(axom::utilities::isNearlyEqual(tetVolumesForCell[ti], 0.0, EPS))
        {
          tetLabel = LabelType::LABEL_OUT;
          continue;
        }

        bool hasIn = false;
        bool hasOut = false;
        for(int vi = 0; vi < TetrahedronType::NUM_VERTS; ++vi)
        {
          const auto& vert = tet[vi];
          double signedDist = plane.signedDistance(vert);
          hasIn |= signedDist > 0;
          hasOut |= signedDist < 0;
        }
        tetLabel = !hasOut ? LabelType::LABEL_IN
          : !hasIn         ? LabelType::LABEL_OUT
                           : LabelType::LABEL_ON;
      }
    });

  return;
}

template <typename ExecSpace>
void Plane3DClipper::specializedClipCellsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                              axom::ArrayView<double> ovlap,
                                              conduit::Node& statistics)
{
  axom::IndexType cellCount = shapeMesh.getCellCount();
  axom::Array<IndexType> cellIds(cellCount, cellCount, shapeMesh.getAllocatorID());
  auto cellIdsView = cellIds.view();
  axom::for_all<ExecSpace>(cellCount, AXOM_LAMBDA(axom::IndexType i) { cellIdsView[i] = i; });
  specializedClipCellsImpl<ExecSpace>(shapeMesh, ovlap, cellIds, statistics);
}

template <typename ExecSpace>
void Plane3DClipper::specializedClipCellsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                              axom::ArrayView<double> ovlap,
                                              const axom::ArrayView<IndexType>& cellIds,
                                              conduit::Node& statistics)
{
  constexpr double EPS = 1e-10;

  auto cellsAsTets = shapeMesh.getCellsAsTets();

  auto plane = m_plane;

  axom::ReduceSum<ExecSpace, std::int64_t> missSum {0};

  axom::for_all<ExecSpace>(
    cellIds.size(),
    AXOM_LAMBDA(axom::IndexType i) {
      axom::IndexType cellId = cellIds[i];
      const TetrahedronType* tetsInHex = cellsAsTets.data() + cellId * NUM_TETS_PER_HEX;
      double vol = 0.0;
      for(int ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        const auto& tet = tetsInHex[ti];
        primal::Polyhedron<double, 3> overlap = primal::clip(tet, plane, EPS);
        if(overlap.numVertices() >= 4)
        {
          auto volume = overlap.volume();
          vol += volume;
        }
        else
        {
          missSum += 1;
        }
      }
      ovlap[cellId] = vol;
    });

  statistics["clipsOn"].set_int64(cellIds.size() * NUM_TETS_PER_HEX);
  statistics["clipsSum"].set_int64(cellIds.size() * NUM_TETS_PER_HEX);
  statistics["missSum"].set_int64(missSum.get());
}

template <typename ExecSpace>
void Plane3DClipper::specializedClipTetsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                             axom::ArrayView<double> ovlap,
                                             const axom::ArrayView<IndexType>& tetIds,
                                             conduit::Node& statistics)
{
  constexpr double EPS = 1e-10;
  using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

  auto meshTets = shapeMesh.getCellsAsTets();
  IndexType tetCount = tetIds.size();
  auto plane = m_plane;

  axom::for_all<ExecSpace>(
    tetCount,
    AXOM_LAMBDA(axom::IndexType ti) {
      axom::IndexType tetId = tetIds[ti];
      axom::IndexType cellId = tetId / NUM_TETS_PER_HEX;
      const auto& tet = meshTets[tetId];
      primal::Polyhedron<double, 3> overlap = primal::clip(tet, plane, EPS);
      double vol = overlap.volume();
      RAJA::atomicAdd<ATOMIC_POL>(ovlap.data() + cellId, vol);
    });

  // Because the tet screening is perfect, all tets in tetIds are on the plane.
  statistics["onSum"].set_int64(tetCount);
  statistics["clipsSum"].set_int64(tetCount);
}

void Plane3DClipper::extractClipperInfo()
{
  const auto normal = m_info.fetch_existing("normal").as_double_array();
  const double offset = m_info.fetch_existing("offset").as_double();
  Vector3DType nVec;
  for(int d = 0; d < 3; ++d)
  {
    nVec[d] = normal[d];
  }
  m_plane = Plane3DType(nVec, offset);
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
