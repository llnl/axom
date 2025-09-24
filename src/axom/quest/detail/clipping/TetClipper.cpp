// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/detail/clipping/TetClipper.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

TetClipper::TetClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Tet") : name)
  , m_transformer(m_extTrans)
{
  extractClipperInfo();

  for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
  {
    m_tet[i] = m_transformer.getTransformed(m_tetBeforeTrans[i]);
  }

  for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
  {
    m_bb.addPoint(m_tet[i]);
  }

  for(int iPlane = 0; iPlane < 4; ++iPlane)
  {
    const Point3DType& a = m_tet[iPlane % 4];
    const Point3DType& b = m_tet[(iPlane + 1) % 4];
    const Point3DType& c = m_tet[(iPlane + 2) % 4];
    m_planes[iPlane] = axom::primal::make_plane(a, b, c);

    // For tet points ordered by right hand rule, odd planes
    // face outside.  Make them face inside.
    if(iPlane % 2 == 1) m_planes[iPlane].flip();

    const Point3DType& apex = m_tet[(iPlane + 3) % 4];
    m_heights[iPlane] = m_planes[iPlane].signedDistance(apex);
    SLIC_ASSERT(m_heights[iPlane] >= 0);
  }
}

bool TetClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh, axom::Array<LabelType>& labels)
{
  AXOM_ANNOTATE_SCOPE("TetClipper::labelCellsInOut");
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelInOutImpl<axom::SEQ_EXEC>(shapeMesh, labels);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelInOutImpl<axom::OMP_EXEC>(shapeMesh, labels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, labels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, labels);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void TetClipper::labelInOutImpl(quest::experimental::ShapeMesh& shapeMesh, axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "TetClipper requires a 3D mesh.");
  /*
   * Compute whether the mesh vertices are above the tet or below as
   * the tet rests on its 4 facets.
   *
   * - For any facet the tet rests on, if all cell vertices are
   *   above the tet or all are below, cell is OUT.
   *
   * - If all cell vertices are on the tet side of all facets,
   *   the cell is IN.
   *
   * - Otherwise, cell is ON.
  */

  int allocId = shapeMesh.getAllocatorID();
  auto vertCount = shapeMesh.getVertexCount();
  auto cellCount = shapeMesh.getCellCount();

  const auto& vertCoords = shapeMesh.getVertexCoords3D();
  const auto& vX = vertCoords[0];
  const auto& vY = vertCoords[1];
  const auto& vZ = vertCoords[2];
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vX.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vY.getAllocatorID()));
  SLIC_ASSERT(axom::execution_space<ExecSpace>::usesAllocId(vZ.getAllocatorID()));

  axom::ArrayView<const axom::IndexType, 2> connView = shapeMesh.getCellNodeConnectivity();
  SLIC_ASSERT(connView.shape() ==
              (axom::StackArray<axom::IndexType, 2> {cellCount, HexahedronType::NUM_HEX_VERTS}));

  axom::Array<bool> below[4];
  axom::Array<bool> above[4];
  axom::ArrayView<bool> belowView[4];
  axom::ArrayView<bool> aboveView[4];
  for(IndexType p = 0; p < 4; ++p)
  {
    below[p] = axom::Array<bool>(ArrayOptions::Uninitialized(), vertCount, vertCount, allocId);
    above[p] = axom::Array<bool>(ArrayOptions::Uninitialized(), vertCount, vertCount, allocId);
    belowView[p] = below[p].view();
    aboveView[p] = above[p].view();
  }

  const auto tetHeights = m_heights;
  auto planes = m_planes;

  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      primal::Point3D vert {vX[vertId], vY[vertId], vZ[vertId]};

      for(IndexType p = 0; p < 4; ++p)
      {
        double vertHeight = planes[p].signedDistance(vert);
        belowView[p][vertId] = vertHeight < 0;
        aboveView[p][vertId] = vertHeight > tetHeights[p];
      }
    });

  if(labels.size() < cellCount || labels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(),
                                    cellCount,
                                    cellCount,
                                    shapeMesh.getAllocatorID());
  }
  auto labelsView = labels.view();

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      LabelType& cellLabel = labelsView[cellId];
      auto cellVertIds = connView[cellId];

      cellLabel = LABEL_ON;
      bool vertsAreOnTetSideOfAllPlanes = true;
      for(IndexType p = 0; p < 4; ++p)
      {
        bool allVertsBelow = true;
        bool allVertsAbove = true;
        for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
        {
          int vertId = cellVertIds[vi];
          auto vertIsBelow = belowView[p][vertId];
          auto vertIsAbove = aboveView[p][vertId];
          allVertsBelow &= vertIsBelow;
          allVertsAbove &= vertIsAbove;
          vertsAreOnTetSideOfAllPlanes &= !vertIsBelow;
        }
        if(allVertsBelow || allVertsAbove)
        {
          cellLabel = LABEL_OUT;
          break;
        }
      }
      if(cellLabel != LABEL_OUT && vertsAreOnTetSideOfAllPlanes)
      {
        cellLabel = LABEL_IN;
      }
    });

  return;
}

bool TetClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                axom::ArrayView<const axom::IndexType> cellsOnBdry,
                                axom::Array<LabelType>& tetLabels)
{
  AXOM_ANNOTATE_SCOPE("TetClipper::labelTetsInOut");
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelTetsInOutImpl<axom::SEQ_EXEC>(shapeMesh, cellsOnBdry, tetLabels);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelTetsInOutImpl<axom::OMP_EXEC>(shapeMesh, cellsOnBdry, tetLabels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelTetsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, cellsOnBdry, tetLabels);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelTetsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, cellsOnBdry, tetLabels);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void TetClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::ArrayView<const axom::IndexType> cellsOnBdry,
                                    axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "TetClipper requires a 3D mesh.");
  /*
   * Compute whether the tets in hexes listed in cellsOnBdry
   * as in, out or on the boundary.
   * the tet rests on its 4 facets.
   *
   * - For any facet the tet rests on, if all cell vertices are
   *   above the tet or all are below, cell is OUT.
   *
   * - If all cell vertices are on the tet side of all facets,
   *   the cell is IN.
   *
   * - Otherwise, cell is ON.
  */

  int allocId = shapeMesh.getAllocatorID();

  auto meshTets = shapeMesh.getCellsAsTets();

  const axom::IndexType cellCount = cellsOnBdry.size();

  if(tetLabels.size() < cellCount * TETS_PER_HEXAHEDRON ||
     tetLabels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(),
                                       cellCount * TETS_PER_HEXAHEDRON,
                                       cellCount * TETS_PER_HEXAHEDRON,
                                       allocId);
  }

  auto tetLabelsView = tetLabels.view();

  const auto geomHeights = m_heights;
  auto planes = m_planes;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
      axom::IndexType cellId = cellsOnBdry[ci];

      const TetrahedronType* tetsForCell = &meshTets[cellId * TETS_PER_HEXAHEDRON];

      for(IndexType ti = 0; ti < TETS_PER_HEXAHEDRON; ++ti)
      {
        const TetrahedronType& tet = tetsForCell[ti];
        LabelType& tetLabel = tetLabelsView[ci * TETS_PER_HEXAHEDRON + ti];

        tetLabel = LABEL_ON;
        bool vertsAreOnTetSideOfAllPlanes = true;

        bool allVertsBelow = true;
        bool allVertsAbove = true;

        for(IndexType p = 0; p < 4; ++p)
        {
          const auto& plane = planes[p];

          for(IndexType vi = 0; vi < 4; ++vi)
          {
            const auto& vert = tet[vi];
            double vertHeight = plane.signedDistance(vert);
            bool vertIsBelow = vertHeight < 0;
            bool vertIsAbove = vertHeight > geomHeights[p];

            allVertsBelow &= vertIsBelow;
            allVertsAbove &= vertIsAbove;
            vertsAreOnTetSideOfAllPlanes &= !vertIsBelow;
          }

          if(allVertsBelow || allVertsAbove)
          {
            tetLabel = LABEL_OUT;
            break;
          }
        }

        if(tetLabel != LABEL_OUT && vertsAreOnTetSideOfAllPlanes)
        {
          tetLabel = LABEL_IN;
        }
      }
    });

  return;
}

bool TetClipper::getGeometryAsTets(quest::experimental::ShapeMesh& shapeMesh, axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_SCOPE("TetClipper::getGeometryAsTets");
  int allocId = shapeMesh.getAllocatorID();
  if(tets.getAllocatorID() != allocId || tets.size() != 1)
  {
    tets = axom::Array<TetrahedronType>(1, 1, allocId);
  }
  axom::copy(tets.data(), &m_tet, sizeof(TetrahedronType));
  return true;
}

void TetClipper::extractClipperInfo()
{
  const auto v0 = m_info.fetch_existing("v0").as_double_array();
  const auto v1 = m_info.fetch_existing("v1").as_double_array();
  const auto v2 = m_info.fetch_existing("v2").as_double_array();
  const auto v3 = m_info.fetch_existing("v3").as_double_array();
  for(int d = 0; d < 3; ++d)
  {
    m_tetBeforeTrans[0][d] = v0[d];
    m_tetBeforeTrans[1][d] = v1[d];
    m_tetBeforeTrans[2][d] = v2[d];
    m_tetBeforeTrans[3][d] = v3[d];
  }
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
