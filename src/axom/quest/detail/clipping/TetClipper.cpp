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
  , m_extTransformer(m_extTrans)
{
  extractClipperInfo();

  for(int i = 0; i < TetrahedronType::NUM_VERTS; ++i)
  {
    m_tet[i] = m_extTransformer.getTransformed(m_tetBeforeTrans[i]);
  }

  /*
   * Compute the transformation from m_tet to a unit tet.  Location of
   * points w.r.t. the tet are done in the unit tet space.  The unit
   * tet has heights 1, 1, 1 and 1/sqrt(3) as it rests on each of its 4
   * faces.  Height w.r.t. the 4th face are scaled by sqrt(3) to make
   * it come out to 1.  So height < 0 means means below the tet and
   * height > 1 means above the tet.
   *
   * Points with any height outside [0,1] cannot possibly be in the tet.
   * Points with all 4 heights in [0,1] must be in the tet.
   */
  Point3DType unitTetPts[] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  m_toUnitTet.setByTerminusPts(&m_tet[0], unitTetPts);
}

bool TetClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                 axom::Array<LabelType>& cellLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "FSorClipper requires a 3D mesh.");

  const int allocId = shapeMesh.getAllocatorID();
  const auto cellCount = shapeMesh.getCellCount();
  if(cellLabels.size() < cellCount || cellLabels.getAllocatorID() != allocId)
  {
    cellLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelCellsInOutImpl<axom::SEQ_EXEC>(shapeMesh, cellLabels.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelCellsInOutImpl<axom::OMP_EXEC>(shapeMesh, cellLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelCellsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, cellLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelCellsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, cellLabels.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void TetClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                     axom::ArrayView<LabelType> cellLabels)
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
  auto meshCellVolumes = shapeMesh.getCellVolumes();

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

  /*
   * Compute whether mesh vertices are above/below the tet.
   */
  axom::Array<bool> below[4];
  axom::Array<bool> above[4];
  axom::ArrayView<bool> belowView[4];
  axom::ArrayView<bool> aboveView[4];
  for(IndexType p = 0; p < 4; ++p)
  {
    below[p] = axom::Array<bool>(ArrayOptions::Uninitialized(), vertCount, 0, allocId);
    above[p] = axom::Array<bool>(ArrayOptions::Uninitialized(), vertCount, 0, allocId);
    belowView[p] = below[p].view();
    aboveView[p] = above[p].view();
  }

  auto toUnitTet = m_toUnitTet;

  axom::for_all<ExecSpace>(
    vertCount,
    AXOM_LAMBDA(axom::IndexType vertId) {
      // vh is the heights of the vertex in the space of the unit tet.
      // See comment on m_toUnitTet.
      axom::NumericArray<double, 4> vh({vX[vertId], vY[vertId], vZ[vertId], 0});
      toUnitTet.transform(vh[0], vh[1], vh[2]);
      vh[3] = 1 - vh[0] - vh[1] - vh[2];

      for(int p = 0; p < 4; ++p)
      {
        belowView[p][vertId] = vh[p] < 0;
        aboveView[p][vertId] = vh[p] > 1;
      }
    });

  constexpr double EPS = 1e-10;

  /*
   * Compute whether mesh cells are above/below the tet.
   */
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      if(axom::utilities::isNearlyEqual(meshCellVolumes[cellId], 0.0, EPS))
      {
        cellLabels[cellId] = LabelType::LABEL_OUT;
        return;
      }

      LabelType& cellLabel = cellLabels[cellId];
      auto cellVertIds = connView[cellId];

      cellLabel = LabelType::LABEL_ON;
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
          cellLabel = LabelType::LABEL_OUT;
          break;
        }
      }
      if(cellLabel != LabelType::LABEL_OUT && vertsAreOnTetSideOfAllPlanes)
      {
        cellLabel = LabelType::LABEL_IN;
      }
    });

  return;
}

bool TetClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                axom::ArrayView<const axom::IndexType> cellIds,
                                axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "TetClipper requires a 3D mesh.");

  const int allocId = shapeMesh.getAllocatorID();
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

template <typename ExecSpace>
void TetClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::ArrayView<const axom::IndexType> cellIds,
                                    axom::ArrayView<LabelType> tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "TetClipper requires a 3D mesh.");
  /*
   * Compute whether the tets in hexes listed in cellIds
   * are in, out or on the boundary.
   *
   * Picture the tet resting on each of of its 4 faces.
   * - For any facet the tet rests on, if all cell vertices are
   *   above the tet or all are below, cell is OUT.
   *
   * - If all cell vertices are on the tet side of all facets,
   *   the cell is IN.
   *
   * - Otherwise, cell is ON.
   */

  auto meshTets = shapeMesh.getCellsAsTets();
  auto tetVolumes = shapeMesh.getTetVolumes();

  const axom::IndexType cellCount = cellIds.size();

  auto toUnitTet = m_toUnitTet;
  constexpr double EPS = 1e-10;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
      axom::IndexType cellId = cellIds[ci];

      const TetrahedronType* tetsForCell = &meshTets[cellId * NUM_TETS_PER_HEX];

      for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        const TetrahedronType& cellTet = tetsForCell[ti];
        LabelType& tetLabel = tetLabels[ci * NUM_TETS_PER_HEX + ti];
        const axom::IndexType tetId = cellId * NUM_TETS_PER_HEX + ti;

        if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
        {
          tetLabel = LabelType::LABEL_OUT;
          continue;
        }

        tetLabel = LabelType::LABEL_ON;

        bool allVertsBelow = true;
        bool allVertsAbove = true;
        bool vertsAreOnTetSideOfAllPlanes = true;

        for(IndexType vi = 0; vi < 4; ++vi)
        {
          const auto& vert = cellTet[vi];

          // vh is the heights of vert in the space of the unit tet.
          // See comment on m_toUnitTet.
          axom::NumericArray<double, 4> vh({vert[0], vert[1], vert[2], 0});
          toUnitTet.transform(vh[0], vh[1], vh[2]);
          vh[3] = 1 - vh[0] - vh[1] - vh[2];

          // Where vertex vi is w.r.t. the tet resting on side pj.
          for(int pj = 0; pj < 4; ++pj)
          {
            bool vertIsBelow = vh[pj] < 0;
            bool vertIsAbove = vh[pj] > 1;

            allVertsBelow &= vertIsBelow;
            allVertsAbove &= vertIsAbove;
            vertsAreOnTetSideOfAllPlanes &= !vertIsBelow;
          }

          if(allVertsBelow || allVertsAbove)
          {
            tetLabel = LabelType::LABEL_OUT;
            break;
          }
        }

        if(tetLabel != LabelType::LABEL_OUT && vertsAreOnTetSideOfAllPlanes)
        {
          tetLabel = LabelType::LABEL_IN;
        }
      }
    });

  return;
}

bool TetClipper::getGeometryAsTets(quest::experimental::ShapeMesh& shapeMesh,
                                   axom::Array<TetrahedronType>& tets)
{
  AXOM_ANNOTATE_SCOPE("TetClipper::getGeometryAsTets");
  int allocId = shapeMesh.getAllocatorID();
  if(tets.getAllocatorID() != allocId || tets.size() != 1)
  {
    tets = axom::Array<TetrahedronType>(1, 1, allocId);
  }
  // Copy tet into tets array, which may be in non-host memory.
  axom::copy(tets.data(), &m_tet, sizeof(TetrahedronType));
  return true;
}

void TetClipper::extractClipperInfo()
{
  const auto v0 = m_info.fetch_existing("v0").as_double_array();
  const auto v1 = m_info.fetch_existing("v1").as_double_array();
  const auto v2 = m_info.fetch_existing("v2").as_double_array();
  const auto v3 = m_info.fetch_existing("v3").as_double_array();
  SLIC_ASSERT(v0.number_of_elements() == 3);
  SLIC_ASSERT(v1.number_of_elements() == 3);
  SLIC_ASSERT(v2.number_of_elements() == 3);
  SLIC_ASSERT(v3.number_of_elements() == 3);
  for(int d = 0; d < 3; ++d)
  {
    m_tetBeforeTrans[0][d] = v0[d];
    m_tetBeforeTrans[1][d] = v1[d];
    m_tetBeforeTrans[2][d] = v2[d];
    m_tetBeforeTrans[3][d] = v3[d];
  }

  bool fixOrientation = false;
  if(m_info.has_child("fixOrientation"))
  {
    fixOrientation = bool(m_info.fetch_existing("fixOrientation").as_int());
  }

  if(fixOrientation)
  {
    m_tetBeforeTrans.checkAndFixOrientation();
  }
  else
  {
    constexpr double EPS = 1e-10;
    double signedVol = m_tetBeforeTrans.signedVolume();
    if(signedVol < -EPS)
    {
      SLIC_ERROR(
        axom::fmt::format("TetClipper tet {} has negative volume {}.:"
                          "  (See TetClipper's 'fixOrientation' flag.)",
                          m_tetBeforeTrans,
                          signedVol));
    }
  }
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
