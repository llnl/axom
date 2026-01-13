// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/detail/clipping/HexClipper.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

HexClipper::HexClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Hex") : name)
  , m_extTransformer(m_extTrans)
{
  extractClipperInfo();

  for(int i = 0; i < HexahedronType::NUM_HEX_VERTS; ++i)
  {
    m_hex[i] = m_extTransformer.getTransformed(m_hexBeforeTrans[i]);
  }

  axom::StackArray<TetrahedronType, ShapeMesh::NUM_TETS_PER_HEX> geomTets;
  ShapeMesh::hexToTets(m_hex, geomTets.data());
  m_tets.reserve(geomTets.size());
  constexpr double EPS = 1e-10;
  for(const auto& tet : geomTets)
  {
    if(!axom::utilities::isNearlyEqual(tet.volume(), 0.0, EPS))
    {
      m_tets.push_back(tet);
    }
  }

  for(int i = 0; i < HexahedronType::NUM_HEX_VERTS; ++i)
  {
    m_hexBb.addPoint(m_hex[i]);
  }

  computeSurface();
}

bool HexClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                 axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "HexClipper requires a 3D mesh.");

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

bool HexClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                axom::ArrayView<const axom::IndexType> cellIds,
                                axom::Array<LabelType>& tetLabels)
{
  const axom::IndexType cellCount = cellIds.size();
  const int allocId = shapeMesh.getAllocatorID();

  if(tetLabels.size() < cellCount * NUM_TETS_PER_HEX ||
     tetLabels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(),
                                       cellCount * NUM_TETS_PER_HEX,
                                       cellCount * NUM_TETS_PER_HEX,
                                       allocId);
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
void HexClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                     axom::ArrayView<LabelType> labels)
{
  const auto cellCount = shapeMesh.getCellCount();
  const int allocId = shapeMesh.getAllocatorID();
  const auto cellBbs = shapeMesh.getCellBoundingBoxes();
  const auto cellsAsHexes = shapeMesh.getCellsAsHexes();
  const auto cellVolumes = shapeMesh.getCellVolumes();
  const auto hexBb = m_hexBb;
  const auto surfaceTriangles = m_surfaceTriangles;
  axom::Array<TetrahedronType> tets(m_tets, allocId);
  axom::ArrayView<const TetrahedronType> tetsView = tets.view();
  constexpr double EPS = 1e-10;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType cellId) {
      auto& cellLabel = labels[cellId];
      if(axom::utilities::isNearlyEqual(cellVolumes[cellId], 0.0, EPS))
      {
        cellLabel = LabelType::LABEL_OUT;
        return;
      }
      auto& cellBb = cellBbs[cellId];
      const auto& cellHex = cellsAsHexes[cellId];
      cellLabel = polyhedronToLabel(cellHex, cellBb, hexBb, tetsView, surfaceTriangles);
    });

  return;
}

template <typename ExecSpace>
void HexClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::ArrayView<const axom::IndexType> cellIds,
                                    axom::ArrayView<LabelType> tetLabels)
{
  const axom::IndexType cellCount = cellIds.size();
  const int allocId = shapeMesh.getAllocatorID();
  auto meshHexes = shapeMesh.getCellsAsHexes();
  auto tetVolumes = shapeMesh.getTetVolumes();
  const auto hexBb = m_hexBb;
  const auto surfaceTriangles = m_surfaceTriangles;
  axom::Array<TetrahedronType> tets(m_tets, allocId);
  axom::ArrayView<const TetrahedronType> tetsView = tets.view();
  constexpr double EPS = 1e-10;

  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType ci) {
      axom::IndexType cellId = cellIds[ci];
      const HexahedronType& hex = meshHexes[cellId];

      TetrahedronType cellTets[NUM_TETS_PER_HEX];
      ShapeMesh::hexToTets(hex, cellTets);

      for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        const TetrahedronType& cellTet = cellTets[ti];
        LabelType& tetLabel = tetLabels[ci * NUM_TETS_PER_HEX + ti];
        axom::IndexType tetId = cellId * NUM_TETS_PER_HEX + ti;
        if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
        {
          tetLabel = LabelType::LABEL_OUT;
          continue;
        }
        BoundingBox3DType cellTetBb {cellTet[0], cellTet[1], cellTet[2], cellTet[3]};
        tetLabel = polyhedronToLabel(cellTet, cellTetBb, hexBb, tetsView, surfaceTriangles);
      }
    });
  return;
}

template <typename Polyhedron>
AXOM_HOST_DEVICE inline MeshClipperStrategy::LabelType HexClipper::polyhedronToLabel(
  const Polyhedron& verts,
  const BoundingBox3DType& vertsBb,
  const BoundingBox3DType& hexBb,
  const axom::ArrayView<const TetrahedronType>& hexTets,
  const axom::StackArray<Triangle3DType, 24>& surfaceTriangles) const
{
  /*
    If vertsBb and hexBb don't intersect, nothing intersects.
    This check is not technically needed, because the checks
    below can catch it, but it is fast and can avoid the more
    expensive surface triangle intersection checks below.
  */
  if(!hexBb.intersectsWith(vertsBb))
  {
    return LabelType::LABEL_OUT;
  }

  // If vertsBb intersects hex surface, there's a high chance cell does too.
  for(int ti = 0; ti < 24; ++ti)
  {
    const auto& surfTri = surfaceTriangles[ti];
    if(axom::primal::intersect(surfTri, vertsBb))
    {
      return LabelType::LABEL_ON;
    }
  }
  /*
    After eliminating possibility that polyhedron is on the surface,
    it's either completely inside or completely out.
    It's IN if any part of it is IN, so check an arbitrary vertex.
    Note: Should the arbitrary vertex be some weird corner case, we could
    use an alternative, like averaging two opposite corners, 0 and 6.
  */
  constexpr double eps = 1e-12;
  const Point3DType& ptInCell(verts[0]);
  for(const auto& tet : hexTets)
  {
    if(tet.contains(ptInCell, eps))
    {
      return LabelType::LABEL_IN;
    }
  }
  return LabelType::LABEL_OUT;
}

bool HexClipper::getGeometryAsTets(quest::experimental::ShapeMesh& shapeMesh,
                                   axom::Array<TetrahedronType>& tets)
{
  int allocId = shapeMesh.getAllocatorID();
  if(tets.getAllocatorID() != allocId || tets.size() != m_tets.size())
  {
    tets = axom::Array<TetrahedronType>(m_tets.size(), m_tets.size(), allocId);
  }
  axom::copy(tets.data(), m_tets.data(), m_tets.size() * sizeof(TetrahedronType));
  return true;
}

void HexClipper::extractClipperInfo()
{
  const auto v0 = m_info.fetch_existing("v0").as_double_array();
  const auto v1 = m_info.fetch_existing("v1").as_double_array();
  const auto v2 = m_info.fetch_existing("v2").as_double_array();
  const auto v3 = m_info.fetch_existing("v3").as_double_array();
  const auto v4 = m_info.fetch_existing("v4").as_double_array();
  const auto v5 = m_info.fetch_existing("v5").as_double_array();
  const auto v6 = m_info.fetch_existing("v6").as_double_array();
  const auto v7 = m_info.fetch_existing("v7").as_double_array();
  for(int d = 0; d < 3; ++d)
  {
    m_hexBeforeTrans[0][d] = v0[d];
    m_hexBeforeTrans[1][d] = v1[d];
    m_hexBeforeTrans[2][d] = v2[d];
    m_hexBeforeTrans[3][d] = v3[d];
    m_hexBeforeTrans[4][d] = v4[d];
    m_hexBeforeTrans[5][d] = v5[d];
    m_hexBeforeTrans[6][d] = v6[d];
    m_hexBeforeTrans[7][d] = v7[d];
  }
}

/*
  Compute the triangulated surface of the hex.  There 4 triangles per hex face.
  Each touches two vertices on the face and the face centroid.
  All are orientated inward.
*/
void HexClipper::computeSurface()
{
  // Hex vertex shorthands
  // See the Hexahedron class documentation, especially the ASCII art.
  const auto& p = m_hex[0];
  const auto& q = m_hex[1];
  const auto& r = m_hex[2];
  const auto& s = m_hex[3];
  const auto& t = m_hex[4];
  const auto& u = m_hex[5];
  const auto& v = m_hex[6];
  const auto& w = m_hex[7];

  // 6 face centroids.
  Point3DType pswt(axom::NumericArray<double, 3> {p.array() + s.array() + w.array() + t.array()} / 4);
  Point3DType quvr(axom::NumericArray<double, 3> {q.array() + u.array() + v.array() + r.array()} / 4);

  Point3DType ptuq(axom::NumericArray<double, 3> {p.array() + t.array() + u.array() + q.array()} / 4);
  Point3DType srvw(axom::NumericArray<double, 3> {s.array() + r.array() + v.array() + w.array()} / 4);

  Point3DType pqrs(axom::NumericArray<double, 3> {p.array() + q.array() + r.array() + s.array()} / 4);
  Point3DType twvu(axom::NumericArray<double, 3> {t.array() + w.array() + v.array() + u.array()} / 4);

  m_surfaceTriangles[0] = Triangle3DType(pswt, p, s);
  m_surfaceTriangles[1] = Triangle3DType(pswt, s, w);
  m_surfaceTriangles[2] = Triangle3DType(pswt, w, t);
  m_surfaceTriangles[3] = Triangle3DType(pswt, t, p);

  m_surfaceTriangles[4] = Triangle3DType(quvr, q, u);
  m_surfaceTriangles[5] = Triangle3DType(quvr, u, v);
  m_surfaceTriangles[6] = Triangle3DType(quvr, v, r);
  m_surfaceTriangles[7] = Triangle3DType(quvr, r, q);

  m_surfaceTriangles[8] = Triangle3DType(ptuq, p, t);
  m_surfaceTriangles[9] = Triangle3DType(ptuq, t, u);
  m_surfaceTriangles[10] = Triangle3DType(ptuq, u, q);
  m_surfaceTriangles[11] = Triangle3DType(ptuq, q, p);

  m_surfaceTriangles[12] = Triangle3DType(srvw, s, r);
  m_surfaceTriangles[13] = Triangle3DType(srvw, r, v);
  m_surfaceTriangles[14] = Triangle3DType(srvw, v, w);
  m_surfaceTriangles[15] = Triangle3DType(srvw, w, s);

  m_surfaceTriangles[16] = Triangle3DType(pqrs, p, q);
  m_surfaceTriangles[17] = Triangle3DType(pqrs, q, r);
  m_surfaceTriangles[18] = Triangle3DType(pqrs, r, s);
  m_surfaceTriangles[19] = Triangle3DType(pqrs, s, p);

  m_surfaceTriangles[20] = Triangle3DType(twvu, t, w);
  m_surfaceTriangles[21] = Triangle3DType(twvu, w, v);
  m_surfaceTriangles[22] = Triangle3DType(twvu, v, u);
  m_surfaceTriangles[23] = Triangle3DType(twvu, u, t);
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
