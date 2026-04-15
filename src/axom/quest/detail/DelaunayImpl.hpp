// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_DETAIL_DELAUNAY_IMPL_HPP_
#define AXOM_QUEST_DETAIL_DELAUNAY_IMPL_HPP_

namespace axom
{
namespace quest
{

template <int DIM>
template <typename FacetSubsetType>
inline typename Delaunay<DIM>::FacetKey Delaunay<DIM>::makeSortedFaceKey(const FacetSubsetType& facet)
{
  FacetKey key {};
  for(int i = 0; i < VERTS_PER_FACET; ++i)
  {
    key[i] = facet[i];
  }
  std::sort(key.begin(), key.end());
  return key;
}

template <int DIM>
inline std::string Delaunay<DIM>::facetKeyString(const FacetKey& facet_key)
{
  return fmt::format("[{}]", fmt::join(facet_key, ", "));
}

template <int DIM>
inline double Delaunay<DIM>::getBoundaryCoordinateTolerance() const
{
  const auto min_pt = m_bounding_box.getMin();
  const auto max_pt = m_bounding_box.getMax();

  double max_extent = 1.;
  for(int dim = 0; dim < DIM; ++dim)
  {
    max_extent = axom::utilities::max(max_extent, max_pt[dim] - min_pt[dim]);
  }

  return 64. * std::numeric_limits<double>::epsilon() * max_extent;
}

template <int DIM>
inline bool Delaunay<DIM>::isFacetOnBoundingBox(const FacetKey& facet_key) const
{
  const auto& min_pt = m_bounding_box.getMin();
  const auto& max_pt = m_bounding_box.getMax();
  const double tol = getBoundaryCoordinateTolerance();

  for(int dim = 0; dim < DIM; ++dim)
  {
    bool on_min_face = true;
    bool on_max_face = true;
    for(const IndexType vertex_idx : facet_key)
    {
      const auto& vertex = m_mesh.getVertexPosition(vertex_idx);
      on_min_face &= axom::utilities::abs(vertex[dim] - min_pt[dim]) <= tol;
      on_max_face &= axom::utilities::abs(vertex[dim] - max_pt[dim]) <= tol;
    }

    if(on_min_face || on_max_face)
    {
      return true;
    }
  }

  return false;
}

template <int DIM>
inline double Delaunay<DIM>::getElementSignedMeasure(IndexType element_idx) const
{
  if constexpr(DIM == 2)
  {
    return getElement(element_idx).signedArea();
  }
  else
  {
    return getElement(element_idx).signedVolume();
  }
}

template <int DIM>
inline double Delaunay<DIM>::getElementMeasureTolerance() const
{
  const double scale = axom::utilities::max(1., m_bounding_box.range().norm());
  return 256. * std::numeric_limits<double>::epsilon() * std::pow(scale, static_cast<double>(DIM));
}

template <int DIM>
inline void Delaunay<DIM>::validateInsertionSeed(IndexType element_idx,
                                                 const PointType& query_pt,
                                                 const BaryCoordType& bary_coord,
                                                 const IndexArray& seed_elements) const
{
  if(m_insertion_validation_mode == InsertionValidationMode::None)
  {
    return;
  }

  fmt::memory_buffer out;
  bool valid = true;

  if(!isSearchableElement(element_idx))
  {
    fmt::format_to(std::back_inserter(out), "\n\tContaining element {} is not searchable", element_idx);
    valid = false;
  }
  else
  {
    const auto verts = m_mesh.boundaryVertices(element_idx);
    for(auto idx : verts.positions())
    {
      if(!m_mesh.isValidVertex(verts[idx]))
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tContaining element {} references invalid vertex {}",
                       element_idx,
                       verts[idx]);
        valid = false;
        break;
      }
    }
  }

  if(!isPointInsideForLocation(element_idx, query_pt, bary_coord))
  {
    fmt::format_to(std::back_inserter(out),
                   "\n\tPoint {} is not inside containing element {}",
                   query_pt,
                   element_idx);
    valid = false;
  }

  if(seed_elements.empty())
  {
    fmt::format_to(std::back_inserter(out), "\n\tNo cavity seed elements were generated");
    valid = false;
  }

  bool found_containing_element = false;
  for(const IndexType seed_element : seed_elements)
  {
    found_containing_element |= (seed_element == element_idx);
    if(!m_mesh.isValidElement(seed_element))
    {
      fmt::format_to(std::back_inserter(out), "\n\tSeed element {} is invalid", seed_element);
      valid = false;
    }
  }

  if(!found_containing_element)
  {
    fmt::format_to(std::back_inserter(out),
                   "\n\tContaining element {} is missing from the seed set",
                   element_idx);
    valid = false;
  }

  SLIC_ERROR_IF(!valid, "Delaunay insertion seed validation failed:" << fmt::to_string(out));
}

template <int DIM>
inline bool Delaunay<DIM>::isSearchableElement(IndexType element_idx) const
{
  return m_mesh.isValidElement(element_idx);
}

template <int DIM>
inline double Delaunay<DIM>::orientationTolerance(const std::array<PointType, 4>& pts)
{
  const double scale = getPointMagnitudeScale(pts);
  if constexpr(DIM == 2)
  {
    return 64. * std::numeric_limits<double>::epsilon() * scale * scale;
  }
  else
  {
    return 64. * std::numeric_limits<double>::epsilon() * scale * scale * scale;
  }
}

template <int DIM>
inline double Delaunay<DIM>::determinant3(const PointType& p0, const PointType& p1, const PointType& p2)
{
  return axom::numerics::determinant(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
}

template <int DIM>
inline double Delaunay<DIM>::orientationDeterminant(const std::array<PointType, 4>& pts)
{
  return axom::numerics::determinant(pts[0][0],
                                     pts[0][1],
                                     pts[0][2],
                                     1.,
                                     pts[1][0],
                                     pts[1][1],
                                     pts[1][2],
                                     1.,
                                     pts[2][0],
                                     pts[2][1],
                                     pts[2][2],
                                     1.,
                                     pts[3][0],
                                     pts[3][1],
                                     pts[3][2],
                                     1.);
}

template <int DIM>
inline int Delaunay<DIM>::symbolicOrientationSign(const std::array<PointType, 4>& pts,
                                                  const std::array<IndexType, 4>& ranks)
{
  const double det = orientationDeterminant(pts);
  const int det_sign = signWithTolerance(det, orientationTolerance(pts));
  if(det_sign != 0)
  {
    return det_sign;
  }

  const std::array<double, 4> cofactors {-determinant3(pts[1], pts[2], pts[3]),
                                         determinant3(pts[0], pts[2], pts[3]),
                                         -determinant3(pts[0], pts[1], pts[3]),
                                         determinant3(pts[0], pts[1], pts[2])};

  std::array<int, 4> order {{0, 1, 2, 3}};
  std::sort(order.begin(), order.end(), [&](int lhs, int rhs) { return ranks[lhs] < ranks[rhs]; });

  const double cofactor_tol = 64. * std::numeric_limits<double>::epsilon() *
    axom::utilities::max(1., orientationTolerance(pts));
  for(const int row : order)
  {
    const int sign = signWithTolerance(cofactors[row], cofactor_tol);
    if(sign != 0)
    {
      return sign;
    }
  }

  return 0;
}

template <int DIM>
inline typename Delaunay<DIM>::CircumsphereEval Delaunay<DIM>::evaluateCircumsphereOnMesh(
  const IAMeshType& mesh,
  IndexType element_idx)
{
  using axom::numerics::determinant;

  const auto verts = mesh.boundaryVertices(element_idx);

  if constexpr(DIM == 2)
  {
    const PointType& p0 = mesh.getVertexPosition(verts[0]);
    const PointType& p1 = mesh.getVertexPosition(verts[1]);
    const PointType& p2 = mesh.getVertexPosition(verts[2]);

    const double vx0 = p1[0] - p0[0];
    const double vx1 = p2[0] - p0[0];
    const double vy0 = p1[1] - p0[1];
    const double vy1 = p2[1] - p0[1];

    const double sq0 = vx0 * vx0 + vy0 * vy0;
    const double sq1 = vx1 * vx1 + vy1 * vy1;

    const double a = determinant(vx0, vx1, vy0, vy1);
    const double eps = (a >= 0.) ? primal::PRIMAL_TINY : -primal::PRIMAL_TINY;
    const double ood = 1. / (2. * a + eps);

    const double center_offset_x = determinant(sq0, sq1, vy0, vy1) * ood;
    const double center_offset_y = -determinant(sq0, sq1, vx0, vx1) * ood;

    return CircumsphereEval(p0, VectorType {center_offset_x, center_offset_y});
  }
  else
  {
    const PointType& p0 = mesh.getVertexPosition(verts[0]);
    const PointType& p1 = mesh.getVertexPosition(verts[1]);
    const PointType& p2 = mesh.getVertexPosition(verts[2]);
    const PointType& p3 = mesh.getVertexPosition(verts[3]);

    const double vx0 = p1[0] - p0[0];
    const double vx1 = p2[0] - p0[0];
    const double vx2 = p3[0] - p0[0];
    const double vy0 = p1[1] - p0[1];
    const double vy1 = p2[1] - p0[1];
    const double vy2 = p3[1] - p0[1];
    const double vz0 = p1[2] - p0[2];
    const double vz1 = p2[2] - p0[2];
    const double vz2 = p3[2] - p0[2];

    const double sq0 = vx0 * vx0 + vy0 * vy0 + vz0 * vz0;
    const double sq1 = vx1 * vx1 + vy1 * vy1 + vz1 * vz1;
    const double sq2 = vx2 * vx2 + vy2 * vy2 + vz2 * vz2;

    const double a = determinant(vx0, vx1, vx2, vy0, vy1, vy2, vz0, vz1, vz2);
    const double eps = (a >= 0.) ? primal::PRIMAL_TINY : -primal::PRIMAL_TINY;
    const double ood = 1. / (2. * a + eps);

    const double center_offset_x = determinant(sq0, sq1, sq2, vy0, vy1, vy2, vz0, vz1, vz2) * ood;
    const double center_offset_y = determinant(sq0, sq1, sq2, vz0, vz1, vz2, vx0, vx1, vx2) * ood;
    const double center_offset_z = determinant(sq0, sq1, sq2, vx0, vx1, vx2, vy0, vy1, vy2) * ood;

    return CircumsphereEval(p0, VectorType {center_offset_x, center_offset_y, center_offset_z});
  }
}

template <int DIM>
inline double Delaunay<DIM>::sphereSquaredDistanceTolerance(const CircumsphereEval& sphere,
                                                            const PointType& x,
                                                            double distance_sq)
{
  double scale = 1.;

  for(int dim = 0; dim < DIM; ++dim)
  {
    scale = axom::utilities::max(scale, axom::utilities::abs(sphere.center[dim]));
    scale = axom::utilities::max(scale, axom::utilities::abs(x[dim]));
  }

  const double local_span =
    axom::utilities::max(1., 2. * std::sqrt(axom::utilities::max(distance_sq, sphere.radius_sq)));

  return 256. * std::numeric_limits<double>::epsilon() * scale * local_span;
}

template <int DIM>
template <typename SphereType>
inline double Delaunay<DIM>::sphereSignedDistanceTolerance(const SphereType& sphere,
                                                           const PointType& x)
{
  const auto& center = sphere.getCenter();
  double scale = axom::utilities::max(1., sphere.getRadius());
  for(int dim = 0; dim < DIM; ++dim)
  {
    scale = axom::utilities::max(scale, axom::utilities::abs(center[dim]));
    scale = axom::utilities::max(scale, axom::utilities::abs(x[dim]));
  }

  return 256. * std::numeric_limits<double>::epsilon() * scale;
}

template <int DIM>
inline int Delaunay<DIM>::inSphereOrientationOnMesh(const IAMeshType& mesh,
                                                    const PointType& q,
                                                    IndexType element_idx)
{
  const auto sphere = evaluateCircumsphereOnMesh(mesh, element_idx);
  const double distance_sq = primal::squared_distance(sphere.center, q);

  const double delta_sq = distance_sq - sphere.radius_sq;
  const double tol = sphereSquaredDistanceTolerance(sphere, q, distance_sq);
  if(axom::utilities::abs(delta_sq) <= tol)
  {
    return primal::ON_BOUNDARY;
  }

  return delta_sq < 0. ? primal::ON_NEGATIVE_SIDE : primal::ON_POSITIVE_SIDE;
}

template <int DIM>
inline double Delaunay<DIM>::inSphereDeterminantOnMesh(const IAMeshType& mesh,
                                                       const PointType& q,
                                                       IndexType element_idx)
{
  const auto verts = mesh.boundaryVertices(element_idx);

  if constexpr(DIM == 2)
  {
    return primal::in_sphere_determinant(q,
                                         mesh.getVertexPosition(verts[0]),
                                         mesh.getVertexPosition(verts[1]),
                                         mesh.getVertexPosition(verts[2]));
  }
  else
  {
    return primal::in_sphere_determinant(q,
                                         mesh.getVertexPosition(verts[0]),
                                         mesh.getVertexPosition(verts[1]),
                                         mesh.getVertexPosition(verts[2]),
                                         mesh.getVertexPosition(verts[3]));
  }
}

template <int DIM>
inline const char* Delaunay<DIM>::orientationResultName(int result)
{
  switch(result)
  {
  case primal::ON_NEGATIVE_SIDE:
    return "inside";
  case primal::ON_BOUNDARY:
    return "boundary";
  case primal::ON_POSITIVE_SIDE:
    return "outside";
  default:
    return "unknown";
  }
}

template <int DIM>
inline bool Delaunay<DIM>::isPointInSphereOnMesh(const IAMeshType& mesh,
                                                 const PointType& q,
                                                 IndexType element_idx,
                                                 bool includeBoundary)
{
  const int res = inSphereOrientationOnMesh(mesh, q, element_idx);
  return includeBoundary ? (res != primal::ON_POSITIVE_SIDE) : (res == primal::ON_NEGATIVE_SIDE);
}

template <int DIM>
inline int Delaunay<DIM>::classifyOrientationDeterminant(double det, double tol)
{
  const int sign = signWithTolerance(det, tol);
  return sign < 0 ? primal::ON_NEGATIVE_SIDE
                  : (sign > 0 ? primal::ON_POSITIVE_SIDE : primal::ON_BOUNDARY);
}

template <int DIM>
inline typename Delaunay<DIM>::OrientationEval Delaunay<DIM>::evaluateElementOrientationDeterminant(
  IndexType element_idx) const
{
  const double scale = (DIM == 2) ? 2. : 6.;
  const double det = scale * getElementSignedMeasure(element_idx);
  const double tol = scale * getElementMeasureTolerance();
  return {det, tol, classifyOrientationDeterminant(det, tol)};
}

template <int DIM>
inline void Delaunay<DIM>::validateInsertionResult() const
{
  if(m_insertion_validation_mode == InsertionValidationMode::None)
  {
    return;
  }

  if(m_insertion_validation_mode == InsertionValidationMode::Local)
  {
    return;
  }

  if(!m_mesh.isValid(false))
  {
    m_mesh.isValid(true);
    SLIC_ERROR("Delaunay insertion produced an invalid IAMesh");
  }

  if(!isConforming(false))
  {
    isConforming(true);
    SLIC_ERROR("Delaunay insertion produced a non-conforming triangulation");
  }

  if(m_insertion_validation_mode == InsertionValidationMode::Full)
  {
    if(!isValid(false))
    {
      isValid(true);
      SLIC_ERROR(
        "Delaunay insertion produced a triangulation that violates the empty-circumsphere "
        "condition");
    }
  }
}

template <int DIM>
inline bool Delaunay<DIM>::shouldCompactMesh() const
{
  constexpr int MIN_REMOVED_ELEMENTS = DIM == 2 ? 512 : 2048;
  constexpr double REMOVED_ELEMENT_FRACTION = DIM == 2 ? 0.25 : 0.35;
  return static_cast<int>(m_deleted_elements.size()) > MIN_REMOVED_ELEMENTS &&
    (static_cast<double>(m_deleted_elements.size()) >
     REMOVED_ELEMENT_FRACTION * static_cast<double>(m_mesh.elements().size()));
}

template <int DIM>
inline void Delaunay<DIM>::compactMesh()
{
  m_mesh.compact();
  m_deleted_elements.clear();
  m_element_finder.recomputeGrid(m_mesh, m_bounding_box);
  if(m_next_regrid_vertex_count > 0)
  {
    while(m_next_regrid_vertex_count <= m_mesh.vertices().size())
    {
      m_next_regrid_vertex_count *= 2;
    }
  }
}

template <int DIM>
inline void Delaunay<DIM>::writeToVTKFile(const std::string& filename)
{
  const auto CELL_TYPE = DIM == 2 ? mint::TRIANGLE : mint::TET;
  mint::UnstructuredMesh<mint::SINGLE_SHAPE> mint_mesh(DIM, CELL_TYPE);

  this->compactMesh();

  for(auto v : m_mesh.vertices().positions())
  {
    mint_mesh.appendNodes(m_mesh.getVertexPosition(v).data(), 1);
  }

  for(auto e : m_mesh.elements().positions())
  {
    mint_mesh.appendCell(&(m_mesh.boundaryVertices(e)[0]), CELL_TYPE);
  }

  mint::write_vtk(&mint_mesh, filename);
}

template <int DIM>
inline void Delaunay<DIM>::removeBoundary()
{
  if(m_has_boundary)
  {
    const int num_boundary_pts = 1 << DIM;

    for(auto e : m_mesh.elements().positions())
    {
      if(m_mesh.isValidElement(e))
      {
        const auto verts = m_mesh.boundaryVertices(e);
        bool touches_boundary = false;
        for(int i = 0; i < VERT_PER_ELEMENT; ++i)
        {
          touches_boundary |= (verts[i] >= 0 && verts[i] < num_boundary_pts);
        }

        if(touches_boundary)
        {
          m_mesh.removeElement(e);
        }
      }
    }

    for(int v = 0; v < num_boundary_pts; ++v)
    {
      m_mesh.removeVertex(v);
    }

    for(auto e : m_mesh.elements().positions())
    {
      if(!m_mesh.isValidElement(e))
      {
        continue;
      }

      const auto verts = m_mesh.boundaryVertices(e);
      bool has_invalid_vertex = false;
      for(int i = 0; i < VERT_PER_ELEMENT; ++i)
      {
        has_invalid_vertex |= !m_mesh.isValidVertex(verts[i]);
      }

      if(has_invalid_vertex)
      {
        m_mesh.removeElement(e);
      }
    }

    this->compactMesh();
    m_has_boundary = false;
  }
}

template <int DIM>
inline bool Delaunay<DIM>::isValid(bool verboseOutput) const
{
  using ImplicitGridType = spin::ImplicitGrid<DIM, axom::SEQ_EXEC, IndexType>;
  using UniformGridType = spin::UniformGrid<IndexType, DIM>;
  using NumericArrayType = NumericArray<DataType, DIM>;
  using axom::numerics::dot_product;

  bool valid = true;

  std::vector<std::pair<IndexType, IndexType>> invalidEntries;
  std::vector<IndexType> invalidElements;

  const IndexType totalVertices = m_mesh.vertices().size();
  const IndexType totalElements = m_mesh.elements().size();
  const IndexType res = axom::utilities::ceil(0.33 * std::pow(totalVertices, 1. / DIM));
  UniformGridType grid(m_bounding_box, NumericArray<int, DIM>(res).data());

  axom::Array<typename ElementType::SphereType> circumspheres(totalElements);

  auto vertexInsideCircumsphere = [&](const PointType& vertex, IndexType element_idx) {
    return isPointInSphereOnMesh(m_mesh, vertex, element_idx, /*includeBoundary=*/false);
  };

  auto circumsphereSignedDistanceTol = [](const typename ElementType::SphereType& sphere,
                                          const PointType& x) {
    const auto& center = sphere.getCenter();
    double scale = axom::utilities::max(1., sphere.getRadius());
    for(int dim = 0; dim < DIM; ++dim)
    {
      scale = axom::utilities::max(scale, axom::utilities::abs(center[dim]));
      scale = axom::utilities::max(scale, axom::utilities::abs(x[dim]));
    }

    return 256. * std::numeric_limits<double>::epsilon() * scale;
  };

  {
    using GridCell = typename ImplicitGridType::GridCell;

    const auto resCell = GridCell(res);
    ImplicitGridType implicitGrid(m_bounding_box, &resCell, totalElements);
    for(auto element_idx : m_mesh.elements().positions())
    {
      if(m_mesh.isValidElement(element_idx))
      {
        const auto verts = m_mesh.boundaryVertices(element_idx);
        bool has_all_vertices = true;
        for(int i = 0; i < VERT_PER_ELEMENT; ++i)
        {
          has_all_vertices &= m_mesh.isValidVertex(verts[i]);
        }

        if(!has_all_vertices)
        {
          valid = false;
          if(verboseOutput)
          {
            invalidElements.push_back(element_idx);
          }
          continue;
        }

        circumspheres[element_idx] = this->getElement(element_idx).circumsphere();
        const auto& sphere = circumspheres[element_idx];
        const auto& center = sphere.getCenter().array();
        const auto offset = NumericArrayType(sphere.getRadius());

        BoundingBox bb;
        bb.addPoint(PointType(center - offset));
        bb.addPoint(PointType(center + offset));

        implicitGrid.insert(bb, element_idx);
      }
    }

    const int kUpper = (DIM == 2) ? 0 : res;
    const IndexType stride[3] = {1, res, (DIM == 2) ? 0 : res * res};
    for(IndexType k = 0; k < kUpper; ++k)
    {
      for(IndexType j = 0; j < res; ++j)
      {
        for(IndexType i = 0; i < res; ++i)
        {
          const IndexType vals[3] = {i, j, k};
          const GridCell cell(vals);
          const auto idx = dot_product(cell.data(), stride, DIM);
          const auto binValues = implicitGrid.getCandidatesAsArray(cell);
          grid.getBinContents(idx).insert(0, binValues.size(), binValues.data());
        }
      }
    }
  }

  for(auto vertex_idx : m_mesh.vertices().positions())
  {
    if(!m_mesh.isValidVertex(vertex_idx))
    {
      continue;
    }

    const auto& vertex = m_mesh.getVertexPosition(vertex_idx);
    for(const auto element_idx : grid.getBinContents(grid.getBinIndex(vertex)))
    {
      if(slam::is_subset(vertex_idx, m_mesh.boundaryVertices(element_idx)))
      {
        continue;
      }

      if(vertexInsideCircumsphere(vertex, element_idx))
      {
        valid = false;

        if(verboseOutput)
        {
          invalidEntries.push_back(std::make_pair(vertex_idx, element_idx));
        }
      }
    }
  }

  if(verboseOutput)
  {
    if(valid)
    {
      SLIC_INFO("Delaunay complex was valid");
    }
    else
    {
      fmt::memory_buffer out;
      if(!invalidElements.empty())
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\t{} valid elements referenced invalid vertices: {}",
                       invalidElements.size(),
                       fmt::join(invalidElements, ", "));
      }
      for(const auto& pr : invalidEntries)
      {
        const auto vertex_idx = pr.first;
        const auto element_idx = pr.second;
        const auto& pos = m_mesh.getVertexPosition(vertex_idx);
        const auto element = this->getElement(element_idx);
        const auto circumsphere = element.circumsphere();
        const double tol = circumsphereSignedDistanceTol(circumsphere, pos);
        fmt::format_to(std::back_inserter(out),
                       "\n\tVertex {} @ {}"
                       "\n\tElement {}: {} w/ circumsphere: {}"
                       "\n\tDistance to circumcenter: {} (tol={})",
                       vertex_idx,
                       pos,
                       element_idx,
                       element,
                       circumsphere,
                       circumsphere.computeSignedDistance(pos),
                       tol);
      }

      SLIC_INFO(
        fmt::format("Delaunay complex was NOT valid. There were {} "
                    "vertices in the circumsphere of an element. {}",
                    invalidEntries.size(),
                    fmt::to_string(out)));
    }
  }

  return valid;
}

template <int DIM>
inline bool Delaunay<DIM>::isConforming(bool verboseOutput) const
{
  fmt::memory_buffer out;

  bool valid = m_mesh.isConforming(verboseOutput);

  for(auto element_idx : m_mesh.elements().positions())
  {
    if(!m_mesh.isValidElement(element_idx))
    {
      continue;
    }

    const auto orient = evaluateElementOrientationDeterminant(element_idx);
    if(orient.orientation != primal::ON_POSITIVE_SIDE)
    {
      if(verboseOutput)
      {
        fmt::format_to(
          std::back_inserter(out),
          "\n\tElement {} has non-positive orientation determinant {:.17g} (tol={:.3g})",
          element_idx,
          orient.det,
          orient.tol);
      }
      valid = false;
    }

    if(m_has_boundary)
    {
      const auto neighbors = m_mesh.adjacentElements(element_idx);
      for(int facet_idx = 0; facet_idx < VERT_PER_ELEMENT; ++facet_idx)
      {
        if(m_mesh.isValidElement(neighbors[facet_idx]))
        {
          continue;
        }

        const FacetKey facet_key = m_mesh.getSortedFacetKey(element_idx, facet_idx);
        if(!isFacetOnBoundingBox(facet_key))
        {
          if(verboseOutput)
          {
            fmt::format_to(std::back_inserter(out),
                           "\n\tBoundary facet {} is not on the initial bounding box",
                           facetKeyString(facet_key));
          }
          valid = false;
        }
      }
    }
  }

  if(verboseOutput)
  {
    if(valid)
    {
      SLIC_INFO("Delaunay mesh was conforming");
    }
    else
    {
      SLIC_INFO("Delaunay mesh was NOT conforming. Summary: " << fmt::to_string(out));
    }
  }

  return valid;
}

template <int DIM>
inline void Delaunay<DIM>::initializeBoundary(const BoundingBox& bb)
{
  std::vector<DataType> points;
  IndexArray elem;

  generateInitialMesh(points, elem, bb);

  m_mesh = IAMeshType(points, elem);
  m_element_finder.recomputeGrid(m_mesh, bb);
  m_next_regrid_vertex_count = 1024;
  m_walk_visited = slam::BitSet(static_cast<int>(m_mesh.elements().size()));
  m_total_removed_elements = 0;
  m_max_removed_elements = 0;
  m_num_insertions = 0;
  m_num_walk_calls = 0;
  m_num_walk_found = 0;
  m_num_walk_outside = 0;
  m_num_walk_failed = 0;
  m_total_walk_steps = 0;
  m_max_walk_steps = 0;
  m_num_linear_fallbacks = 0;
  m_num_empty_seed_fallbacks = 0;

  m_candidate_elements_scratch.clear();
  m_candidate_elements_scratch.reserve(QUERY_CANDIDATE_LIMIT);
  m_walked_elements_scratch.clear();
  m_walked_elements_scratch.reserve(256);
  m_walk_local_elements_scratch.clear();
  m_walk_local_elements_scratch.reserve(256);
  m_initial_vertices_scratch.clear();
  m_initial_vertices_scratch.reserve(8);
  m_fallback_vertices_scratch.clear();
  m_fallback_vertices_scratch.reserve(QUERY_CANDIDATE_LIMIT);
  m_deleted_elements.clear();
  m_deleted_elements.reserve(DIM == 2 ? 128 : 512);

  if(!m_insertion_helper)
  {
    m_insertion_helper = std::make_unique<InsertionHelper>(m_mesh);
  }

  m_bounding_box = bb;
  m_has_boundary = true;
}

template <int DIM>
inline void Delaunay<DIM>::insertPoint(const PointType& new_pt)
{
  SLIC_ASSERT_MSG(m_has_boundary, "Error: Need a predefined boundary box prior to adding points.");
  SLIC_ASSERT_MSG(m_insertion_helper != nullptr,
                  "Error: Insertion helper was not initialized. "
                  "Delaunay::initializeBoundary() needs to be called first.");
  SLIC_ASSERT_MSG(m_bounding_box.contains(new_pt),
                  "Error: new point is outside of the boundary box.");

  IndexType element_i = findContainingElement(new_pt);

  if(element_i == INVALID_INDEX)
  {
    SLIC_WARNING(
      fmt::format("Could not insert point {} into Delaunay triangulation: "
                  "Element containing that point was not found",
                  new_pt));
    return;
  }

  const BaryCoordType bary_coord = getBaryCoords(element_i, new_pt);
  const IndexArray seed_elements = getSeedElements(element_i, bary_coord);
  validateInsertionSeed(element_i, new_pt, bary_coord, seed_elements);

  auto& insertionHelper = *m_insertion_helper;
  insertionHelper.reset();
  insertionHelper.containing_element = element_i;
  insertionHelper.containing_bary = bary_coord;
  insertionHelper.seed_elements_debug = seed_elements;
  insertionHelper.findCavityElements(
    new_pt,
    seed_elements,
    [&](const PointType& query_pt, IndexType element_idx) {
      return isPointInSphereOnMesh(m_mesh, query_pt, element_idx, /*includeBoundary=*/true);
    });
  ++m_num_insertions;
  m_total_removed_elements += static_cast<std::uint64_t>(insertionHelper.numRemovedElements());
  m_max_removed_elements =
    axom::utilities::max(m_max_removed_elements,
                         static_cast<std::uint64_t>(insertionHelper.numRemovedElements()));
  validateCavityBoundary(insertionHelper);
  insertionHelper.createCavity(m_deleted_elements);
  IndexType new_pt_i = m_mesh.addVertex(new_pt);
  insertionHelper.delaunayBall(new_pt_i, m_deleted_elements);
  validateInsertedBall(new_pt_i, insertionHelper);
  validateInsertionResult();

  m_element_finder.updateBin(new_pt, new_pt_i);
  if(m_next_regrid_vertex_count > 0 && m_mesh.vertices().size() >= m_next_regrid_vertex_count)
  {
    m_element_finder.recomputeGrid(m_mesh, m_bounding_box);
    while(m_next_regrid_vertex_count <= m_mesh.vertices().size())
    {
      m_next_regrid_vertex_count *= 2;
    }
  }

  if(shouldCompactMesh())
  {
    this->compactMesh();
  }
}

template <int DIM>
inline void Delaunay<DIM>::validateCavityBoundary(const InsertionHelper& insertion_helper) const
{
  if(m_insertion_validation_mode == InsertionValidationMode::None)
  {
    return;
  }

  struct FacetInfo
  {
    IndexType neighbor_idx {INVALID_INDEX};
    bool matched {false};
  };

  fmt::memory_buffer out;
  bool valid = true;
  std::map<FacetKey, FacetInfo> cavity_boundary;

  for(const auto& facet : insertion_helper.boundary_facets)
  {
    const FacetKey facet_key = makeSortedFaceKey(facet.vertices);
    const auto insert_status = cavity_boundary.insert({facet_key, {facet.neighbor, false}});
    if(!insert_status.second)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tCavity boundary facet {} was recorded more than once",
                     facetKeyString(facet_key));
      valid = false;
    }
  }

  for(const IndexType cavity_element : insertion_helper.cavity_elems)
  {
    const auto neighbors = m_mesh.adjacentElements(cavity_element);
    for(int facet_idx = 0; facet_idx < VERT_PER_ELEMENT; ++facet_idx)
    {
      const IndexType neighbor_idx = neighbors[facet_idx];
      if(m_mesh.isValidElement(neighbor_idx) && insertion_helper.containsCavityElement(neighbor_idx))
      {
        continue;
      }

      const FacetKey facet_key = m_mesh.getSortedFacetKey(cavity_element, facet_idx);
      auto facet_it = cavity_boundary.find(facet_key);
      if(facet_it == cavity_boundary.end())
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCavity facet {} on element {} face {} is missing from the "
                       "boundary facet set",
                       facetKeyString(facet_key),
                       cavity_element,
                       facet_idx);
        valid = false;
        continue;
      }

      if(facet_it->second.neighbor_idx != neighbor_idx)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCavity face {} expects neighbor {} but facet relation stores {}",
                       facetKeyString(facet_key),
                       neighbor_idx,
                       facet_it->second.neighbor_idx);
        valid = false;
      }

      if(!m_mesh.isValidElement(neighbor_idx) && m_has_boundary && !isFacetOnBoundingBox(facet_key))
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCavity boundary face {} exits the mesh away from the bounding box",
                       facetKeyString(facet_key));
        valid = false;
      }

      facet_it->second.matched = true;
    }
  }

  for(const auto& facet_entry : cavity_boundary)
  {
    if(!facet_entry.second.matched)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tFacet {} does not correspond to a cavity boundary face",
                     facetKeyString(facet_entry.first));
      valid = false;
    }
  }

  SLIC_ERROR_IF(!valid,
                "Delaunay cavity validation failed after insertion " << m_num_insertions << ":"
                                                                     << fmt::to_string(out));
}

template <int DIM>
inline void Delaunay<DIM>::validateInsertedBall(IndexType new_pt_i,
                                                const InsertionHelper& insertion_helper) const
{
  if(m_insertion_validation_mode == InsertionValidationMode::None)
  {
    return;
  }

  struct BoundaryFacetInfo
  {
    IndexType neighbor_idx {INVALID_INDEX};
    bool matched {false};
  };

  fmt::memory_buffer out;
  bool valid = true;
  std::map<FacetKey, BoundaryFacetInfo> cavity_boundary;
  std::map<FacetKey, std::vector<FacetRecord>> inserted_faces;
  auto findOppositeVertex = [&](IndexType element_idx, const FacetKey& facet_key) {
    const auto verts = m_mesh.boundaryVertices(element_idx);
    for(int i = 0; i < VERT_PER_ELEMENT; ++i)
    {
      bool on_face = false;
      for(int j = 0; j < VERTS_PER_FACET; ++j)
      {
        on_face |= (verts[i] == facet_key[j]);
      }
      if(!on_face)
      {
        return verts[i];
      }
    }

    return INVALID_INDEX;
  };

  if(insertion_helper.inserted_elems.size() != insertion_helper.boundary_facets.size())
  {
    fmt::format_to(
      std::back_inserter(out),
      "\n\tInserted ball element count {} does not match cavity boundary facet count {}",
      insertion_helper.inserted_elems.size(),
      insertion_helper.boundary_facets.size());
    valid = false;
  }

  if(insertion_helper.containing_element != INVALID_INDEX)
  {
    fmt::format_to(std::back_inserter(out),
                   "\n\tContaining element {} barycentric coordinates {}",
                   insertion_helper.containing_element,
                   insertion_helper.containing_bary);
    if(!insertion_helper.seed_elements_debug.empty())
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tInsertion seeds: [{}]",
                     fmt::join(insertion_helper.seed_elements_debug, ", "));
    }
  }

  for(const auto& facet : insertion_helper.boundary_facets)
  {
    cavity_boundary.insert({makeSortedFaceKey(facet.vertices), {facet.neighbor, false}});
  }

  for(const IndexType element_idx : insertion_helper.inserted_elems)
  {
    if(!m_mesh.isValidElement(element_idx))
    {
      fmt::format_to(std::back_inserter(out), "\n\tInserted element {} is invalid", element_idx);
      valid = false;
      continue;
    }

    const auto verts = m_mesh.boundaryVertices(element_idx);
    if(!slam::is_subset(new_pt_i, verts))
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tInserted element {} does not contain the new vertex {}",
                     element_idx,
                     new_pt_i);
      valid = false;
    }

    const auto orient = evaluateElementOrientationDeterminant(element_idx);
    if(orient.orientation != primal::ON_POSITIVE_SIDE)
    {
      fmt::format_to(
        std::back_inserter(out),
        "\n\tInserted element {} has non-positive orientation determinant {:.17g} (tol={:.3g})",
        element_idx,
        orient.det,
        orient.tol);
      valid = false;
    }

    const auto neighbors = m_mesh.adjacentElements(element_idx);
    for(int facet_idx = 0; facet_idx < VERT_PER_ELEMENT; ++facet_idx)
    {
      inserted_faces[m_mesh.getSortedFacetKey(element_idx, facet_idx)].push_back(
        {element_idx, facet_idx, neighbors[facet_idx]});
    }
  }

  for(auto& face_entry : inserted_faces)
  {
    auto cavity_it = cavity_boundary.find(face_entry.first);
    auto& records = face_entry.second;
    if(cavity_it != cavity_boundary.end())
    {
      if(records.size() != 1)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tBoundary face {} of the inserted ball is used by {} new elements",
                       facetKeyString(face_entry.first),
                       records.size());
        valid = false;
        continue;
      }

      if(records.front().neighbor_idx != cavity_it->second.neighbor_idx)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tInserted boundary face {} points to neighbor {} instead of {}",
                       facetKeyString(face_entry.first),
                       records.front().neighbor_idx,
                       cavity_it->second.neighbor_idx);
        valid = false;
      }

      if(m_mesh.isValidElement(cavity_it->second.neighbor_idx))
      {
        const IndexType inserted_opposite =
          findOppositeVertex(records.front().element_idx, face_entry.first);
        const IndexType neighbor_opposite =
          findOppositeVertex(cavity_it->second.neighbor_idx, face_entry.first);

        if(inserted_opposite == INVALID_INDEX || neighbor_opposite == INVALID_INDEX)
        {
          fmt::format_to(
            std::back_inserter(out),
            "\n\tCould not identify opposite vertices across inserted boundary face {}",
            facetKeyString(face_entry.first));
          valid = false;
        }
        else
        {
          const PointType& inserted_point = m_mesh.getVertexPosition(inserted_opposite);
          const PointType& neighbor_point = m_mesh.getVertexPosition(neighbor_opposite);

          if(isPointInSphereOnMesh(m_mesh,
                                   inserted_point,
                                   cavity_it->second.neighbor_idx,
                                   /*includeBoundary=*/false))
          {
            const int query_in_neighbor =
              inSphereOrientationOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx);
            const int inserted_in_neighbor =
              inSphereOrientationOnMesh(m_mesh, inserted_point, cavity_it->second.neighbor_idx);
            fmt::format_to(
              std::back_inserter(out),
              "\n\tInserted boundary face {} leaves new vertex {} inside neighbor {} "
              "circumsphere (query {}={}, det={:.17g}; opposite {}={}, det={:.17g})",
              facetKeyString(face_entry.first),
              inserted_opposite,
              cavity_it->second.neighbor_idx,
              new_pt_i,
              orientationResultName(query_in_neighbor),
              inSphereDeterminantOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx),
              inserted_opposite,
              orientationResultName(inserted_in_neighbor),
              inSphereDeterminantOnMesh(m_mesh, inserted_point, cavity_it->second.neighbor_idx));
            valid = false;
          }

          if(isPointInSphereOnMesh(m_mesh,
                                   neighbor_point,
                                   records.front().element_idx,
                                   /*includeBoundary=*/false))
          {
            const int query_in_neighbor =
              inSphereOrientationOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx);
            const int neighbor_in_inserted =
              inSphereOrientationOnMesh(m_mesh, neighbor_point, records.front().element_idx);
            fmt::format_to(
              std::back_inserter(out),
              "\n\tInserted boundary face {} leaves neighbor vertex {} inside new element {} "
              "circumsphere (query {} in neighbor {} is {}, det={:.17g}; opposite {} in new "
              "element is {}, det={:.17g})",
              facetKeyString(face_entry.first),
              neighbor_opposite,
              records.front().element_idx,
              new_pt_i,
              cavity_it->second.neighbor_idx,
              orientationResultName(query_in_neighbor),
              inSphereDeterminantOnMesh(m_mesh,
                                        m_mesh.getVertexPosition(new_pt_i),
                                        cavity_it->second.neighbor_idx),
              neighbor_opposite,
              orientationResultName(neighbor_in_inserted),
              inSphereDeterminantOnMesh(m_mesh, neighbor_point, records.front().element_idx));
            valid = false;
          }
        }
      }

      cavity_it->second.matched = true;
    }
    else if(records.size() != 2 || records[0].neighbor_idx != records[1].element_idx ||
            records[1].neighbor_idx != records[0].element_idx)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tInternal face {} of the inserted ball has inconsistent adjacency",
                     facetKeyString(face_entry.first));
      valid = false;
    }
    else
    {
      const IndexType lhs_opposite = findOppositeVertex(records[0].element_idx, face_entry.first);
      const IndexType rhs_opposite = findOppositeVertex(records[1].element_idx, face_entry.first);

      if(lhs_opposite == INVALID_INDEX || rhs_opposite == INVALID_INDEX)
      {
        fmt::format_to(std::back_inserter(out),
                       "\n\tCould not identify opposite vertices across inserted internal face {}",
                       facetKeyString(face_entry.first));
        valid = false;
      }
      else
      {
        const PointType& lhs_point = m_mesh.getVertexPosition(lhs_opposite);
        const PointType& rhs_point = m_mesh.getVertexPosition(rhs_opposite);

        if(isPointInSphereOnMesh(m_mesh, lhs_point, records[1].element_idx, /*includeBoundary=*/false))
        {
          const int query_in_rhs = inSphereOrientationOnMesh(m_mesh,
                                                             m_mesh.getVertexPosition(new_pt_i),
                                                             records[1].element_idx);
          const int lhs_in_rhs = inSphereOrientationOnMesh(m_mesh, lhs_point, records[1].element_idx);
          fmt::format_to(
            std::back_inserter(out),
            "\n\tInserted internal face {} leaves vertex {} inside adjacent element {} "
            "circumsphere (query {}={}, det={:.17g}; opposite {}={}, det={:.17g})",
            facetKeyString(face_entry.first),
            lhs_opposite,
            records[1].element_idx,
            new_pt_i,
            orientationResultName(query_in_rhs),
            inSphereDeterminantOnMesh(m_mesh,
                                      m_mesh.getVertexPosition(new_pt_i),
                                      records[1].element_idx),
            lhs_opposite,
            orientationResultName(lhs_in_rhs),
            inSphereDeterminantOnMesh(m_mesh, lhs_point, records[1].element_idx));
          valid = false;
        }

        if(isPointInSphereOnMesh(m_mesh, rhs_point, records[0].element_idx, /*includeBoundary=*/false))
        {
          const int query_in_lhs = inSphereOrientationOnMesh(m_mesh,
                                                             m_mesh.getVertexPosition(new_pt_i),
                                                             records[0].element_idx);
          const int rhs_in_lhs = inSphereOrientationOnMesh(m_mesh, rhs_point, records[0].element_idx);
          fmt::format_to(
            std::back_inserter(out),
            "\n\tInserted internal face {} leaves vertex {} inside adjacent element {} "
            "circumsphere (query {}={}, det={:.17g}; opposite {}={}, det={:.17g})",
            facetKeyString(face_entry.first),
            rhs_opposite,
            records[0].element_idx,
            new_pt_i,
            orientationResultName(query_in_lhs),
            inSphereDeterminantOnMesh(m_mesh,
                                      m_mesh.getVertexPosition(new_pt_i),
                                      records[0].element_idx),
            rhs_opposite,
            orientationResultName(rhs_in_lhs),
            inSphereDeterminantOnMesh(m_mesh, rhs_point, records[0].element_idx));
          valid = false;
        }
      }
    }
  }

  for(const auto& facet_entry : cavity_boundary)
  {
    if(!facet_entry.second.matched)
    {
      fmt::format_to(std::back_inserter(out),
                     "\n\tCavity boundary face {} was not covered by the inserted ball",
                     facetKeyString(facet_entry.first));
      valid = false;
    }
  }

  SLIC_ERROR_IF(!valid,
                "Delaunay ball validation failed after insertion " << m_num_insertions << ":"
                                                                   << fmt::to_string(out));
}

template <>
inline void Delaunay<2>::generateInitialMesh(std::vector<DataType>& points,
                                             std::vector<IndexType>& elem,
                                             const BoundingBox& bb)
{
  const PointType& mins = bb.getMin();
  const PointType& maxs = bb.getMax();

  std::vector<DataType> pt {mins[0], mins[1], mins[0], maxs[1], maxs[0], mins[1], maxs[0], maxs[1]};

  std::vector<IndexType> el {0, 2, 1, 3, 1, 2};

  points.swap(pt);
  elem.swap(el);
}

template <>
inline void Delaunay<3>::generateInitialMesh(std::vector<DataType>& points,
                                             std::vector<IndexType>& elem,
                                             const BoundingBox& bb)
{
  const PointType& mins = bb.getMin();
  const PointType& maxs = bb.getMax();

  std::vector<DataType> pt {mins[0], mins[1], mins[2], mins[0], mins[1], maxs[2], mins[0], maxs[1],
                            mins[2], mins[0], maxs[1], maxs[2], maxs[0], mins[1], mins[2], maxs[0],
                            mins[1], maxs[2], maxs[0], maxs[1], mins[2], maxs[0], maxs[1], maxs[2]};

  std::vector<IndexType> el {3, 2, 4, 0, 3, 4, 1, 0, 3, 2, 6, 4,
                             3, 6, 7, 4, 3, 5, 1, 4, 3, 7, 5, 4};

  points.swap(pt);
  elem.swap(el);
}

}  // namespace quest
}  // namespace axom

#endif
