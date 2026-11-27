// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file GregoryMesh.hpp
 * \brief Named collections of Gregory facets built from a Blueprint surface mesh.
 */

#include "axom/core/Array.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/GregoryPatch.hpp"
#include "axom/primal/geometry/GregoryTriangle.hpp"

#ifndef AXOM_USE_CONDUIT
  #error "GregoryMesh usage requires Conduit"
#endif

#include "conduit_node.hpp"

#include <array>
#include <string>
#include <variant>

namespace axom::quest
{

enum class GregoryFacetShape
{
  Tri,
  Quad
};

/*!
* \brief A wrapper for a Gregory patch (quadrilateral) or triangle.
*/
struct GregoryFacet
{
  using GregoryPatch = primal::GregoryPatch<double>;
  using GregoryTriangle = primal::GregoryTriangle<double>;

  // Get types from the patch (although they should match the triangle too)
  using PointType = GregoryPatch::PointType;
  using VectorType = GregoryPatch::VectorType;
  using BoundingBoxType = GregoryPatch::BoundingBoxType;
  using GeometryType = std::variant<GregoryTriangle, GregoryPatch>;

  // Store a list of coordinate IDs to simplify identification of adjacent facets.
  std::array<int, 4> node_ids {};

  // Store a variant for the type of gregory facet
  GregoryFacetShape shape {GregoryFacetShape::Tri};
  GeometryType geometry;

  int numVertices() const { return shape == GregoryFacetShape::Tri ? 3 : 4; }

  // Visit-based methods which dispatch to the respective geometry type
  PointType evaluate(double u, double v) const
  {
    return std::visit([&](const auto& facet) { return facet.evaluate(u, v); }, geometry);
  }

  BoundingBoxType boundingBox() const
  {
    return std::visit([](const auto& facet) { return facet.boundingBox(); }, geometry);
  }
};

/*!
 * \brief A named collection of Gregory facets (quads or tris).
 *
 * Each surface is derived from a single Blueprint topology, which may
 *  contain multiple disconnected pieces.
 */
struct GregorySurface
{
  std::string name;
  axom::Array<GregoryFacet> facets;
};

/*!
 * \brief A collection of Gregory surfaces constructed from a Blueprint mesh.
 */
class GregoryMesh
{
public:
  using PointType = GregoryFacet::PointType;
  using VectorType = GregoryFacet::VectorType;

  /*!
   * \brief Build named Gregory surfaces from a Blueprint mesh.
   *
   * Assumes the surfaces share one explicit 3D coordinate set with possibly multiple
   * unstructured tri, quad, or mixed tri/quad topologies. If multiple domains
   * are supplied, only the first is used. Unsupported or degenerate topologies
   * are skipped with a warning. Blueprint-invalid domains and unusable shared
   * coordinates are reported through SLIC_ERROR.
   *
   * \param[in] mesh The input Conduit Blueprint mesh.
   */
  explicit GregoryMesh(const conduit::Node& mesh);

  /// List of gregory surfaces in the order they're read in
  const axom::Array<GregorySurface>& getSurfaces() const { return m_surfaces; }

  /*!
   * \brief Write a triangulated view of the Gregory surfaces to a VTK file.
   *
   * Each facet is sampled independently. The cell centered `surface_id` field
   * contains the zero based index into getSurfaces(), which identifies the
   * original Blueprint topology name. The cell centered `facet_id` field
   * contains a zero based facet index across surfaces in that order.
   *
   * \param[in] file_path Output VTK file path.
   * \param[in] subdivisions Number of subdivisions along each parameter edge;
   *                         must be positive.
   */
  void writeVtkFile(const std::string& file_path, int subdivisions = 8) const;

private:
  axom::Array<GregorySurface> m_surfaces;
};

}  // namespace axom::quest
