// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SPHERECLIPPER_HPP
#define AXOM_QUEST_SPHERECLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Geometry clipping operations for sphere geometries.
 */
class SphereClipper : public MeshClipperStrategy
{
public:
  /*!
   * @brief Constructor.
   *
   * @param [in] kGeom Describes the shape to place
   *   into the mesh.
   * @param [in] name To override the default strategy name
   *
   * \c kGeom.asHierarchy() must contain the following data:
   * - center: 3D coordinates of the center
   * - radius: radius of the sphere
   * - levelOfRefinement: number of refinement levels used
   *   to approximate the sphere with octahedra.  The number
   *   of octs grows exponentially with this level.
   *   @see discretize(const SphereType& sphere, int levels, axom::Array<OctType>& out, int& octcount).
   *   In practice, keep this number to 11 or less.
   */
  SphereClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~SphereClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelCellsInOut(quest::experimental::ShapeMesh& shappeMesh,
                       axom::Array<LabelType>& label) override;

  bool labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                      axom::ArrayView<const axom::IndexType> cellIds,
                      axom::Array<LabelType>& tetLabels) override;

  /*!
   * @brief Compute overlap volumes for tetrahedra labeled as boundary tets.
   *
   * This is a sphere-specific fast path used by MeshClipper after the
   * tet screening step has already removed interior and exterior tets.
   * The implementation avoids the generic octahedron/BVH clipping path by
   * working directly against the analytic sphere.
   */
  bool specializedClipTets(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<double> ovlap,
                           const axom::ArrayView<IndexType>& tetIds,
                           conduit::Node& statistics) override;

  bool getGeometryAsOcts(quest::experimental::ShapeMesh& shappeMesh,
                         axom::Array<axom::primal::Octahedron<double, 3>>& octs) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Sphere before external transformations.
  SphereType m_sphereBeforeTrans;

  //!@brief Sphere after external transformations.
  SphereType m_sphere;

  //!@brief External transformations.
  axom::primal::experimental::CoordinateTransformer<double> m_transformer;

  int m_levelOfRefinement;

  template <typename ExecSpace>
  void labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<LabelType> label);

  template <typename ExecSpace>
  void labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                          axom::ArrayView<const axom::IndexType> cellIds,
                          axom::ArrayView<LabelType> tetLabels);

  /*!
   * @brief Execution-space implementation of the sphere-specific boundary-tet clip.
   *
   * Each boundary tet is adaptively subdivided into smaller tets. Fully interior
   * and exterior subtets are accumulated cheaply, while only the smallest
   * unresolved subtets use a local linearization of the sphere boundary.
   */
  template <typename ExecSpace>
  void specializedClipTetsImpl(quest::experimental::ShapeMesh& shapeMesh,
                               axom::ArrayView<double> ovlap,
                               const axom::ArrayView<IndexType>& tetIds,
                               conduit::Node& statistics);

  //!@brief Compute LabelType for a polyhedron (hex or tet in our case).
  template <typename Polyhedron>
  AXOM_HOST_DEVICE static inline MeshClipperStrategy::LabelType polyhedronToLabel(
    const Polyhedron& verts,
    const SphereType& sphere);

  void extractClipperInfo();

  void transformSphere();
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SPHERECLIPPER_HPP
