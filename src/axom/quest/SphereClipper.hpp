// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SPHERECLIPPER_HPP
#define AXOM_QUEST_SPHERECLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief GeometryClipper specialized for sphere geometries.
*/
class SphereClipper : public GeometryClipperStrategy
{
public:
  /*!
    @brief Constructor.

    @param [in] kGeom Describes the shape to place
      into the mesh.
    @param [in] name To override the default strategy name
  */
  SphereClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~SphereClipper() = default;

  const std::string& name() const override { return m_name; }

  /*!
    The algorithm isn't perfect and may mislabel a cell as outside
    when all its vertices are outside.  We check whether segments of
    the cell's tetrahedral representation touches the sphere, but we
    don't check whether their faces do.

    This problem can be fixed but it's significantly more complexity
    to reduce an error that is probably O(h^3).  The larger error
    would likely be from discretizing the sphere (which can be reduced
    with more refinement).
  */
  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getGeometryAsOcts(quest::ShapeeMesh& shappeMesh,
                         axom::Array<axom::primal::Octahedron<double, 3>>& octs) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //!@brief Sphere before external transformations.
  axom::primal::Sphere<double, 3> m_sphereBeforeTrans;

  //!@brief Sphere after external transformations from m_transMat.
  axom::primal::Sphere<double, 3> m_sphere;

  //!@brief External transformations, equivalent to m_transMat.
  axom::primal::CoordinateTransformer<double> m_transformer;

  int m_levelOfRefinement;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  void extractClipperInfo();

  void transformSphere();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SPHERECLIPPER_HPP
