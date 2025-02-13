// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_PLANE3DCLIPPER_HPP
#define AXOM_QUEST_PLANE3DCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief GeometryClipper specialized for plane geometries.
*/
class Plane3DClipper : public GeometryClipperStrategy
{
public:
  /*!
    @brief Constructor.

    @param [in] kGeom Describes the shape to place
      into the mesh.
    @param [in] name To override the default strategy name

    Clipping operations for a semi-infinite half-space
    on the positive normal direction of a plane.
  */
  Plane3DClipper(const klee::Geometry& kGeom,
                 const std::string& name="");

  virtual ~Plane3DClipper() = default;

  std::string name() const override { return m_name; }

  bool labelInOut(quest::ShapeeMesh& shappeMesh,
                  axom::Array<char>& label) override;

  bool specializedClip(quest::ShapeeMesh& shappeMesh,
                       axom::ArrayView<double> ovlap,
                       const axom::ArrayView<IndexType>& cellIds) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  axom::primal::Plane<double, 3> m_plane;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  template <typename ExecSpace>
  void specializedClipImpl(quest::ShapeeMesh& shapeeMesh,
                           axom::ArrayView<double>& ovlap,
                           const axom::ArrayView<IndexType>& cellIds);
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_PLANE3DCLIPPER_HPP
