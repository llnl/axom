// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_PLANE3DCLIPPER_HPP
#define AXOM_QUEST_PLANE3DCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/MeshClipperStrategy.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

/*!
 * @brief Geometry clipping operations for plane geometries.
*/
class Plane3DClipper : public MeshClipperStrategy
{
public:
  /*!
   * @brief Constructor.
   *
   * @param [in] kGeom Describes the shape to place
   *   into the mesh.
   * @param [in] name To override the default strategy name
   *
   * Clipping operations for a semi-infinite half-space
   * on the positive normal direction of a plane.
   *
   * @internal Because this class provides screening via the
   * labelCellsInOut and labelTetsInOut methods, the
   * specializedClipCells methods below are not essential.  They are
   * implemented only to avoid crashing when the MeshClipper's screen
   * level is overridden (for performance comparisons).  The
   * specializedClipTets method is much faster and is the one used for
   * the default screen level.  If the screen level override is remove,
   * the non-essential methods can also be removed.
   */
  Plane3DClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~Plane3DClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelCellsInOut(quest::experimental::ShapeMesh& shappeMesh,
                       axom::Array<LabelType>& label) override;

  bool labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                      axom::ArrayView<const axom::IndexType> cellIds,
                      axom::Array<LabelType>& tetLabels) override;

  bool specializedClipCells(quest::experimental::ShapeMesh& shappeMesh,
                            axom::ArrayView<double> ovlap,
                            conduit::Node& statistics) override;

  bool specializedClipCells(quest::experimental::ShapeMesh& shappeMesh,
                            axom::ArrayView<double> ovlap,
                            const axom::ArrayView<IndexType>& cellIds,
                            conduit::Node& statistics) override;

  bool specializedClipTets(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<double> ovlap,
                           const axom::ArrayView<IndexType>& tetIds,
                           conduit::Node& statistics) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  axom::primal::Plane<double, 3> m_plane;

  template <typename ExecSpace>
  void labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                           axom::ArrayView<LabelType> label);

  template <typename ExecSpace>
  void labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                          axom::ArrayView<const axom::IndexType> cellIds,
                          axom::ArrayView<LabelType> label);

  template <typename ExecSpace>
  void specializedClipCellsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                axom::ArrayView<double> ovlap,
                                conduit::Node& statistics);

  template <typename ExecSpace>
  void specializedClipCellsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                axom::ArrayView<double> ovlap,
                                const axom::ArrayView<IndexType>& cellIds,
                                conduit::Node& statistics);

  template <typename ExecSpace>
  void specializedClipTetsImpl(quest::experimental::ShapeMesh& shapeMesh,
                               axom::ArrayView<double> ovlap,
                               const axom::ArrayView<IndexType>& tetIds,
                               conduit::Node& statistics);

  void extractClipperInfo();
};

}  // namespace experimental
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_PLANE3DCLIPPER_HPP
