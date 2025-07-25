// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SORCLIPPER_HPP
#define AXOM_QUEST_SORCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/quest/FSorClipper.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Geometry clipping operations for 3D
  surface-of-revolution geometries.

  This implementation allows the SOR function to have non-monotonic
  axial coordinates.  For SOR curves that don't, the less complex
  FSorClipper is sufficient.

  The SOR specification may include rotation and translation
  internally, in addition to any external transformation.
*/
class SorClipper : public GeometryClipperStrategy
{
public:
  /*!
    @brief Constructor.

    @param [in] kGeom Describes the shape to place
      into the mesh.
    @param [in] name To override the default strategy name
  */
  SorClipper(const klee::Geometry& kGeom, const std::string& name = "");

  virtual ~SorClipper() = default;

  const std::string& name() const override { return m_name; }

  bool specializedClip(quest::ShapeeMesh& shapeeMesh,
                       axom::ArrayView<double> ovlap) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  axom::Array<std::shared_ptr<FSorClipper>> m_fsorStrategies;

  /*!
    @brief The discrete r(z) curve as an array of y(x) points.

    This data is before internal or external transformations.
    It may include points on each end to connect the curve to
    the axis of rotation.
  */
  axom::Array<Point2DType> m_sorCurve;

  //! @brief Bounding box of points in m_sorCurve;
  BoundingBox2DType m_curveBb;

  //! @brief Maximum radius of the SOR.
  double m_maxRadius;

  //! @brief Minimum radius of the SOR.
  double m_minRadius;

  //!@brief The point corresponding to z=0 on the SOR axis.
  Point3DType m_sorOrigin;

  //!@brief SOR axis in 3D space, in the direction of increasing z.
  Vector3DType m_sorDirection;

  //!@brief Level of refinement for discretizing curved
  // analytical shapes and surfaces of revolutions.
  axom::IndexType m_levelOfRefinement = 0;

  //!@brief Array implementation of a += b.
  void accumulateData(axom::ArrayView<double> a,
                      axom::ArrayView<const double> b,
                      double scale,
                      axom::runtime_policy::Policy runtimePolicy);

  // Extract clipper info from GeometryClipperStrategy::m_info.
  void extractClipperInfo();

  void splitIntoMonotonicSections(axom::ArrayView<const Point2DType> pts,
                                  axom::Array<axom::Array<Point2DType>>& sections);

  void initializeFSorClippers();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SORCLIPPER_HPP
