// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_FSORCLIPPER_HPP
#define AXOM_QUEST_FSORCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief Geometry clipping operations for 3D
  surface-of-revolution geometries.

  This implementation requires the SOR curve to be a function
  (has monotonically changing axial coordinates).  For SOR curves
  that are not functions, use SorClipper.

  The SOR specification may include rotation and translation
  internally, in addition to any external transformation.
*/
class FSorClipper : public GeometryClipperStrategy
{
public:
  /*!
    @brief Constructor.

    @param [in] kGeom Describes the shape to place
      into the mesh.
    @param [in] name To override the default strategy name
  */
  FSorClipper(const klee::Geometry& kGeom, const std::string& name = "");

  /*!
    @brief Construct from geometric specifications.
  */
  FSorClipper(const klee::Geometry& kGeom,
              const std::string& name,
              axom::ArrayView<const Point2DType> sorCurve,
              const Point3DType& sorOrigin,
              const Vector3DType& sorDirection,
              axom::IndexType levelOfRefinement);

  virtual ~FSorClipper() = default;

  const std::string& name() const override { return m_name; }

  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getGeometryAsOcts(quest::ShapeeMesh& shappeMesh,
                         axom::Array<OctahedronType>& octs) override;

  axom::ArrayView<const Point2DType> getSorCurve() const
    { return m_sorCurve.view(); }

  //@{
  //! @name Utilities shared with SorClipper for handling SOR.
  /*!
    @brief Find division points between curve sections where z (x)
    changes directions.

    @param sorCurve Set of at least 2 2D points describing a curve
      in r-z space (in host array).

    @return Indices of switchbacks, plus the first and last indices.
  */
  static axom::Array<axom::IndexType> findZSwitchbacks(
    axom::ArrayView<const Point2DType> pts);

  /*
    @brief Combine consecutive radial segments of the curve into a
    single segment.

    This step is necessary because some other steps assume there are
    no consecutive radial segments.
  */
  static void combineRadialSegments(axom::Array<Point2DType>& sorCurve);
  //@}

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

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

  /*!
    @brief Boxes (in rz space) on the curve.

    The curve lies completely in these boxes and includes the planes
    of the base and top.  Points in these boxes are require more
    computation to determine their signed distance.
  */
  axom::Array<BoundingBox2DType> m_bbOn;

  /*!
    @brief Boxes (in rz space) completely under the curve.

    These boxes lie completely under the curve.
  */
  axom::Array<BoundingBox2DType> m_bbUnder;

  //!@brief Internal and external transforms (includes m_sorDirection and m_sorOrigin).
  axom::primal::CoordinateTransformer<double> m_transformer;

  /*!
    @brief Inverse of m_transformer.

    Axom supports vector scaling.  @see axom::klee::Scale.  This means
    a SOR may be transformed into a shape that we cannot represent.
    Therefore, we don't transform the shape until after it's discretized.
    When needed, we will inverse-transform the mesh.
  */
  axom::primal::CoordinateTransformer<double> m_inverseTransformer;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  // Extract clipper info from GeometryClipperStrategy::m_info.
  void extractClipperInfo();

  void clusterSorFunction();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_FSORCLIPPER_HPP
