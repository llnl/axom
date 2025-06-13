// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SORCLIPPER_HPP
#define AXOM_QUEST_SORCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"
#include "axom/primal/geometry/CoordinateTransformer.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief GeometryClipper specialized for 3D surface-of-revolution geometries.

  The SOR internally support rotation and translation, in addition to any
  external transformation.
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

  /*!
    TODO: This method can be implemented but has not been.
  */
  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getGeometryAsOcts(quest::ShapeeMesh& shappeMesh, axom::Array<OctahedronType>& octs) override;
  bool getCurveWithAxisPoints(axom::Array<Point2DType>& curveWithAxisPoints);

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  /*!
    @brief The discrete r(z) curve, as an Nx2 array, if used.

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

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  // Extract clipper info from GeometryClipperStrategy::m_info.
  void extractClipperInfo();

  // Compute blocking of areas on and under m_sorCurve.
  void computeRoughBlockings();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SORCLIPPER_HPP
