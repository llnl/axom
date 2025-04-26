// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_SORCLIPPER_HPP
#define AXOM_QUEST_SORCLIPPER_HPP

#include "axom/klee/Geometry.hpp"
#include "axom/quest/GeometryClipperStrategy.hpp"

namespace axom
{
namespace quest
{

/*!
  @brief GeometryClipper specialized for 3D surface-of-revolution geometries.
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

  bool labelInOut(quest::ShapeeMesh& shappeMesh, axom::Array<char>& label) override;

  bool getShapeAsOcts(quest::ShapeeMesh& shappeMesh, axom::Array<OctahedronType>& octs) override;

#if !defined(__CUDACC__)
private:
#endif
  std::string m_name;

  //! @brief The discrete r(z) function, as an Nx2 array, if used.
  axom::Array<double, 2> m_discreteFcn;

  //!@brief The point corresponding to z=0 on the SOR axis.
  Point3D m_sorBase;

  //!@brief SOR axis in 3D space, in the direction of increasing z.
  Vector3DType m_sorDirection;

  //!@brief Level of refinement for discretizing curved
  // analytical shapes and surfaces of revolutions.
  axom::IndexType m_levelOfRefinement = 0;

  //! @brief Number of octahedral cells in the discrete SOR.
  axom::IndexType m_cellCount;

  template <typename ExecSpace>
  void labelInOutImpl(quest::ShapeeMesh& shapeeMesh, axom::Array<char>& label);

  // Return a 3x3 matrix that rotates coordinates from the x-axis to the given direction.
  numerics::Matrix<double> sorAxisRotMatrix(const Vector3DType& dir);

  // Extract clipper info from GeometryClipperStrategy::m_info.
  void extractClipperInfo();
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SORCLIPPER_HPP
