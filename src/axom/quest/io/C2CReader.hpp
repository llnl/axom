// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_C2CREADER_HPP_
#define QUEST_C2CREADER_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_C2C
  #error C2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "c2c/C2C.hpp"

#include <string>
#include <vector>

namespace axom
{
namespace quest
{
/*
 * \class C2CReader
 *
 * \brief A class to help with reading C2C contour files and converting them to 1D meshes
 *
 * We treat all contours as NURBS curves and sample them to create segment meshes.
 */
class C2CReader
{
public:
  using NURBSCurve = axom::primal::NURBSCurve<double, 2>;

public:
  C2CReader() = default;

  virtual ~C2CReader() = default;

  /// Sets the name of the contour file to load. Must be called before \a read()
  void setFileName(const std::string &fileName) { m_fileName = fileName; }

  /// Sets the length unit. All lengths will be converted to this unit when reading the mesh
  void setLengthUnit(c2c::LengthUnit lengthUnit) { m_lengthUnit = lengthUnit; }

  /// Sets the threshold for welding vertices of adjacent Pieces of curves
  void setVertexWeldingThreshold(double thresh) { m_vertexWeldThreshold = thresh; }

  /// Clears data associated with this reader
  void clear();

  /*!
   * \brief Read the contour file provided by \a setFileName()
   * 
   * \return 0 for a successful read; non-zero otherwise
   */
  virtual int read();

  /// \brief Utility function to log details about the read in file
  virtual void log();

  /*!
   * \brief Projects high-order NURBS contours onto a linear mesh using \a segmentsPerPiece 
   * linear segments per knot span of the contour
   * 
   * Knot spans are the sub-intervals within a spline
   */
  void getLinearMeshUniform(mint::UnstructuredMesh<mint::SINGLE_SHAPE> *mesh,
                            int segmentsPerKnotSpan);

  /*!
   * \brief Projects high-order NURBS contours onto a linear mesh using \a percentError 
   *        to decide when to stop refinement.
   * 
   * \param[in] mesh The mesh object that will contain the linearized line segments.
   * \param[in] percentError A percent of error that is acceptable to stop refinement.
   */
  void getLinearMeshNonUniform(mint::UnstructuredMesh<mint::SINGLE_SHAPE> *mesh, double percentError);

  /*!
   * \brief Compute the revolved volume of the shape using quadrature.
   *
   * \param transform A 4x4 matrix transform to apply to the shape before computing
   *                  the revolved volume.
   *
   * \note We compute revolved volume on the actual shapes so we can get a
   *       real revolved volume computed using the curve functions rather than
   *       relying on a linearized curve. The revolved volume is the volume
   *       enclosed by the surface of revolution when the shape is revolved
   *       about axis of revolution.
   *
   * \return The revolved volume.
   */
  double getRevolvedVolume(const numerics::Matrix<double> &transform) const;

protected:
  int readContour();

protected:
  std::string m_fileName;
  c2c::LengthUnit m_lengthUnit {c2c::LengthUnit::cm};
  double m_vertexWeldThreshold {1E-9};

  std::vector<NURBSCurve> m_nurbsData;
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_C2CREADER_HPP_
