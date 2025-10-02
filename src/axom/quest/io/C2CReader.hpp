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

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
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
 * \brief A class to help with reading C2C contour files.
 *
 * We treat all contours as NURBS curves.
 */
class C2CReader
{
public:
  using NURBSCurve = axom::primal::NURBSCurve<double, 2>;
  using CurveArray = axom::Array<NURBSCurve>;
  using CurveArrayView = axom::ArrayView<NURBSCurve>;

public:
  C2CReader() = default;

  virtual ~C2CReader() = default;

  /// Sets the name of the contour file to load. Must be called before \a read()
  void setFileName(const std::string &fileName) { m_fileName = fileName; }

  /// Sets the length unit. All lengths will be converted to this unit when reading the mesh
  void setLengthUnit(c2c::LengthUnit lengthUnit) { m_lengthUnit = lengthUnit; }

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
   * \brief Get a view that contains the curves.
   * 
   * \return A view that contains the curves.
   */
  CurveArrayView getCurvesView() { return m_nurbsData.view(); }

protected:
  int readContour();

protected:
  std::string m_fileName;
  c2c::LengthUnit m_lengthUnit {c2c::LengthUnit::cm};
  CurveArray m_nurbsData;
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_C2CREADER_HPP_
