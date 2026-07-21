// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"

#ifndef AXOM_USE_C2C
  #error C2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/utilities/Units.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"

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

  enum class ResultType
  {
    Success,
    Failure
  };
public:
  C2CReader() = default;

  virtual ~C2CReader() = default;

  /// Sets the name of the contour file to load. Must be called before \a read()
  void setFileName(const std::string &fileName) { m_fileName = fileName; }

  /// Sets the length unit. All lengths will be converted to this unit when reading the mesh
  void setLengthUnit(utilities::LengthUnit lengthUnit);

  /// Clears data associated with this reader
  void clear();

  /// Returns true if the file has a recognized c2c extension.
  static bool hasValidExtension(const std::string &filename);

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
  /*!
   * \brief Internal helper for reading files.
   *
   * \param filename The name of the file to read.
   * \param[inout] inputCurves The array of curves to append.
   *
   * \return Success on success, Failure otherwise.
   */
  ResultType readInternal(const std::string &filename, CurveArray &inputCurves);

  /*!
   * \brief Internal helper for reading a contour file.
   *
   * \param filename The name of the file to read.
   * \param[inout] inputCurves The array of curves to append.
   *
   * \return Success on success, Failure otherwise.
   */
  ResultType readContour(const std::string &filename, CurveArray &inputCurves);

  /*!
   * \brief Internal helper for reading an assembly file.
   *
   * \param filename The name of the file to read.
   * \param[inout] inputCurves The array of curves to append.
   *
   * \return Success on success, Failure otherwise.
   */
  ResultType readAssembly(const std::string &filename, CurveArray &inputCurves);

protected:
  std::string m_fileName;
  utilities::LengthUnit m_lengthUnit {utilities::LengthUnit::cm};
  CurveArray m_nurbsData;
};

}  // namespace quest
}  // namespace axom
