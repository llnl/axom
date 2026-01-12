// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_MFEMREADER_HPP_
#define QUEST_MFEMREADER_HPP_

#include "axom/config.hpp"

#if !defined(AXOM_USE_MFEM) || !defined(AXOM_USE_SIDRE)
  #error MFEMReader should only be included when Axom is configured with MFEM, SIDRE (and MFEM_SIDRE_DATACOLLECTION)
#endif

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"

#include <string>
#include <vector>

namespace axom
{
namespace quest
{
/*
 * \class MFEMReader
 *
 * \brief A class to help with reading MFEM files that contain contours.
 */
class MFEMReader
{
public:
  using NURBSCurve = axom::primal::NURBSCurve<double, 2>;
  using CurveArray = axom::Array<NURBSCurve>;

  using CurvedPolygon = axom::primal::CurvedPolygon<NURBSCurve>;
  using CurvedPolygonArray = axom::Array<CurvedPolygon>;

  static constexpr int READ_FAILED = 1;
  static constexpr int READ_SUCCESS = 0;

public:
  /// Sets the name of the contour file to load. Must be called before \a read()
  void setFileName(const std::string &fileName) { m_fileName = fileName; }

  /*!
   * \brief Read the contour file provided by \a setFileName()
   *
   * \param[out] curves The curve array that will contain curves read from the MFEM file.
   *
   * \return READ_SUCCESS for a successful read; READ_FAILED (non-zero) otherwise
   */
  int read(CurveArray &curves);

  /*!
   * \brief Read the contour file provided by \a setFileName()
   *
   * \param[out] curves The curve array that will contain curves read from the MFEM file.
   * \param[out] attributes The MFEM attribute value associated with each curve.
   *
   * \note The i-th entry in \a attributes contains the MFEM attribute for the i-th curve.
   *
   * \return READ_SUCCESS for a successful read; READ_FAILED (non-zero) otherwise
   */
  int read(CurveArray &curves, axom::Array<int> &attributes);

  /*!
   * \brief Read the contour file provided by \a setFileName()
   * 
   * \param[out] curvedPolygons The curved polygon array that will contain curved polygons created from reading
   *                            the MFEM file.
   * \note The returned polygons are stored in a 0-based, contiguous array. When MFEM attributes are not
   *       contiguous, the polygon index will not match the MFEM attribute value.
   *
   * \return READ_SUCCESS for a successful read; READ_FAILED (non-zero) otherwise
   */
  int read(CurvedPolygonArray &curvedPolygons);

  /*!
   * \brief Read the contour file provided by \a setFileName()
   *
   * \param[out] curvedPolygons The curved polygon array that will contain curved polygons created from reading
   *                            the MFEM file.
   * \param[out] attributes The MFEM attribute value associated with each curved polygon.
   *
   * \note The returned polygons are stored in a 0-based, contiguous array. The i-th entry in \a attributes
   *       contains the MFEM attribute for the i-th polygon.
   *
   * \return READ_SUCCESS for a successful read; READ_FAILED (non-zero) otherwise
   */
  int read(CurvedPolygonArray &curvedPolygons, axom::Array<int> &attributes);

protected:
  std::string m_fileName;
};

}  // namespace quest
}  // namespace axom

#endif
