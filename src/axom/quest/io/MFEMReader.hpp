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

public:
  /// Sets the name of the contour file to load. Must be called before \a read()
  void setFileName(const std::string &fileName) { m_fileName = fileName; }

  /*!
   * \brief Read the contour file provided by \a setFileName()
   *
   * \param[out] curves The curve array that will contain curves read from the MFEM file.
   *
   * \return 0 for a successful read; non-zero otherwise
   */
  int read(CurveArray &curves);

  /*!
   * \brief Read the contour file provided by \a setFileName()
   * 
   * \param[out] curvedPolygons The curved polygon array that will contain curved polygons created from reading
   *                            the MFEM file.
   *
   * \return 0 for a successful read; non-zero otherwise
   */
  int read(CurvedPolygonArray &curvedPolygons);

protected:
  std::string m_fileName;
};

}  // namespace quest
}  // namespace axom

#endif
