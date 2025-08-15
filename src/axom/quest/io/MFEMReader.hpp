// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_MFEMREADER_HPP_
#define QUEST_MFEMREADER_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_MFEM
  #error MFEMReader should only be included when Axom is configured with MFEM
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
 * \brief A class to help with reading C2C contour files.
 *
 * We treat all contours as NURBS curves.
 */
class MFEMReader
{
public:
  using NURBSCurve = axom::primal::NURBSCurve<double, 2>;
  using CurveArray = axom::Array<NURBSCurve>;

  using CurvedPolygon = axom::primal::CurvedPolygon<double, 2>;
  using CurvedPolygonArray = axom::Array<CurvedPolygon>;
public:
  MFEMReader() = default;

  ~MFEMReader() = default;

  /// Sets the name of the contour file to load. Must be called before \a read()
  void setFileName(const std::string &fileName) { m_fileName = fileName; }

  /*!
   * \brief Read the contour file provided by \a setFileName()
   * 
   * \return 0 for a successful read; non-zero otherwise
   */
  int read(CurveArray &curves);

  /*!
   * \brief Read the contour file provided by \a setFileName()
   * 
   * \return 0 for a successful read; non-zero otherwise
   */
  int read(CurvedPolygonArray &curves);

protected:
  std::string m_fileName;
};

}  // namespace quest
}  // namespace axom

#endif
