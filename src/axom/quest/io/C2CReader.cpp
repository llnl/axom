// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/C2CReader.hpp"

#ifndef AXOM_USE_C2C
  #error C2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/fmt.hpp"

#include <fstream>
#include <string>

namespace axom
{
namespace quest
{

void C2CReader::clear() { m_nurbsData.clear(); }

int C2CReader::read()
{
  SLIC_WARNING_IF(m_fileName.empty(), "Missing a filename in C2CReader::read()");

  using axom::utilities::string::endsWith;

  int ret = 1;

  if(endsWith(m_fileName, ".contour"))
  {
    ret = readContour();
  }
  else if(endsWith(m_fileName, ".assembly"))
  {
    SLIC_WARNING("Input was an assembly! This is not currently supported");
  }
  else
  {
    SLIC_WARNING("Not a valid c2c file");
  }

  return ret;
}

int C2CReader::readContour()
{
  using PointType = primal::Point<double, 2>;

  c2c::Contour contour = c2c::parseContour(m_fileName);

  SLIC_INFO(fmt::format("Loading contour with {} pieces", contour.getPieces().size()));

  m_nurbsData.clear();
  m_nurbsData.reserve(contour.getPieces().size());

  int piece_index = 0;
  for(auto* piece : contour.getPieces())
  {
    const auto nurbsData = c2c::toNurbs(*piece, m_lengthUnit);

    // Load control points
    axom::Array<PointType> controlPoints;
    controlPoints.reserve(nurbsData.controlPoints.size());
    for(const auto& pt : nurbsData.controlPoints)
    {
      controlPoints.emplace_back(PointType {pt.getZ().getValue(), pt.getR().getValue()});
    }
    const int npts = static_cast<int>(controlPoints.size());

    // Load and check knot vector; check degree first then knots
    const auto nkts = static_cast<axom::IndexType>(nurbsData.knots.size());
    const int degree = static_cast<int>(nkts - npts - 1);
    if(degree < 0)
    {
      SLIC_WARNING(
        fmt::format("Invalid contour file '{}': computed negative NURBS degree for piece "
                    "{} (npts={}, nkts={})",
                    m_fileName,
                    piece_index,
                    npts,
                    nkts));
      return 1;
    }

    const axom::ArrayView<const double> knots_view(nurbsData.knots.data(), nkts);
    primal::KnotVector<double> knotvec(knots_view,
                                       degree,
                                       primal::KnotVector<double>::SkipValidityChecks {});

    if(!knotvec.isValid(true))
    {
      SLIC_WARNING(
        fmt::format("Invalid contour file '{}': piece {} converted to an invalid NURBS knot vector "
                    "(degree={}).",
                    m_fileName,
                    piece_index,
                    degree));
      return 1;
    }

    // Load and check weights; count must either be 0 or match control points
    axom::Array<double> weights;
    if(!nurbsData.weights.empty())
    {
      if(static_cast<axom::IndexType>(nurbsData.weights.size()) != controlPoints.size())
      {
        SLIC_WARNING(
          fmt::format("Invalid contour file '{}': piece {} has {} weights for {} control points",
                      m_fileName,
                      piece_index,
                      nurbsData.weights.size(),
                      npts));
        return 1;
      }

      // Check if weights are non-trivial (present and not all equal to 1)
      bool has_non_trivial_weights = false;
      for(const double& wt : nurbsData.weights)
      {
        if(wt != 1.0)
        {
          has_non_trivial_weights = true;
          break;
        }
      }

      // Only copy weights if they are non-trivial
      if(has_non_trivial_weights)
      {
        weights.reserve(nurbsData.weights.size());
        for(const double& wt : nurbsData.weights)
        {
          weights.push_back(wt);
        }
      }
    }

    // Construct NURBSCurve using Array constructors to avoid use-after-free
    if(weights.empty())
    {
      m_nurbsData.emplace_back(controlPoints, knotvec);
    }
    else
    {
      m_nurbsData.emplace_back(controlPoints, weights, knotvec);
    }

    ++piece_index;
  }

  return 0;
}

void C2CReader::log()
{
  std::stringstream sstr;

  sstr << fmt::format("The contour has {} pieces\n", m_nurbsData.size());

  int index = 0;
  for(const auto& curve : m_nurbsData)
  {
    sstr << fmt::format("\tCurve {}: {}\n", index, curve);
    ++index;
  }

  SLIC_INFO(sstr.str());
}

}  // end namespace quest
}  // end namespace axom
