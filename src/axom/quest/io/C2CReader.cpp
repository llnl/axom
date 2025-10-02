// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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

  for(auto* piece : contour.getPieces())
  {
    const auto nurbsData = c2c::toNurbs(*piece, m_lengthUnit);

    axom::Array<PointType> controlPoints;
    for(const auto& pt : nurbsData.controlPoints)
    {
      controlPoints.emplace_back(PointType {pt.getZ().getValue(), pt.getR().getValue()});
    }

    m_nurbsData.emplace_back(controlPoints.data(),
                             nurbsData.weights.data(),
                             controlPoints.size(),
                             nurbsData.knots.data(),
                             nurbsData.knots.size());
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
