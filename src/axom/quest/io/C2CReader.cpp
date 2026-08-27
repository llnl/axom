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

#include "c2c/C2C.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace axom::quest
{
namespace
{
constexpr c2c::LengthUnit toC2CLengthUnit(utilities::LengthUnit unit)
{
  switch(unit)
  {
  case utilities::LengthUnit::km:
    return c2c::LengthUnit::km;
  case utilities::LengthUnit::m:
    return c2c::LengthUnit::m;
  case utilities::LengthUnit::cm:
    return c2c::LengthUnit::cm;
  case utilities::LengthUnit::mm:
    return c2c::LengthUnit::mm;
  case utilities::LengthUnit::um:
    return c2c::LengthUnit::um;
  case utilities::LengthUnit::miles:
    return c2c::LengthUnit::miles;
  case utilities::LengthUnit::feet:
    return c2c::LengthUnit::ft;
  case utilities::LengthUnit::inches:
    return c2c::LengthUnit::in;
  case utilities::LengthUnit::mils:
    return c2c::LengthUnit::mils;
  case utilities::LengthUnit::dm:
  case utilities::LengthUnit::hm:
  case utilities::LengthUnit::dam:
  case utilities::LengthUnit::am:
  case utilities::LengthUnit::fm:
  case utilities::LengthUnit::pm:
  case utilities::LengthUnit::nm:
  case utilities::LengthUnit::angstrom:
  case utilities::LengthUnit::unspecified:
    throw std::invalid_argument("Length unit is not supported by c2c");
  }

  throw std::invalid_argument("Unknown length unit");
}
}  // namespace

void C2CReader::clear() { m_nurbsData.clear(); }

void C2CReader::setLengthUnit(utilities::LengthUnit lengthUnit)
{
  toC2CLengthUnit(lengthUnit);
  m_lengthUnit = lengthUnit;
}

bool C2CReader::hasValidExtension(const std::string& filename)
{
  return utilities::string::endsWith(filename, ".contour") ||
    utilities::string::endsWith(filename, ".assembly");
}

int C2CReader::read()
{
  SLIC_WARNING_ROOT_IF(m_fileName.empty(), "Missing a filename in C2CReader::read()");

  // Always clear prior results so callers never observe stale curves after a failed read
  this->clear();

  CurveArray inputCurves;
  const auto result = readInternal(m_fileName, inputCurves);
  const int ret = (result == ResultType::Success) ? 0 : 1;
  if(result == ResultType::Success)
  {
    m_nurbsData = std::move(inputCurves);
  }

  return ret;
}

C2CReader::ResultType C2CReader::readInternal(const std::string& filename, CurveArray& inputCurves)
{
  try
  {
    if(!hasValidExtension(filename))
    {
      SLIC_WARNING_ROOT(axom::fmt::format("{} is not a valid c2c file", filename));
      return ResultType::Failure;
    }

    if(utilities::string::endsWith(filename, ".contour"))
    {
      return readContour(filename, inputCurves);
    }
    else if(utilities::string::endsWith(filename, ".assembly"))
    {
      return readAssembly(filename, inputCurves);
    }
  }
  catch(const std::exception& e)
  {
    SLIC_WARNING_ROOT(axom::fmt::format("Failed to read c2c file '{}': {}", filename, e.what()));
  }
  catch(...)
  {
    SLIC_WARNING_ROOT(axom::fmt::format("Failed to read c2c file '{}'", filename));
  }

  return ResultType::Failure;
}

C2CReader::ResultType C2CReader::readAssembly(const std::string& filename, CurveArray& inputCurves)
{
  const c2c::Assembly assembly = c2c::parseAssembly(filename);
  std::string assemblyDir;
  utilities::filesystem::getDirName(assemblyDir, filename);

  SLIC_INFO_ROOT(fmt::format("Loading assembly with {} pieces", assembly.getNumEntries()));

  // Make an initial guess at the number of curves we may need.
  constexpr int contours_per_file_guess = 6;
  CurveArray assemblyCurves;
  assemblyCurves.reserve(assembly.getNumEntries() * contours_per_file_guess);

  ResultType ret = ResultType::Success;
  for(auto it = assembly.begin(); it != assembly.end(); it++)
  {
    if(ret == ResultType::Success)
    {
      const auto entryPath = utilities::filesystem::prefixRelativePath(it->getPath(), assemblyDir);
      ret = readInternal(entryPath, assemblyCurves);
    }
  }

  if(ret == ResultType::Success)
  {
    // Move the curves out to inputCurves.
    inputCurves.reserve(inputCurves.size() + assemblyCurves.size());
    for(auto& curve : assemblyCurves)
    {
      inputCurves.emplace_back(std::move(curve));
    }
  }

  return ret;
}

C2CReader::ResultType C2CReader::readContour(const std::string& filename, CurveArray& inputCurves)
{
  using PointType = primal::Point<double, 2>;

  c2c::Contour contour = c2c::parseContour(filename);
  const c2c::LengthUnit c2cLengthUnit = toC2CLengthUnit(m_lengthUnit);

  SLIC_INFO_ROOT(fmt::format("Loading contour with {} pieces", contour.getPieces().size()));

  inputCurves.reserve(inputCurves.size() + contour.getPieces().size());

  int piece_index = 0;
  for(auto* piece : contour.getPieces())
  {
    const auto nurbsData = c2c::toNurbs(*piece, c2cLengthUnit);

    // Load control points
    axom::Array<PointType> controlPoints;
    controlPoints.reserve(nurbsData.controlPoints.size());
    for(const auto& pt : nurbsData.controlPoints)
    {
      controlPoints.emplace_back(PointType {pt.getZ().getValue(), pt.getR().getValue()});
    }
    const int npts = static_cast<int>(controlPoints.size());
    if(npts <= 0)
    {
      continue;
    }

    // Load and check knot vector; check degree first then knots
    const auto nkts = static_cast<axom::IndexType>(nurbsData.knots.size());
    const int degree = static_cast<int>(nkts - npts - 1);
    if(degree < 0)
    {
      SLIC_WARNING_ROOT(
        fmt::format("Invalid contour file '{}': computed negative NURBS degree for piece "
                    "{} (npts={}, nkts={})",
                    filename,
                    piece_index,
                    npts,
                    nkts));
      return ResultType::Failure;
    }

    if(npts <= degree)
    {
      SLIC_WARNING_ROOT(
        fmt::format("Invalid contour file '{}': piece {} has too few control points for degree "
                    "(degree={}, npts={}, nkts={})",
                    filename,
                    piece_index,
                    degree,
                    npts,
                    nkts));
      return ResultType::Failure;
    }

    const axom::ArrayView<const double> knots_view(nurbsData.knots.data(), nkts);
    primal::KnotVector<double> knotvec(knots_view,
                                       degree,
                                       primal::KnotVector<double>::SkipValidityChecks {});

    if(!knotvec.isValid(true))
    {
      SLIC_WARNING_ROOT(
        fmt::format("Invalid contour file '{}': piece {} converted to an invalid NURBS knot vector "
                    "(degree={}).",
                    filename,
                    piece_index,
                    degree));
      return ResultType::Failure;
    }

    // Load and check weights; count must either be 0 or match control points
    axom::Array<double> weights;
    if(!nurbsData.weights.empty())
    {
      if(static_cast<axom::IndexType>(nurbsData.weights.size()) != controlPoints.size())
      {
        SLIC_WARNING_ROOT(
          fmt::format("Invalid contour file '{}': piece {} has {} weights for {} control points",
                      filename,
                      piece_index,
                      nurbsData.weights.size(),
                      npts));
        return ResultType::Failure;
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
      inputCurves.emplace_back(controlPoints, knotvec);
    }
    else
    {
      inputCurves.emplace_back(controlPoints, weights, knotvec);
    }

    ++piece_index;
  }

  return ResultType::Success;
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

  SLIC_INFO_ROOT(sstr.str());
}

}  // namespace axom::quest
