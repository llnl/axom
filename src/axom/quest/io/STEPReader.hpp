// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_STEPREADER_HPP_
#define QUEST_STEPREADER_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_OPENCASCADE
  #error STEPReader should only be included when Axom is configured with opencascade
#endif

#include "axom/core.hpp"
#include "axom/primal.hpp"

#include <string>
#include <memory>
#include <map>

class TopoDS_Shape;

namespace axom
{
namespace quest
{

namespace internal
{
class StepFileProcessor;

/// Struct to hold data associated with each surface patch of the mesh
struct PatchData
{
  int patchIndex {-1};
  bool wasOriginallyPeriodic_u {false};
  bool wasOriginallyPeriodic_v {false};
  axom::primal::BoundingBox<double, 2> parametricBBox;
  axom::primal::BoundingBox<double, 3> physicalBBox;
  axom::Array<bool> trimmingCurves_originallyPeriodic;
};

using PatchDataMap = std::map<int, PatchData>;

}  // namespace internal

/*
 * \class STEPReader
 *
 * \brief A class to help with reading a STEP file containing a parametric BRep (Boundary Representation)
 * consisting of trimmed NURBS patches.
 */
class STEPReader
{
public:
  using NURBSPatch = axom::primal::NURBSPatch<double, 3>;
  using PatchArray = axom::Array<NURBSPatch>;

  STEPReader() = default;
  virtual ~STEPReader();

public:
  /// Sets the name of the step file to load. Must be called before \a read()
  void setFileName(const std::string& fileName) { m_fileName = fileName; }

  void setVerbosity(bool verbosity) { m_verbosity = verbosity; }

  /*!
   * \brief Read the contour file provided by \a setFileName()
   *
   * \param[out] patches An array of NURBS patches that will contain the trimmed patches read from the STEP file
   *
   * \return 0 for a successful read; non-zero otherwise
   */
  int read();

  std::string getFileUnits() const;
  const internal::PatchDataMap& getPatchDataMap() const;
  const TopoDS_Shape& getShape() const;

  PatchArray& getPatchArray() { return m_patches; }
  const PatchArray& getPatchArray() const { return m_patches; }

  /// Logs some information about the loaded BRep
  void printBRepStats() const;

protected:
  // open cascade does not appear to offer a direct way to get the number of patches
  int numPatchesInFile() const;

protected:
  std::string m_fileName;
  bool m_verbosity {false};
  internal::StepFileProcessor* m_stepProcessor {nullptr};
  PatchArray m_patches;
};

}  // namespace quest
}  // namespace axom

#endif
