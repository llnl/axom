// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_STEPREADER_HPP_
#define QUEST_STEPREADER_HPP_

#include "axom/config.hpp"
#include "axom/mint.hpp"

#ifndef AXOM_USE_OPENCASCADE
  #error STEPReader should only be included when Axom is configured with opencascade
#endif

#include "axom/core.hpp"
#include "axom/primal.hpp"

#include <string>
#include <memory>
#include <map>

namespace axom
{
namespace quest
{
namespace internal
{
class StepFileProcessor;
}

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

  using NURBSCurve = axom::primal::NURBSCurve<double, 2>;

  STEPReader() = default;
  virtual ~STEPReader();

public:
  /// Sets the name of the step file to load. Must be called before \a read()
  void setFileName(const std::string& fileName) { m_fileName = fileName; }

  void setVerbosity(bool verbosity) { m_verbosity = verbosity; }

  /*!
   * \brief Read the contour file provided by \a setFileName()
   *
   * \param[in] validate Adds validation tests on the model, when true
   * \return 0 for a successful read; non-zero otherwise
   */
  virtual int read(bool validate);

  std::string getFileUnits() const;

  PatchArray& getPatchArray() { return m_patches; }
  const PatchArray& getPatchArray() const { return m_patches; }

  /// Returns some information about the loaded BRep
  std::string getBRepStats() const;

  /*!
   * \brief Generates a triangulated representation of the STEP file as a Mint unstructured triangle mesh.
   *
   * \param[inout] mesh Pointer to a Mint unstructured mesh that will be populated
   *            with triangular elements approximating the STEP geometry.
   * \param[in] linear_deflection Maximum allowed deviation between the
   *            original geometry and the triangulated approximation.
   * \param[in] angular_deflection Maximum allowed angular deviation (in radians)
   *            between normals of adjacent triangles.
   * \param[in] is_relative When false (default), linear deflection is in mesh units. When true,
                linear deflection is relative to the local edge length of the triangles.
   * \param[in] trimmed If true (default), the triangulation respects trimming curves.
   *            otherwise, we triangulate the untrimmed patches. The latter is mostly to aid 
   *            in understanding the model's patches and is not generally useful.
   */
  virtual int getTriangleMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* mesh,
                              double linear_deflection = 0.1,
                              double angular_deflection = 0.5,
                              bool is_relative = false,
                              bool trimmed = true);

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
