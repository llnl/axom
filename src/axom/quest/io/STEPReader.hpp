// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/config.hpp"
#include "axom/mint.hpp"

#ifndef AXOM_USE_OPENCASCADE
  #error STEPReader should only be included when Axom is configured with Open Cascade
#endif

#include "axom/core.hpp"
#include "axom/primal.hpp"

#include <string>
#include <memory>

namespace axom::quest
{
namespace internal
{
class StepFileProcessor;
}

/*!
 * \class STEPReader
 *
 * \brief A class to help with reading a STEP file containing a parametric 
 *  BRep (Boundary Representation) consisting of trimmed NURBS patches.
 */
class STEPReader
{
public:
  using NURBSPatch = axom::primal::NURBSPatch<double, 3>;
  using PatchArray = axom::Array<NURBSPatch>;

  using NURBSCurve = axom::primal::NURBSCurve<double, 2>;
  using IndexArray = axom::Array<int>;

  STEPReader();
  virtual ~STEPReader();

public:
  /// Sets the name of the step file to load. Must be called before \a read()
  void setFileName(const std::string& fileName) { m_fileName = fileName; }

  void setVerbosity(bool verbosity) { m_verbosity = verbosity; }

  /*!
   * \brief Read the STEP file provided by \a setFileName()
   *
   * \param[in] validate Logs validation checks, sampled UV containment,
   *  and trimming curve gaps, when true.
   *
   * \note UV checks use an absolute tolerance of 1e-7 on the exported patches.
   *
   * Validation reports problems without changing the geometry or read status.
   * A successful read can still contain out-of-domain trims or open wires.
   * Each read replaces the previous model; a failed read leaves no patches.
   * \return 0 for a successful read; non-zero otherwise
   */
  [[nodiscard]] virtual int read(bool validate);

  /// Returns the file units, or an empty string if no model is loaded on this rank.
  std::string getFileUnits() const;

  PatchArray& getPatchArray() { return m_patches; }
  const PatchArray& getPatchArray() const { return m_patches; }

  /// Get the number of patches in the read file
  int numPatches() const { return m_patches.size(); }

  /*!
   * \brief Returns a 0-based patch id per entry in \a getPatchArray()
   *
   * This id corresponds to the face index in the input STEP file enumeration.
   * If consumers skip patches, these ids will be non-contiguous.
   */
  axom::ArrayView<const int> getPatchIds() const { return m_patchIds.view(); }

  /*!
   * \brief Returns the 0-based wire index for each extracted trimming curve
   *
   * The i-th entry corresponds to `getPatchArray()[patchArrayIndex].getTrimmingCurves()[i]`
   * and stores the 0-based wire index from the input face enumeration.
   */
  axom::ArrayView<const int> getTrimmingCurveWireIds(int patchArrayIndex) const;

  /// \brief Returns whether the input STEP surface was originally periodic in u for this patch
  bool patchWasOriginallyPeriodic_u(int patchArrayIndex) const;

  /// \brief Returns whether the input STEP surface was originally periodic in v for this patch
  bool patchWasOriginallyPeriodic_v(int patchArrayIndex) const;

  /// Returns some information about the loaded BRep
  std::string getBRepStats() const;

  /// Returns an AABB for the loaded BRep as evaluated by Open Cascade
  axom::primal::BoundingBox<double, 3> getBRepBoundingBox(bool useTriangulation = false) const;

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
  [[nodiscard]] virtual int getTriangleMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* mesh,
                                            double linear_deflection = 0.1,
                                            double angular_deflection = 0.5,
                                            bool is_relative = false,
                                            bool trimmed = true);

protected:
  /// Releases the loaded model and exported patches, retaining reader settings.
  void clear();

  // Open Cascade does not appear to offer a direct way to get the number of patches
  int numPatchesInFile() const;

protected:
  std::string m_fileName;
  bool m_verbosity {false};
  std::unique_ptr<internal::StepFileProcessor> m_stepProcessor;
  PatchArray m_patches;
  IndexArray m_patchIds;
  axom::Array<IndexArray> m_trimmingCurveWireIds;
  IndexArray m_patchOriginallyPeriodic_u;
  IndexArray m_patchOriginallyPeriodic_v;

private:
  DISABLE_COPY_AND_ASSIGNMENT(STEPReader);
  DISABLE_MOVE_AND_ASSIGNMENT(STEPReader);
};

}  // namespace axom::quest
