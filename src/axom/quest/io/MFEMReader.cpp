// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/MFEMReader.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#ifndef AXOM_USE_MFEM
  #error MFEMReader should only be included when Axom is configured with MFEM
#endif

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include <mfem.hpp>

#include <map>
#include <string>
#include <memory>

namespace axom
{
namespace quest
{
namespace internal
{
/*!
 * \brief Read the MFEM file and build the desired type of geometry from it using a supplied function.
 *        The MFEM files must contain only 1D curves in 2D space.
 *
 * \tparam BuildGeometry A function/lambda that will be used to build geometry from the MFEM mesh.
 *
 * \param fileName The name of the file to read.
 * \param build The function/lambda that generates the geometry using the map of zones to curves.
 *              The function must take 2 parameters, an MFEM mesh pointer, and a reference to
 *              std::map<int, axom::Array<int>>. The latter map contains contourId:zoneIdList mapping, which
 *              can be used to group related MFEM zones/contours as edges in a shape.
 *
 * \return 0 on success; non-zero on failure.
 */
template <typename BuildGeometry>
int read_mfem(const std::string &fileName, BuildGeometry &&build)
{
  constexpr int READ_FAILED = 1;
  constexpr int READ_SUCCESS = 0;

  int retval = READ_FAILED;

  // Load the MFEM file
  std::unique_ptr<mfem::Mesh> mesh;
  try
  {
    mesh = std::make_unique<mfem::Mesh>(fileName, 1, 1, true);

    // This code is only supporting 1D meshes in 2D space.
    if(mesh->Dimension() == 1 && mesh->SpaceDimension() == 2)
    {
      // Examine the mesh attributes and group all of the related zones that are
      // edges of the same contour.
      std::map<int, axom::Array<int>> contourZones;
      for(int zoneId = 0; zoneId < mesh->GetNE(); zoneId++)
      {
        // Get element attribute and make it zero-origin.
        const int contourId = mesh->GetAttribute(zoneId) - 1;
        contourZones[contourId].push_back(zoneId);
      }

      // Use the map to build the geometry.
      build(mesh.get(), contourZones);

      retval = READ_SUCCESS;
    }
    else
    {
      SLIC_WARNING(
        axom::fmt::format("Mesh must have dimension 1 and spatial dimension 2. The supplied mesh "
                          "is dimension {} with spatial dimension {}.",
                          mesh->Dimension(),
                          mesh->SpaceDimension()));
      retval = READ_FAILED;
    }
  }
  catch(std::exception &e)
  {
    retval = READ_FAILED;
  }
  return retval;
}

}  // end namespace internal

int MFEMReader::read(CurveArray &curves)
{
  SLIC_WARNING_IF(m_fileName.empty(), "Missing a filename in MFEMReader::read()");

  curves.clear();
  return internal::read_mfem(
    m_fileName,
    [&](mfem::Mesh *mesh, const std::map<int, axom::Array<int>> &contourZones) {
      // Build NURBSCurves from the MFEM zones.
      for(auto &[contourId, zoneIds] : contourZones)
      {
        for(int zoneId : zoneIds)
        {
          curves.push_back(axom::quest::internal::segment_to_nurbs(mesh, zoneId));
        }
      }
    });
}

int MFEMReader::read(CurvedPolygonArray &curvedPolygons)
{
  SLIC_WARNING_IF(m_fileName.empty(), "Missing a filename in MFEMReader::read()");

  return internal::read_mfem(
    m_fileName,
    [&](mfem::Mesh *mesh, const std::map<int, axom::Array<int>> &contourZones) {
      // Build CurvedPolygons from the MFEM zones.

      // Resize the array.
      curvedPolygons.clear();
      curvedPolygons.resize(contourZones.size());

      for(auto &[contourId, zoneIds] : contourZones)
      {
        auto &poly = curvedPolygons[contourId];
        for(int zoneId : zoneIds)
        {
          auto curve = axom::quest::internal::segment_to_nurbs(mesh, zoneId);
          poly.addEdge(std::move(curve));
        }
      }
    });
}

}  // end namespace quest
}  // end namespace axom
