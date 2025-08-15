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

namespace axom
{
namespace quest
{
namespace internal
{
/*!
 * \brief Read the MFEM file and build the desired type of geometry from it using a supplied function.
 *
 * \param func The function/lambda that generates the geometry using the map of zones to curves.
 *
 * \return 0 on success; non-zero on failure.
 */
template <typename Func>
int read_mfem(const std::string &fileName, Func &&func)
{
  constexpr int READ_FAILED = 1;
  constexpr int READ_SUCCESS = 0;

  int retval = READ_FAILED;

  // Load the MFEM file
  mfem::Mesh* mesh = nullptr;
  try
  {
    mesh = new mfem::Mesh(fileName, 1, 1, true);

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
      func(mesh, contourZones);

      retval = READ_SUCCESS;
    }
    else
    {
      SLIC_INFO("Mesh must have dimension 1 and spatial dimension 2.");
    }

    delete mesh;
  }
  catch(std::exception& e)
  {
    delete mesh;
  }
  return retval;
}

} // end namespace internal

int MFEMReader::read(CurveArray &curves)
{
  SLIC_WARNING_IF(m_fileName.empty(), "Missing a filename in MFEMReader::read()");

  return internal::read_mfem(m_fileName,
    [&](mfem::Mesh *mesh, const std::map<int, axom::Array<int>> &contourZones)
  {
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

  return internal::read_mfem(m_fileName,
    [&](mfem::Mesh *mesh, const std::map<int, axom::Array<int>> &contourZones)
  {
    // Build CurvedPolygons from the MFEM zones.

    // Resize the array.
    curvedPolygons.clear();
    curvedPolygons.resize(contourZones.size());

    for(auto &[contourId, zoneIds] : contourZones)
    {
      auto& poly = curvedPolygons[contourId];
      for(int zoneId : zoneIds)
      {
        auto curve = axom::quest::internal::segment_to_curve(mesh, zoneId);
        poly.addEdge(std::move(curve));
      }
    }
  });
}

}  // end namespace quest
}  // end namespace axom
