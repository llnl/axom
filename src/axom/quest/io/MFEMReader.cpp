// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/MFEMReader.hpp"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#ifndef AXOM_USE_MFEM
  #error MFEMReader should only be included when Axom is configured with MFEM
#endif
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
 * \param fileName The name of the file to read.
 * \param curvemap Output map from MFEM attribute value to the associated curves.
 *
 * \return 0 on success; non-zero on failure.
 */
int read_mfem(const std::string &fileName,
              std::map<int, axom::Array<primal::NURBSCurve<double, 2>>> &curvemap)
{
  if(!axom::utilities::filesystem::pathExists(fileName))
  {
    SLIC_WARNING(axom::fmt::format("Cannot open the provided MFEM mesh file '{}'", fileName));
    return MFEMReader::READ_FAILED;
  }

  // Load the MFEM file
  constexpr int generate_edges = 1;
  constexpr int refine = 1;
  constexpr bool fix_orientation = true;
  auto mesh = std::make_unique<mfem::Mesh>(fileName, generate_edges, refine, fix_orientation);

  if(mesh->Dimension() != 1 || mesh->SpaceDimension() != 2)
  {
    SLIC_WARNING(
      axom::fmt::format("Mesh must have dimension 1 and spatial dimension 2. The supplied mesh "
                        "is dimension {} with spatial dimension {}.",
                        mesh->Dimension(),
                        mesh->SpaceDimension()));
    return MFEMReader::READ_FAILED;
  }

  const auto *nodes = mesh->GetNodes();
  const auto *fes = nodes != nullptr ? nodes->FESpace() : nullptr;
  const auto *fec = fes != nullptr ? fes->FEColl() : nullptr;
  if(nodes == nullptr || fes == nullptr || fec == nullptr)
  {
    SLIC_WARNING("Mesh does not have a valid nodes grid function");
    return MFEMReader::READ_FAILED;
  }

  // lambda to extract the knot vector associated with curve idx. Converts from mfem::KnotVector to primal::KnotVector
  auto get_knots = [&mesh](int idx) -> primal::KnotVector<double> {
    mfem::Array<const mfem::KnotVector *> kvs;
    mesh->NURBSext->GetPatchKnotVectors(idx, kvs);
    const mfem::KnotVector &kv = *kvs[0];

    axom::ArrayView<const double> knots_view(&kv[0], kv.Size());
    return primal::KnotVector<double>(knots_view, kv.GetOrder());
  };

  // lambda to extract the control points for curve idx from the mfem mesh as an array of primal::Point
  using ControlPoint = primal::Point<double, 2>;
  auto get_controlpoints = [nodes, fes](int idx) -> axom::Array<ControlPoint> {
    mfem::Array<int> vdofs;
    mfem::Vector v;
    fes->GetElementVDofs(idx, vdofs);
    nodes->GetSubVector(vdofs, v);
    const auto ord = fes->GetOrdering();

    const int ncp = v.Size() / 2;
    axom::Array<ControlPoint> cp(0, ncp);
    for(int i = 0; i < ncp; ++i)
    {
      ord == mfem::Ordering::byVDIM ? cp.push_back({v[i], v[i + ncp]})
                                    : cp.push_back({v[2 * i], v[2 * i + 1]});
    }

    return cp;
  };

  // lambda to extract the weights for curve idx from the mfem mesh
  auto get_weights = [&mesh, fes](int idx) -> axom::Array<double> {
    mfem::Array<int> dofs;
    fes->GetElementDofs(idx, dofs);

    const int NW = dofs.Size();
    axom::Array<double> w(NW, NW);

    // wrap our array's buffer w/ an mfem::Vector for GetSubVector
    mfem::Vector mfem_vec_weights(w.data(), NW);
    mesh->NURBSext->GetWeights().GetSubVector(dofs, mfem_vec_weights);

    return w;
  };

  // lambda to check if the weights correspond to a rational curve. If they are all equal it is not rational
  auto is_rational = [](const axom::Array<double> &weights) -> bool {
    const int sz = weights.size();
    if(sz == 0)
    {
      return false;
    }

    const double first = weights[0];
    for(int i = 1; i < sz; ++i)
    {
      if(weights[i] != first)  // strict equality is fine for this
      {
        return true;
      }
    }

    return false;
  };

  // Examine the mesh attributes and group all of the related curves w/ same attribute
  // Assumption is that they're part of the same contour
  const bool isNURBS = dynamic_cast<const mfem::NURBSFECollection *>(fec) != nullptr;
  if(isNURBS)
  {
    const int num_patches = fes->GetNURBSext()->GetNP();
    for(int patchId = 0; patchId < num_patches; ++patchId)
    {
      const int attribute = mesh->GetPatchAttribute(patchId);

      const auto kv = get_knots(patchId);
      const auto cp = get_controlpoints(patchId);
      const auto w = get_weights(patchId);

      is_rational(w) ? curvemap[attribute].push_back({cp, w, kv})
                     : curvemap[attribute].push_back({cp, kv});
    }
  }
  else
  {
    const bool is_bernstein = dynamic_cast<const mfem::H1Pos_FECollection *>(fec) != nullptr;
    if(!is_bernstein)
    {
      SLIC_WARNING(axom::fmt::format(
        "Non-NURBS meshes must define their nodes in the positive Bernstein basis "
        "(mfem::H1Pos_FECollection). Got FECollection '{}'.",
        fec->Name()));
      return MFEMReader::READ_FAILED;
    }

    for(int zoneId = 0; zoneId < mesh->GetNE(); zoneId++)
    {
      const int attribute = mesh->GetAttribute(zoneId);
      auto control_points = get_controlpoints(zoneId);
      const int order = fes->GetOrder(zoneId);

      // mfem::H1Pos_FECollection stores 1D element DOFs grouped by entity (vertices, then edges).
      // For a segment, this yields local ordering: [v0, v1, e0, e1, ...].
      // Convert to Bezier control point ordering: [v0, e0, e1, ..., v1].
      const int num_cp = control_points.size();
      if(order > 1 && num_cp == order + 1)
      {
        axom::Array<ControlPoint> reordered(0, num_cp);
        reordered.push_back(control_points[0]);
        for(int i = 2; i < num_cp; ++i)
        {
          reordered.push_back(control_points[i]);
        }
        reordered.push_back(control_points[1]);
        control_points = std::move(reordered);
      }

      curvemap[attribute].push_back({control_points, order});
    }
  }

  return MFEMReader::READ_SUCCESS;
}

}  // end namespace internal

int MFEMReader::read(CurveArray &curves)
{
  axom::Array<int> attributes;
  return read(curves, attributes);
}

int MFEMReader::read(CurveArray &curves, axom::Array<int> &attributes)
{
  SLIC_WARNING_IF(m_fileName.empty(), "Missing a filename in MFEMReader::read()");

  curves.clear();
  attributes.clear();
  std::map<int, CurveArray> curvemap;
  const int ret = internal::read_mfem(m_fileName, curvemap);
  if(ret == READ_SUCCESS)
  {
    for(auto &[attribute, nurbs] : curvemap)
    {
      for(const auto &curve : nurbs)
      {
        curves.push_back(curve);
        attributes.push_back(attribute);
      }
    }
  }

  return ret;
}

int MFEMReader::read(CurvedPolygonArray &curvedPolygons)
{
  axom::Array<int> attributes;
  return read(curvedPolygons, attributes);
}

int MFEMReader::read(CurvedPolygonArray &curvedPolygons, axom::Array<int> &attributes)
{
  SLIC_WARNING_IF(m_fileName.empty(), "Missing a filename in MFEMReader::read()");

  curvedPolygons.clear();
  attributes.clear();
  std::map<int, CurveArray> curvemap;
  const int ret = internal::read_mfem(m_fileName, curvemap);
  if(ret == READ_SUCCESS)
  {
    curvedPolygons.resize(curvemap.size());
    attributes.resize(curvemap.size());

    int polygon_index = 0;
    for(auto &[attribute, nurbs] : curvemap)
    {
      attributes[polygon_index] = attribute;
      auto &poly = curvedPolygons[polygon_index];
      for(auto &cur : nurbs)
      {
        poly.addEdge(cur);
      }
      ++polygon_index;
    }
  }
  return ret;
}

}  // end namespace quest
}  // end namespace axom
