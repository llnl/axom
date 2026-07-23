// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

// MFEM does not support reading patch-based 1D NURBS meshes until after the v4.9
// release. Prefer patch-based extraction only when the MFEM version is new
// enough, otherwise fall back to element-based extraction.
#ifndef MFEM_VERSION
  #define MFEM_VERSION 0
#endif

#define AXOM_MFEM_MIN_VERSION_PATCH_BASED_1D_NURBS 40901

namespace axom
{
namespace quest
{
namespace internal
{
namespace
{
enum class PatchReadResult
{
  NotPatchFormat,
  Success,
  Failure
};

// MFEM 4.9.1 patch-format 1D NURBS reconstruction is intermittently corrupted
// in CUDA builds before Axom sees the patch data, so for actual "patches"
// files we re-read the patch blocks with MFEM's own stream parser instead.
PatchReadResult read_nurbs_patches_1d(const std::string &fileName,
                                      const mfem::Mesh &mesh,
                                      std::map<int, axom::Array<primal::NURBSCurve<double, 2>>> &curvemap)
{
  using ControlPoint = primal::Point<double, 2>;
  auto is_rational = [](const axom::Array<double> &weights) -> bool {
    if(weights.empty())
    {
      return false;
    }

    const double first = weights[0];
    for(int i = 1; i < weights.size(); ++i)
    {
      if(weights[i] != first)
      {
        return true;
      }
    }
    return false;
  };

  std::ifstream input(fileName);
  if(!input)
  {
    SLIC_WARNING(axom::fmt::format("Cannot open the provided MFEM mesh file '{}'", fileName));
    return PatchReadResult::Failure;
  }

  // MFEM's NURBSPatch stream constructor does not skip inline/comment-only
  // lines inside a patch block. Strip comments once up front, then let MFEM
  // parse the cleaned patch payload.
  std::ostringstream cleaned;
  std::string line;
  while(std::getline(input, line))
  {
    const auto comment_pos = line.find('#');
    if(comment_pos != std::string::npos)
    {
      line.erase(comment_pos);
    }
    cleaned << line << '\n';
  }

  std::istringstream sanitized_input(cleaned.str());
  std::string ident;
  bool found_patches = false;
  while(sanitized_input)
  {
    if(!(sanitized_input >> std::ws >> ident))
    {
      break;
    }
    if(ident == "patches")
    {
      found_patches = true;
      break;
    }
  }

  if(!found_patches)
  {
    return PatchReadResult::NotPatchFormat;
  }

  const auto *nurbsExt = mesh.NURBSext;
  if(nurbsExt == nullptr)
  {
    SLIC_WARNING("MFEM mesh is missing its NURBS extension");
    return PatchReadResult::Failure;
  }

  const int num_patches = nurbsExt->GetNP();
  for(int patchId = 0; patchId < num_patches; ++patchId)
  {
    // Keep patch attributes/topology from the loaded MFEM mesh, but source the
    // patch control net from the sanitized file stream to avoid the corrupted
    // in-memory CUDA reconstruction path.
    const int attribute = mesh.GetPatchAttribute(patchId);
    mfem::NURBSPatch patch(sanitized_input);
    if(patch.GetNKV() != 1 || patch.GetNC() != 3)
    {
      SLIC_WARNING(axom::fmt::format(
        "MFEM patch {} has unsupported dimensions for a 1D curve in 2D space (nkv={}, nc={}).",
        patchId,
        patch.GetNKV(),
        patch.GetNC()));
      return PatchReadResult::Failure;
    }

    mfem::KnotVector &kv0 = *patch.GetKV(0);
    const int ncp = kv0.GetNCP();
    axom::Array<ControlPoint> cp(0, ncp);
    axom::Array<double> weights(0, ncp);
    for(int i = 0; i < ncp; ++i)
    {
      const double wt = patch(i, 2);
      if(wt <= 0.0)
      {
        SLIC_WARNING(
          axom::fmt::format("MFEM patch {} has non-positive weight {} at control point {}.",
                            patchId,
                            wt,
                            i));
        return PatchReadResult::Failure;
      }

      cp.push_back({patch(i, 0) / wt, patch(i, 1) / wt});
      weights.push_back(wt);
    }

    const int degree = kv0.Size() - ncp - 1;
    if(degree < 0)
    {
      SLIC_WARNING(axom::fmt::format(
        "MFEM patch {} has inconsistent NURBS data: {} knots and {} control points.",
        patchId,
        kv0.Size(),
        ncp));
      return PatchReadResult::Failure;
    }

    axom::ArrayView<const double> knots_view(&kv0[0], kv0.Size());
    const primal::KnotVector<double> kv(knots_view, degree);

    is_rational(weights) ? curvemap[attribute].push_back({cp, weights, kv})
                         : curvemap[attribute].push_back({cp, kv});
  }

  return PatchReadResult::Success;
}
}  // namespace

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

  // lambda to extract the control points for element idx from the mfem mesh as an array of primal::Point
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

#if MFEM_VERSION < AXOM_MFEM_MIN_VERSION_PATCH_BASED_1D_NURBS
  auto get_element_weights = [fes, &mesh](int elemId) -> axom::Array<double> {
    mfem::Array<int> dofs;
    fes->GetElementDofs(elemId, dofs);

    const int nw = dofs.Size();
    axom::Array<double> w(nw, nw);

    mfem::Vector mfem_vec_weights(w.data(), nw);
    mesh->NURBSext->GetWeights().GetSubVector(dofs, mfem_vec_weights);

    return w;
  };

  auto get_element_degree = [fes, &mesh](int elemId) -> int {
    const mfem::Array<int> &orders = mesh->NURBSext->GetOrders();
    // MFEM calls this "order"; in Axom terminology this is the polynomial degree.
    return (elemId < orders.Size()) ? orders[elemId] : fes->GetOrder(elemId);
  };
#endif

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
#if MFEM_VERSION >= AXOM_MFEM_MIN_VERSION_PATCH_BASED_1D_NURBS
    // Parse patch-format 1D NURBS directly from the mesh file using MFEM's
    // NURBSPatch stream parser, bypassing the intermittently corrupted CUDA
    // in-memory NURBS reconstruction path.
    switch(read_nurbs_patches_1d(fileName, *mesh, curvemap))
    {
    case PatchReadResult::Success:
      return MFEMReader::READ_SUCCESS;
    case PatchReadResult::Failure:
      return MFEMReader::READ_FAILED;
    case PatchReadResult::NotPatchFormat:
      // Older MFEM NURBS files use the knotvectors-based format and should
      // continue through the normal instantiated-patch path below.
      break;
    }

    mfem::Array<mfem::NURBSPatch *> patches;
    mesh->GetNURBSPatches(patches);

    const int num_patches = patches.Size();
    for(int patchId = 0; patchId < num_patches; ++patchId)
    {
      const int attribute = mesh->GetPatchAttribute(patchId);
      mfem::NURBSPatch *patch = patches[patchId];
      if(patch == nullptr)
      {
        SLIC_WARNING(
          axom::fmt::format("MFEM patch {} could not be instantiated; cannot extract NURBS curve.",
                            patchId));
        for(int i = 0; i < patches.Size(); ++i)
        {
          delete patches[i];
        }
        return MFEMReader::READ_FAILED;
      }

      if(patch->GetNKV() != 1 || patch->GetNC() != 3)
      {
        SLIC_WARNING(axom::fmt::format(
          "MFEM patch {} has unsupported dimensions for a 1D curve in 2D space (nkv={}, nc={}).",
          patchId,
          patch->GetNKV(),
          patch->GetNC()));
        for(int i = 0; i < patches.Size(); ++i)
        {
          delete patches[i];
        }
        return MFEMReader::READ_FAILED;
      }

      const mfem::KnotVector &kv0 = *patch->GetKV(0);
      if(kv0.Size() <= 0)
      {
        SLIC_WARNING(
          axom::fmt::format("MFEM patch {} has an empty knot vector; cannot extract NURBS curve.",
                            patchId));
        for(int i = 0; i < patches.Size(); ++i)
        {
          delete patches[i];
        }
        return MFEMReader::READ_FAILED;
      }

      const int ncp = kv0.GetNCP();
      axom::Array<ControlPoint> cp(0, ncp);
      axom::Array<double> w(0, ncp);
      for(int i = 0; i < ncp; ++i)
      {
        const double wt = (*patch)(i, 2);
        if(wt <= 0.0)
        {
          SLIC_WARNING(
            axom::fmt::format("MFEM patch {} has non-positive weight {} at control point {}.",
                              patchId,
                              wt,
                              i));
          for(int j = 0; j < patches.Size(); ++j)
          {
            delete patches[j];
          }
          return MFEMReader::READ_FAILED;
        }

        cp.push_back({(*patch)(i, 0) / wt, (*patch)(i, 1) / wt});
        w.push_back(wt);
      }

      const int degree = kv0.Size() - ncp - 1;
      if(degree < 0)
      {
        SLIC_WARNING(axom::fmt::format(
          "MFEM patch {} has inconsistent NURBS data: {} knots and {} control points.",
          patchId,
          kv0.Size(),
          ncp));
        for(int i = 0; i < patches.Size(); ++i)
        {
          delete patches[i];
        }
        return MFEMReader::READ_FAILED;
      }

      axom::ArrayView<const double> knots_view(&kv0[0], kv0.Size());
      const primal::KnotVector<double> kv(knots_view, degree);

      is_rational(w) ? curvemap[attribute].push_back({cp, w, kv})
                     : curvemap[attribute].push_back({cp, kv});
    }

    for(int i = 0; i < patches.Size(); ++i)
    {
      delete patches[i];
    }
#else
    {
      // MFEM versions prior to AXOM_MFEM_MIN_VERSION_PATCH_BASED_1D_NURBS do not
      // support reading patch-based 1D NURBS meshes. In that case, treat each
      // MFEM element as a single (rational) Bezier span.
      for(int zoneId = 0; zoneId < mesh->GetNE(); ++zoneId)
      {
        const int attribute = mesh->GetAttribute(zoneId);
        const int degree = get_element_degree(zoneId);

        const auto cp = get_controlpoints(zoneId);
        const auto w = get_element_weights(zoneId);

        is_rational(w) ? curvemap[attribute].push_back({cp, w, degree})
                       : curvemap[attribute].push_back({cp, degree});
      }
    }
#endif
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
