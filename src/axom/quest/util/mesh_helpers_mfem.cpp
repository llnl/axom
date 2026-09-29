// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "mesh_helpers.hpp"

namespace axom::quest::util
{
mfem::Mesh* make_cartesian_mfem_mesh_2D(const primal::BoundingBox<double, 2>& bbox,
                                        const NumericArray<int, 2>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling)
{
  constexpr int DIM = 2;
  const auto range = bbox.range();

  auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian2D(res[0],
                                                          res[1],
                                                          mfem::Element::QUADRILATERAL,
                                                          true,
                                                          range[0],
                                                          range[1],
                                                          reorder_space_filling));

  // Offset the mesh to lie w/in the bounding box
  const int NV = mesh->GetNV();
  for(int i = 0; i < NV; ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < DIM; ++d)
    {
      v[d] += bbox.getMin()[d];
    }
  }

  // Ensure that mesh has high order nodes
  mesh->SetCurvature(polynomial_order);

  return mesh;
}

mfem::Mesh* make_cartesian_mfem_mesh_3D(const primal::BoundingBox<double, 3>& bbox,
                                        const NumericArray<int, 3>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling)
{
  constexpr int DIM = 3;
  const auto range = bbox.range();

  auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian3D(res[0],
                                                          res[1],
                                                          res[2],
                                                          mfem::Element::HEXAHEDRON,
                                                          range[0],
                                                          range[1],
                                                          range[2],
                                                          reorder_space_filling));

  // Offset the mesh to lie w/in the bounding box
  const int NV = mesh->GetNV();
  for(int i = 0; i < NV; ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < DIM; ++d)
    {
      v[d] += bbox.getMin()[d];
    }
  }

  // Ensure that mesh has high order nodes
  mesh->SetCurvature(polynomial_order);

  return mesh;
}

}  // namespace axom::quest::util
