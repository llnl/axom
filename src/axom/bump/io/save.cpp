// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/bump/io/save.hpp"
#include "axom/bump/views/Shapes.hpp"
#include "axom/bump/views/dispatch_coordset.hpp"

#include <string>
#include <cstdio>

namespace views = axom::bump::views;

namespace axom
{
namespace bump
{
namespace io
{

/**
 * \brief Turns a ShapeID to a VTK cell type value.
 *
 * \param shape_value A ShapeID.
 * \return A corresponding VTK cell type.
 */
static int ShapeID_to_vtk_cell(int shape_value)
{
  int vtktype = 0;
  if(shape_value == views::Tri_ShapeID)
    vtktype = 5;  // VTK_TRIANGLE
  else if(shape_value == views::Quad_ShapeID)
    vtktype = 9;  // VTK_QUAD
  else if(shape_value == views::Polygon_ShapeID)
    vtktype = 7;  // VTK_POLYGON
  else if(shape_value == views::Tet_ShapeID)
    vtktype = 10;  // VTK_TETRA
  else if(shape_value == views::Pyramid_ShapeID)
    vtktype = 14;  // VTK_PYRAMID
  else if(shape_value == views::Wedge_ShapeID)
    vtktype = 13;  // VTK_WEDGE
  else if(shape_value == views::Hex_ShapeID)
    vtktype = 12;  // VTK_HEXAHEDRON
  else if(shape_value == views::Polyhedron_ShapeID)
    vtktype = 42;  // VTK_POLYHEDRON

  return vtktype;
}

static void save_unstructured_vtk(const conduit::Node &mesh, const std::string &path)
{
  FILE *file = fopen(path.c_str(), "wt");
  if(file == nullptr)
  {
    SLIC_ERROR(fmt::format("The file {} could not be created.", path));
    return;
  }

  // Write the VTK file header
  fprintf(file, "# vtk DataFile Version 3.0\n");
  fprintf(file, "Unstructured Grid Example\n");
  fprintf(file, "ASCII\n");
  fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

  // Write the points
  const conduit::Node &coordset = mesh["coordsets"][0];
  const conduit::Node &points = coordset["values"];
  const auto x = points["x"].as_double_accessor();
  const auto y = points["y"].as_double_accessor();
  size_t num_points = 0;
  views::dispatch_coordset(coordset, [&](auto coordsetView) {
    num_points = coordsetView.size();
    fprintf(file, "POINTS %zu float\n", num_points);

    if(coordsetView.dimension() == 3)
    {
      const auto z = points["z"].as_double_accessor();
      for(size_t i = 0; i < num_points; ++i)
      {
        const auto p = coordsetView[i];
        fprintf(file,
                "%f %f %f\n",
                static_cast<float>(p[0]),
                static_cast<float>(p[1]),
                static_cast<float>(p[2]));
      }
    }
    else
    {
      for(size_t i = 0; i < num_points; ++i)
      {
        const auto p = coordsetView[i];
        fprintf(file, "%f %f 0\n", static_cast<float>(p[0]), static_cast<float>(p[1]));
      }
    }
  });

  // Write the cells
  const conduit::Node &topologies = mesh["topologies"];
  const conduit::Node &topo = topologies[0];
  const conduit::Node &elements = topo["elements"];
  const conduit::Node &connectivity = elements["connectivity"];
  size_t num_cells = elements["sizes"].dtype().number_of_elements();
  size_t total_num_indices = connectivity.dtype().number_of_elements();

  fprintf(file, "CELLS %zu %zu\n", num_cells, total_num_indices + num_cells);
  if(elements.has_path("offsets"))
  {
    for(size_t i = 0; i < num_cells; ++i)
    {
      size_t cell_size = elements["sizes"].as_int32_array()[i];
      size_t index = elements["offsets"].as_int32_array()[i];
      fprintf(file, "%zu", cell_size);
      for(size_t j = 0; j < cell_size; ++j)
      {
        fprintf(file, " %d", connectivity.as_int32_array()[index++]);
      }
      fprintf(file, "\n");
    }
  }
  else
  {
    size_t index = 0;
    for(size_t i = 0; i < num_cells; ++i)
    {
      size_t cell_size = elements["sizes"].as_int32_array()[i];
      fprintf(file, "%zu", cell_size);
      for(size_t j = 0; j < cell_size; ++j)
      {
        fprintf(file, " %d", connectivity.as_int32_array()[index++]);
      }
      fprintf(file, "\n");
    }
  }

  // Write the cell types
  fprintf(file, "CELL_TYPES %zu\n", num_cells);
  if(elements.has_child("shapes"))
  {
    const conduit::Node &shapes = elements["shapes"];
    for(size_t i = 0; i < num_cells; ++i)
    {
      const auto type = ShapeID_to_vtk_cell(shapes.as_int32_array()[i]);
      fprintf(file, "%d\n", type);
    }
  }
  else
  {
    const auto type = ShapeID_to_vtk_cell(views::shapeNameToID(elements["shape"].as_string()));
    for(size_t i = 0; i < num_cells; ++i)
    {
      fprintf(file, "%d\n", type);
    }
  }

  // TODO: write fields.

  // Close the file
  fclose(file);
}

void save_vtk(const conduit::Node &mesh, const std::string &path)
{
  const conduit::Node &n_topologies = mesh.fetch_existing("topologies");
  if(n_topologies.number_of_children() != 1)
  {
    SLIC_ERROR("The mesh must have a single topology.");
    return;
  }

  // For now.
  if(n_topologies[0].fetch_existing("type").as_string() != "unstructured")
  {
    SLIC_ERROR("The mesh must have a single unstructured topology.");
    return;
  }

  save_unstructured_vtk(mesh, path);
}

}  // end namespace io
}  // end namespace bump
}  // end namespace axom
