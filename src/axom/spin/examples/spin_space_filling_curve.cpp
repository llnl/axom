// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *  \file spin_space_filling_curve.cpp
 *  \brief Visualize Morton and Hilbert ordering of random points and a BVH.
 *
 *  Sample uniform points in the unit square or cube, and use a RectangularLattice
 *  to maps them to integer cells for Morton or Hilbert encoding.
 *  Point keys use 32 bits, with 16 bits per axis in 2D and 10 in 3D.
 *
 *  A separate display grid has 2^level cells per axis.
 *  Level 0 chooses its resolution from the point count and changing the level
 *  does not affect point keys, sort order, or BVH construction.
 *
 *  The example uses SLAM to encode the coordinates, fields, and line connectivity
 *  and outputs them as Blueprint meshes:
 *
 *  - _points.root: Points with their original IDs, keys, and sort ranks.
 *  - _grid.root: Display grid with a key for each cell.
 *  - _curve.root: Curve joining all display cell centers in key order.
 *  - _point_segments.root: Segments joining sampled points in key order.
 *  - _bvh.vtk: BVH boxes with a "level" cell field for depth, starting at zero.
 *
 *  The BVH uses the selected ordering and one box with zero extent per point.
 *  It normalizes coordinates to the input bounds. Point mesh keys use the unit square or cube.
 */

#include "axom/config.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/spin.hpp"
#include "axom/slam.hpp"
#include "axom/slic.hpp"
#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#if !defined(AXOM_USE_CONDUIT)
  #error "The space filling curves example requires Conduit"
#endif

#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay_io.hpp"
#include "conduit_relay_io_blueprint.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

namespace
{
using KeyType = std::uint64_t;
using SortKeyType = std::uint32_t;
using IndexType = axom::IndexType;
using axom::spin::SpaceFillingCurve;

using EntitySet = axom::slam::PositionSet<IndexType, IndexType>;
template <typename T>
using Field = axom::slam::Map<T, EntitySet>;
using SegmentToVertex = axom::slam::ConstantRelation<EntitySet, EntitySet, 2>;

template <int DIM>
using Point = axom::primal::Point<double, DIM>;
template <int DIM, typename CoordType>
using Lattice = axom::spin::RectangularLattice<DIM, double, CoordType>;

// Options and grid sizes.
struct Input
{
  SpaceFillingCurve order {SpaceFillingCurve::Hilbert};
  int dimension {2};
  int level {0};
  IndexType num_points {1000};
  unsigned int seed {42};
  KeyType max_cells {1 << 20};
  std::string output {"space_filling_curve"};
  std::string protocol {"json"};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    const std::map<std::string, SpaceFillingCurve> orders {{"morton", SpaceFillingCurve::Morton},
                                                           {"hilbert", SpaceFillingCurve::Hilbert}};
    app.add_option("--order", order, "Space-filling curve")
      ->transform(axom::CLI::CheckedTransformer(orders))
      ->default_str("hilbert");
    app.add_option("-d,--dimension,--dim", dimension, "Spatial dimension")
      ->check(axom::CLI::IsMember({2, 3}))
      ->capture_default_str();
    app.add_option("-l,--level", level, "Display grid depth; 0 chooses from the point count")
      ->check(axom::CLI::Range(0, 16))
      ->capture_default_str();
    app.add_option("-n,--num-points", num_points)
      ->description("Number of uniform points in the unit square or cube")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
    app.add_option("--seed", seed, "Random seed. Reuse it to compare curves on the same points")
      ->check(axom::CLI::NonNegativeNumber)
      ->capture_default_str();
    app.add_option("--max-cells", max_cells, "Maximum cells to create in the grid and curve")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();
    app.add_option("-o,--output", output)
      ->description("Output prefix for the four Blueprint meshes and _bvh.vtk")
      ->capture_default_str();
    app.add_option("--protocol", protocol, "Blueprint output format")
      ->check(axom::CLI::IsMember({"json", "hdf5"}))
      ->capture_default_str();
    app.parse(argc, argv);

    // Choose the deepest grid with at most num_points cells, with a minimum level of one.
    if(level == 0)
    {
      // floor(log2(num_points) / dimension), with at least two cells per axis.
      for(auto count = num_points; count >= (IndexType {1} << dimension); count >>= dimension)
      {
        ++level;
      }
      level = std::max(1, level);
    }
    if(level > sortLevel())
    {
      throw axom::CLI::ValidationError("--level", "Display grid exceeds 32-bit sorting resolution");
    }

    if(output.empty())
    {
      throw axom::CLI::ValidationError("--output", "The output prefix must not be empty");
    }

    if(numCells() > max_cells)
    {
      throw axom::CLI::ValidationError(
        "--level",
        axom::fmt::format(
          "{}D level {} requires {} cells; reduce --level or raise --max-cells ({})",
          dimension,
          level,
          numCells(),
          max_cells));
    }
    // The VTK writer uses 32-bit counts for 2*N-1 boxes and their corner connectivity.
    const auto max_bvh_nodes = std::numeric_limits<std::int32_t>::max() / ((1 << dimension) + 1);
    if(num_points > (max_bvh_nodes + 1) / 2)
    {
      throw axom::CLI::ValidationError("--num-points", "BVH export exceeds VTK index capacity");
    }
    // Each segment needs two indices. Check array sizes and Conduit's byte offsets.
    const auto point_bytes = dimension * sizeof(double);
    const auto segment_bytes = 2 * sizeof(IndexType);
    const auto max_entries = std::min<KeyType>(
      {std::numeric_limits<IndexType>::max() / 2 + 1,
       std::numeric_limits<std::size_t>::max() / std::max(2 * point_bytes, segment_bytes),
       std::numeric_limits<conduit::index_t>::max() / std::max(point_bytes, segment_bytes)});
    if(numCells() > max_entries || static_cast<KeyType>(num_points) > max_entries)
    {
      throw axom::CLI::ValidationError("--level/--num-points",
                                       "Requested arrays exceed index or storage capacity");
    }

    conduit::Node relay_info;
    conduit::relay::io::about(relay_info);
    if(relay_info["protocols/" + protocol].as_string() != "enabled")
    {
      throw axom::CLI::ValidationError("--protocol",
                                       protocol + " is unavailable in this Conduit build");
    }
  }

  // Number of display cells along each axis.
  KeyType cellsPerAxis() const { return KeyType {1} << level; }

  // Number of display cells and curve vertices.
  KeyType numCells() const { return KeyType {1} << (dimension * level); }

  // Bits per axis in a point key.
  int sortLevel() const { return std::numeric_limits<SortKeyType>::digits / dimension; }

  std::string orderName() const
  {
    return order == SpaceFillingCurve::Hilbert ? "hilbert" : "morton";
  }

  // Shared name for the key field in all four meshes.
  std::string indexField() const { return orderName() + "_index"; }
};

template <int DIM>
using GridPoint = axom::primal::Point<KeyType, DIM>;

// Encode a cell with the selected codec. Hilbert also needs the curve's level.
template <typename CoordType, int DIM>
CoordType encode(const axom::primal::Point<CoordType, DIM>& cell, SpaceFillingCurve order, int level)
{
  return order == SpaceFillingCurve::Hilbert
    ? axom::spin::convertPointToHilbert<CoordType>(cell, level)
    : axom::spin::convertPointToMorton<CoordType>(cell);
}

// Recover integer cell coordinates from a key at the display grid's level.
template <int DIM>
GridPoint<DIM> decode(KeyType index, const Input& params)
{
  return params.order == SpaceFillingCurve::Hilbert
    ? axom::spin::convertHilbertToPoint<KeyType, DIM>(index, params.level)
    : axom::spin::convertMortonToPoint<KeyType, DIM>(index);
}

// Draw a point inside the box from the seeded random sequence.
template <int DIM>
Point<DIM> getRandomPoint(const axom::primal::BoundingBox<double, DIM>& bbox, unsigned int seed)
{
  Point<DIM> point;
  for(int d = 0; d < DIM; ++d)
  {
    point[d] = axom::utilities::random_real(bbox.getMin()[d], bbox.getMax()[d], seed);
  }
  return point;
}

// Sample points and connect them in key order. Keep coordinates in input order.
template <int DIM>
void orderPoints(const Input& params,
                 const axom::primal::BoundingBox<double, DIM>& domain,
                 const Lattice<DIM, SortKeyType>& lattice,
                 Field<Point<DIM>>& positions,
                 Field<SortKeyType>& keys,
                 Field<IndexType>& ranks,
                 SegmentToVertex& segment_to_point)
{
  const int bits = params.sortLevel();
  // Inset the sampling box so gridCell() stays below the lattice's upper boundary.
  constexpr double epsilon = 1e-6;
  auto sampling_domain = domain;
  sampling_domain.expand(-epsilon);

  // Encode points on the sorting lattice, independently of the display level.
  const auto& points = *positions.set();
  for(const auto i : points)
  {
    positions[i] = getRandomPoint(sampling_domain, params.seed);
    const auto cell = lattice.gridCell(positions[i]);
    keys[i] = encode(cell, params.order, bits);
  }

  // Stable sorting preserves input order for equal keys.
  Field<IndexType> permutation(&points);
  std::iota(permutation.begin(), permutation.end(), IndexType {0});
  std::stable_sort(permutation.begin(), permutation.end(), [&](IndexType a, IndexType b) {
    return keys[a] < keys[b];
  });

  // Map each original ID to its sort rank.
  for(const auto rank : points)
  {
    ranks[permutation[rank]] = rank;
  }

  // Endpoint IDs refer to the original point array. A single point has no segments.
  for(const auto segment : *segment_to_point.fromSet())
  {
    auto endpoints = segment_to_point[segment];
    endpoints[0] = permutation[segment];
    endpoints[1] = permutation[segment + 1];
  }
  SLIC_ASSERT(segment_to_point.isValid());
}

// Connect display cell centers in key order.
// Grid keys use Cartesian cell IDs; curve vertices use their keys as IDs.
template <int DIM>
void makeSpaceFillingCurve(const Input& params,
                           const Lattice<DIM, KeyType>& lattice,
                           Field<KeyType>& grid_keys,
                           Field<Point<DIM>>& centers,
                           SegmentToVertex& segment_to_vertex)
{
  const auto side = params.cellsPerAxis();
  for(const auto key : *centers.set())
  {
    const auto cell = decode<DIM>(static_cast<KeyType>(key), params);
    centers[key] = lattice.cellBounds(cell).getCentroid();

    // The grid's cell IDs run with x varying fastest, independently of curve order.
    KeyType ordinal = 0;
    for(int d = DIM - 1; d >= 0; --d)
    {
      ordinal = ordinal * side + cell[d];
    }
    grid_keys[static_cast<IndexType>(ordinal)] = static_cast<KeyType>(key);
  }

  for(const auto segment : *segment_to_vertex.fromSet())
  {
    auto endpoints = segment_to_vertex[segment];
    endpoints[0] = segment;
    endpoints[1] = segment + 1;
  }
  SLIC_ASSERT(segment_to_vertex.isValid());
}

// Build a BVH over point boxes in input order.
template <int DIM>
void buildBVH(const Input& params, const Field<Point<DIM>>& positions, axom::spin::BVH<DIM>& bvh)
{
  using Box = typename axom::spin::BVH<DIM>::BoxType;
  Field<Box> boxes(positions.set());
  for(const auto i : *positions.set())
  {
    boxes[i] = Box(positions[i]);
  }
  bvh.setScaleFactor(1.);
  if(bvh.initialize(boxes.data().view(), boxes.size(), {.ordering = params.order}) !=
     axom::spin::BVH_BUILD_OK)
  {
    throw std::runtime_error("Failed to build the point BVH");
  }
}

constexpr const char* AXES[] = {"x", "y", "z"};
constexpr const char* GRID_AXES[] = {"i", "j", "k"};
constexpr const char* SPACING_AXES[] = {"dx", "dy", "dz"};

// Share the point buffer with Conduit using one strided array per axis.
template <int DIM>
void writeCoordinates(conduit::Node& mesh, Field<Point<DIM>>& positions)
{
  auto& coords = mesh["coordsets/coords"];
  coords["type"] = "explicit";
  for(int d = 0; d < DIM; ++d)
  {
    // Conduit offsets and strides are in bytes, not coordinate counts.
    coords["values"][AXES[d]].set_external(positions[0].data(),
                                           positions.size(),
                                           d * sizeof(double),
                                           sizeof(Point<DIM>));
  }
}

// Share a field's buffer with Conduit.
template <typename T>
void writeField(conduit::Node& mesh,
                const std::string& name,
                const std::string& association,
                axom::ArrayView<T> values)
{
  auto& field = mesh["fields"][name];
  field["association"] = association;
  field["topology"] = "mesh";
  field["values"].set_external(values.data(), values.size());
}

// Add metadata, verify Blueprint, and write the mesh.
void saveMesh(conduit::Node& mesh, const Input& params, const std::string& suffix, int key_level)
{
  mesh["state/domain_id"] = 0;
  mesh["state/cycle"] = 0;
  mesh["state/order"] = params.orderName();
  mesh["state/level"] = params.level;
  mesh["state/key_level"] = key_level;
  mesh["state/seed"] = params.seed;
  conduit::Node info;
  if(!conduit::blueprint::verify("mesh", mesh, info))
  {
    throw std::runtime_error("Invalid " + suffix + " Blueprint mesh:\n" + info.to_yaml());
  }
  conduit::Node options;
  options["file_style"] = "root_only";
  options["suffix"] = "none";
  const std::string path = params.output + "_" + suffix;
  conduit::relay::io::blueprint::save_mesh(mesh, path, params.protocol, options);
  SLIC_INFO(axom::fmt::format("Wrote {}.root ({})", path, params.protocol));
}

// Write the four Blueprint meshes and the BVH.
template <int DIM>
void writeMeshes(Input params,
                 const Lattice<DIM, KeyType>& grid,
                 Field<Point<DIM>>& positions,
                 Field<SortKeyType>& keys,
                 Field<IndexType>& ranks,
                 SegmentToVertex& segment_to_point,
                 Field<KeyType>& grid_keys,
                 Field<Point<DIM>>& centers,
                 SegmentToVertex& segment_to_vertex,
                 const axom::spin::BVH<DIM>& bvh)
{
  // Relay treats ':' as a subpath separator, so save from the output directory
  // to avoid passing a Windows drive letter to Relay.
  const auto output_prefix = std::filesystem::absolute(params.output);
  std::filesystem::current_path(output_prefix.parent_path());
  params.output = output_prefix.filename().string();

  // Share one ID array across the output fields. Keep all buffers alive until saving finishes.
  axom::Array<KeyType> ids(std::max(positions.size(), centers.size()));
  std::iota(ids.begin(), ids.end(), KeyType {0});

  // Points and sorted segments share coordinates and vertex fields.
  {
    conduit::Node mesh;
    writeCoordinates(mesh, positions);
    mesh["topologies/mesh/type"] = "points";
    mesh["topologies/mesh/coordset"] = "coords";
    writeField(mesh, params.indexField(), "vertex", keys.data().view());
    writeField(mesh, "point_id", "vertex", ids.view().subspan(0, positions.size()));
    writeField(mesh, "sort_rank", "vertex", ranks.data().view());
    saveMesh(mesh, params, "points", params.sortLevel());

    auto& topology = mesh["topologies/mesh"];
    topology["type"] = "unstructured";
    topology["elements/shape"] = "line";
    auto& connectivity = segment_to_point.relationData();
    topology["elements/connectivity"].set_external(connectivity.data(), connectivity.size());
    writeField(mesh,
               "segment_index",
               "element",
               ids.view().subspan(0, segment_to_point.fromSet()->size()));
    saveMesh(mesh, params, "point_segments", params.sortLevel());
  }

  // Display grid with a key for each cell.
  {
    conduit::Node mesh;
    auto& coords = mesh["coordsets/coords"];
    coords["type"] = "uniform";
    for(int d = 0; d < DIM; ++d)
    {
      // Blueprint dimensions count vertices, not cells.
      coords["dims"][GRID_AXES[d]] = params.cellsPerAxis() + 1;
      coords["origin"][AXES[d]] = grid.origin()[d];
      coords["spacing"][SPACING_AXES[d]] = grid.spacing()[d];
    }
    mesh["topologies/mesh/type"] = "uniform";
    mesh["topologies/mesh/coordset"] = "coords";
    writeField(mesh, params.indexField(), "element", grid_keys.data().view());
    saveMesh(mesh, params, "grid", params.level);
  }

  // Curve through all display cell centers.
  {
    conduit::Node mesh;
    writeCoordinates(mesh, centers);
    auto& topology = mesh["topologies/mesh"];
    topology["type"] = "unstructured";
    topology["coordset"] = "coords";
    topology["elements/shape"] = "line";
    auto& connectivity = segment_to_vertex.relationData();
    topology["elements/connectivity"].set_external(connectivity.data(), connectivity.size());
    writeField(mesh, params.indexField(), "vertex", ids.view().subspan(0, centers.size()));
    writeField(mesh,
               "segment_index",
               "element",
               ids.view().subspan(0, segment_to_vertex.fromSet()->size()));
    saveMesh(mesh, params, "curve", params.level);
  }

  const auto bvh_path = params.output + "_bvh.vtk";
  bvh.writeVtkFile(bvh_path);
  SLIC_INFO(axom::fmt::format("Wrote {}", bvh_path));
}

// Generate the points, display curve, and BVH, then write the output files.
template <int DIM>
void run(const Input& params)
{
  SLIC_INFO(
    axom::fmt::format("{} order, {}D display level {}, point-key depth {}: {} points, {} cells",
                      params.orderName(),
                      DIM,
                      params.level,
                      params.sortLevel(),
                      params.num_points,
                      params.numCells()));

  const axom::primal::BoundingBox<double, DIM> domain(Point<DIM>(0.), Point<DIM>(1.));
  const auto sorting_lattice = axom::spin::rectangular_lattice_from_bounding_box(
    domain,
    axom::NumericArray<SortKeyType, DIM>(SortKeyType {1} << params.sortLevel()));
  const auto display_lattice = axom::spin::rectangular_lattice_from_bounding_box(
    domain,
    axom::NumericArray<KeyType, DIM>(params.cellsPerAxis()));

  EntitySet points(params.num_points);
  Field<Point<DIM>> positions(&points);
  Field<SortKeyType> keys(&points);
  Field<IndexType> ranks(&points);
  EntitySet point_segments(points.size() - 1);
  axom::Array<IndexType> point_connectivity(2 * point_segments.size());
  auto segment_to_point =
    axom::slam::make_constant_relation_ct<2>(&point_segments, &points, point_connectivity);
  orderPoints(params, domain, sorting_lattice, positions, keys, ranks, segment_to_point);

  axom::spin::BVH<DIM> bvh;
  buildBVH(params, positions, bvh);

  EntitySet cells(static_cast<IndexType>(params.numCells()));
  EntitySet vertices(cells.size());
  EntitySet segments(vertices.size() - 1);
  Field<KeyType> grid_keys(&cells);
  Field<Point<DIM>> centers(&vertices);
  axom::Array<IndexType> connectivity(2 * segments.size());
  auto segment_to_vertex =
    axom::slam::make_constant_relation_ct<2>(&segments, &vertices, connectivity);
  makeSpaceFillingCurve(params, display_lattice, grid_keys, centers, segment_to_vertex);

  writeMeshes(params,
              display_lattice,
              positions,
              keys,
              ranks,
              segment_to_point,
              grid_keys,
              centers,
              segment_to_vertex,
              bvh);
}
}  // namespace

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;
  Input params;
  axom::CLI::App app {"Write Morton or Hilbert Blueprint meshes and BVH boxes in VTK format"};
  try
  {
    params.parse(argc, argv, app);
    if(params.dimension == 2)
    {
      run<2>(params);
    }
    else
    {
      run<3>(params);
    }
  }
  catch(const axom::CLI::ParseError& error)
  {
    return app.exit(error);
  }
  catch(const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
