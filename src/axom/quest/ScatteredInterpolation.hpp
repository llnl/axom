// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_SCATTERED_INTERPOLATION_H_
#define QUEST_SCATTERED_INTERPOLATION_H_

#include "axom/core.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"
#include "axom/spin.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"

#include "axom/fmt.hpp"

#include "conduit.hpp"
#include "conduit_blueprint.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>
#include <vector>

namespace
{
/// Helper function to extract the dimension from the coordinate values group
/// of a mesh blueprint coordset
inline int extractDimension(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node.has_child("z") ? 3 : (values_node.has_child("y") ? 2 : 1);
}

/// Helper function to extract the number of points from the coordinate values group
/// of a mesh blueprint coordset
inline int extractSize(const conduit::Node& values_node)
{
  SLIC_ASSERT(values_node.has_child("x"));
  return values_node["x"].dtype().number_of_elements();
}

inline bool getScatteredInterpSeed(std::uint64_t& seed)
{
  const char* env = std::getenv("AXOM_SCATTERED_INTERP_SEED");
  if(env == nullptr || env[0] == '\0')
  {
    return false;
  }

  char* end = nullptr;
  const auto parsed = std::strtoull(env, &end, 10);
  if(end == env)
  {
    return false;
  }

  seed = static_cast<std::uint64_t>(parsed);
  return true;
}

/**
 * \brief Utility function to create an axom::ArrayView over the array
 * of native types stored by a conduit::Node
 */
template <typename T>
inline axom::ArrayView<T> ArrayView_from_Node(conduit::Node& node, int sz)
{
  T* ptr = node.value();
  return axom::ArrayView<T>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,2>
 *
 * \warning Assumes the underlying data is an MCArray with stride 2 access
 */
template <>
inline axom::ArrayView<axom::primal::Point<double, 2>> ArrayView_from_Node(conduit::Node& node, int sz)
{
  using PointType = axom::primal::Point<double, 2>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

/**
 * \brief Template specialization of ArrayView_from_Node for Point<double,3>
 *
 * \warning Assumes the underlying data is an MCArray with stride 3 access
 */
template <>
inline axom::ArrayView<axom::primal::Point<double, 3>> ArrayView_from_Node(conduit::Node& node, int sz)
{
  using PointType = axom::primal::Point<double, 3>;

  PointType* ptr = static_cast<PointType*>(node.data_ptr());
  return axom::ArrayView<PointType>(ptr, sz);
}

/**
 * \brief Returns an ArrayView to the desired field from a mesh blueprint, which must be present
 *
 * \param [in] mesh_node The root conduit node of a valid mesh blueprint
 * \param [in] field_name The name of the field to return; Can either be either the name of the field,
 *  or the path to the field, e.g. fields/<field_name>/values
 *
 * \pre The field named \a field_name must be present in \a mesh_node or code will error out
 * \note Only currently supports scalar fields
 */
template <typename T>
inline axom::ArrayView<T> getField(conduit::Node& mesh_node, const std::string& field_name)
{
  using axom::utilities::string::startsWith;

  const std::string field_path =
    startsWith(field_name, "field") ? field_name : axom::fmt::format("fields/{}/values", field_name);
  SLIC_ERROR_IF(!mesh_node.has_path(field_path),
                axom::fmt::format("Mesh blueprint is missing required field '{}'", field_name));

  auto& values = mesh_node[field_path];
  const auto sz = values.dtype().number_of_elements();

  return ::ArrayView_from_Node<T>(values, sz);
}

/// Check validity of blueprint group
inline bool isValidBlueprint(const conduit::Node& mesh_node)
{
  bool success = true;
  conduit::Node info;
  if(!conduit::blueprint::verify("mesh", mesh_node, info))
  {
    SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
    success = false;
  }

  return success;
}

}  // namespace

namespace axom
{
namespace quest
{
namespace detail
{
/**
 * \brief Utility class to enable processing an array of points whose layout
 * is either interleaved or separated strided arrays
 */
template <typename T, int NDIMS>
struct InterleavedOrStridedPoints
{
public:
  static constexpr int DIM = NDIMS;
  using CoordType = T;
  using PointType = primal::Point<CoordType, NDIMS>;
  using StridedPoints = primal::detail::ZipBase<PointType>;
  using InterleavedPoints = axom::ArrayView<PointType>;

  /// Constructor from a multi-component array Conduit node
  explicit InterleavedOrStridedPoints(conduit::Node& values)
  {
    SLIC_ASSERT(isMultiComponentArray(values));

    m_is_interleaved = conduit::blueprint::mcarray::is_interleaved(values);

    m_npts = ::extractSize(values);

    if(m_is_interleaved)
    {
      m_interleaved = ::ArrayView_from_Node<PointType>(values["x"], m_npts);
    }
    else
    {
      const int dim = ::extractDimension(values);
      SLIC_ASSERT(dim == NDIMS);

      m_strided =
        StridedPoints {{static_cast<CoordType*>(values["x"].data_ptr()),
                        dim >= 2 ? static_cast<CoordType*>(values["y"].data_ptr()) : nullptr,
                        dim >= 3 ? static_cast<CoordType*>(values["z"].data_ptr()) : nullptr}};
    }
  }

  /// Constructor from a multi-component array Sidre group node
  explicit InterleavedOrStridedPoints(const sidre::Group* values)
  {
    conduit::Node vals;
    SLIC_ASSERT(values != nullptr);
    values->createNativeLayout(vals);

    SLIC_ASSERT(isMultiComponentArray(vals));
    m_is_interleaved = conduit::blueprint::mcarray::is_interleaved(vals);
    m_npts = ::extractSize(vals);

    if(m_is_interleaved)
    {
      m_interleaved = ::ArrayView_from_Node<PointType>(vals["x"], m_npts);
    }
    else
    {
      const int dim = ::extractDimension(vals);
      SLIC_ASSERT(dim == NDIMS);

      m_strided =
        StridedPoints {{static_cast<CoordType*>(vals["x"].data_ptr()),
                        dim >= 2 ? static_cast<CoordType*>(vals["y"].data_ptr()) : nullptr,
                        dim >= 3 ? static_cast<CoordType*>(vals["z"].data_ptr()) : nullptr}};
    }
  }

  /// Returns the number of points in the array
  int size() const { return m_npts; }

  /// Access the point at index \a idx in the array
  PointType operator[](int idx) const
  {
    return m_is_interleaved ? m_interleaved[idx] : m_strided[idx];
  }

private:
  /// Predicate to check that a Conduit node is a valid mcarray
  /// and print some debug information if it is not
  bool isMultiComponentArray(conduit::Node& node) const
  {
    conduit::Node info;
    if(!conduit::blueprint::verify("mcarray", node, info))
    {
      SLIC_INFO("Input was not a valid multicomponent array: " << info.to_yaml());
      return false;
    }
    return true;
  }

private:
  StridedPoints m_strided;
  InterleavedPoints m_interleaved;
  bool m_is_interleaved;
  int m_npts;
};

}  // namespace detail

/**
 * \brief A class to perform scattered data interpolation at arbitrary points
 * over an input point set
 *
 * The class uses linear interpolation over a Delaunay triangulation of the point set.
 */
template <int NDIMS = 3>
class ScatteredInterpolation
{
public:
  static constexpr int DIM = NDIMS;
  using DelaunayTriangulation = Delaunay<DIM>;
  using PointType = typename DelaunayTriangulation::PointType;
  using BoundingBoxType = typename DelaunayTriangulation::BoundingBox;
  using CoordType = typename PointType::CoordType;

private:
  using MortonIndexType = std::uint64_t;

  using VertexSet = typename DelaunayTriangulation::IAMeshType::VertexSet;
  using VertexIndirectionSet =
    slam::ArrayIndirectionSet<typename VertexSet::PositionType, axom::IndexType>;

private:
  /**
   *  \brief Helper struct for sorting input points in the Biased Randomized Incremental Order (BRIO)
   *
   *  BRIO helps improve worst-case performance on poorly ordered point sets.
   *  It was introduced in the following paper:
   *    N. Amenta, S. Choi, and G. Rote. "Incremental constructions con BRIO."
   *    Proceedings of the 19th annual symposium on Computational geometry, 2003.
   */
  struct BrioComparator
  {
    axom::IndexType m_index;
    int m_level;
    MortonIndexType m_morton;

    BrioComparator(IndexType index, int level, MortonIndexType morton)
      : m_index(index)
      , m_level(level)
      , m_morton(morton)
    { }

    friend bool operator<(const BrioComparator& lhs, const BrioComparator& rhs)
    {
      return (lhs.m_level == rhs.m_level) ? lhs.m_morton < rhs.m_morton : lhs.m_level < rhs.m_level;
    }
  };

  /**
   * \brief Generates a permutation of [0, pts.size()) following BRIO 
   * 
   * \sa BrioComparator
   */
  template <typename PointArray>
  axom::Array<axom::IndexType> computeInsertionOrder(const PointArray& pts, const BoundingBoxType& bb)
  {
    // This function will compute a permutation of pts following BRIO.
    // Each point gets a level from the computeLevel() lambda
    // and a quantized Morton index from a rectangular lattice over the bounding box

    const int npts = pts.size();
    const int nlevels = axom::utilities::ceil(axom::utilities::log2<CoordType>(npts));

    // Each point has a 50% chance of being at the max level; of the remaining points
    // from the previous level, there's a 50% chance of being at the current level.
    // Any remaining points are at level 0.
    //
    // Implementation note:
    // The original implementation used repeated calls to `random_real()` to
    // simulate coin flips. At large N, this becomes costly (millions of calls
    // into `std::uniform_real_distribution`). We can generate the same level
    // distribution using a single 64-bit random integer per point:
    // - Take the top `nlevels` bits.
    // - The BRIO level is the position of the highest set bit (1..nlevels), or 0 if all are 0.
    std::mt19937_64 mt;
    std::uint64_t seed = 0;
    if(::getScatteredInterpSeed(seed))
    {
      // Use a dimension/size-specific offset so the BRIO ordering remains
      // reproducible without sharing the example's point-generation stream.
      mt.seed(seed ^
              (0x9e3779b97f4a7c15ULL + static_cast<std::uint64_t>(DIM) +
               (static_cast<std::uint64_t>(npts) << 8)));
    }
    else
    {
      std::random_device rd;
      mt.seed(rd());
    }

    auto computeLevel = [&mt, nlevels]() -> int {
      constexpr int RNG_BITS = 64;
      AXOM_STATIC_ASSERT_MSG(std::numeric_limits<std::uint64_t>::digits == RNG_BITS,
                             "Expected 64-bit RNG output");

      const int used_bits = axom::utilities::min(nlevels, RNG_BITS);
      std::uint64_t bits = mt();
      if(used_bits < RNG_BITS)
      {
        bits >>= (RNG_BITS - used_bits);
      }

      if(bits == 0)
      {
        return 0;
      }

      // level in [1, used_bits]
      const int level = RNG_BITS - __builtin_clzll(bits);
      return axom::utilities::min(level, nlevels);
    };

    // We use a Morton index, quantized over the mesh bounding box to
    // order the points on each level
    using QuantizedCoordType = std::uint32_t;
    using MortonizerType = spin::Mortonizer<QuantizedCoordType, MortonIndexType, DIM>;

    // Fit as many bits as possible per dimension into an int64, i.e. floor(63/DIM)
    constexpr QuantizedCoordType shift_bits = (DIM == 2) ? 31 : 21;
    NumericArray<QuantizedCoordType, DIM> res(static_cast<QuantizedCoordType>(1) << shift_bits, DIM);
    auto quantizer =
      spin::rectangular_lattice_from_bounding_box<DIM, CoordType, QuantizedCoordType>(bb, res);

    // Phase 1: compute keys and count points per level (for a counting-sort by level).
    std::vector<std::uint8_t> levels(static_cast<std::size_t>(npts));
    std::vector<std::pair<MortonIndexType, axom::IndexType>> keys(static_cast<std::size_t>(npts));
    std::vector<axom::IndexType> level_counts(static_cast<std::size_t>(nlevels + 1), 0);

    for(int idx = 0; idx < npts; ++idx)
    {
      const int level = computeLevel();
      levels[static_cast<std::size_t>(idx)] = static_cast<std::uint8_t>(level);
      ++level_counts[static_cast<std::size_t>(level)];

      keys[static_cast<std::size_t>(idx)] = {MortonizerType::mortonize(quantizer.gridCell(pts[idx])),
                                             idx};
    }

    // Phase 2: counting-sort keys by level into a single contiguous array in level order.
    std::vector<axom::IndexType> level_offsets(static_cast<std::size_t>(nlevels + 1), 0);
    for(int level = 1; level <= nlevels; ++level)
    {
      level_offsets[static_cast<std::size_t>(level)] =
        level_offsets[static_cast<std::size_t>(level - 1)] +
        level_counts[static_cast<std::size_t>(level - 1)];
    }
    std::vector<axom::IndexType> level_write = level_offsets;

    std::vector<std::pair<MortonIndexType, axom::IndexType>> bucketed(static_cast<std::size_t>(npts));
    for(int idx = 0; idx < npts; ++idx)
    {
      const int level = static_cast<int>(levels[static_cast<std::size_t>(idx)]);
      const axom::IndexType out_pos = level_write[static_cast<std::size_t>(level)]++;
      bucketed[static_cast<std::size_t>(out_pos)] = keys[static_cast<std::size_t>(idx)];
    }

    // Phase 3: sort within each level by Morton index.
    for(int level = 0; level <= nlevels; ++level)
    {
      const axom::IndexType begin = level_offsets[static_cast<std::size_t>(level)];
      const axom::IndexType end = begin + level_counts[static_cast<std::size_t>(level)];
      if(end - begin <= 1)
      {
        continue;
      }

      auto first = bucketed.begin() + begin;
      auto last = bucketed.begin() + end;
      std::sort(first, last, [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });
    }

    // Extract and return the reordered point indices.
    axom::Array<axom::IndexType> reordered(npts, npts);
    for(int idx = 0; idx < npts; ++idx)
    {
      reordered[idx] = bucketed[static_cast<std::size_t>(idx)].second;
    }

    return reordered;
  }

public:
  /**
   * \brief Builds a Delaunay triangulation over the point set from \a mesh_node
   *
   * \param [in] mesh_node Conduit node for the input mesh
   * \param [in] coordset The name of the coordinate set for the input mesh
   */
  void buildTriangulation(conduit::Node& mesh_node, const std::string& coordset)
  {
    // Perform some simple error checking
    SLIC_ASSERT(::isValidBlueprint(mesh_node));

    // Extract coordinates as ArrayView of PointType
    const auto valuesPath = fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(mesh_node.has_path(valuesPath));
    auto coords = detail::InterleavedOrStridedPoints<CoordType, DIM>(mesh_node[valuesPath]);
    const int npts = coords.size();

    // Compute the bounding box
    m_bounding_box.clear();
    for(int i = 0; i < npts; ++i)
    {
      m_bounding_box.addPoint(coords[i]);
    }

    // Reorder the points according to the Biased Random Insertion Order (BRIO) algorithm
    // and store the mapping since we'll need to apply it during interpolation
    const bool report_timing = std::getenv("AXOM_SCATTERED_INTERP_TIMING") != nullptr;
    m_delaunay.setCollectPointLocationStats(report_timing);
    axom::utilities::Timer phase_timer(false);
    if(report_timing)
    {
      phase_timer.start();
    }

    m_brio_data = computeInsertionOrder(coords, m_bounding_box);

    if(report_timing)
    {
      phase_timer.stop();
      SLIC_INFO(axom::fmt::format("ScatteredInterpolation BRIO ordering: {:.6f} sec",
                                  phase_timer.elapsedTimeInSec()));
    }
    m_brio =
      VertexIndirectionSet(typename VertexIndirectionSet::SetBuilder().size(npts).data(&m_brio_data));

    // Scale the Delaunay bounding box to ensure that all input points are contained
    BoundingBoxType bb = m_bounding_box;
    bb.scale(1.5);

    m_delaunay.initializeBoundary(bb);
    m_delaunay.reserveForPointCount(npts);

    if(report_timing)
    {
      phase_timer.reset();
      phase_timer.start();
    }
    for(int i = 0; i < npts; ++i)
    {
      m_delaunay.insertPoint(coords[m_brio[i]]);
    }
    if(report_timing)
    {
      phase_timer.stop();
      SLIC_INFO(axom::fmt::format("ScatteredInterpolation Delaunay insertion: {:.6f} sec",
                                  phase_timer.elapsedTimeInSec()));
      const auto stats = m_delaunay.getInsertionStats();
      SLIC_INFO(axom::fmt::format(
        "ScatteredInterpolation Delaunay cavity removals: mean {:.2f}, max {}, over {} insertions",
        stats.mean_removed(),
        stats.max_removed,
        stats.insertions));
      const auto location_stats = m_delaunay.getPointLocationStats();
      SLIC_INFO(axom::fmt::format(
        "ScatteredInterpolation point location: walks {}, mean steps {:.2f}, max steps {}, "
        "linear fallbacks {}, empty-seed fallbacks {}",
        location_stats.walk_calls,
        location_stats.mean_walk_steps(),
        location_stats.max_walk_steps,
        location_stats.linear_fallbacks,
        location_stats.empty_seed_fallbacks));
    }

    m_delaunay.removeBoundary();
  }

  /**
   * \brief Locates cell from Delaunay complex containing each point in \a query_mesh
   *
   * \param [inout] query_node Conduit node for the query points in mesh Blueprint format;
   * results will be placed into the `cell_idx` field
   * \param [in] coordset The name of the coordinate set for the query mesh
   *
   * \pre query_mesh is the root of a valid mesh blueprint with an unstructured
   * coordinate set \a coordset and a scalar field named `cell_idx` to store the results
   * \note Uses `Delaunay::INVALID_INDEX` for points that cannot be located within the mesh
   */
  void locatePoints(conduit::Node& query_mesh, const std::string& coordset)
  {
    // Perform some simple error checking
    SLIC_ASSERT(::isValidBlueprint(query_mesh));

    const auto valuesPath = fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(query_mesh.has_path(valuesPath));
    auto coords = detail::InterleavedOrStridedPoints<CoordType, DIM>(query_mesh[valuesPath]);
    const int npts = coords.size();

    SLIC_ERROR_IF(!query_mesh.has_path("fields/cell_idx/values"),
                  "Query mesh for ScatteredInterpolation::locatePoints() is "
                  "missing required 'cell_idx' field");

    auto cell_idx =
      ::ArrayView_from_Node<axom::IndexType>(query_mesh["fields/cell_idx/values"], npts);

    // we expect that some points will be outside the mesh
    constexpr bool warnOnInvalid = false;
    for(int idx = 0; idx < npts; ++idx)
    {
      cell_idx[idx] = m_delaunay.findContainingElement(coords[idx], warnOnInvalid);
    }
  }

  /**
   * \brief Given a location in space, find the associated indices
   * and interpolation weights with respect to the input mesh points
   *
   * \param [in]  query_pt The point at which we want to interpolate
   * \param [out] indices The indices of the points from the input mesh in the support of \a query_pt
   * \param [out] weights The interpolation weights associated with each input point in \a indices
   *
   * \returns true if \a query_pt is found within a cell of the Delaunay complex; false otherwise.
   *  If true, the associated indices from points in the input mesh are returned in \a indices
   *  and the interpolation weights for each point are returned in \a weights
   */
  bool getInterpolationWeights(const PointType& query_pt,
                               primal::Point<axom::IndexType, NDIMS + 1>& indices,
                               primal::Point<CoordType, NDIMS + 1>& weights) const
  {
    constexpr bool warnOnInvalid = false;
    constexpr auto INVALID_INDEX = DelaunayTriangulation::INVALID_INDEX;

    const auto cell_id = m_delaunay.findContainingElement(query_pt, warnOnInvalid);

    if(cell_id != INVALID_INDEX)
    {
      // apply BRIO mapping to input vertex indices to match Delaunay insertion order
      const auto verts = m_delaunay.getMeshData()->boundaryVertices(cell_id);
      for(auto idx : verts.positions())
      {
        indices[idx] = m_brio[verts[idx]];
      }
      weights = m_delaunay.getBaryCoords(cell_id, query_pt);

      return true;
    }

    // Cell not found
    return false;
  }

  /**
   * \brief Interpolates a field from an \a input_mesh to one on a \a query_mesh
   *
   * \param [inout] query_mesh Root node of mesh (in blueprint format) containing field to generate
   * \param [in] coordset The name of the coords for query points in \a query_mesh
   * \param [in] input_mesh Root node of mesh (in blueprint format) containing the field to interpolate
   * \param [in] input_field_name Name of field on \a input_mesh
   * \param [in] output_field_name Name of field on \a query_mesh
   * \param [in] INVALID_VALUE Value to use for points that are not in the \a input_mesh
   *
   * \pre \a input_mesh and \a query_mesh must conform to the mesh blueprint schema for point meshes
   * \pre \a input_mesh must contain a nodal scalar field named \a input_field_name
   * \pre \a query_mesh must contain a nodal scalar field named \a output_field_name
   * \pre \a query_mesh must contain a nodal scalar field named \a cell_idx whose values
   * are the index of the cell from the Delaunay triangulation containing the query points.
   * These can be computed in the \a locatePoints() function
   * \post Scalar field \a output_field_name on \a query_mesh will contain the interpolated
   * values of \a input_field_name from \a input_mesh for all query points that have a valid
   * \a cell_idx to a cell in the \a input_mesh.
   */
  void interpolateField(conduit::Node& query_mesh,
                        const std::string& coordset,
                        conduit::Node& input_mesh,
                        const std::string& input_field_name,
                        const std::string& output_field_name,
                        const double INVALID_VALUE = axom::numeric_limits<double>::quiet_NaN())
  {
    constexpr auto INVALID_INDEX = DelaunayTriangulation::INVALID_INDEX;

    SLIC_ASSERT(::isValidBlueprint(query_mesh));
    SLIC_ASSERT(::isValidBlueprint(input_mesh));

    // Extract the required fields from the input and query meshes
    auto in_fld = ::getField<double>(input_mesh, input_field_name);
    auto out_fld = ::getField<double>(query_mesh, output_field_name);
    auto containing_cell = ::getField<axom::IndexType>(query_mesh, "cell_idx");

    const auto valuesPath = fmt::format("coordsets/{}/values", coordset);
    SLIC_ASSERT(query_mesh.has_path(valuesPath));
    auto coords = detail::InterleavedOrStridedPoints<CoordType, DIM>(query_mesh[valuesPath]);

    // Interpolate field at query points
    const int npts = coords.size();
    for(int idx = 0; idx < npts; ++idx)
    {
      const auto cell_id = containing_cell[idx];
      if(cell_id == INVALID_INDEX)
      {
        out_fld[idx] = INVALID_VALUE;
      }
      else
      {
        double res = 0.;
        const auto baryCoords = m_delaunay.getBaryCoords(cell_id, coords[idx]);
        const auto verts = m_delaunay.getMeshData()->boundaryVertices(cell_id);
        for(auto it = verts.begin(); it < verts.end(); ++it)
        {
          // apply BRIO mapping to input vertex indices to match Delaunay insertion order
          res += in_fld[m_brio[*it]] * baryCoords[it.index()];
        }
        out_fld[idx] = res;
      }
    }
  }

  /**
   * \brief Exports the Delaunay complex with scalar fields as a vtk file
   * 
   * \param [in] mesh_node Conduit node for the input mesh
   * \param [in] filename The name of the output file
   * 
   * \note Currently only includes scalar fields of type float64
   */
  void exportDelaunayComplex(conduit::Node& mesh_node, std::string&& filename) const
  {
    constexpr auto CELL_TYPE = DIM == 2 ? mint::TRIANGLE : mint::TET;
    mint::UnstructuredMesh<mint::SINGLE_SHAPE> mint_mesh(DIM, CELL_TYPE);

    const auto* iaMesh = m_delaunay.getMeshData();

    // Add the mesh vertices
    for(auto v : iaMesh->vertices().positions())
    {
      mint_mesh.appendNodes(iaMesh->getVertexPosition(v).data(), 1);
    }

    // Add the mesh cells
    for(auto e : iaMesh->elements().positions())
    {
      mint_mesh.appendCell(&(iaMesh->boundaryVertices(e)[0]), CELL_TYPE);
    }

    // Add each of the scalar fields from mesh_node to the Delaunay complex
    // Note: Currently only adds fields that have type float64
    const int num_verts = iaMesh->vertices().size();
    auto itr = mesh_node["fields"].children();
    while(itr.has_next())
    {
      auto& node = itr.next();
      const auto& fieldName = node.name();

      if(node.has_child("values") && node["values"].dtype().number_of_elements() == num_verts)
      {
        SLIC_DEBUG(fmt::format("Processing field '{}' of type '{}' for Delaunay vtk file",
                               fieldName,
                               conduit::DataType::id_to_name(node["values"].dtype().id())));

        if(node["values"].dtype().is_float64())
        {
          double* vals = node["values"].as_float64_ptr();
          auto* fld = mint_mesh.createField<double>(fieldName, mint::NODE_CENTERED);

          for(auto idx : m_brio.positions())
          {
            // apply BRIO mapping to input vertex indices to match Delaunay insertion order
            fld[idx] = vals[m_brio[idx]];
          }
        }
      }
    }

    mint::write_vtk(&mint_mesh, filename);
  }

  /// Returns the number of vertices in the underlying Deluanay complex
  int numVertices() const { return m_delaunay.getMeshData()->vertices().size(); }

  /// Returns the number of simplices (triangles/tetrahedra) in the underlying Deluanay complex
  int numSimplices() const { return m_delaunay.getMeshData()->elements().size(); }

  /// Returns the bounding box of the input data points
  const BoundingBoxType& boundingBox() const { return m_bounding_box; }

private:
  DelaunayTriangulation m_delaunay;

  axom::Array<axom::IndexType> m_brio_data;
  VertexIndirectionSet m_brio;
  BoundingBoxType m_bounding_box;
};

template <int NDIMS>
constexpr int ScatteredInterpolation<NDIMS>::DIM;

}  // namespace quest
}  // namespace axom

#endif  // QUEST_SCATTERED_INTERPOLATION_H_
