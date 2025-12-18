// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_TOPOLOGY_MAPPER_HPP_
#define AXOM_BUMP_TOPOLOGY_MAPPER_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/spin.hpp"

#include "axom/bump/PrimalAdaptor.hpp"
#include "axom/bump/VariableShape.hpp"
#include "axom/bump/utilities/utilities.hpp"

#include <conduit.hpp>
#include <conduit_relay.hpp>

#include <iostream>

// Uncomment to emit debugging messages
//#define AXOM_DEBUG_TOPOLOGY_MAPPER

namespace axom
{
namespace bump
{
namespace detail
{

//------------------------------------------------------------------------------
/*!
 * \brief Return area where 2 polygons overlap.
 * \param shape1 The subject polygon.
 * \param shape2 The clip polygon.
 * \return The area common to 2 polygons.
 */
AXOM_SUPPRESS_HD_WARN
template <typename T, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> &shape1,
                                     const axom::primal::Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> &shape2,
                                     double eps = 1.e-10)
{
  constexpr bool tryFixOrientation = false;
  const auto p = axom::primal::clip(shape1, shape2, eps, tryFixOrientation);
  return p.area();
}

// We define various shapeOverlap methods to handle
// Tetrahedron, Hexahedron, Polyhedron, PolyhedralFaces shapes.

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
///@{

// Tetrahedron first
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Tetrahedron<T, 3> &shape1,
                                     const axom::primal::Tetrahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Tetrahedron<T, 3> &shape1,
                                     const axom::primal::Hexahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Tetrahedron<T, 3> &shape1,
                                     const axom::primal::Polyhedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Tetrahedron<T, 3> &shape1,
                                     const axom::bump::PolyhedralFaces<T> &shape2,
                                     double eps = 1.e-10)
{
  const bool tryFixOrientation = false;
  auto clipped = axom::primal::Polyhedron<T, 3>::from_primitive(shape1, tryFixOrientation);
  axom::primal::detail::clipPolyhedron(clipped, shape2.getFaces(), eps);
  return clipped.volume();
}

// Hexahedron first
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Hexahedron<T, 3> &shape1,
                                     const axom::primal::Tetrahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Hexahedron<T, 3> &shape1,
                                     const axom::primal::Hexahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Hexahedron<T, 3> &shape1,
                                     const axom::primal::Polyhedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Hexahedron<T, 3> &shape1,
                                     const axom::bump::PolyhedralFaces<T> &shape2,
                                     double eps = 1.e-10)
{
  const bool tryFixOrientation = false;
  auto clipped = axom::primal::Polyhedron<T, 3>::from_primitive(shape1, tryFixOrientation);
  axom::primal::detail::clipPolyhedron(clipped, shape2.getFaces(), eps);
  return clipped.volume();
}

// Polyhedron first
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Polyhedron<T, 3> &shape1,
                                     const axom::primal::Tetrahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Polyhedron<T, 3> &shape1,
                                     const axom::primal::Hexahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Polyhedron<T, 3> &shape1,
                                     const axom::primal::Polyhedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Polyhedron<T, 3> &shape1,
                                     const axom::bump::PolyhedralFaces<T> &shape2,
                                     double eps = 1.e-10)
{
  auto clipped = shape1;
  axom::primal::detail::clipPolyhedron(clipped, shape2.getFaces(), eps);
  return clipped.volume();
}

// PolyhedralFaces first
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::bump::PolyhedralFaces<T> &shape1,
                                     const axom::primal::Tetrahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  return shapeOverlap(shape2, shape1, eps);
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::bump::PolyhedralFaces<T> &shape1,
                                     const axom::primal::Hexahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  return shapeOverlap(shape2, shape1, eps);
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::bump::PolyhedralFaces<T> &shape1,
                                     const axom::primal::Polyhedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  return shapeOverlap(shape2, shape1, eps);
}

template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::bump::PolyhedralFaces<T> &shape1,
                                     const axom::bump::PolyhedralFaces<T> &shape2,
                                     double eps = 1.e-10)
{
  using PointType = axom::primal::Point<T, 3>;
  // Find largest plane offset.
  T maxOffset {};
  for(const auto &plane : shape1.getFaces())
  {
    maxOffset = axom::utilities::max(maxOffset, axom::utilities::abs(plane.getOffset()));
  }
  for(const auto &plane : shape2.getFaces())
  {
    maxOffset = axom::utilities::max(maxOffset, axom::utilities::abs(plane.getOffset()));
  }
  // maxOffset is a radius from the origin. Make it bigger so we can use it as the bounds of a hex.
  maxOffset *= T {2};

  // Make a hex to clip
  axom::primal::Hexahedron<T, 3> hex(PointType {-maxOffset, -maxOffset, -maxOffset},
                                     PointType {maxOffset, -maxOffset, -maxOffset},
                                     PointType {maxOffset, maxOffset, -maxOffset},
                                     PointType {-maxOffset, maxOffset, -maxOffset},
                                     PointType {-maxOffset, -maxOffset, maxOffset},
                                     PointType {maxOffset, -maxOffset, maxOffset},
                                     PointType {maxOffset, maxOffset, maxOffset},
                                     PointType {-maxOffset, maxOffset, maxOffset});
  // Turn the large hex polyhedral and then clip by all planes.
  auto clipped = axom::primal::Polyhedron<T, 3>::from_primitive(hex);
  axom::primal::detail::clipPolyhedron(clipped, shape1.getFaces(), eps);
  axom::primal::detail::clipPolyhedron(clipped, shape2.getFaces(), eps);
  return clipped.volume();
}

///@}

// VariableShape helpers

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T, typename Shape2Type>
AXOM_HOST_DEVICE double shapeOverlap(const VariableShape<T, 3> &shape1,
                                     const Shape2Type &shape2,
                                     double eps = 1.e-10)
{
  const int id = shape1.id();
  double retval = 0.;
  if(id == axom::bump::views::Tet_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tet(shape1[0], shape1[1], shape1[2], shape1[3]);
    retval = shapeOverlap(tet, shape2, eps);
  }
  else if(id == axom::bump::views::Pyramid_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[2];
    shape1.splitPyramid(tets);
    for(int i = 0; i < 2; i++)
    {
      retval += shapeOverlap(tets[i], shape2, eps);
    }
  }
  else if(id == axom::bump::views::Wedge_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[3];
    shape1.splitWedge(tets);
    for(int i = 0; i < 3; i++)
    {
      const auto vol = shapeOverlap(tets[i], shape2, eps);
      retval += vol;
    }
  }
  else if(id == axom::bump::views::Hex_ShapeID)
  {
    axom::primal::Hexahedron<T, 3>
      hex(shape1[0], shape1[1], shape1[2], shape1[3], shape1[4], shape1[5], shape1[6], shape1[7]);
    retval = shapeOverlap(hex, shape2, eps);
  }
  else
  {
    SLIC_ASSERT("Unsupported shape type");
  }
  return retval;
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T, typename Shape1Type>
AXOM_HOST_DEVICE double shapeOverlap(const Shape1Type &shape1,
                                     const VariableShape<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const int id = shape2.id();
  double retval = 0.;
  if(id == axom::bump::views::Tet_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tet(shape2[0], shape2[1], shape2[2], shape2[3]);
    retval = shapeOverlap(shape1, tet, eps);
  }
  else if(id == axom::bump::views::Pyramid_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[2];
    shape2.splitPyramid(tets);
    for(int i = 0; i < 2; i++)
    {
      retval += shapeOverlap(shape1, tets[i], eps);
    }
  }
  else if(id == axom::bump::views::Wedge_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[3];
    shape2.splitWedge(tets);
    for(int i = 0; i < 3; i++)
    {
      const auto vol = shapeOverlap(shape1, tets[i], eps);
      retval += vol;
    }
  }
  else if(id == axom::bump::views::Hex_ShapeID)
  {
    axom::primal::Hexahedron<T, 3>
      hex(shape2[0], shape2[1], shape2[2], shape2[3], shape2[4], shape2[5], shape2[6], shape2[7]);
    retval = shapeOverlap(shape1, hex, eps);
  }
  else
  {
    SLIC_ASSERT("Unsupported shape type");
  }
  return retval;
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const VariableShape<T, 3> &shape1,
                                     const VariableShape<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  int id = shape1.id();
  double retval = 0.;
  if(id == axom::bump::views::Tet_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tet(shape2[0], shape2[1], shape2[2], shape2[3]);
    retval = shapeOverlap(shape1, tet, eps);
  }
  else if(id == axom::bump::views::Pyramid_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[2];
    shape2.splitPyramid(tets);
    for(int i = 0; i < 2; i++)
    {
      retval += shapeOverlap(shape1, tets[i], eps);
    }
  }
  else if(id == axom::bump::views::Wedge_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[3];
    shape2.splitWedge(tets);
    for(int i = 0; i < 3; i++)
    {
      const auto vol = shapeOverlap(shape1, tets[i], eps);
      retval += vol;
    }
  }
  else if(id == axom::bump::views::Hex_ShapeID)
  {
    axom::primal::Hexahedron<T, 3>
      hex(shape2[0], shape2[1], shape2[2], shape2[3], shape2[4], shape2[5], shape2[6], shape2[7]);
    retval = shapeOverlap(shape1, hex, eps);
  }
  else
  {
    SLIC_ASSERT("Unsupported shape type");
  }
  return retval;
}

}  // end namespace detail

//------------------------------------------------------------------------------
/**
 * \brief Take a source topology with clean matset and a target topology and
 *        intersect the source and target zones to build up a new matset on
 *        the target mesh that represents how the source mesh overlaps the
 *        target mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will execute.
 * \tparam SrcTopologyView The view type for the source topology.
 * \tparam SrcCoordsetView The view type for the source coordset.
 * \tparam SrcMatsetView The view type for the source matset.
 * \tparam TargetTopologyView The view type for the target topology.
 * \tparam TargetCoordsetView The view type for the target coordset.
 * \tparam MAX_VERTS_2D The maximum number of vertices allowed in a polygon. This
 *                      value is used only in 2D.
 * \tparam makeFaces Make faces instead of proper Polyhedron zones when polyhedra
 *                   are involved. This enables faster conversion between Blueprint
 *                   and Axom since making planes is less complicated than Axom's
 *                   Polyhedron.
 *
 * \note The use of topology and coordset views as template parameters allows
 *       this class to be instantiated for use with various topology and coordset
 *       types. We also use the source matsetview so we can make an output matset
 *       with the same types.
 */
template <typename ExecSpace,
          typename SrcTopologyView,
          typename SrcCoordsetView,
          typename SrcMatsetView,
          typename TargetTopologyView,
          typename TargetCoordsetView,
          int MAX_VERTS_2D = 12,
          bool makeFaces = true>
class TopologyMapper
{
public:
  static_assert(SrcCoordsetView::dimension() == TargetCoordsetView::dimension(),
                "coordset dimension mismatch");

  using SrcShapeView = PrimalAdaptor<SrcTopologyView, SrcCoordsetView, MAX_VERTS_2D, makeFaces>;
  using TargetShapeView = PrimalAdaptor<TargetTopologyView, TargetCoordsetView, MAX_VERTS_2D>;

  /**
   * \brief Constructor
   *
   * \param srcTopoView The source topology view.
   * \param srcCoordsetView The source coordset view.
   * \param srcMatsetView The source matset view.
   * \param targetTopoView The target topology view.
   * \param targetCoordsetView The target coordset view.
   */
  TopologyMapper(const SrcTopologyView &srcTopoView,
                 const SrcCoordsetView &srcCoordsetView,
                 const SrcMatsetView &srcMatsetView,
                 const TargetTopologyView &targetTopoView,
                 const TargetCoordsetView &targetCoordsetView)
    : m_srcView({srcTopoView, srcCoordsetView})
    , m_srcMatsetView(srcMatsetView)
    , m_targetView({targetTopoView, targetCoordsetView})
  { }

  /**
   * \brief Intersect the source and target topologies and map the source
   *        material onto the target mesh.
   *
   * \param n_srcMesh The Conduit node that contains the coordset, topology, matset for the source mesh.
   * \param n_options A Conduit node that contains the algorithm options. A copy will be made in the appropriate memory space.
   * \param n_targetMesh The node that contains the coordset and topology for the source mesh.
   *
   * \verbatim
   * The n_options node must contain the following keys.
   *
   * source/matsetName: "name of source mesh matset"
   * target/topologyName: "name of target mesh topology"
   * target/matsetName "name of new target mesh matset"
   *
   * The n_options node may contain these optional parameters.
   *
   * source/selectedZones: [0,1,2,...]  # List of selected zone ids in source mesh.
   * target/selectedZones: [0,1,2,...]  # List of selected zone ids in target mesh.
   *
   * \endverbatim
   *
   * \note After executing, the n_targetMesh node will contain a new matset containing
   *       the results of the intersections with the src/target meshes.
   */
  void execute(const conduit::Node &n_srcMesh,
               const conduit::Node &n_options,
               conduit::Node &n_targetMesh) const
  {
    AXOM_ANNOTATE_SCOPE("TopologyMapper::execute");
    namespace utils = axom::bump::utilities;

    // Pick output matset types (use input types)
    using MatIntType = typename SrcMatsetView::IndexType;
    using MatFloatType = typename SrcMatsetView::FloatType;

    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    const char *SRC_MATSET_NAME = "source/matsetName";
    const char *SRC_SELECTED_ZONES = "source/selectedZones";
    const char *TARGET_TOPOLOGY_NAME = "target/topologyName";
    const char *TARGET_MATSET_NAME = "target/matsetName";
    const char *TARGET_SELECTED_ZONES = "target/selectedZones";

    // Make sure options are in the right memory space in case we are given lists of
    // selected zone ids.
    conduit::Node n_options_copy;
    utils::copy<ExecSpace>(n_options_copy, n_options);

    // Ensure required options exist.
    const char *required[] = {SRC_MATSET_NAME, TARGET_TOPOLOGY_NAME, TARGET_MATSET_NAME};
    for(const auto &key : required)
    {
      if(!n_options_copy.has_path(key))
      {
        SLIC_ERROR(axom::fmt::format("Key \"{}\" missing from options.", key));
        return;
      }
    }
    const std::string srcMatsetName = n_options_copy[SRC_MATSET_NAME].as_string();
    const std::string targetTopologyName = n_options_copy[TARGET_TOPOLOGY_NAME].as_string();
    const std::string targetMatsetName = n_options_copy[TARGET_MATSET_NAME].as_string();

    // Look at the source mesh's matset. Count the number of materials.
    const conduit::Node &n_matset = n_srcMesh.fetch_existing("matsets/" + srcMatsetName);
    const conduit::Node &n_materialMap = n_matset.fetch_existing("material_map");
    const auto nmats = n_materialMap.number_of_children();
    const auto numMaterialSlots = nmats + 1;  // leave space for empty material.

    // Look through the material map to get the largest matId in use.
    int maxMatId = 0;
    for(conduit::index_t i = 0; i < nmats; i++)
    {
      int matId = n_materialMap[i].to_int();
      if(i == 0 || matId > maxMatId)
      {
        maxMatId = matId;
      }
    }
    const auto MaterialEmpty = static_cast<MatIntType>(maxMatId + 1);

    // -------------------------------------------------------------------------
    // Build up BVH that contains the src polygon bounding boxes.
    using SrcBoundingBox = typename SrcShapeView::BoundingBox;
    using src_value_type = typename SrcCoordsetView::value_type;
    AXOM_ANNOTATE_BEGIN("bbox");
    const auto srcView = m_srcView;
    SelectedZones<ExecSpace> srcSelection(srcView.numberOfZones(), n_options_copy, SRC_SELECTED_ZONES);
    srcSelection.setSorted(false);
    const auto srcSelectionView = srcSelection.view();
    const axom::IndexType nSrcZones = srcSelectionView.size();
    axom::Array<SrcBoundingBox> srcBoundingBoxes(nSrcZones, nSrcZones, allocatorID);
    auto srcBoundingBoxesView = srcBoundingBoxes.view();
    axom::for_all<ExecSpace>(
      nSrcZones,
      AXOM_LAMBDA(axom::IndexType index) {
        const auto zi = srcSelectionView[index];
        srcBoundingBoxesView[index] = srcView.getBoundingBox(zi);
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
        std::cout << "source zone " << zi << ": bbox=" << srcBoundingBoxesView[index] << std::endl;
#endif
      });
    AXOM_ANNOTATE_END("bbox");

    // -------------------------------------------------------------------------
    // Build the list of selected zones in the target.
    AXOM_ANNOTATE_BEGIN("target");
    const auto targetView = m_targetView;
    const auto nTargetZones = targetView.numberOfZones();
    SelectedZones<ExecSpace> targetSelection(nTargetZones, n_options_copy, TARGET_SELECTED_ZONES);
    targetSelection.setSorted(false);
    const auto targetSelectionView = targetSelection.view();
    AXOM_ANNOTATE_END("target");

    // -------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("build");
    axom::spin::BVH<SrcCoordsetView::dimension(), ExecSpace, src_value_type> bvh;
    bvh.setAllocatorID(allocatorID);
    bvh.initialize(srcBoundingBoxesView, srcBoundingBoxesView.size());
    AXOM_ANNOTATE_END("build");

    // -------------------------------------------------------------------------
    // Set up storage for a new matset.
    AXOM_ANNOTATE_BEGIN("allocation");
    utils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Make target matset.
    conduit::Node &n_targetMatset = n_targetMesh["matsets/" + targetMatsetName];
    n_targetMatset["material_map"].set(n_materialMap);
    n_targetMatset["topology"].set(targetTopologyName);

    conduit::Node &n_volume_fractions = n_targetMatset["volume_fractions"];
    conduit::Node &n_material_ids = n_targetMatset["material_ids"];
    conduit::Node &n_indices = n_targetMatset["indices"];
    conduit::Node &n_sizes = n_targetMatset["sizes"];
    conduit::Node &n_offsets = n_targetMatset["offsets"];

    // Allocate memory for the output matset.
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set_allocator(c2a.getConduitAllocatorID());

    n_volume_fractions.set(
      conduit::DataType(utils::cpp2conduit<MatFloatType>::id, numMaterialSlots * nTargetZones));
    n_material_ids.set(
      conduit::DataType(utils::cpp2conduit<MatIntType>::id, numMaterialSlots * nTargetZones));
    n_sizes.set(conduit::DataType(utils::cpp2conduit<MatIntType>::id, nTargetZones));
    n_offsets.set(conduit::DataType(utils::cpp2conduit<MatIntType>::id, nTargetZones));
    // n_indices are allocated later

    // Wrap the output matset data in some array views.
    auto material_ids = utils::make_array_view<MatIntType>(n_material_ids);
    auto volume_fractions = utils::make_array_view<MatFloatType>(n_volume_fractions);
    auto sizes = utils::make_array_view<MatIntType>(n_sizes);
    auto offsets = utils::make_array_view<MatIntType>(n_offsets);
    volume_fractions.fill(MatFloatType(0.));
    material_ids.fill(MaterialEmpty);
    sizes.fill(MatIntType(0));
    AXOM_ANNOTATE_END("allocation");

    // -------------------------------------------------------------------------
    // Iterate over the target zones and intersect them with source zones.
    AXOM_ANNOTATE_BEGIN("intersection");
    const SrcMatsetView srcMatsetView(m_srcMatsetView);
    const auto bvh_device = bvh.getTraverser();
    axom::for_all<ExecSpace>(
      targetSelectionView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        SLIC_ASSERT(index >= 0 && index < targetSelectionView.size());

        // Get the target zone as a primal shape.
        const axom::IndexType zi = targetSelectionView[index];
        SLIC_ASSERT(zi >= 0 && zi < targetView.numberOfZones());

        const auto targetBBox = targetView.getBoundingBox(zi);
        const auto targetShape = targetView.getShape(zi);
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
        std::cout << "-------------------------------\ntarget zone " << zi << ": " << targetShape
                  << ", bbox=" << targetBBox << std::endl;
#endif
        // Get the area or volume of the target shape (depends on the dimension).
        double targetAmount =
          utils::ComputeShapeAmount<TargetCoordsetView::dimension()>::execute(targetShape);

        // Handle intersection in-depth of the bounding boxes intersected.
        auto handleIntersection = [&](std::int32_t currentNode, const std::int32_t *leafNodes) {
          const auto srcBboxIndex = leafNodes[currentNode];
          SLIC_ASSERT(srcBboxIndex >= 0 && srcBboxIndex < srcSelectionView.size());

          const auto srcZone = srcSelectionView[srcBboxIndex];
          SLIC_ASSERT(srcZone >= 0 && srcZone < srcView.numberOfZones());
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
          std::cout << "handleIntersection: targetZone=" << zi << ", srcZone=" << srcZone
                    << std::endl;
#endif
          // Get the current zone as a primal shape.
          const auto srcShape = srcView.getShape(srcZone);

          // Determine the overlap of the src and target shapes.
          constexpr double eps = 1.e-6;
          const double srcOverlapsTarget = detail::shapeOverlap(srcShape, targetShape, eps);

          if(srcOverlapsTarget > 0.)
          {
            MatFloatType vf = srcOverlapsTarget / targetAmount;

            // Get the src material - there should just be one because we assume
            // that a clean matset is being mapped.
            auto zoneMat = srcMatsetView.beginZone(srcZone);
            SLIC_ASSERT(zoneMat.size() == 1);
            const auto mat = zoneMat.material_id();

#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
            std::cout << "\tintersection:" << std::endl
                      << "\t\ttargetShape=" << targetShape << std::endl
                      << "\t\tsrcShape=" << srcShape << std::endl
                      << "\t\tmat=" << mat << std::endl
                      << "\t\tsrcOverlapsTarget=" << srcOverlapsTarget << std::endl
                      << "\t\ttargetAmount=" << targetAmount << std::endl
                      << "\t\tvf=" << vf << std::endl;
#endif

            // Add the src material contribution into the target material.
            MatIntType *matids = material_ids.data() + zi * numMaterialSlots;
            MatFloatType *vfs = volume_fractions.data() + zi * numMaterialSlots;
            for(int m = 0; m < nmats; m++)
            {
              if(matids[m] == mat)
              {
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
                std::cout << "\t\tAdded " << vf << " to slot " << m << std::endl;
#endif
                vfs[m] += vf;
                break;
              }
              else if(matids[m] == MaterialEmpty)
              {
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
                std::cout << "\t\tAdded new slot " << m << ", mat=" << mat << ", vf=" << vf
                          << std::endl;
#endif
                matids[m] = mat;
                vfs[m] = vf;
                sizes[zi]++;
                break;
              }
            }
          }
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
          else
          {
            std::cout << "\tno intersection" << std::endl;
          }
#endif
        };

        // This predicate determines whether 2 bboxes intersect.
        auto bbIsect = [](const SrcBoundingBox &queryBbox, const SrcBoundingBox &bvhBbox) -> bool {
          bool rv = queryBbox.intersectsWith(bvhBbox);
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
          std::cout << "bbIsect: rv=" << rv << ", q=" << queryBbox << ", bvh=" << bvhBbox
                    << std::endl;
#endif
          return rv;
        };
        // Traverse BVH, looking for bboxes that intersect the current target bbox.
        bvh_device.traverse_tree(targetBBox, handleIntersection, bbIsect);
      });  // axom::for_all
    AXOM_ANNOTATE_END("intersection");

    // -------------------------------------------------------------------------
    // Make a pass through the new material to see if any zones have no materials
    // of if they are partially covered. If there are such zones, make sure the
    // empty material is used.
    axom::IndexType totalSize = 0;
    {
      AXOM_ANNOTATE_SCOPE("sizes");
      axom::ReduceSum<ExecSpace, int> reduceSize(0), emptyCount(0);
      axom::for_all<ExecSpace>(
        nTargetZones,
        AXOM_LAMBDA(axom::IndexType index) {
          // Sum the material within the zone.
          MatFloatType *vfs = volume_fractions.data() + index * numMaterialSlots;
          MatFloatType vfSum(0);
          for(MatIntType m = 0; m < sizes[index]; m++)
          {
            vfSum += vfs[m];
          }
          // If the zone was not completely covered by other materials, increment
          // its size to include the empty material and set its VF.
          constexpr MatFloatType MatTolerance = 1.e-6;
          if(sizes[index] == 0 || (1. - vfSum) > MatTolerance)
          {
            vfs[sizes[index]] = 1. - vfSum;
            sizes[index]++;
            emptyCount += 1;
          }
          reduceSize += sizes[index];
        });
      if(emptyCount.get() > 0)
      {
        // Add an empty material entry to the material map in case some of the slots did
        // not get covered.
        n_targetMatset["material_map"]["empty"] = MaterialEmpty;
      }
      totalSize = reduceSize.get();
    }

    // -------------------------------------------------------------------------
    // All the contributions have been added to the target matset. Finish building it.
    AXOM_ANNOTATE_BEGIN("finish");
    axom::exclusive_scan<ExecSpace>(sizes, offsets);
    n_indices.set(conduit::DataType(utils::cpp2conduit<MatIntType>::id, totalSize));
    auto indices = utils::make_array_view<MatIntType>(n_indices);

    // The volume_fractions and material_ids arrays contain gaps that we can compress out.
    conduit::Node n_new_volume_fractions, n_new_material_ids;
    n_new_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_new_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_new_volume_fractions.set(conduit::DataType(utils::cpp2conduit<MatFloatType>::id, totalSize));
    n_new_material_ids.set(conduit::DataType(utils::cpp2conduit<MatIntType>::id, totalSize));
    auto new_volume_fractions = utils::make_array_view<MatFloatType>(n_new_volume_fractions);
    auto new_material_ids = utils::make_array_view<MatIntType>(n_new_material_ids);
    axom::for_all<ExecSpace>(
      nTargetZones,
      AXOM_LAMBDA(axom::IndexType index) {
        const auto destOffset = offsets[index];
        for(MatIntType m = 0; m < sizes[index]; m++)
        {
          const auto destIndex = destOffset + m;
          const auto srcIndex = index * numMaterialSlots + m;
          new_volume_fractions[destIndex] = volume_fractions[srcIndex];
          new_material_ids[destIndex] = material_ids[srcIndex];
          indices[destIndex] = destIndex;
        }
      });
    // Move the reorganized data into the output.
    n_volume_fractions.move(n_new_volume_fractions);
    n_material_ids.move(n_new_material_ids);

    AXOM_ANNOTATE_END("finish");
  }

  SrcShapeView m_srcView;
  SrcMatsetView m_srcMatsetView;
  TargetShapeView m_targetView;
};

}  // namespace bump
}  // namespace axom

#endif
