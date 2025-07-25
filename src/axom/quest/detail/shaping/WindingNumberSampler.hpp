// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_WINDING_NUMBER_SAMPLER__HPP_
#define AXOM_QUEST_WINDING_NUMBER_SAMPLER__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"
#include "axom/quest/InOutOctree.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"

#include "axom/fmt.hpp"

#include "mfem.hpp"

namespace axom
{
namespace quest
{
namespace shaping
{

using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;

template <int NDIMS>
class WindingNumberSampler
{
public:
  static constexpr int DIM = NDIMS;
#if 0
  using InOutOctreeType = quest::InOutOctree<DIM>;

  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;
  using SpacePt = typename InOutOctreeType::SpacePt;
  using SpaceVector = typename InOutOctreeType::SpaceVector;
  using GridPt = typename InOutOctreeType::GridPt;
  using BlockIndex = typename InOutOctreeType::BlockIndex;
#endif
  using GeometryView = typename axom::ArrayView<axom::primal::CurvedPolygon<double, DIM>>;
  using PointType = primal::Point<double, DIM>;
  using GeometricBoundingBox = axom::primal::BoundingBox<double, DIM>;

public:
  /**
   * \brief Constructor for a WindingNumberSampler
   *
   * \param shapeName The name of the shape; will be used for the field for the associated samples
   * \param geomView A view that contains the shapes being queried.
   *
   */
  WindingNumberSampler(const std::string& shapeName, GeometryView geomView)
    : m_shapeName(shapeName)
    , m_geometryView(geomView)
  { }

  ~WindingNumberSampler() = default;

  /// Computes the bounding box of the surface mesh
  void computeBounds()
  {
    AXOM_ANNOTATE_SCOPE("compute bounding box");

    m_bbox.clear();
    PointType pt;

    for(axom::IndexType i = 0; i < m_geometryView.size(); ++i)
    {
      const auto bbox = m_geometryView[i].boundingBox();
      m_bbox.addPoint(bbox.getMin());
      m_bbox.addPoint(bbox.getMax());
    }

    SLIC_ASSERT(m_bbox.isValid());

    SLIC_INFO_ROOT("Mesh bounding box: " << m_bbox);
  }

  void initSpatialIndex(double AXOM_UNUSED_PARAM(vertexWeldThreshold))
  {
#if 0 // FOR NOW
    AXOM_ANNOTATE_SCOPE("generate InOutOctree");
    // Create octree over mesh's bounding box
    m_octree = new InOutOctreeType(m_bbox, m_surfaceMesh);
    m_octree->setVertexWeldThreshold(vertexWeldThreshold);
    m_octree->generateIndex();
#endif
  }

  /**
   * \brief Samples the inout field over the indexed geometry, possibly using a
   * callback function to project the input points (from the computational mesh)
   * to query points on the spatial index
   * 
   * \tparam FromDim The dimension of points from the input mesh
   * \tparam ToDim The dimension of points on the indexed shape
   * \param [in] dc The data collection containing the mesh and associated query points
   * \param [inout] inoutQFuncs A collection of quadrature functions for the shape and material
   * inout samples
   * \param [in] sampleRes The quadrature order at which to sample the inout field
   * \param [in] projector A callback function to apply to points from the input mesh
   * before querying them on the spatial index
   * 
   * \note A projector callback must be supplied when \a FromDim is not equal 
   * to \a ToDim, the projector
   * \note \a ToDim must be equal to \a DIM, the dimension of the spatial index
   */
  template <int FromDim, int ToDim = DIM>
  std::enable_if_t<ToDim == DIM, void> sampleInOutField(mfem::DataCollection* dc,
                                                        shaping::QFunctionCollection& inoutQFuncs,
                                                        int sampleRes,
                                                        PointProjector<FromDim, ToDim> projector = {})
  {
    const auto geometryView = m_geometryView;
    auto checkInside = [=](const PointType &pt) -> bool
    {
      // TODO: figure out curved polygons that might contain point from index.

      // Check each candidate 
      bool inside = false;
      for(axom::IndexType i = 0; i < geometryView.size() && !inside; i++)
      {
        const auto &shapeGeom = geometryView[i];
        double wn {};
        for(int c = 0; c < shapeGeom.numEdges(); c++)
        {
          wn += axom::primal::winding_number(pt, shapeGeom[c]);
        }
        // A point inside the polygon should have non-zero winding number.
        inside |= (std::round(wn) != 0);
      }
      return inside;
    };
    shaping::sampleInOutField<FromDim, ToDim>(m_shapeName, dc, inoutQFuncs, sampleRes, checkInside, projector);
  }

  /** 
   * \warning Do not call this overload with \a ToDim != \a DIM. The compiler needs it to be
   * defined to support various callback specializations for the \a PointProjector.
   */
  template <int FromDim, int ToDim>
  std::enable_if_t<ToDim != DIM, void> sampleInOutField(mfem::DataCollection*,
                                                        shaping::QFunctionCollection&,
                                                        int,
                                                        PointProjector<FromDim, ToDim>)
  {
    static_assert(ToDim != DIM,
                  "Do not call this function -- it only exists to appease the compiler!"
                  "Projector's return dimension (ToDim), must match class dimension (DIM)");
  }

  /**
   * Compute "baseline" volume fractions by sampling at grid function degrees of freedom
   * (instead of at quadrature points)
  */
  void computeVolumeFractionsBaseline(mfem::DataCollection* dc,
                                      int sampleRes,
                                      int outputOrder)
  {
    const auto geometryView = m_geometryView;
    auto checkInside = [=](const PointType &pt) -> bool
    {
      // TODO: figure out curved polygons that might contain point from index.

      // Check each candidate 
      bool inside = false;
      for(axom::IndexType i = 0; i < geometryView.size() && !inside; i++)
      {
        const auto &shapeGeom = geometryView[i];
        double wn {};
        for(int c = 0; c < shapeGeom.numEdges(); c++)
        {
          wn += axom::primal::winding_number(pt, shapeGeom[c]);
        }
        // A point inside the polygon should have non-zero winding number.
        inside |= (std::round(wn) != 0);
      }
      return inside;
    };
    shaping::computeVolumeFractionsBaseline<DIM>(m_shapeName, dc, sampleRes, outputOrder, checkInside);
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(WindingNumberSampler);
  DISABLE_MOVE_AND_ASSIGNMENT(WindingNumberSampler);

  std::string m_shapeName;
  GeometricBoundingBox m_bbox;
  GeometryView m_geometryView;
};

}  // namespace shaping
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_WINDING_NUMBER_SAMPLER__HPP_
