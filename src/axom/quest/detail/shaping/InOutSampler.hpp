// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file InOutSampler.hpp
 *
 * \brief Helper class for sampling-based shaping queries using the InOutOctree
 */

#ifndef AXOM_QUEST_INOUT_SAMPLER__HPP_
#define AXOM_QUEST_INOUT_SAMPLER__HPP_

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
class InOutSampler
{
public:
  static constexpr int DIM = NDIMS;
  using InOutOctreeType = quest::InOutOctree<DIM>;

  using GeometricBoundingBox = typename InOutOctreeType::GeometricBoundingBox;
  using SpacePt = typename InOutOctreeType::SpacePt;
  using SpaceVector = typename InOutOctreeType::SpaceVector;
  using GridPt = typename InOutOctreeType::GridPt;
  using BlockIndex = typename InOutOctreeType::BlockIndex;

public:
  /*!
   * \brief Constructor for a InOutSampler
   *
   * \param shapeName The name of the shape; will be used for the field for the associated samples
   * \param surfaceMesh Pointer to the surface mesh
   *
   * \note Does not take ownership of the surface mesh
   */
  InOutSampler(const std::string& shapeName, std::shared_ptr<mint::Mesh> surfaceMesh)
    : m_shapeName(shapeName)
    , m_surfaceMesh(surfaceMesh)
  { }

  ~InOutSampler() { delete m_octree; }

  std::shared_ptr<mint::Mesh> getSurfaceMesh() const { return m_surfaceMesh; }

  /// Computes the bounding box of the surface mesh
  void computeBounds()
  {
    AXOM_ANNOTATE_SCOPE("compute bounding box");
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    m_bbox.clear();
    SpacePt pt;

    for(int i = 0; i < m_surfaceMesh->getNumberOfNodes(); ++i)
    {
      m_surfaceMesh->getNode(i, pt.data());
      m_bbox.addPoint(pt);
    }

    SLIC_ASSERT(m_bbox.isValid());

    SLIC_INFO_ROOT("Mesh bounding box: " << m_bbox);
  }

  void initSpatialIndex(double vertexWeldThreshold)
  {
    AXOM_ANNOTATE_SCOPE("generate InOutOctree");
    // Create octree over mesh's bounding box
    m_octree = new InOutOctreeType(m_bbox, m_surfaceMesh);
    m_octree->setVertexWeldThreshold(vertexWeldThreshold);
    m_octree->generateIndex();
  }

  /*!
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
    using PointType = primal::Point<double, DIM>;

    const InOutOctreeType* octree = m_octree;
    auto checkInside = [=](const PointType& pt) -> bool { return octree->within(pt); };
    shaping::sampleInOutField<FromDim, ToDim>(m_shapeName,
                                              dc,
                                              inoutQFuncs,
                                              sampleRes,
                                              checkInside,
                                              projector);
  }

  /*!
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

  /*!
   * Compute "baseline" volume fractions by sampling at grid function degrees of freedom
   * (instead of at quadrature points)
   */
  template <int FromDim, int ToDim = DIM>
  std::enable_if_t<ToDim == DIM, void> computeVolumeFractionsBaseline(
    mfem::DataCollection* dc,
    int sampleRes,
    int outputOrder,
    PointProjector<FromDim, ToDim> projector = {})
  {
    using PointType = primal::Point<double, DIM>;
    const InOutOctreeType* octree = m_octree;
    auto checkInside = [=](const PointType& pt) -> bool { return octree->within(pt); };
    shaping::computeVolumeFractionsBaseline<FromDim, ToDim>(m_shapeName,
                                                            dc,
                                                            sampleRes,
                                                            outputOrder,
                                                            checkInside,
                                                            projector);
  }

  /*!
   * \warning Do not call this overload with \a ToDim != \a DIM. The compiler needs it to be
   * defined to support various callback specializations for the \a PointProjector.
   */
  template <int FromDim, int ToDim>
  std::enable_if_t<ToDim != DIM, void> computeVolumeFractionsBaseline(
    mfem::DataCollection* AXOM_UNUSED_PARAM(dc),
    int AXOM_UNUSED_PARAM(sampleRes),
    int AXOM_UNUSED_PARAM(outputOrder),
    PointProjector<FromDim, ToDim> AXOM_UNUSED_PARAM(projector))
  {
    static_assert(ToDim != DIM,
                  "Do not call this function -- it only exists to appease the compiler!"
                  "Projector's return dimension (ToDim), must match class dimension (DIM)");
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(InOutSampler);
  DISABLE_MOVE_AND_ASSIGNMENT(InOutSampler);

  std::string m_shapeName;

  GeometricBoundingBox m_bbox;
  std::shared_ptr<mint::Mesh> m_surfaceMesh {nullptr};
  InOutOctreeType* m_octree {nullptr};
};

}  // namespace shaping
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_INOUT_SAMPLER__HPP_
