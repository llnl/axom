// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file PrimitiveSampler.hpp
 *
 * \brief Helper class for sampling-based shaping queries using primal geometric primitives
 */

#ifndef AXOM_QUEST_PRIMITIVE_SAMPLER__HPP_
#define AXOM_QUEST_PRIMITIVE_SAMPLER__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"
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

template <int NDIMS, typename ExecSpace>
class PrimitiveSampler
{
public:
  static constexpr int DIM = NDIMS;
  using BVHType = spin::BVH<DIM, ExecSpace, double>;
  using GeometricBoundingBox = primal::BoundingBox<double, DIM>;
  using SpacePt = primal::Point<double, DIM>;

  using SimplexType =
    std::conditional_t<NDIMS == 3, primal::Tetrahedron<double, 3>, primal::Triangle<double, 2>>;

  static_assert(NDIMS == 3 || NDIMS == 2, "PrimitiveSampler needs to be 2 or 3");

public:
  /**
    * \brief Constructor for a PrimitiveSampler over a collection of triangles (in 2D) or tetrahedra (in 3D)
    *
    * \param shapeName The name of the shape; will be used for the field for the associated samples
    * \param surfaceMesh Pointer to the surface mesh
    *
    * \note Does not take ownership of the surface mesh
    */
  PrimitiveSampler(const std::string& shapeName, std::shared_ptr<mint::Mesh> surfaceMesh)
    : m_shapeName(shapeName)
    , m_surfaceMesh(surfaceMesh)
  { }

  ~PrimitiveSampler() { }

  std::shared_ptr<mint::Mesh> getSurfaceMesh() const { return m_surfaceMesh; }

  /// Computes the bounding box of the surface mesh
  void computeBounds()
  {
    AXOM_ANNOTATE_SCOPE("compute bounding box");
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    auto* pmesh = dynamic_cast<mint::UnstructuredMesh<mint::SINGLE_SHAPE>*>(m_surfaceMesh.get());
    if(pmesh == nullptr)
    {
      SLIC_ERROR_ROOT(
        axom::fmt::format("Expected a {} mesh", NDIMS == 2 ? "triangle" : "tetrahedral"));
      return;
    }

    m_bbox.clear();

    const double* x = pmesh->getCoordinateArray(mint::X_COORDINATE);
    const double* y = pmesh->getCoordinateArray(mint::Y_COORDINATE);
    const double* z = DIM == 3 ? pmesh->getCoordinateArray(mint::Z_COORDINATE) : nullptr;
    primal::ZipIndexable<SpacePt> verts = {{x, y, z}};

    // compute the overall bounding box
    const int num_nodes = pmesh->getNumberOfNodes();
    for(int i = 0; i < num_nodes; ++i)
    {
      m_bbox.addPoint(verts[i]);
    }
    SLIC_ASSERT(m_bbox.isValid());

    constexpr double EPS = 1e-8;

    // Optimization opportunity -- only consider primitives within the bounding box of this domain

    // extract the primitives and their bounding boxes
    const int num_cells = pmesh->getNumberOfCells();
    m_aabbs.resize(num_cells);
    m_primitives.resize(num_cells);
    for(int i = 0; i < num_cells; ++i)
    {
      const axom::IndexType* connec = pmesh->getCellNodeIDs(i);
      auto& simplex = m_primitives[i];
      for(int j = 0; j < NDIMS + 1; ++j)
      {
        simplex[j] = verts[connec[j]];
      }

      // ensure positive volumes
      // note: this effectively ignores the orientation of simplices
      bool is_degenerate = false;
      const double vol = simplex.signedVolume();
      if(axom::utilities::isNearlyEqual(vol, 0., EPS))
      {
        is_degenerate = true;
      }
      else if(vol < 0)
      {
        axom::utilities::swap(simplex[0], simplex[1]);
      }

      // we will skip degenerate volumes
      // TODO: WE should only consider simplices in the bounding box of the current domain
      if(!is_degenerate)
      {
        m_aabbs[i] = primal::compute_bounding_box(simplex);
      }
    }

    SLIC_INFO_ROOT("Mesh bounding box: " << m_bbox);

#if defined(AXOM_USE_RAJA)
    // Print out the total volume of all the tetrahedra
    auto prim_view = m_primitives.view();
    using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
    RAJA::ReduceSum<REDUCE_POL, double> total_tet_vol(0.0);
    axom::for_all<ExecSpace>(
      num_cells,
      AXOM_LAMBDA(axom::IndexType i) { total_tet_vol += prim_view[i].volume(); });

    SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                     "Total volume of all generated tetrahedra is {:.2Lf}",
                                     total_tet_vol.get()));
#endif  //defined(AXOM_USE_RAJA)
  }

  void initSpatialIndex()
  {
    AXOM_ANNOTATE_SCOPE("generate BVH tree");

    m_bvh.initialize(m_aabbs.view(), m_aabbs.size());
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
    using FromPoint = primal::Point<double, FromDim>;
    using ToPoint = primal::Point<double, ToDim>;
    AXOM_ANNOTATE_SCOPE("sample containment");

    SLIC_ERROR_IF(FromDim != ToDim && !projector,
                  "A projector callback function is required when FromDim != ToDim");

    auto* mesh = dc->GetMesh();
    SLIC_ASSERT(mesh != nullptr);
    //const int NE = mesh->GetNE();
    //const int dim = mesh->Dimension();

    // Generate a Quadrature Function with the geometric positions, if not already available
    if(!inoutQFuncs.Has("positions"))
    {
      shaping::generatePositionsQFunction(mesh, inoutQFuncs, sampleRes);
    }

    // Access the positions QFunc and associated QuadratureSpace
    mfem::QuadratureFunction* pos_coef = inoutQFuncs.Get("positions");
    auto* sp = pos_coef->GetSpace();
    const int nq = sp->GetSize();

    // Sample the in/out field at each point
    // store in QField which we register with the QFunc collection
    const std::string inoutName = axom::fmt::format("inout_{}", m_shapeName);
    const int vdim = 1;
    auto* inout = new mfem::QuadratureFunction(sp, vdim);
    inoutQFuncs.Register(inoutName, inout, true);

    // TODO: For debugging it'd also be useful to see how many candidates each query point has
    // optimization: We might want to first check the bounding boxes of the elements
    //               and only do the tetrahedral containment for the elements w/ candidates

    axom::utilities::Timer timer(true);

    SLIC_INFO_ROOT(axom::fmt::format("{:-^80}", " Finding shape candidates for each mesh element "));

    // Get the positions of the query points, project them if needed
    axom::ArrayView<FromPoint> orig_qpts_v(reinterpret_cast<FromPoint*>(pos_coef->HostReadWrite()),
                                           nq);
    axom::Array<ToPoint> projected_qpts(0);
    if(projector)
    {
      AXOM_ANNOTATE_SCOPE("project query points");
      projected_qpts.resize(nq);
      auto proj_pts_v = projected_qpts.view();
      axom::for_all<ExecSpace>(
        nq,
        AXOM_LAMBDA(axom::IndexType i) { proj_pts_v[i] = projector(orig_qpts_v[i]); });
    }
    // We need to reinterpret_cast since the compiler can't rule out that FromPoint is a different type from ToPoint
    // in the else case, despite our SLIC_ERROR above that checks for this.
    // This should look a lot cleaner w/ `if constexpr` when we move to C++17
    const auto query_view = projector
      ? projected_qpts.view()
      : axom::ArrayView<ToPoint>(reinterpret_cast<ToPoint*>(pos_coef->HostReadWrite()), nq);

    axom::ArrayView<double> inout_view(const_cast<double*>(inout->HostRead()), nq);
    axom::for_all<ExecSpace>(nq, AXOM_LAMBDA(axom::IndexType i) { inout_view[i] = 0.; });

    axom::Array<IndexType> offsets(nq, nq);
    axom::Array<IndexType> counts(nq, nq);
    axom::Array<IndexType> candidates;

    m_bvh.findPoints(offsets.view(), counts.view(), candidates, nq, query_view);

    auto counts_view = counts.view();
    auto offsets_view = offsets.view();
    auto candidates_view = candidates.view();
    auto aabbs_view = m_aabbs.view();
    auto prims_view = m_primitives.view();
    AXOM_UNUSED_VAR(aabbs_view);

    AXOM_ANNOTATE_BEGIN("checking containment");
    axom::for_all<ExecSpace>(
      nq,
      AXOM_LAMBDA(axom::IndexType i) {
        for(int j = 0; j < counts_view[i]; j++)
        {
          const auto shapeIdx = candidates_view[offsets_view[i] + j];

          SLIC_ASSERT(aabbs_view[shapeIdx].scale(1.05).contains(query_view[i]));

          if(prims_view[shapeIdx].contains(query_view[i]))
          {
            inout_view[i] = 1.;
          }
        }
      });
    AXOM_ANNOTATE_END("checking containment");

    timer.stop();

    // print stats for root rank
    SLIC_INFO_ROOT(axom::fmt::format(
      axom::utilities::locale(),
      "\t Sampling inout field '{}' took {:.3Lf} seconds (@ {:L} queries per second)",
      inoutName,
      timer.elapsed(),
      static_cast<int>(nq / timer.elapsed())));
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
   * \warning Not yet implemented
   */
  template <int FromDim, int ToDim = DIM>
  void computeVolumeFractionsBaseline(mfem::DataCollection* AXOM_UNUSED_PARAM(dc),
                                      int AXOM_UNUSED_PARAM(sampleRes),
                                      int AXOM_UNUSED_PARAM(outputOrder),
                                      PointProjector<FromDim, ToDim> AXOM_UNUSED_PARAM(projector))
  {
    AXOM_ANNOTATE_SCOPE("computeVolumeFractionsBaseline");
    SLIC_WARNING_ROOT("computeVolumeFractionsBaseline() not implemented yet");
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(PrimitiveSampler);
  DISABLE_MOVE_AND_ASSIGNMENT(PrimitiveSampler);

  std::string m_shapeName;

  BVHType m_bvh;

  axom::Array<GeometricBoundingBox> m_aabbs;
  axom::Array<SimplexType> m_primitives;

  GeometricBoundingBox m_bbox;
  std::shared_ptr<mint::Mesh> m_surfaceMesh {nullptr};
};

}  // namespace shaping
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_PRIMITIVE_SAMPLER__HPP_
