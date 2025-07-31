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

namespace detail
{

/*!
 * \brief Check whether point \a pt is inside \a shape.
 *
 * \param shape The shape being checked for inside/outside.
 * \param pt The point being checked against the shape.
 *
 * \return True if pt is inside shape; False otherwise.
 */
template <typename ShapeType, typename PointType>
bool AXOM_HOST_DEVICE checkInside(const ShapeType &shape, const PointType &pt)
{
  bool inside = false;
  double wn {};
  for(int c = 0; c < shape.numEdges(); c++)
  {
    wn += axom::primal::winding_number(pt, shape[c]);
  }
  // A point inside the polygon should have non-zero winding number.
  inside |= (std::round(wn) != 0);

  return inside;
}

} // end namespace detail

/*!
 * \brief This class samples a geometry view of shapes against quad points in a supplied MFEM mesh.
 */
template <int NDIMS>
class WindingNumberSampler
{
public:
  static constexpr int DIM = NDIMS;

  // For now.
  using ExecSpace = axom::SEQ_EXEC;

  using GeometryView = typename axom::ArrayView<axom::primal::CurvedPolygon<double, DIM>>;
  using PointType = primal::Point<double, DIM>;
  using GeometricBoundingBox = axom::primal::BoundingBox<double, DIM>;
  using BVH = typename axom::spin::BVH<NDIMS, ExecSpace, double>;
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

  /// Computes the bounding box of the surface mesh. This version does nothing.
  void computeBounds()
  {
    // no-op - We do it in initSpatialIndex.
  }

  /*!
   * \brief Initialize the BVH that is used for queries from the geometry view.
   */
  void initSpatialIndex(double AXOM_UNUSED_PARAM(vertexWeldThreshold))
  {
    AXOM_ANNOTATE_SCOPE("Initialize spatial index");

    // Figure out bounding boxes for each geometric object.
    const axom::IndexType geometrySize = m_geometryView.size();
    axom::Array<GeometricBoundingBox> aabbs(geometrySize, geometrySize, axom::execution_space<ExecSpace>::allocatorID());
    auto aabbsView = aabbs.view();
    const auto geometryView = m_geometryView;
    axom::for_all<ExecSpace>(geometrySize, AXOM_LAMBDA(axom::IndexType i)
    {
      aabbsView[i] = geometryView[i].boundingBox();
    });

    // Initialize the BVH using the bounding boxes.
    m_bvh.initialize(aabbs, aabbs.size());
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
    static_assert(axom::execution_space<ExecSpace>::onDevice() == false,
      "This sampler does not work on GPU yet due to some MFEM usage and lack of support for GPU in CurvedPolygon.");

  using FromPoint = primal::Point<double, FromDim>;
  using ToPoint = primal::Point<double, ToDim>;
  AXOM_ANNOTATE_SCOPE("sampleInOutField");

  SLIC_ERROR_IF(FromDim != ToDim && !projector,
                "A projector callback function is required when FromDim != ToDim");

  auto* mesh = dc->GetMesh();
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  // Generate a Quadrature Function with the geometric positions, if not already available
  if(!inoutQFuncs.Has("positions"))
  {
    shaping::generatePositionsQFunction(mesh, inoutQFuncs, sampleRes);
  }

  // Access the positions QFunc and associated QuadratureSpace
  mfem::QuadratureFunction* pos_coef = inoutQFuncs.Get("positions");
  auto* sp = pos_coef->GetSpace();
  const int nq = sp->GetIntRule(0).GetNPoints();

  // Sample the in/out field at each point
  // store in QField which we register with the QFunc collection
  const std::string inoutName = axom::fmt::format("inout_{}", m_shapeName);
  const int vdim = 1;
  auto* inout = new mfem::QuadratureFunction(sp, vdim);
  inoutQFuncs.Register(inoutName, inout, true);

  // Build an array of query points at the quad points.
  axom::utilities::Timer timer(true);
  AXOM_ANNOTATE_BEGIN("Create query points");
  const int numQueryPoints = NE * nq;
  const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<ToPoint> queryPoints(numQueryPoints, numQueryPoints, allocatorID);
  auto queryPointsView = queryPoints.view();
  axom::for_all<ExecSpace>(NE, AXOM_LAMBDA(axom::IndexType i)
  {
    // FIXME: GPU portability MFEM usage.
    mfem::DenseMatrix m;
    pos_coef->GetValues(i, m);

    int qpi = i * nq;
    if(projector)
    {
      for(int p = 0; p < nq; ++p)
      {
        queryPointsView[qpi++] = projector(FromPoint(m.GetColumn(p), dim));
      }
    }
    else
    {
      for(int p = 0; p < nq; ++p)
      {
        queryPointsView[qpi++] = ToPoint(m.GetColumn(p), dim);
      }
    }
  });
  AXOM_ANNOTATE_END("Create query points");

  // Look up all of the query points. This will allocate the candidates array.
  AXOM_ANNOTATE_BEGIN("findPoints");
  axom::Array<axom::IndexType> offsets(numQueryPoints, numQueryPoints, allocatorID);
  axom::Array<axom::IndexType> sizes(numQueryPoints, numQueryPoints, allocatorID);
  axom::Array<axom::IndexType> candidates;
  auto offsetsView = offsets.view();
  auto sizesView = sizes.view();
  m_bvh.findPoints(offsetsView,
                   sizesView,
                   candidates,
                   numQueryPoints,
                   queryPointsView);
  AXOM_ANNOTATE_END("findPoints");

  // Check each element's quad points for in/out.
  AXOM_ANNOTATE_BEGIN("InOut tests");
  axom::Array<bool> inOutResult(numQueryPoints, numQueryPoints, allocatorID);
  auto inOutResultView = inOutResult.view();
  const auto candidatesView = candidates.view();
  const auto geometryView = m_geometryView;
  axom::for_all<ExecSpace>(numQueryPoints, AXOM_LAMBDA(axom::IndexType qpi)
  {
    // Check whether the current query point is inside candidate shapes.
    bool in = false;
    const auto numCandidates = sizesView[qpi];
    const auto &queryPoint = queryPointsView[qpi];
    for(axom::IndexType ci = 0; ci < numCandidates && in == false; ci++)
    {
      const auto candidateIndex = candidatesView[offsetsView[qpi] + ci];
      in |= detail::checkInside(geometryView[candidateIndex], queryPoint);
    }
    inOutResultView[qpi] = in;
  });

  // Store the results back into the MFEM quad function.
  axom::for_all<ExecSpace>(NE, AXOM_LAMBDA(axom::IndexType i)
  {
    // FIXME: GPU portability MFEM usage.
    mfem::Vector res;
    inout->GetValues(i, res);

    // Check this element i's quad points against candidate shapes.
    axom::IndexType qpi = i * nq;
    for(int p = 0; p < nq; p++, qpi++)
    {
      res(p) = inOutResultView[qpi] ? 1. : 0.;
    }
  });
  AXOM_ANNOTATE_END("InOut tests");
  timer.stop();

  // print stats for rank 0
  SLIC_INFO_ROOT(axom::fmt::format(
    axom::utilities::locale(),
    "\t Sampling inout field '{}' took {:.3Lf} seconds (@ {:L} queries per second)",
    inoutName,
    timer.elapsed(),
    static_cast<int>((NE * nq) / timer.elapsed())));
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
  void computeVolumeFractionsBaseline(mfem::DataCollection* dc, int sampleRes, int outputOrder)
  {
    AXOM_ANNOTATE_SCOPE("computeVolumeFractionsBaseline");
    const auto geometryView = m_geometryView;
    auto checkInside = [=](const PointType& pt) -> bool {
      // TODO: figure out curved polygons that might contain point from index.

      // Check each candidate
      bool inside = false;
      for(axom::IndexType i = 0; i < geometryView.size() && !inside; i++)
      {
        const auto& shapeGeom = geometryView[i];
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
  GeometricBoundingBox m_bbox {};
  GeometryView m_geometryView;
  BVH m_bvh {};
};

}  // namespace shaping
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_WINDING_NUMBER_SAMPLER__HPP_
