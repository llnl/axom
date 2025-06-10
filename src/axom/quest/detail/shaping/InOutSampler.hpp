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
  /**
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
    AXOM_ANNOTATE_SCOPE("sample InOutOctree");

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

    mfem::DenseMatrix m;
    mfem::Vector res;

    axom::utilities::Timer timer(true);
    for(int i = 0; i < NE; ++i)
    {
      pos_coef->GetValues(i, m);
      inout->GetValues(i, res);

      if(projector)
      {
        for(int p = 0; p < nq; ++p)
        {
          const ToPoint pt = projector(FromPoint(m.GetColumn(p), dim));
          const bool in = m_octree->within(pt);
          res(p) = in ? 1. : 0.;
        }
      }
      else
      {
        for(int p = 0; p < nq; ++p)
        {
          const ToPoint pt(m.GetColumn(p), dim);
          const bool in = m_octree->within(pt);
          res(p) = in ? 1. : 0.;
        }
      }
    }
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
  void computeVolumeFractionsBaseline(mfem::DataCollection* dc,
                                      int AXOM_UNUSED_PARAM(sampleRes),
                                      int outputOrder)
  {
    AXOM_ANNOTATE_SCOPE("computeVolumeFractionsBaseline");

    // Step 1 -- generate a QField w/ the spatial coordinates
    mfem::Mesh* mesh = dc->GetMesh();
    const int NE = mesh->GetNE();
    const int dim = mesh->Dimension();

    if(NE < 1)
    {
      SLIC_WARNING("Mesh has no elements!");
      return;
    }

    const auto volFracName = axom::fmt::format("vol_frac_{}", m_shapeName);
    mfem::GridFunction* volFrac =
      shaping::getOrAllocateL2GridFunction(dc, volFracName, outputOrder, dim, mfem::BasisType::Positive);
    const mfem::FiniteElementSpace* fes = volFrac->FESpace();

    auto* fe = fes->GetFE(0);
    auto& ir = fe->GetNodes();

    // Assume all elements have the same integration rule
    const int nq = ir.GetNPoints();
    const auto* geomFactors = mesh->GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

    mfem::DenseTensor pos_coef(dim, nq, NE);

    // Rearrange positions into quadrature function
    {
      for(int i = 0; i < NE; ++i)
      {
        for(int j = 0; j < dim; ++j)
        {
          for(int k = 0; k < nq; ++k)
          {
            pos_coef(j, k, i) = geomFactors->X((i * nq * dim) + (j * nq) + k);
          }
        }
      }
    }

    // Step 2 -- sample the in/out field at each point -- store directly in volFrac grid function
    mfem::Vector res(nq);
    mfem::Array<int> dofs;
    for(int i = 0; i < NE; ++i)
    {
      mfem::DenseMatrix& m = pos_coef(i);
      for(int p = 0; p < nq; ++p)
      {
        const SpacePt pt(m.GetColumn(p), dim);
        const bool in = m_octree->within(pt);
        res(p) = in ? 1. : 0.;
      }

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, res);
    }
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
