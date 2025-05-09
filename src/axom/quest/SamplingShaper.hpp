// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file SamplingShaper.hpp
 *
 * \brief Helper class for sampling-based shaping queries
 */

#ifndef AXOM_QUEST_SAMPLING_SHAPER__HPP_
#define AXOM_QUEST_SAMPLING_SHAPER__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"
#include "axom/klee.hpp"

#ifndef AXOM_USE_MFEM
  #error Shaping functionality requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif

#include "axom/quest/Shaper.hpp"
#include "axom/quest/InOutOctree.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"

#include "mfem.hpp"

#include "axom/fmt.hpp"

#include <functional>

namespace axom
{
namespace quest
{
namespace shaping
{
using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;

enum class VolFracSampling : int
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

/// Alias to function pointer that projects a \a FromDim dimensional input point to
/// a \a ToDim dimensional query point when sampling the InOut field
template <int FromDim, int ToDim>
using PointProjector = std::function<primal::Point<double, ToDim>(primal::Point<double, FromDim>)>;

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
   * \brief Constructor for a SamplingShaper
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
    SLIC_ASSERT(m_surfaceMesh != nullptr);

    m_bbox.clear();
    SpacePt pt;

    for(int i = 0; i < m_surfaceMesh->getNumberOfNodes(); ++i)
    {
      m_surfaceMesh->getNode(i, pt.data());
      m_bbox.addPoint(pt);
    }

    SLIC_ASSERT(m_bbox.isValid());

    SLIC_INFO("Mesh bounding box: " << m_bbox);
  }

  void initSpatialIndex(double vertexWeldThreshold)
  {
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

    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "\t Sampling inout field '{}' took {:.3Lf} seconds "
                                "(@ {:L} queries per second)",
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
  * Compute volume fractions function for shape on a grid of resolution \a gridRes
  * in region defined by bounding box \a queryBounds
  */
  void computeVolumeFractionsBaseline(mfem::DataCollection* dc,
                                      int AXOM_UNUSED_PARAM(sampleRes),
                                      int outputOrder)
  {
    // Step 1 -- generate a QField w/ the spatial coordinates
    mfem::Mesh* mesh = dc->GetMesh();
    const int NE = mesh->GetNE();
    const int dim = mesh->Dimension();

    if(NE < 1)
    {
      SLIC_WARNING("Mesh has no elements!");
      return;
    }

    mfem::L2_FECollection* coll =
      new mfem::L2_FECollection(outputOrder, dim, mfem::BasisType::Positive);
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, coll);
    mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
    volFrac->MakeOwner(coll);
    auto volFracName = axom::fmt::format("vol_frac_{}", m_shapeName);
    dc->RegisterField(volFracName, volFrac);

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
};  // class InOutSampler

}  // end namespace shaping

/*!
 * \brief Concrete class for sample based shaping
 */
class SamplingShaper : public Shaper
{
public:
  SamplingShaper(const klee::ShapeSet& shapeSet, sidre::MFEMSidreDataCollection* dc)
    : Shaper(Shaper::RuntimePolicy::seq,
             axom::policyToDefaultAllocatorID(axom::runtime_policy::Policy::seq),
             shapeSet,
             dc)
  { }

  ~SamplingShaper()
  {
    m_inoutShapeQFuncs.DeleteData(true);
    m_inoutShapeQFuncs.clear();

    m_inoutMaterialQFuncs.DeleteData(true);
    m_inoutMaterialQFuncs.clear();

    m_inoutDofs.DeleteData(true);
    m_inoutDofs.clear();
  }

  //@{
  //!  @name Functions to get and set shaping parameters related to sampling; supplements parameters in base class

  void setSamplingType(shaping::VolFracSampling vfSampling) { m_vfSampling = vfSampling; }

  void setQuadratureOrder(int quadratureOrder) { m_quadratureOrder = quadratureOrder; }

  void setVolumeFractionOrder(int volfracOrder) { m_volfracOrder = volfracOrder; }

  /// Registers a function to project from 2D input points to 2D query points
  void setPointProjector(shaping::PointProjector<2, 2> projector) { m_projector22 = projector; }

  /// Registers a function to project from 3D input points to 2D query points
  void setPointProjector(shaping::PointProjector<3, 2> projector) { m_projector32 = projector; }

  /// Registers a function to project from 2D input points to 3D query points
  void setPointProjector(shaping::PointProjector<2, 3> projector) { m_projector23 = projector; }

  /// Registers a function to project from 3D input points to 3D query points
  void setPointProjector(shaping::PointProjector<3, 3> projector) { m_projector33 = projector; }

  //@}

private:
  klee::Dimensions getShapeDimension() const
  {
    const bool has2D = (m_inoutSampler2D != nullptr);
    const bool has3D = (m_inoutSampler3D != nullptr);
    SLIC_ERROR_IF(!(has2D || has3D), "Shape not initialized");
    SLIC_ERROR_IF(has2D && has3D, "Cannot have concurrent 2D and 3D shapes");

    return has2D ? klee::Dimensions::Two : klee::Dimensions::Three;
  }

public:
  //@{
  //!  @name Functions related to the stages for a given shape

  /// Initializes the spatial index for shaping
  void prepareShapeQuery(klee::Dimensions shapeDimension, const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("prepareShapeQuery");

    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    if(!shape.getGeometry().hasGeometry())
    {
      return;
    }

    SLIC_INFO(axom::fmt::format("{:-^80}", " Generating the octree "));

    const auto& shapeName = shape.getName();

    switch(shapeDimension)
    {
    case klee::Dimensions::Two:
      m_inoutSampler2D = new shaping::InOutSampler<2>(shapeName, m_surfaceMesh);
      m_inoutSampler2D->computeBounds();
      m_inoutSampler2D->initSpatialIndex(this->m_vertexWeldThreshold);
      break;

    case klee::Dimensions::Three:
      m_inoutSampler3D = new shaping::InOutSampler<3>(shapeName, m_surfaceMesh);
      m_inoutSampler3D->computeBounds();
      m_inoutSampler3D->initSpatialIndex(this->m_vertexWeldThreshold);
      break;

    default:
      SLIC_ERROR("Shaping dimension must be 2 or 3, but requested dimension was "
                 << static_cast<int>(shapeDimension));
      break;
    }

    // Check that one of sampling shapers (2D or 3D) is null and the other is not
    SLIC_ASSERT((m_inoutSampler2D == nullptr && m_inoutSampler3D != nullptr) ||
                (m_inoutSampler3D == nullptr && m_inoutSampler2D != nullptr));

    // Output some logging info and dump the mesh
    if(this->isVerbose() && this->getRank() == 0)
    {
      const int nVerts = m_surfaceMesh->getNumberOfNodes();
      const int nCells = m_surfaceMesh->getNumberOfCells();

      SLIC_INFO(axom::fmt::format("After welding, surface mesh has {} vertices  and {} triangles.",
                                  nVerts,
                                  nCells));
      mint::write_vtk(m_surfaceMesh.get(), axom::fmt::format("meldedTriMesh_{}.vtk", shapeName));
    }
  }

  void runShapeQuery(const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("runShapeQuery");

    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    if(!shape.getGeometry().hasGeometry())
    {
      return;
    }

    SLIC_INFO(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format(" Querying the octree for shape '{}'", shape.getName())));

    switch(getShapeDimension())
    {
    case klee::Dimensions::Two:
      runShapeQueryImpl(m_inoutSampler2D);
      break;
    case klee::Dimensions::Three:
      runShapeQueryImpl(m_inoutSampler3D);
      break;
    }
  }

  void applyReplacementRules(const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("applyReplacementRules");

    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    const auto& shapeName = shape.getName();
    const auto& thisMatName = shape.getMaterial();

    SLIC_INFO(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Applying replacement rules for shape '{}'", shapeName)));

    mfem::QuadratureFunction* shapeQFunc = nullptr;

    if(shape.getGeometry().hasGeometry())
    {
      // Get inout qfunc for this shape
      shapeQFunc = m_inoutShapeQFuncs.Get(axom::fmt::format("inout_{}", shapeName));

      SLIC_ASSERT_MSG(shapeQFunc != nullptr,
                      axom::fmt::format("Missing inout samples for shape '{}'", shapeName));
    }
    else
    {
      // No input geometry for the shape, get inout qfunc for associated material
      shapeQFunc = m_inoutMaterialQFuncs.Get(axom::fmt::format("mat_inout_{}", thisMatName));

      SLIC_ASSERT_MSG(shapeQFunc != nullptr,
                      axom::fmt::format("Missing inout samples for material '{}'", thisMatName));
    }

    // Create a copy of the inout samples for this shape
    // Replacements will be applied to this and then copied into our shape's material
    auto* shapeQFuncCopy = new mfem::QuadratureFunction(*shapeQFunc);

    // apply replacement rules to all other materials
    for(auto& otherMatName : m_knownMaterials)
    {
      // We'll handle the current shape's material at the end
      if(otherMatName == thisMatName)
      {
        continue;
      }

      const bool shouldReplace = shape.replaces(otherMatName);
      SLIC_INFO(
        axom::fmt::format("Should we replace material '{}' with shape '{}' of material '{}'? {}",
                          otherMatName,
                          shapeName,
                          thisMatName,
                          shouldReplace ? "yes" : "no"));

      auto* otherMatQFunc =
        m_inoutMaterialQFuncs.Get(axom::fmt::format("mat_inout_{}", otherMatName));
      SLIC_ASSERT_MSG(otherMatQFunc != nullptr,
                      axom::fmt::format("Missing inout samples for material '{}'", otherMatName));

      quest::shaping::replaceMaterial(shapeQFuncCopy, otherMatQFunc, shouldReplace);
    }

    // Get inout qfunc for the current material
    const std::string materialQFuncName = axom::fmt::format("mat_inout_{}", thisMatName);
    if(!m_inoutMaterialQFuncs.Has(materialQFuncName))
    {
      // initialize material from shape inout, the QFunc registry takes ownership
      m_inoutMaterialQFuncs.Register(materialQFuncName, shapeQFuncCopy, true);
    }
    else
    {
      // copy shape data into current material and delete the copy
      auto* matQFunc = m_inoutMaterialQFuncs.Get(materialQFuncName);
      SLIC_ASSERT_MSG(matQFunc != nullptr,
                      axom::fmt::format("Missing inout samples for material '{}'", thisMatName));

      const bool reuseExisting = shape.getGeometry().hasGeometry();
      quest::shaping::copyShapeIntoMaterial(shapeQFuncCopy, matQFunc, reuseExisting);

      delete shapeQFuncCopy;
      shapeQFuncCopy = nullptr;
    }

    m_knownMaterials.insert(thisMatName);
  }

  void finalizeShapeQuery() override
  {
    AXOM_ANNOTATE_SCOPE("finalizeShapeQuery");

    delete m_inoutSampler2D;
    m_inoutSampler2D = nullptr;

    delete m_inoutSampler3D;
    m_inoutSampler3D = nullptr;

    m_surfaceMesh.reset();
  }

  //@}

public:
  /**
   * \brief Import an initial set of material volume fractions before shaping
   *
   * \param [in] initialGridFuncions The input data as a map from material names to grid functions
   * 
   * The imported grid functions are interpolated at quadrature points and registered
   * with the supplied names as material-based quadrature fields
   */
  void importInitialVolumeFractions(const std::map<std::string, mfem::GridFunction*>& initialGridFunctions)
  {
    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    auto* mesh = m_dc->GetMesh();

    // ensure we have a starting quadrature field for the positions
    if(!m_inoutShapeQFuncs.Has("positions"))
    {
      shaping::generatePositionsQFunction(mesh, m_inoutShapeQFuncs, m_quadratureOrder);
    }
    auto* positionsQSpace = m_inoutShapeQFuncs.Get("positions")->GetSpace();

    // Interpolate grid functions at quadrature points & register material quad functions
    // assume all elements have same integration rule
    for(auto& entry : initialGridFunctions)
    {
      const auto& name = entry.first;
      auto* gf = entry.second;

      SLIC_INFO(axom::fmt::format("Importing volume fraction field for '{}' material", name));

      if(gf == nullptr)
      {
        SLIC_WARNING(
          axom::fmt::format("Skipping missing volume fraction field for material '{}'", name));
        continue;
      }

      auto* matQFunc = new mfem::QuadratureFunction(*positionsQSpace);
      const auto& ir = matQFunc->GetSpace()->GetIntRule(0);
      const auto* interp = gf->FESpace()->GetQuadratureInterpolator(ir);
      interp->Values(*gf, *matQFunc);

      const auto matName = axom::fmt::format("mat_inout_{}", name);
      m_inoutMaterialQFuncs.Register(matName, matQFunc, true);
    }
  }

  void adjustVolumeFractions() override
  {
    AXOM_ANNOTATE_SCOPE("adjustVolumeFractions");

    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    for(auto& mat : m_inoutMaterialQFuncs)
    {
      const std::string matName = mat.first;
      SLIC_INFO(axom::fmt::format("Generating volume fraction fields for '{}' material", matName));

      // Sample the InOut field at the mesh quadrature points
      switch(m_vfSampling)
      {
      case shaping::VolFracSampling::SAMPLE_AT_QPTS:
        quest::shaping::computeVolumeFractions(matName, m_dc, m_inoutMaterialQFuncs, m_volfracOrder);
        break;
      case shaping::VolFracSampling::SAMPLE_AT_DOFS:
        /* no-op for now */
        break;
      }
    }
  }

  /// Prints out the names of the registered fields related to shapes and materials
  /// This function is intended to help with debugging
  void printRegisteredFieldNames(const std::string& initialMessage)
  {
    // helper lambda to extract the keys of a map<string,*> as a vector of strings
    auto extractKeys = [](const auto& map) {
      std::vector<std::string> keys;
      for(const auto& kv : map)
      {
        keys.push_back(kv.first);
      }
      return keys;
    };

    axom::fmt::memory_buffer out;

    axom::fmt::format_to(std::back_inserter(out),
                         "List of registered fields in the SamplingShaper {}"
                         "\n\t* Data collection grid funcs: {}"
                         "\n\t* Data collection qfuncs: {}"
                         "\n\t* Known materials: {}",
                         initialMessage,
                         axom::fmt::join(extractKeys(m_dc->GetFieldMap()), ", "),
                         axom::fmt::join(extractKeys(m_dc->GetQFieldMap()), ", "),
                         axom::fmt::join(m_knownMaterials, ", "));

    if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
    {
      axom::fmt::format_to(std::back_inserter(out),
                           "\n\t* Shape qfuncs: {}"
                           "\n\t* Mat qfuncs: {}",
                           axom::fmt::join(extractKeys(m_inoutShapeQFuncs), ", "),
                           axom::fmt::join(extractKeys(m_inoutMaterialQFuncs), ", "));
    }
    else if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_DOFS)
    {
      axom::fmt::format_to(std::back_inserter(out),
                           "\n\t* Shape samples at DOFs: {}",
                           axom::fmt::join(extractKeys(m_inoutDofs), ", "));
    }
    SLIC_INFO(axom::fmt::to_string(out));
  }

private:
  // Handles 2D or 3D shaping, based on the template and associated parameter
  template <typename InOutSamplerType>
  void runShapeQueryImpl(InOutSamplerType* shaper)
  {
    // Sample the InOut field at the mesh quadrature points
    const int meshDim = m_dc->GetMesh()->Dimension();
    switch(m_vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      switch(InOutSamplerType::DIM)
      {
      case 2:
        if(meshDim == 2)
        {
          m_projector22
            ? shaper->template sampleInOutField<2>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_quadratureOrder,
                                                   m_projector22)
            : shaper->template sampleInOutField<2>(m_dc, m_inoutShapeQFuncs, m_quadratureOrder);
        }
        else if(meshDim == 3)
        {
          m_projector32
            ? shaper->template sampleInOutField<3>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_quadratureOrder,
                                                   m_projector32)
            : shaper->template sampleInOutField<3>(m_dc, m_inoutShapeQFuncs, m_quadratureOrder);
        }
        break;
      case 3:
        if(meshDim == 2)
        {
          m_projector23
            ? shaper->template sampleInOutField<2>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_quadratureOrder,
                                                   m_projector23)
            : shaper->template sampleInOutField<2>(m_dc, m_inoutShapeQFuncs, m_quadratureOrder);
        }
        else if(meshDim == 3)
        {
          m_projector33
            ? shaper->template sampleInOutField<3>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_quadratureOrder,
                                                   m_projector33)
            : shaper->template sampleInOutField<3>(m_dc, m_inoutShapeQFuncs, m_quadratureOrder);
        }
        break;
      }
      break;
    case shaping::VolFracSampling::SAMPLE_AT_DOFS:
      shaper->computeVolumeFractionsBaseline(m_dc, m_quadratureOrder, m_volfracOrder);
      break;
    }
  }

private:
  shaping::QFunctionCollection m_inoutShapeQFuncs;
  shaping::QFunctionCollection m_inoutMaterialQFuncs;
  shaping::DenseTensorCollection m_inoutDofs;

  shaping::InOutSampler<2>* m_inoutSampler2D {nullptr};
  shaping::InOutSampler<3>* m_inoutSampler3D {nullptr};

  std::set<std::string> m_knownMaterials;

  shaping::PointProjector<2, 2> m_projector22 {};
  shaping::PointProjector<3, 2> m_projector32 {};
  shaping::PointProjector<2, 3> m_projector23 {};
  shaping::PointProjector<3, 3> m_projector33 {};

  shaping::VolFracSampling m_vfSampling {shaping::VolFracSampling::SAMPLE_AT_QPTS};
  int m_quadratureOrder {5};
  int m_volfracOrder {2};
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SAMPLING_SHAPER__HPP_
