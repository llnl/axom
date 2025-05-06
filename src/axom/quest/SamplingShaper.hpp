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
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"
#include "axom/quest/detail/shaping/InOutSampler.hpp"
#include "axom/quest/detail/shaping/PrimitiveSampler.hpp"

#include "mfem.hpp"
#include "mfem/linalg/dtensor.hpp"

#include "axom/fmt.hpp"

#include <functional>

namespace axom
{
namespace quest
{

/// \brief Concrete class for sample based shaping
class SamplingShaper : public Shaper
{
public:
  SamplingShaper(RuntimePolicy execPolicy,
                 int allocatorId,
                 const klee::ShapeSet& shapeSet,
                 sidre::MFEMSidreDataCollection* dc)
    : Shaper(execPolicy, allocatorId, shapeSet, dc)
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
  void setPointProjector22(shaping::PointProjector<2, 2> projector) { m_projector22 = projector; }

  /// Registers a function to project from 3D input points to 2D query points
  void setPointProjector32(shaping::PointProjector<3, 2> projector) { m_projector32 = projector; }

  /// Registers a function to project from 2D input points to 3D query points
  void setPointProjector23(shaping::PointProjector<2, 3> projector) { m_projector23 = projector; }

  /// Registers a function to project from 3D input points to 3D query points
  void setPointProjector33(shaping::PointProjector<3, 3> projector) { m_projector33 = projector; }

  //@}

private:
  int numSamplersInitialized(int dim) const
  {
    int count = 0;
    switch(dim)
    {
    case 2:
      count += m_inoutSampler2D ? 1 : 0;
      break;
    case 3:
      count += m_inoutSampler3D ? 1 : 0;
      count += m_primitiveSampler3D_seq ? 1 : 0;
      count += m_primitiveSampler3D_omp ? 1 : 0;
      count += m_primitiveSampler3D_cuda ? 1 : 0;
      count += m_primitiveSampler3D_hip ? 1 : 0;
      break;
    default:
      SLIC_ERROR("Invalid dimension " << dim);
      break;
    }
    return count;
  }

  klee::Dimensions getShapeDimension() const
  {
    const int count2D = numSamplersInitialized(2);
    const int count3D = numSamplersInitialized(3);
    SLIC_ERROR_IF(count2D + count3D < 1, "Shape not initialized");
    SLIC_ERROR_IF(count2D > 0 && count3D > 0, "Cannot have concurrent 2D and 3D shapes");
    SLIC_ERROR_IF(count2D > 1, "Cannot have more than one 2D");
    SLIC_ERROR_IF(count3D > 1, "Cannot have more than one 3D");

    return count2D > 0 ? klee::Dimensions::Two : klee::Dimensions::Three;
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

    SLIC_INFO(axom::fmt::format("{:-^80}", " Generating the spatial index "));

    const auto& shapeName = shape.getName();

    // note: ignoring the global shapeDimension for now since it's causing problems
    // reading c2c when the dimension is Three
    AXOM_UNUSED_VAR(shapeDimension);
    if(this->shapeFormat(shape) == "c2c")
    {
      m_inoutSampler2D = std::make_unique<shaping::InOutSampler<2>>(shapeName, m_surfaceMesh);
      m_inoutSampler2D->computeBounds();
      m_inoutSampler2D->initSpatialIndex(this->m_vertexWeldThreshold);
    }
    else if(this->shapeFormat(shape) == "stl")
    {
      m_inoutSampler3D = std::make_unique<shaping::InOutSampler<3>>(shapeName, m_surfaceMesh);
      m_inoutSampler3D->computeBounds();
      m_inoutSampler3D->initSpatialIndex(this->m_vertexWeldThreshold);
    }
    else if(this->shapeFormat(shape) == "proe")
    {
      switch(this->getExecutionPolicy())
      {
      case runtime_policy::Policy::seq:
        m_primitiveSampler3D_seq =
          std::make_unique<shaping::PrimitiveSampler<3, seq_exec>>(shapeName, m_surfaceMesh);
        m_primitiveSampler3D_seq->computeBounds();
        m_primitiveSampler3D_seq->initSpatialIndex();
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case runtime_policy::Policy::omp:
        m_primitiveSampler3D_omp =
          std::make_unique<shaping::PrimitiveSampler<3, omp_exec>>(shapeName, m_surfaceMesh);
        m_primitiveSampler3D_omp->computeBounds();
        m_primitiveSampler3D_omp->initSpatialIndex();
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case runtime_policy::Policy::cuda:
        m_primitiveSampler3D_cuda =
          std::make_unique<shaping::PrimitiveSampler<3, cuda_exec>>(shapeName, m_surfaceMesh);
        m_primitiveSampler3D_cuda->computeBounds();
        m_primitiveSampler3D_cuda->initSpatialIndex();
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case runtime_policy::Policy::hip:
        m_primitiveSampler3D_hip =
          std::make_unique<shaping::PrimitiveSampler<3, hip_exec>>(shapeName, m_surfaceMesh);
        m_primitiveSampler3D_hip->computeBounds();
        m_primitiveSampler3D_hip->initSpatialIndex();
        break;
#endif
      default:
        SLIC_ERROR("Unsupported execution policy for PrimitiveSampler3D");
        break;
      }
    }

    // Check that one of sampling shapers (2D or 3D) is null and the other is not
    SLIC_ASSERT((numSamplersInitialized(2) + numSamplersInitialized(3)) == 1);

    // Output some logging info and dump the mesh
    if(this->isVerbose() && this->getRank() == 0)
    {
      const int nVerts = m_surfaceMesh->getNumberOfNodes();
      const int nCells = m_surfaceMesh->getNumberOfCells();
      SLIC_INFO(axom::fmt::format("After welding, surface mesh has {} vertices  and {} elements.",
                                  nVerts,
                                  nCells));
      mint::write_vtk(m_surfaceMesh.get(), axom::fmt::format("melded_shape_mesh_{}.vtk", shapeName));
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
      axom::fmt::format("{:-^80}", axom::fmt::format(" Querying for shape '{}'", shape.getName())));

    switch(getShapeDimension())
    {
    case klee::Dimensions::Two:
      runShapeQueryImpl(m_inoutSampler2D.get());
      break;
    case klee::Dimensions::Three:
      if(this->shapeFormat(shape) == "stl")
      {
        runShapeQueryImpl(m_inoutSampler3D.get());
      }
      else if(this->shapeFormat(shape) == "proe")
      {
        switch(this->getExecutionPolicy())
        {
        case runtime_policy::Policy::seq:
          runShapeQueryImpl(m_primitiveSampler3D_seq.get());
          break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
        case runtime_policy::Policy::omp:
          runShapeQueryImpl(m_primitiveSampler3D_omp.get());
          break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
        case runtime_policy::Policy::cuda:
          runShapeQueryImpl(m_primitiveSampler3D_cuda.get());
          break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
        case runtime_policy::Policy::hip:
          runShapeQueryImpl(m_primitiveSampler3D_hip.get());
          break;
#endif
        default:
          SLIC_ERROR("Unsupported execution policy for PrimitiveSampler3D");
          break;
        }
      }
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

    m_inoutSampler2D.reset();
    m_inoutSampler3D.reset();
    m_primitiveSampler3D_seq.reset();
    m_primitiveSampler3D_omp.reset();
    m_primitiveSampler3D_cuda.reset();
    m_primitiveSampler3D_hip.reset();

    SLIC_WARNING_IF(
      m_surfaceMesh.use_count() > 1,
      axom::fmt::format(
        "in finalizeShapeQuery -- Surface mesh pointer has {} references -- should be at most 1",
        m_surfaceMesh.use_count()));
    slic::flushStreams();

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
        this->computeVolumeFractionsForMaterial(matName);
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
  // Handles 2D or 3D shaping for InOutSampler, based on the template and associated parameter
  template <int DIM>
  void runShapeQueryImpl(shaping::InOutSampler<DIM>* shaper)
  {
    // Sample the InOut field at the mesh quadrature points
    const int meshDim = m_dc->GetMesh()->Dimension();
    switch(m_vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      switch(DIM)
      {
      case 2:
        if(meshDim == 2)
        {
          shaper->template sampleInOutField<2>(m_dc,
                                               m_inoutShapeQFuncs,
                                               m_quadratureOrder,
                                               m_projector22);
        }
        else if(meshDim == 3)
        {
          shaper->template sampleInOutField<3>(m_dc,
                                               m_inoutShapeQFuncs,
                                               m_quadratureOrder,
                                               m_projector32);
        }
        break;
      case 3:
        if(meshDim == 2)
        {
          shaper->template sampleInOutField<2>(m_dc,
                                               m_inoutShapeQFuncs,
                                               m_quadratureOrder,
                                               m_projector23);
        }
        else if(meshDim == 3)
        {
          shaper->template sampleInOutField<3>(m_dc,
                                               m_inoutShapeQFuncs,
                                               m_quadratureOrder,
                                               m_projector33);
        }
        break;
      }
      break;
    case shaping::VolFracSampling::SAMPLE_AT_DOFS:
      shaper->computeVolumeFractionsBaseline(m_dc, m_quadratureOrder, m_volfracOrder);
      break;
    }
  }

  // Handles 2D or 3D shaping for PrimitiveSampler, based on the template and associated parameter
  template <int DIM, typename ExecSpace>
  void runShapeQueryImpl(shaping::PrimitiveSampler<DIM, ExecSpace>* shaper)
  {
    // Sample the InOut field at the mesh quadrature points
    const int meshDim = m_dc->GetMesh()->Dimension();
    switch(m_vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      switch(DIM)
      {
      case 2:
        SLIC_ERROR("Not implemented yet!");
        break;
      case 3:
        if(meshDim == 2)
        {
          shaper->template sampleInOutField<2>(m_dc,
                                               m_inoutShapeQFuncs,
                                               m_quadratureOrder,
                                               m_projector23);
        }
        else if(meshDim == 3)
        {
          shaper->template sampleInOutField<3>(m_dc,
                                               m_inoutShapeQFuncs,
                                               m_quadratureOrder,
                                               m_projector33);
        }
        break;
      }
      break;
    case shaping::VolFracSampling::SAMPLE_AT_DOFS:
      SLIC_ERROR("Not implemented yet!");
      break;
    }
  }

  /**
   * \brief Compute volume fractions for a given material using its associated quadrature function.
   * 
   * The generated grid function will be registered in the data collection and prefixed by `vol_frac_`
   *
   * \param [in] matField The name of the material
   */
  void computeVolumeFractionsForMaterial(const std::string& matField)
  {
    AXOM_ANNOTATE_SCOPE("computeVolumeFractionsForMaterial");

    // Retrieve the inout samples QFunc
    SLIC_ASSERT(axom::utilities::string::startsWith(matField, "mat_inout_"));
    mfem::QuadratureFunction* inout = m_inoutMaterialQFuncs.Get(matField);

    const int sampleOrder = inout->GetSpace()->GetIntRule(0).GetOrder();
    const int sampleNQ = inout->GetSpace()->GetIntRule(0).GetNPoints();
    const int sampleSZ = inout->GetSpace()->GetSize();
    // print info about sampling on rank 0
    // TODO: mpi reduce this for stats on all ranks
    SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                     "In computeVolumeFractions(): sample order {} | "
                                     "sample num qpts {} |  total samples {:L}",
                                     sampleOrder,
                                     sampleNQ,
                                     sampleSZ));

    // extract some properties from computational mesh
    mfem::Mesh* mesh = m_dc->GetMesh();
    const int dim = mesh->Dimension();
    const int NE = mesh->GetNE();
    SLIC_INFO_ROOT(
      axom::fmt::format(axom::utilities::locale(), "Mesh has dim {} and {:L} elements", dim, NE));

    // Access or create a registered volume fraction grid function from the data collection
    const auto volFracName = axom::fmt::format("vol_frac_{}", matField.substr(10));
    mfem::GridFunction* volFrac = shaping::getOrAllocateL2GridFunction(m_dc,
                                                                       volFracName,
                                                                       m_volfracOrder,
                                                                       dim,
                                                                       mfem::BasisType::Positive);
    const mfem::FiniteElementSpace* fes = volFrac->FESpace();

    // Project QField onto volume fractions field using flux corrected transport (FCT)
    // to keep the range of values between 0 and 1
    axom::utilities::Timer timer(true);
    {
      const auto& ir = inout->GetSpace()->GetIntRule(0);  // assume all elements are the same

      const int dofs = fes->GetFE(0)->GetDof();
      mfem::DenseTensor mass_mat(dofs, dofs, NE);
      mass_mat = 0.;
      {
        AXOM_ANNOTATE_SCOPE("mass integrator assemble");

        const int sz = mass_mat.TotalSize();
        mfem::ConstantCoefficient one_coef(1.0);
        mfem::MassIntegrator mass_integrator(one_coef, &ir);

        // wrap mass_mat data as vector for AssembleEA call
        // note: AssembleEA expects the transpose, but it's ok since mass matrices are symmetric
        mfem::Vector mass_vec;
        mfem::Swap(mass_mat.GetMemory(), mass_vec.GetMemory());
        mass_vec.SetSize(sz);
        mass_integrator.AssembleEA(*fes, mass_vec, false);
        mfem::Swap(mass_mat.GetMemory(), mass_vec.GetMemory());
      }

      // Perform batched LU factorization on the mass tensor
      // Note: This was written against mfem@4.7 and should be updated for mfem@4.8 when we upgrade
      mfem::Array<int> pivots(NE * dofs);
      mfem::DenseTensor mInv = mass_mat;
      {
        AXOM_ANNOTATE_SCOPE("batch lu factor");
        mfem::BatchLUFactor(mInv, pivots);
      }

      mfem::Vector b(fes->GetVSize());
      {
        AXOM_ANNOTATE_SCOPE("domain lf integrator assemble");

        mfem::QuadratureFunctionCoefficient qfc(*inout);
        mfem::DomainLFIntegrator rhs(qfc);
        rhs.SetIntRule(&ir);

        mfem::Array<int> elem_marker(fes->GetNE());
        elem_marker = 1;
        b = 0.;
        rhs.AssembleDevice(*fes, elem_marker, b);
      }

      mfem::Vector one(dofs);
      one = 1.;

      constexpr double minY = 0.;
      constexpr double maxY = 1.;

      mfem::Vector fct_mat(dofs * dofs * NE);
      fct_mat = 0.;

      // Reshape returns an indexable view of a multidimensional array
      const auto m_d = mfem::Reshape(mass_mat.HostRead(), dofs, dofs, NE);
      const auto mInv_d = mfem::Reshape(mInv.HostRead(), dofs, dofs, NE);
      const auto P_d = mfem::Reshape(pivots.HostRead(), dofs, NE);
      const auto b_d = mfem::Reshape(b.HostRead(), dofs, NE);
      const auto one_d = mfem::Reshape(one.HostRead(), dofs);
      auto fct_mat_d = mfem::Reshape(fct_mat.HostReadWrite(), dofs, dofs, NE);
      auto x_d = mfem::Reshape(volFrac->HostWrite(), dofs, NE);

      AXOM_ANNOTATE_BEGIN("fct project");
      axom::for_all<axom::SEQ_EXEC>(0, NE, [&](int i) {
        shaping::FCT_project(&m_d(0, 0, i),
                             &mInv_d(0, 0, i),
                             &P_d(0, i),
                             dofs,
                             &b_d(0, i),
                             &one_d(0),
                             minY,
                             maxY,
                             &x_d(0, i),
                             &fct_mat_d(0, 0, i));
      });
      AXOM_ANNOTATE_END("fct project");
    }
    timer.stop();

    // print stats for root rank
    SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                     "\t Generating volume fractions '{}' took {:.3f} seconds (@ "
                                     "{:L} dofs processed per second)",
                                     volFracName,
                                     timer.elapsed(),
                                     static_cast<int>(fes->GetNDofs() / timer.elapsed())));
  }

private:
  shaping::QFunctionCollection m_inoutShapeQFuncs;
  shaping::QFunctionCollection m_inoutMaterialQFuncs;
  shaping::DenseTensorCollection m_inoutDofs;

  // add pointers to all possible samplers for the various dimensions and execution spaces
  // Note: the omp, cuda and hip pointers can only be instantiated with appropriate axom congigs
  // Note: only one of these can be instantiated within a SamplingShaper
  // TODO: This will be a lot cleaner with a std::variant
  std::unique_ptr<shaping::InOutSampler<2>> m_inoutSampler2D;
  std::unique_ptr<shaping::InOutSampler<3>> m_inoutSampler3D;
  std::unique_ptr<shaping::PrimitiveSampler<3, seq_exec>> m_primitiveSampler3D_seq;
  std::unique_ptr<shaping::PrimitiveSampler<3, omp_exec>> m_primitiveSampler3D_omp;
  std::unique_ptr<shaping::PrimitiveSampler<3, cuda_exec>> m_primitiveSampler3D_cuda;
  std::unique_ptr<shaping::PrimitiveSampler<3, hip_exec>> m_primitiveSampler3D_hip;

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
