// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file SamplingShaper.hpp
 *
 * \brief Helper class for sampling-based shaping queries
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/klee.hpp"

#if (!defined(AXOM_USE_MFEM) && !(defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP))) || \
  !defined(AXOM_USE_SIDRE)
  #error SamplingShaper requires Axom to be configured with Sidre and either MFEM or Conduit+Bump
#endif

#include "axom/quest/Shaper.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"
#include "axom/quest/detail/shaping/InOutSampler.hpp"
#include "axom/quest/detail/shaping/PrimitiveSampler.hpp"
#include "axom/quest/detail/shaping/WindingNumberSampler.hpp"
#if defined(AXOM_USE_MFEM)
  #include "axom/quest/io/MFEMReader.hpp"
  #include "mfem.hpp"
  #include "mfem/linalg/dtensor.hpp"
#endif

#include "axom/fmt.hpp"

#include <functional>
#include <numeric>

namespace axom
{
namespace quest
{
/// \brief Concrete class for sample based shaping on MFEM or Blueprint meshes.
class SamplingShaper : public Shaper
{
public:
  /// Struct to help choose sampler method: InOut or WindingNumber.
  enum class SamplingMethod : int
  {
    InOut,
    WindingNumber
  };

private:
  using InOutSampler2D = shaping::InOutSampler<2>;
  using InOutSampler3D = shaping::InOutSampler<3>;
  using PrimitiveSampler3D_seq = shaping::PrimitiveSampler<3, seq_exec>;
  using PrimitiveSampler3D_omp = shaping::PrimitiveSampler<3, omp_exec>;
  using PrimitiveSampler3D_cuda = shaping::PrimitiveSampler<3, cuda_exec>;
  using PrimitiveSampler3D_hip = shaping::PrimitiveSampler<3, hip_exec>;
  using WindingNumberSampler2D = shaping::WindingNumberSampler<2>;

  // Type trait for any InOutSampler type
  template <typename T>
  struct is_inoutsampler
    : std::bool_constant<std::is_same_v<T, InOutSampler2D> || std::is_same_v<T, InOutSampler3D>>
  { };

  template <typename T>
  inline static constexpr bool is_inoutsampler_v = is_inoutsampler<T>::value;

  // Type trait for any WindingNumberSampler type
  template <typename T>
  struct is_wnsampler : std::bool_constant<std::is_same_v<T, WindingNumberSampler2D>>
  { };

  template <typename T>
  inline static constexpr bool is_wnsampler_v = is_wnsampler<T>::value;

  // Type trait for any PrimitiveSampler type
  template <typename T>
  struct is_primitivesampler
    : std::bool_constant<
        std::is_same_v<T, PrimitiveSampler3D_seq> || std::is_same_v<T, PrimitiveSampler3D_omp> ||
        std::is_same_v<T, PrimitiveSampler3D_cuda> || std::is_same_v<T, PrimitiveSampler3D_hip>>
  { };

  template <typename T>
  inline static constexpr bool is_primitivesampler_v = is_primitivesampler<T>::value;

  // Type trait to get the dimension of a sampler
  template <typename T>
  struct sampler_dimension
    : std::integral_constant<int,
                             std::is_same_v<T, InOutSampler2D>              ? 2
                               : std::is_same_v<T, InOutSampler3D>          ? 3
                               : std::is_same_v<T, WindingNumberSampler2D>  ? 2
                               : std::is_same_v<T, PrimitiveSampler3D_seq>  ? 3
                               : std::is_same_v<T, PrimitiveSampler3D_omp>  ? 3
                               : std::is_same_v<T, PrimitiveSampler3D_cuda> ? 3
                               : std::is_same_v<T, PrimitiveSampler3D_hip>  ? 3
                                                                            : 0>
  { };

  template <typename T>
  inline static constexpr int sampler_dimension_v = sampler_dimension<T>::value;

  using SamplerVariant = std::variant<std::monostate,
                                      std::unique_ptr<InOutSampler2D>,
                                      std::unique_ptr<InOutSampler3D>,
                                      std::unique_ptr<WindingNumberSampler2D>,
                                      std::unique_ptr<PrimitiveSampler3D_seq>
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
                                      ,
                                      std::unique_ptr<PrimitiveSampler3D_omp>
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
                                      ,
                                      std::unique_ptr<PrimitiveSampler3D_cuda>
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
                                      ,
                                      std::unique_ptr<PrimitiveSampler3D_hip>
#endif
                                      >;

public:
#if defined(AXOM_USE_MFEM)
  /// MFEM-compatible constructor
  SamplingShaper(RuntimePolicy execPolicy,
                 int allocatorId,
                 const klee::ShapeSet& shapeSet,
                 sidre::MFEMSidreDataCollection* dc)
    : Shaper(execPolicy, allocatorId, shapeSet, dc)
  {
    initializeSamplingResolution();
  }
#endif

#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
  /// Sidre-compatible constructor
  SamplingShaper(RuntimePolicy execPolicy,
                 int allocatorId,
                 const klee::ShapeSet& shapeSet,
                 sidre::Group* bpMesh,
                 const std::string& topo = "")
    : Shaper(execPolicy, allocatorId, shapeSet, bpMesh, topo)
  {
    initializeSamplingResolution();
    m_volfracOrder = 1;
  }

  /// Blueprint-compatible constructor
  SamplingShaper(RuntimePolicy execPolicy,
                 int allocatorId,
                 const klee::ShapeSet& shapeSet,
                 conduit::Node& bpNode,
                 const std::string& topo = "")
    : Shaper(execPolicy, allocatorId, shapeSet, bpNode, topo)
  {
    initializeSamplingResolution();
    m_volfracOrder = 1;
  }
#endif

  ~SamplingShaper() override = default;

  ///@{
  //!  @name Functions to get and set shaping parameters related to sampling; supplements parameters in base class

  void setSamplingType(shaping::VolFracSampling vfSampling) { m_vfSampling = vfSampling; }

  void setSamplingMethod(SamplingMethod samplingMethod) { m_samplingMethod = samplingMethod; }

  /// \brief Controls whether InOutOctree VTK visualization dumps are written during sampling.
  void setInOutOctreeVtkOutputEnabled(bool enabled) { m_inoutOctreeVtkOutputEnabled = enabled; }

  /// \brief Sets the directory for InOutOctree VTK visualization dumps during sampling.
  void setInOutOctreeVtkOutputDirectory(const std::string& directory)
  {
    m_inoutOctreeVtkOutputDirectory = directory;
  }

  /*!
   * \brief Sets the 1D quadrature family used to generate custom sample points.
   *
   * Passing `axom::numerics::QuadratureType::Invalid` selects the default
   * quadrature behavior. Other values request a specific quadrature family.
   * For uniform point sampling over the full zone, including the element
   * edges, `axom::numerics::QuadratureType::ClosedUniform` is often a good
   * choice. Users
   * can experiment with other quadrature families when different sample point
   * patterns are desired.
   *
   * \param [in] qtype Quadrature family selection.
   */
  void setQuadratureType(axom::numerics::QuadratureType qtype);

  /*!
   * \brief Sets an isotropic sampling resolution for custom quadrature.
   *
   * The same positive sample count is used in each logical mesh direction.
   * For custom quadrature families, these values specify the per-direction
   * sample counts directly, which in turn determine the quadrature rule used
   * in each logical direction.
   *
   * \param [in] sampleRes Number of sample points to use per logical
   *                       direction.
   */
  void setSamplingResolution(int sampleRes);

  /*!
   * \brief Sets an anisotropic sampling resolution for custom quadrature.
   *
   * The entries correspond to the logical `I`, `J`, and `K` directions of the
   * reference element. Each entry must be positive. For custom quadrature
   * families, these values specify the per-direction sample counts directly,
   * which in turn determine the quadrature rule used in each logical
   * direction.
   *
   * \param [in] sampleRes ArrayView containing the sample count per logical
   *                       direction. The size needs to match the number of
   *                       mesh dimensions.
   */
  void setSamplingResolution(axom::ArrayView<int> sampleRes);

  /// Deprecated backward compatibility method
  [[deprecated]] void setQuadratureOrder(int order) { setSamplingResolution(order); }

  /*!
   * \brief Set the order for the output volume fractions. This function has no
   *        effect for Blueprint meshes.
   *
   * \param volfracOrder The order for the output volume fractions.
   */
  void setVolumeFractionOrder(int volfracOrder);

  /// Registers a function to project from 2D input points to 2D query points
  void setPointProjector22(shaping::PointProjector<2, 2> projector) { m_projector22 = projector; }

  /// Registers a function to project from 3D input points to 2D query points
  void setPointProjector32(shaping::PointProjector<3, 2> projector) { m_projector32 = projector; }

  /// Registers a function to project from 2D input points to 3D query points
  void setPointProjector23(shaping::PointProjector<2, 3> projector) { m_projector23 = projector; }

  /// Registers a function to project from 3D input points to 3D query points
  void setPointProjector33(shaping::PointProjector<3, 3> projector) { m_projector33 = projector; }

  ///@}

#if defined(AXOM_USE_MFEM)
  // NOTE: These methods are used in tests.

  /// Returns a pointer to the quadrature function associated with shape \a name if it exists, else nullptr
  mfem::QuadratureFunction* getShapeQFunction(const std::string& name) const
  {
    SLIC_ASSERT(m_mfem_state != nullptr);
    return m_mfem_state->shapeQFuncs().Get(name);
  }
  /// Returns a pointer to the quadrature function associated with material \a name if it exists, else nullptr
  mfem::QuadratureFunction* getMaterialQFunction(const std::string& name) const
  {
    SLIC_ASSERT(m_mfem_state != nullptr);
    return m_mfem_state->materialQFuncs().Get(name);
  }
#endif
protected:
  /// Initializes the sampling resolution array based on the mesh dimension.
  void initializeSamplingResolution();

  /*!
   * \brief Verifies the input mesh.
   *
   * \param[out] whyBad A string containing the reason the mesh was bad.
   *
   * \return True if the mesh is ok, false if it is bad. The \a whyBad string is set when false.
   */
  bool verifyInputMeshImpl(std::string& whyBad) const override;

#if defined(AXOM_USE_CONDUIT)
  /*!
   * \brief Save a Blueprint file.
   *
   * \param n_mesh The Blueprint mesh to save.
   * \param filename The name of the file to save.
   */
  void saveBlueprintFile(const conduit::Node& n_mesh, const std::string& filename) const;
#endif

  /*!
   * \brief Saves the sampling quadrature points as a Blueprint point mesh.
   *
   * For MFEM-backed sampling, this converts the `"positions"` quadrature
   * function to a temporary Blueprint point mesh before saving. For
   * Blueprint-backed sampling, this saves the generated quadrature-point
   * topology and any fields associated with it.
   */
  void saveQuadraturePoints(const std::string& filename) const;

  bool hasValidSampler() const { return !std::holds_alternative<std::monostate>(m_sampler); }

  klee::Dimensions getShapeDimension() const
  {
    return std::visit(
      [](const auto& s) {
        using T = std::decay_t<decltype(s)>;
        if constexpr(std::is_same_v<T, std::monostate>)
        {
          return klee::Dimensions::Unspecified;
        }
        else if constexpr(sampler_dimension_v<typename T::element_type> == 2)
        {
          return klee::Dimensions::Two;
        }
        else if constexpr(sampler_dimension_v<typename T::element_type> == 3)
        {
          return klee::Dimensions::Three;
        }
        else
        {
          SLIC_ERROR("Unreachable code reached in getShapeDimension().");
          return klee::Dimensions::Unspecified;
        }
      },
      m_sampler);
  }

  /// Determine whether it is appropriate to use the winding number sampler.
  bool useWindingNumberSampler(const klee::Shape& shape) const
  {
    return this->shapeFormat(shape) == "mfem" &&
      this->m_samplingMethod == SamplingMethod::WindingNumber;
  }

public:
  ///@{
  //!  @name Functions related to the stages for a given shape

  /*!
   * \brief Load the shape geometry. For MFEM files, geometry is loaded into m_contours.
   *        Other formats make discrete geometry and load it into m_surface in the Shaper
   *        base class.
   *
   * \param shape The shape to load.
   */
  void loadShape(const klee::Shape& shape) override;

  /// Initializes the spatial index for shaping
  void prepareShapeQuery(klee::Dimensions shapeDimension, const klee::Shape& shape) override;

  void runShapeQuery(const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("runShapeQuery");

    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    if(!shape.getGeometry().hasGeometry())
    {
      return;
    }

    SLIC_INFO_ROOT(
      axom::fmt::format("{:-^80}", axom::fmt::format(" Querying for shape '{}'", shape.getName())));

    // Impl function allows us to handle different capabilities of each shaper
    std::visit(
      [this](auto& sampler) {
        using T = std::decay_t<decltype(sampler)>;
        if constexpr(!std::is_same_v<T, std::monostate>)
        {
          this->runShapeQueryImpl(sampler.get());
        }
      },
      m_sampler);
  }

  void applyReplacementRules(const klee::Shape& shape) override
  {
    AXOM_ANNOTATE_SCOPE("applyReplacementRules");

    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      applyReplacementRulesImpl(*m_mfem_state, shape);
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
    if(m_bp_state != nullptr)
    {
      applyReplacementRulesImpl(*m_bp_state, shape);
      return;
    }
#endif
    SLIC_ERROR("No mesh state is available for SamplingShaper.");
  }

  void finalizeShapeQuery() override
  {
    AXOM_ANNOTATE_SCOPE("finalizeShapeQuery");

    m_sampler = std::monostate();  // frees memory associated w/ the sampler

    SLIC_WARNING_IF(
      m_surfaceMesh.use_count() > 1,
      axom::fmt::format(
        "in finalizeShapeQuery -- Surface mesh pointer has {} references -- should be at most 1",
        m_surfaceMesh.use_count()));
    slic::flushStreams();

    m_surfaceMesh.reset();
  }

  ///@}

public:
#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
  /**
   * \brief Import an initial set of material volume fractions before shaping
   *
   * \param [in] initialVolumeFractions The input data as a map from material names to fields
   * 
   * The imported fields are interpolated at quadrature points and registered
   * with the supplied names as material-based quadrature fields
   */
  void importInitialVolumeFractions(const std::map<std::string, conduit::Node*>& initialVolumeFractions);
#endif
#if defined(AXOM_USE_MFEM)
  /**
   * \brief Import an initial set of material volume fractions before shaping
   *
   * \param [in] initialGridFuncions The input data as a map from material names to grid functions
   * 
   * The imported grid functions are interpolated at quadrature points and registered
   * with the supplied names as material-based quadrature fields
   */
  void importInitialVolumeFractions(
    const std::map<std::string, mfem::GridFunction*>& initialGridFunctions);
#endif

  /*!
   * \brief Turn the in/out samples into material in/out fields
   */
  void adjustVolumeFractions() override;

  /// Prints out the names of the registered fields related to shapes and materials
  /// This function is intended to help with debugging
  void printRegisteredFieldNames(const std::string& initialMessage);

  /*!
   * \brief Save the shaping results to disk.
   *
   * \param extra Save extra data when available.
   */
  virtual void saveResults(bool extra) override;

private:
  /// Return the mesh dimension.
  int meshDimension() const;

  // Handles 2D or 3D shaping for compatible samplers, based on the template and associated parameter
  template <typename MeshState, typename SamplerType>
  void runShapeQueryImplSampler(MeshState& meshState, SamplerType* sampler)
  {
    // Sample the InOut field at the mesh quadrature points
    if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
    {
      shaping::generateSamplingPositions(meshState, m_samplingResolution.view(), m_quadratureType);
    }

    const int meshDim = meshDimension();
    switch(m_vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      switch(SamplerType::DIM)
      {
      case 2:
        if(meshDim == 2)
        {
          sampler->template sampleInOutField<2, 2>(meshState, m_projector22);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 2>(meshState, m_projector32);
        }
        break;
      case 3:
        if(meshDim == 2)
        {
          sampler->template sampleInOutField<2, 3>(meshState, m_projector23);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 3>(meshState, m_projector33);
        }
        break;
      }
      break;
    case shaping::VolFracSampling::SAMPLE_AT_DOFS:
      switch(SamplerType::DIM)
      {
      case 2:
        if(meshDim == 2)
        {
          sampler->template computeVolumeFractionsBaseline<2, 2>(meshState,
                                                                 m_volfracOrder,
                                                                 m_projector22);
        }
        else if(meshDim == 3)
        {
          sampler->template computeVolumeFractionsBaseline<3, 2>(meshState,
                                                                 m_volfracOrder,
                                                                 m_projector32);
        }
        break;
      case 3:
        if(meshDim == 2)
        {
          sampler->template computeVolumeFractionsBaseline<2, 3>(meshState,
                                                                 m_volfracOrder,
                                                                 m_projector23);
        }
        else if(meshDim == 3)
        {
          sampler->template computeVolumeFractionsBaseline<3, 3>(meshState,
                                                                 m_volfracOrder,
                                                                 m_projector33);
        }
        break;
      }
      break;
    }
  }

  // Handles 2D or 3D shaping for InOutSampler, based on the template and associated parameter
  template <int DIM>
  void runShapeQueryImpl(shaping::InOutSampler<DIM>* sampler)
  {
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      runShapeQueryImplSampler(*m_mfem_state, sampler);
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
    if(m_bp_state != nullptr)
    {
      runShapeQueryImplSampler(*m_bp_state, sampler);
      return;
    }
#endif
    SLIC_ERROR("No mesh state is available for SamplingShaper.");
  }

  // Handles 2D or 3D shaping for WindingNumberSampler, based on the template and associated parameter
  template <int DIM>
  void runShapeQueryImpl(shaping::WindingNumberSampler<DIM>* sampler)
  {
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      runShapeQueryImplSampler(*m_mfem_state, sampler);
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
    if(m_bp_state != nullptr)
    {
      runShapeQueryImplSampler(*m_bp_state, sampler);
      return;
    }
#endif
    SLIC_ERROR("No mesh state is available for SamplingShaper.");
  }

  // Handles 2D or 3D shaping for PrimitiveSampler, based on the template and associated parameter
  template <int DIM, typename ExecSpace>
  void runShapeQueryImpl(shaping::PrimitiveSampler<DIM, ExecSpace>* sampler)
  {
    auto runImpl = [this, sampler](auto& meshState) {
      int meshDim = meshDimension();
      if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
      {
        shaping::generateSamplingPositions(meshState, m_samplingResolution.view(), m_quadratureType);
      }

      // Sample the InOut field at the mesh quadrature points
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
            sampler->template sampleInOutField<2, 3>(meshState, m_projector23);
          }
          else if(meshDim == 3)
          {
            sampler->template sampleInOutField<3, 3>(meshState, m_projector33);
          }
          break;
        }
        break;
      case shaping::VolFracSampling::SAMPLE_AT_DOFS:
        SLIC_ERROR("Not implemented yet!");
        break;
      }
    };

#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      runImpl(*m_mfem_state);
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
    if(m_bp_state != nullptr)
    {
      runImpl(*m_bp_state);
      return;
    }
#endif
    SLIC_ERROR("No mesh state is available for SamplingShaper.");
  }

  /*!
   * \brief Apply replacement rules using for the supplied shape, adjusting functions in \a meshState.
   *
   * \param meshState The object that contains the mesh and fields.
   * \param shape The shape being considered.
   */
  template <typename MeshState>
  void applyReplacementRulesImpl(MeshState& meshState, const klee::Shape& shape)
  {
    const auto& shapeName = shape.getName();
    const auto& thisMatName = shape.getMaterial();

    SLIC_INFO_ROOT(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Applying replacement rules for shape '{}'", shapeName)));

    auto* shapeFunc = shape.getGeometry().hasGeometry()
      ? meshState.getShapeFunction(shaping::shapeInOutFieldName(shapeName))
      : meshState.getMaterialFunction(shaping::materialInOutFieldName(thisMatName));

    if(shape.getGeometry().hasGeometry())
    {
      SLIC_ERROR_IF(shapeFunc == nullptr,
                    axom::fmt::format("Missing inout samples for shape '{}'. "
                                      "This indicates the shape query did not produce a "
                                      "quadrature field before replacement rules were applied.",
                                      shapeName));
    }
    else
    {
      SLIC_ERROR_IF(shapeFunc == nullptr,
                    axom::fmt::format("Missing inout samples for material '{}' while applying "
                                      "replacement rules for shape '{}', which has no input "
                                      "geometry. Initialize that material before shaping, e.g. "
                                      "pass '--background-material {}' in the shaping driver or "
                                      "import initial volume fractions for it.",
                                      thisMatName,
                                      shapeName,
                                      thisMatName));
    }

    auto* shapeFuncCopy = quest::shaping::cloneInOutFunction(shapeFunc);

    for(auto& otherMatName : m_knownMaterials)
    {
      if(otherMatName == thisMatName)
      {
        continue;
      }

      const bool shouldReplace = shape.replaces(otherMatName);
      SLIC_INFO_ROOT(
        axom::fmt::format("Should we replace material '{}' with shape '{}' of material '{}'? {}",
                          otherMatName,
                          shapeName,
                          thisMatName,
                          shouldReplace ? "yes" : "no"));

      auto* otherMatFunc =
        meshState.getMaterialFunction(shaping::materialInOutFieldName(otherMatName));
      SLIC_ERROR_IF(otherMatFunc == nullptr,
                    axom::fmt::format("Missing inout samples for material '{}' while applying "
                                      "replacement rules for shape '{}'.",
                                      otherMatName,
                                      shapeName));

      quest::shaping::replaceMaterial(shapeFuncCopy, otherMatFunc, shouldReplace);
    }

    const std::string materialFunctionName = shaping::materialInOutFieldName(thisMatName);
    auto* materialFunc = meshState.getMaterialFunction(materialFunctionName);
    const bool hadExistingMaterial = (materialFunc != nullptr);

    if(!hadExistingMaterial)
    {
      materialFunc = meshState.createMaterialFunction(materialFunctionName);
    }

    SLIC_ERROR_IF(materialFunc == nullptr,
                  axom::fmt::format("Missing inout samples for material '{}' while updating "
                                    "the material field for shape '{}'.",
                                    thisMatName,
                                    shapeName));

    const bool reuseExisting = hadExistingMaterial && shape.getGeometry().hasGeometry();
    quest::shaping::copyShapeIntoMaterial(shapeFuncCopy, materialFunc, reuseExisting);
    if(shape.getGeometry().hasGeometry())
    {
      meshState.deleteShapeFunction(shaping::shapeInOutFieldName(shapeName));
    }

    delete shapeFuncCopy;
    shapeFuncCopy = nullptr;

    m_knownMaterials.insert(thisMatName);
  }

  /**
   * \brief Compute volume fractions for a given material using its associated quadrature function.
   *
   * The generated field uses shaping::volumeFractionFieldName() for its name.
   *
   * \param [in] matField The name of the material
   */
  void computeVolumeFractionsForMaterial(const std::string& matField);

private:
  // Holds an instance of the 2D or 3D sampler; only one can be active at a time
  SamplerVariant m_sampler;
  axom::Array<axom::primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> m_contours;

  std::set<std::string> m_knownMaterials;

  shaping::PointProjector<2, 2> m_projector22 {};
  shaping::PointProjector<3, 2> m_projector32 {};
  shaping::PointProjector<2, 3> m_projector23 {};
  shaping::PointProjector<3, 3> m_projector33 {};

  shaping::VolFracSampling m_vfSampling {shaping::VolFracSampling::SAMPLE_AT_QPTS};
  axom::numerics::QuadratureType m_quadratureType {axom::numerics::QuadratureType::Invalid};
  axom::Array<int> m_samplingResolution {};
  int m_volfracOrder {2};
  SamplingMethod m_samplingMethod {SamplingMethod::InOut};
  bool m_inoutOctreeVtkOutputEnabled {false};
  std::string m_inoutOctreeVtkOutputDirectory;
};

}  // namespace quest
}  // namespace axom
