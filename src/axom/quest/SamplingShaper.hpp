// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
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
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/klee.hpp"

#if !defined(AXOM_USE_MFEM) || !defined(AXOM_USE_SIDRE)
  #error SamplingShaper requires Axom to be configured with MFEM and Sidre
#endif

#include "axom/quest/Shaper.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"
#include "axom/quest/detail/shaping/InOutSampler.hpp"
#include "axom/quest/detail/shaping/PrimitiveSampler.hpp"
#include "axom/quest/detail/shaping/WindingNumberSampler.hpp"
#include "axom/quest/io/MFEMReader.hpp"

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

    m_inoutTensors.DeleteData(true);
    m_inoutTensors.clear();

    m_inoutArrays.DeleteData(true);
    m_inoutArrays.clear();
  }

  ///@{
  //!  @name Functions to get and set shaping parameters related to sampling; supplements parameters in base class

  void setSamplingType(shaping::VolFracSampling vfSampling) { m_vfSampling = vfSampling; }

  void setSamplingMethod(SamplingMethod samplingMethod) { m_samplingMethod = samplingMethod; }

  /*!
   * \brief Sets the 1D quadrature family used to generate custom sample points.
   *
   * Passing `mfem::Quadrature1D::Invalid` selects Axom's default MFEM quadrature
   * behavior. Any other accepted value must correspond to a valid
   * `mfem::Quadrature1D` enum in the inclusive range
   * `[mfem::Quadrature1D::Invalid, mfem::Quadrature1D::ClosedGL]`.
   * For uniform point sampling over the full zone, including the element
   * edges, `mfem::Quadrature1D::ClosedUniform` is often a good choice. Users
   * can experiment with other quadrature families when different sample point
   * patterns are desired.
   *
   * \param [in] qtype Integer value corresponding to an `mfem::Quadrature1D`
   *                   enum entry.
   */
  void setQuadratureType(int qtype)
  {
    if(qtype >= static_cast<int>(mfem::Quadrature1D::Invalid) &&
       qtype <= static_cast<int>(mfem::Quadrature1D::ClosedGL))
    {
      m_quadratureType = qtype;
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Invalid quadrature type value {}", qtype));
    }
  }

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
  void setSamplingResolution(int sampleRes)
  {
    SLIC_ASSERT(sampleRes > 0);
    m_sampleResolution[0] = sampleRes;
    m_sampleResolution[1] = sampleRes;
    m_sampleResolution[2] = sampleRes;
  }

  /*!
   * \brief Sets an anisotropic sampling resolution for custom quadrature.
   *
   * The entries correspond to the logical `I`, `J`, and `K` directions of the
   * reference element. Each entry must be positive. For custom quadrature
   * families, these values specify the per-direction sample counts directly,
   * which in turn determine the quadrature rule used in each logical
   * direction.
   *
   * \param [in] sampleRes Array containing the sample count per logical
   *                       direction.
   */
  void setSamplingResolution(int sampleRes[3])
  {
    SLIC_ASSERT(sampleRes[0] > 0);
    SLIC_ASSERT(sampleRes[1] > 0);
    SLIC_ASSERT(sampleRes[2] > 0);
    m_sampleResolution[0] = sampleRes[0];
    m_sampleResolution[1] = sampleRes[1];
    m_sampleResolution[2] = sampleRes[2];
  }

  // Deprecated backward compatibility method
  [[deprecated]] void setQuadratureOrder(int order) { setSamplingResolution(order); }

  void setVolumeFractionOrder(int volfracOrder) { m_volfracOrder = volfracOrder; }

  /// Registers a function to project from 2D input points to 2D query points
  void setPointProjector22(shaping::PointProjector<2, 2> projector) { m_projector22 = projector; }

  /// Registers a function to project from 3D input points to 2D query points
  void setPointProjector32(shaping::PointProjector<3, 2> projector) { m_projector32 = projector; }

  /// Registers a function to project from 2D input points to 3D query points
  void setPointProjector23(shaping::PointProjector<2, 3> projector) { m_projector23 = projector; }

  /// Registers a function to project from 3D input points to 3D query points
  void setPointProjector33(shaping::PointProjector<3, 3> projector) { m_projector33 = projector; }

  ///@}

  /// Returns a pointer to the quadrature function associated with shape \a name if it exists, else nullptr
  mfem::QuadratureFunction* getShapeQFunction(const std::string& name) const
  {
    return m_inoutShapeQFuncs.Get(name);
  }
  /// Returns a pointer to the quadrature function associated with material \a name if it exists, else nullptr
  mfem::QuadratureFunction* getMaterialQFunction(const std::string& name) const
  {
    return m_inoutMaterialQFuncs.Get(name);
  }

private:
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
  void loadShape(const klee::Shape& shape) override
  {
    if(useWindingNumberSampler(shape))
    {
      const std::string shapePath =
        axom::utilities::filesystem::prefixRelativePath(shape.getGeometry().getPath(), m_prefixPath);
      SLIC_INFO_ROOT("Reading file: " << shapePath << "...");
      // Read the MFEM file as curved polygon contours for winding number intersection.
      quest::MFEMReader reader;
      reader.setFileName(shapePath);
      const int rc = reader.read(m_contours);

      SLIC_ERROR_IF(rc != quest::MFEMReader::READ_SUCCESS,
                    axom::fmt::format("Failed to read MFEM shape '{}' from file '{}'.",
                                      shape.getName(),
                                      shapePath));
    }
    else
    {
      Shaper::loadShape(shape);
    }
  }

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

    SLIC_INFO_ROOT(axom::fmt::format("{:-^80}", " Generating the spatial index "));

    const auto& shapeName = shape.getName();

    // Initialize the sampler based on shape format
    // note: ignoring the global shapeDimension for now since it's causing problems
    // reading c2c when the dimension is Three
    AXOM_UNUSED_VAR(shapeDimension);
    const auto format = this->shapeFormat(shape);
    if(useWindingNumberSampler(shape))
    {
      m_sampler = std::make_unique<WindingNumberSampler2D>(shapeName, m_contours.view());
    }
    else if(format == "c2c" || format == "mfem")
    {
      m_sampler = std::make_unique<InOutSampler2D>(shapeName, m_surfaceMesh);
    }
    else if(format == "stl")
    {
      m_sampler = std::make_unique<InOutSampler3D>(shapeName, m_surfaceMesh);
    }
    else if(format == "proe")
    {
      using Policy = runtime_policy::Policy;
      switch(this->getExecutionPolicy())
      {
      case Policy::seq:
        m_sampler = std::make_unique<PrimitiveSampler3D_seq>(shapeName, m_surfaceMesh);
        break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
      case Policy::omp:
        m_sampler = std::make_unique<PrimitiveSampler3D_omp>(shapeName, m_surfaceMesh);
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
      case Policy::cuda:
        m_sampler = std::make_unique<PrimitiveSampler3D_cuda>(shapeName, m_surfaceMesh);
        break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
      case Policy::hip:
        m_sampler = std::make_unique<PrimitiveSampler3D_hip>(shapeName, m_surfaceMesh);
        break;
#endif
      default:
        SLIC_ERROR("Unsupported execution policy for PrimitiveSampler3D");
        break;
      }
    }

    SLIC_ASSERT(hasValidSampler());

    // Use visitor to initialize the sampler
    std::visit(
      [this](auto& sampler) {
        using T = std::decay_t<decltype(sampler)>;
        if constexpr(std::is_same_v<T, std::monostate>)
        {
          // no op -- monostate
        }
        else if constexpr(is_wnsampler_v<typename T::element_type>)
        {
          sampler->computeBounds();
          sampler->initSpatialIndex(this->m_vertexWeldThreshold);
        }
        else if constexpr(is_inoutsampler_v<typename T::element_type>)
        {
          sampler->computeBounds();
          sampler->initSpatialIndex(this->m_vertexWeldThreshold);
        }
        else if constexpr(is_primitivesampler_v<typename T::element_type>)
        {
          sampler->computeBounds();
          sampler->initSpatialIndex();
        }
      },
      m_sampler);

    // Output some logging info and dump the mesh
    if(this->isVerbose() && this->getRank() == 0)
    {
      if(m_surfaceMesh != nullptr)
      {
        const int nVerts = m_surfaceMesh->getNumberOfNodes();
        const int nCells = m_surfaceMesh->getNumberOfCells();
        SLIC_INFO(axom::fmt::format("After welding, surface mesh has {} vertices  and {} elements.",
                                    nVerts,
                                    nCells));
        mint::write_vtk(m_surfaceMesh.get(),
                        axom::fmt::format("melded_shape_mesh_{}.vtk", shapeName));
      }
      else if(!m_contours.empty())
      {
        SLIC_INFO(axom::fmt::format("Contours contain {} curved polygons.", m_contours.size()));
      }
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

    const auto& shapeName = shape.getName();
    const auto& thisMatName = shape.getMaterial();

    SLIC_INFO_ROOT(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Applying replacement rules for shape '{}'", shapeName)));

    mfem::QuadratureFunction* shapeQFunc = nullptr;

    if(shape.getGeometry().hasGeometry())
    {
      // Get inout qfunc for this shape
      shapeQFunc = m_inoutShapeQFuncs.Get(axom::fmt::format("inout_{}", shapeName));

      SLIC_ERROR_IF(shapeQFunc == nullptr,
                    axom::fmt::format("Missing inout samples for shape '{}'. "
                                      "This indicates the shape query did not produce a "
                                      "quadrature field before replacement rules were applied.",
                                      shapeName));
    }
    else
    {
      // No input geometry for the shape, get inout qfunc for associated material
      shapeQFunc = m_inoutMaterialQFuncs.Get(axom::fmt::format("mat_inout_{}", thisMatName));

      SLIC_ERROR_IF(shapeQFunc == nullptr,
                    axom::fmt::format("Missing inout samples for material '{}' while applying "
                                      "replacement rules for shape '{}', which has no input "
                                      "geometry. Initialize that material before shaping, e.g. "
                                      "pass '--background-material {}' in the shaping driver or "
                                      "import initial volume fractions for it.",
                                      thisMatName,
                                      shapeName,
                                      thisMatName));
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
      SLIC_INFO_ROOT(
        axom::fmt::format("Should we replace material '{}' with shape '{}' of material '{}'? {}",
                          otherMatName,
                          shapeName,
                          thisMatName,
                          shouldReplace ? "yes" : "no"));

      auto* otherMatQFunc =
        m_inoutMaterialQFuncs.Get(axom::fmt::format("mat_inout_{}", otherMatName));
      SLIC_ERROR_IF(otherMatQFunc == nullptr,
                    axom::fmt::format("Missing inout samples for material '{}' while applying "
                                      "replacement rules for shape '{}'.",
                                      otherMatName,
                                      shapeName));

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
      SLIC_ERROR_IF(matQFunc == nullptr,
                    axom::fmt::format("Missing inout samples for material '{}' while updating "
                                      "the material field for shape '{}'.",
                                      thisMatName,
                                      shapeName));

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
      shaping::generatePositionsQFunction(mesh,
                                          m_inoutShapeQFuncs,
                                          m_sampleResolution,
                                          m_quadratureType);
    }
    auto* positionsQSpace = m_inoutShapeQFuncs.Get("positions")->GetSpace();

    // Interpolate grid functions at quadrature points & register material quad functions
    // assume all elements have same integration rule
    for(auto& entry : initialGridFunctions)
    {
      const auto& name = entry.first;
      auto* gf = entry.second;

      SLIC_INFO_ROOT(axom::fmt::format("Importing volume fraction field for '{}' material", name));

      if(gf == nullptr)
      {
        SLIC_WARNING(
          axom::fmt::format("Skipping missing volume fraction field for material '{}'", name));
        continue;
      }

      auto* matQFunc = new mfem::QuadratureFunction(*positionsQSpace);
      const auto& ir = matQFunc->GetSpace()->GetIntRule(0);

      if(usesAnisotropicCustomTensorQuadrature(*mesh))
      {
        // Avoid MFEM's tensor quadrature interpolation path only for
        // anisotropic custom quad/hex rules. MFEM infers a single q1d from
        // ir.GetNPoints(), which cannot represent per-direction sample counts
        // such as 3 x 5 or 3 x 5 x 2.
        mfem::Vector elemValues;
        mfem::Vector qfuncValues;
        for(int elem = 0; elem < mesh->GetNE(); ++elem)
        {
          gf->GetValues(elem, ir, elemValues);
          matQFunc->GetValues(elem, qfuncValues);
          qfuncValues = elemValues;
        }
      }
      else
      {
        const auto* interp = gf->FESpace()->GetQuadratureInterpolator(ir);
        SLIC_ERROR_IF(interp == nullptr,
                      axom::fmt::format("Could not create a quadrature interpolator while "
                                        "importing volume fractions for '{}'.",
                                        name));
        interp->Values(*gf, *matQFunc);
      }

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
      SLIC_INFO_ROOT(
        axom::fmt::format("Generating volume fraction fields for '{}' material", matName));

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
                           "\n\t* Shaping tensors: {}",
                           axom::fmt::join(extractKeys(m_inoutTensors), ", "));
    }
    SLIC_INFO_ROOT(axom::fmt::to_string(out));
  }

private:
  // Handles 2D or 3D shaping for compatible samplers, based on the template and associated parameter
  template <typename SamplerType>
  void runShapeQueryImplSampler(SamplerType* sampler)
  {
    // Sample the InOut field at the mesh quadrature points
    const int meshDim = m_dc->GetMesh()->Dimension();
    switch(m_vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      switch(SamplerType::DIM)
      {
      case 2:
        if(meshDim == 2)
        {
          sampler->template sampleInOutField<2, 2>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_sampleResolution,
                                                   m_quadratureType,
                                                   m_projector22);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 2>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_sampleResolution,
                                                   m_quadratureType,
                                                   m_projector32);
        }
        break;
      case 3:
        if(meshDim == 2)
        {
          sampler->template sampleInOutField<2, 3>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_sampleResolution,
                                                   m_quadratureType,
                                                   m_projector23);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 3>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_sampleResolution,
                                                   m_quadratureType,
                                                   m_projector33);
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
          sampler->template computeVolumeFractionsBaseline<2, 2>(m_dc, m_volfracOrder, m_projector22);
        }
        else if(meshDim == 3)
        {
          sampler->template computeVolumeFractionsBaseline<3, 2>(m_dc, m_volfracOrder, m_projector32);
        }
        break;
      case 3:
        if(meshDim == 2)
        {
          sampler->template computeVolumeFractionsBaseline<2, 3>(m_dc, m_volfracOrder, m_projector23);
        }
        else if(meshDim == 3)
        {
          sampler->template computeVolumeFractionsBaseline<3, 3>(m_dc, m_volfracOrder, m_projector33);
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
    runShapeQueryImplSampler(sampler);
  }

  // Handles 2D or 3D shaping for InOutSampler, based on the template and associated parameter
  template <int DIM>
  void runShapeQueryImpl(shaping::WindingNumberSampler<DIM>* sampler)
  {
    runShapeQueryImplSampler(sampler);
  }

  // Handles 2D or 3D shaping for PrimitiveSampler, based on the template and associated parameter
  template <int DIM, typename ExecSpace>
  void runShapeQueryImpl(shaping::PrimitiveSampler<DIM, ExecSpace>* sampler)
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
          sampler->template sampleInOutField<2, 3>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_sampleResolution,
                                                   m_quadratureType,
                                                   m_projector23);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 3>(m_dc,
                                                   m_inoutShapeQFuncs,
                                                   m_sampleResolution,
                                                   m_quadratureType,
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

    const auto& sampleIR = inout->GetSpace()->GetIntRule(0);  // assume all elements are the same
    const int sampleOrder = sampleIR.GetOrder();
    const int sampleNQ = sampleIR.GetNPoints();
    const int sampleSZ = inout->GetSpace()->GetSize();

    // extract some properties from computational mesh
    mfem::Mesh* mesh = m_dc->GetMesh();
    const int dim = mesh->Dimension();
    const int NE = mesh->GetNE();
    const auto geom = mesh->GetTypicalElementGeometry();

    auto samples_per_dim = [=](int sampleRes[3], mfem::Geometry::Type geom) -> std::string {
      switch(geom)
      {
      case mfem::Geometry::SQUARE:
        return axom::fmt::format(" ({} * {})", sampleRes[0], sampleRes[1]);
      case mfem::Geometry::CUBE:
        return axom::fmt::format(" ({} * {} * {})", sampleRes[0], sampleRes[1], sampleRes[2]);
      default:
        return std::string();
      }
    };

    // print info about sampling on rank 0
    // TODO: mpi reduce this for stats on all ranks
    SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                     "In computeVolumeFractions(): num samples per element {}{} | "
                                     "sample polynomial order {} | total samples {:L}",
                                     sampleNQ,
                                     samples_per_dim(m_sampleResolution, geom),
                                     sampleOrder,
                                     sampleSZ));

    SLIC_INFO_ROOT(
      axom::fmt::format(axom::utilities::locale(), "Mesh has dim {} and {:L} elements", dim, NE));

    // Access or create a registered volume fraction grid function from the data collection
    const auto vf_name = axom::fmt::format("vol_frac_{}", matField.substr(10));
    mfem::GridFunction* vf = shaping::getOrAllocateL2GridFunction(m_dc,
                                                                  vf_name,
                                                                  m_volfracOrder,
                                                                  dim,
                                                                  mfem::BasisType::Positive);
    const mfem::FiniteElementSpace* fes = vf->FESpace();
    const int dofs = fes->GetTypicalFE()->GetDof();

    // access or compute the mass matrix
    mfem::DenseTensor* mass_mat {nullptr};
    const std::string mass_matrix_name = "shaping_mass_matrix";
    if(this->m_inoutTensors.Has(mass_matrix_name))
    {
      mass_mat = m_inoutTensors.Get(mass_matrix_name);
    }
    else
    {
      AXOM_ANNOTATE_SCOPE("mass integrator assemble");

      mass_mat = new mfem::DenseTensor(dofs, dofs, NE);
      mass_mat->HostWrite();
      (*mass_mat) = 0.;
      mass_mat->ReadWrite();

      mfem::ConstantCoefficient one_coef(1.0);
      mfem::MassIntegrator mass_integrator(one_coef, &sampleIR);

      if(usesAnisotropicCustomTensorQuadrature(*fes->GetMesh()))
      {
        mfem::DenseMatrix elemMat;
        mass_mat->HostWrite();
        for(int elem = 0; elem < NE; ++elem)
        {
          mass_integrator.AssembleElementMatrix(*fes->GetFE(elem),
                                                *fes->GetElementTransformation(elem),
                                                elemMat);
          for(int j = 0; j < dofs; ++j)
          {
            for(int i = 0; i < dofs; ++i)
            {
              (*mass_mat)(i, j, elem) = elemMat(i, j);
            }
          }
        }
      }
      else
      {
        const int sz = mass_mat->TotalSize();

        // wrap mass_mat data as vector for AssembleEA call
        // note: AssembleEA expects the transpose, but it's ok since mass matrices are symmetric
        mfem::Vector mass_vec;
        mfem::Swap(mass_mat->GetMemory(), mass_vec.GetMemory());
        mass_vec.SetSize(sz);
        mass_integrator.AssembleEA(*fes, mass_vec, false);
        mfem::Swap(mass_mat->GetMemory(), mass_vec.GetMemory());
      }

      m_inoutTensors.Register(mass_matrix_name, mass_mat, true);
    }
    SLIC_ASSERT(mass_mat->SizeI() == dofs);
    SLIC_ASSERT(mass_mat->SizeJ() == dofs);
    SLIC_ASSERT(mass_mat->SizeK() == NE);

    // access or compute LU factorization of the mass matrix
    mfem::DenseTensor* mass_mat_inv {nullptr};
    mfem::Array<int>* mass_mat_pivots {nullptr};
    const std::string minv_name = "shaping_mass_matrix_inv";
    const std::string pivots_name = "shaping_mass_matrix_pivots";
    if(this->m_inoutTensors.Has(minv_name) && this->m_inoutArrays.Has(pivots_name))
    {
      mass_mat_inv = this->m_inoutTensors.Get(minv_name);
      mass_mat_pivots = this->m_inoutArrays.Get(pivots_name);
    }
    else
    {
      AXOM_ANNOTATE_SCOPE("batch lu factor");

      // Perform batched LU factorization on the mass tensor
      mass_mat->ReadWrite();
      mass_mat_inv = new mfem::DenseTensor(*mass_mat);
      mass_mat_pivots = new mfem::Array<int>(dofs * NE);

      mass_mat_inv->ReadWrite();
      mass_mat_pivots->Write();
      mfem::BatchLUFactor(*mass_mat_inv, *mass_mat_pivots);

      m_inoutTensors.Register(minv_name, mass_mat_inv, true);
      m_inoutArrays.Register(pivots_name, mass_mat_pivots, true);
    }
    SLIC_ASSERT(mass_mat_inv->SizeJ() == dofs);
    SLIC_ASSERT(mass_mat_inv->SizeI() == dofs);
    SLIC_ASSERT(mass_mat_inv->SizeK() == NE);
    SLIC_ASSERT(mass_mat_pivots->Size() == dofs * NE);

    mfem::DenseTensor* shaping_scratch_buffer {nullptr};
    const std::string scratch_buffer_name = "shaping_scratch_buffer";
    if(this->m_inoutTensors.Has(scratch_buffer_name))
    {
      shaping_scratch_buffer = this->m_inoutTensors.Get(scratch_buffer_name);
    }
    else
    {
      shaping_scratch_buffer = new mfem::DenseTensor(dofs, dofs, NE);
      // TODO -- we only need this buffer to be Write
      // and only in the space that FCT_project is called
      shaping_scratch_buffer->HostWrite();
      (*shaping_scratch_buffer) = 0.;

      m_inoutTensors.Register(scratch_buffer_name, shaping_scratch_buffer, true);
    }
    SLIC_ASSERT(shaping_scratch_buffer->SizeJ() == dofs);
    SLIC_ASSERT(shaping_scratch_buffer->SizeI() == dofs);
    SLIC_ASSERT(shaping_scratch_buffer->SizeK() == NE);

    // Project QField onto volume fractions field using flux corrected transport (FCT)
    // to keep the range of values between 0 and 1
    axom::utilities::Timer timer(true);
    {
      // assemble the right hand side integral, incorporating the inout samples
      mfem::Vector b(fes->GetVSize());
      SLIC_ASSERT(b.Size() == dofs * NE);
      {
        AXOM_ANNOTATE_SCOPE("domain lf integrator assemble");

        inout->ReadWrite();

        b.HostWrite();
        b = 0.;
        b.ReadWrite();

        this->assembleVolumeFractionRHS(*fes, *inout, sampleIR, b);
      }
      inout->HostReadWrite();

      {
        AXOM_ANNOTATE_SCOPE("batch lu solve");

        mass_mat_inv->Read();
        mass_mat_pivots->Read();

        vf->HostReadWrite();
        (*vf) = b;
        vf->ReadWrite();
        mfem::BatchLUSolve(*mass_mat_inv, *mass_mat_pivots, *vf);
      }
      mass_mat_inv->HostReadWrite();
      mass_mat_pivots->HostReadWrite();

      constexpr double minY = 0.;
      constexpr double maxY = 1.;

      // Reshape returns an indexable view of a multidimensional array
      auto m_d = mfem::Reshape(mass_mat->HostReadWrite(), dofs, dofs, NE);
      auto b_d = mfem::Reshape(b.HostReadWrite(), dofs, NE);
      auto vf_d = mfem::Reshape(vf->HostReadWrite(), dofs, NE);
      auto fct_mat_d = mfem::Reshape(shaping_scratch_buffer->HostReadWrite(), dofs, dofs, NE);

      AXOM_ANNOTATE_BEGIN("fct project");
      axom::for_all<axom::SEQ_EXEC>(0, NE, [=](int i) {
        shaping::FCT_correct(&m_d(0, 0, i),
                             dofs,
                             &b_d(0, i),
                             minY,
                             maxY,
                             &vf_d(0, i),
                             &fct_mat_d(0, 0, i));
      });
      AXOM_ANNOTATE_END("fct project");
    }
    timer.stop();

    // print stats for root rank
    SLIC_INFO_ROOT(axom::fmt::format(axom::utilities::locale(),
                                     "\t Generating volume fractions '{}' took {:.3f} seconds (@ "
                                     "{:L} dofs processed per second)",
                                     vf_name,
                                     timer.elapsed(),
                                     static_cast<int>(fes->GetNDofs() / timer.elapsed())));

    vf->HostReadWrite();
  }

  bool usesAnisotropicCustomTensorQuadrature(const mfem::Mesh& mesh) const
  {
    if(m_quadratureType == static_cast<int>(mfem::Quadrature1D::Invalid))
    {
      return false;
    }

    switch(mesh.GetTypicalElementGeometry())
    {
    case mfem::Geometry::SQUARE:
      return m_sampleResolution[0] != m_sampleResolution[1];
    case mfem::Geometry::CUBE:
      return m_sampleResolution[0] != m_sampleResolution[1] ||
        m_sampleResolution[0] != m_sampleResolution[2];
    default:
      return false;
    }
  }

  void assembleVolumeFractionRHS(const mfem::FiniteElementSpace& fes,
                                 mfem::QuadratureFunction& inout,
                                 const mfem::IntegrationRule& sampleIR,
                                 mfem::Vector& b) const
  {
    mfem::QuadratureFunctionCoefficient qfc(inout);
    mfem::DomainLFIntegrator rhs(qfc, &sampleIR);

    if(usesAnisotropicCustomTensorQuadrature(*fes.GetMesh()))
    {
      mfem::Vector elemVec;
      mfem::Array<int> elemVDofs;

      for(int elem = 0; elem < fes.GetNE(); ++elem)
      {
        rhs.AssembleRHSElementVect(*fes.GetFE(elem), *fes.GetElementTransformation(elem), elemVec);
        fes.GetElementVDofs(elem, elemVDofs);
        b.AddElementVector(elemVDofs, elemVec);
      }
    }
    else
    {
      mfem::Array<int> elem_marker(fes.GetNE());
      elem_marker.HostWrite();
      elem_marker = 1;
      elem_marker.ReadWrite();
      rhs.AssembleDevice(fes, elem_marker, b);
    }
  }

private:
  shaping::QFunctionCollection m_inoutShapeQFuncs;
  shaping::QFunctionCollection m_inoutMaterialQFuncs;
  shaping::DenseTensorCollection m_inoutTensors;
  shaping::MFEMArrayCollection m_inoutArrays;

  // Holds an instance of the 2D or 3D sampler; only one can be active at a time
  SamplerVariant m_sampler;
  axom::Array<axom::primal::CurvedPolygon<axom::primal::NURBSCurve<double, 2>>> m_contours;

  std::set<std::string> m_knownMaterials;

  shaping::PointProjector<2, 2> m_projector22 {};
  shaping::PointProjector<3, 2> m_projector32 {};
  shaping::PointProjector<2, 3> m_projector23 {};
  shaping::PointProjector<3, 3> m_projector33 {};

  shaping::VolFracSampling m_vfSampling {shaping::VolFracSampling::SAMPLE_AT_QPTS};
  int m_quadratureType {static_cast<int>(mfem::Quadrature1D::Invalid)};
  int m_sampleResolution[3] = {5, 5, 5};
  int m_volfracOrder {2};
  SamplingMethod m_samplingMethod {SamplingMethod::InOut};
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SAMPLING_SHAPER__HPP_
