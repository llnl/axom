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

#ifdef CONDUIT_RELAY_IO_HDF5_ENABLED
  #ifdef CONDUIT_RELAY_MPI_ENABLED
    #include "conduit_relay_mpi_io_blueprint.hpp"
  #else
    #include "conduit_relay_io_blueprint.hpp"
  #endif
#endif

#include "axom/fmt.hpp"

#include <functional>
#include <numeric>

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
  {
    initializeSamplingMFEMState();
  }

#if defined(AXOM_USE_CONDUIT)
  SamplingShaper(RuntimePolicy execPolicy,
                 int allocatorId,
                 const klee::ShapeSet& shapeSet,
                 sidre::Group* bpMesh,
                 const std::string& topo = "")
    : Shaper(execPolicy, allocatorId, shapeSet, bpMesh, topo)
  { }

  SamplingShaper(RuntimePolicy execPolicy,
                 int allocatorId,
                 const klee::ShapeSet& shapeSet,
                 conduit::Node& bpNode,
                 const std::string& topo = "")
    : Shaper(execPolicy, allocatorId, shapeSet, bpNode, topo)
  { }
#endif

  ~SamplingShaper() override = default;

  ///@{
  //!  @name Functions to get and set shaping parameters related to sampling; supplements parameters in base class

  void setSamplingType(shaping::VolFracSampling vfSampling) { m_vfSampling = vfSampling; }

  void setSamplingMethod(SamplingMethod samplingMethod) { m_samplingMethod = samplingMethod; }

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
  void setQuadratureType(axom::numerics::QuadratureType qtype)
  {
    if(axom::numerics::is_valid_quadrature_type(static_cast<int>(qtype)))
    {
      m_quadratureType = qtype;
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Invalid quadrature type value {}", static_cast<int>(qtype)));
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

protected:
  bool verifyInputMeshImpl(std::string& whyBad) const override
  {
    bool rval = true;

#if defined(AXOM_USE_CONDUIT)
    if(m_bp_state != nullptr)
    {
      rval = verifyBlueprintMeshIsStructuredOrUnstructuredQuadHex(whyBad);
    }
#endif

#if defined(AXOM_USE_MFEM)
    if(getDC() != nullptr)
    {
      rval = verifyMFEMInputMesh(whyBad);
    }
#endif

    return rval;
  }

public:
  /// Returns a pointer to the quadrature function associated with shape \a name if it exists, else nullptr
  mfem::QuadratureFunction* getShapeQFunction(const std::string& name) const
  {
    return shapeQFuncs().Get(name);
  }
  /// Returns a pointer to the quadrature function associated with material \a name if it exists, else nullptr
  mfem::QuadratureFunction* getMaterialQFunction(const std::string& name) const
  {
    return materialQFuncs().Get(name);
  }

  /*!
   * \brief Saves the sampling quadrature points as a Blueprint point mesh.
   *
   * For MFEM-backed sampling, this converts the `"positions"` quadrature
   * function to a temporary Blueprint point mesh before saving. For
   * Blueprint-backed sampling, this saves the generated quadrature-point
   * topology and any fields associated with it.
   */
  void saveQuadraturePoints(const std::string& filename) const
  {
#ifdef CONDUIT_RELAY_IO_HDF5_ENABLED
    conduit::Node n_mesh;

#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      auto* positions = getShapeQFunction("positions");
      if(positions == nullptr)
      {
        SLIC_WARNING("No MFEM quadrature positions are available to save.");
        return;
      }

      const int dim = positions->GetSpace()->GetMesh()->Dimension();
      mfem::real_t* X = const_cast<mfem::real_t*>(positions->GetData());
      const int npts = positions->Size() / positions->GetVDim();
      const conduit::index_t stride = dim * sizeof(mfem::real_t);
      n_mesh["coordsets/coords/type"] = "explicit";
      n_mesh["coordsets/coords/values/x"].set_external(X, npts, 0, stride);
      n_mesh["coordsets/coords/values/y"].set_external(X, npts, sizeof(mfem::real_t), stride);
      if(dim > 2)
      {
        n_mesh["coordsets/coords/values/z"].set_external(X, npts, 2 * sizeof(mfem::real_t), stride);
      }
      n_mesh["topologies/points/type"] = "unstructured";
      n_mesh["topologies/points/coordset"] = "coords";
      n_mesh["topologies/points/elements/shape"] = "point";
      std::vector<int> tmp(npts);
      std::iota(tmp.begin(), tmp.end(), 0);
      n_mesh["topologies/points/elements/connectivity"].set(tmp);
      n_mesh["topologies/points/elements/offsets"].set(tmp);
      std::fill(tmp.begin(), tmp.end(), 1);
      n_mesh["topologies/points/elements/sizes"].set(tmp);

  #ifdef CONDUIT_RELAY_MPI_ENABLED
      conduit::relay::mpi::io::blueprint::save_mesh(n_mesh, filename, "hdf5", MPI_COMM_WORLD);
  #else
      conduit::relay::io::blueprint::save_mesh(n_mesh, filename, "hdf5");
  #endif
      SLIC_INFO_ROOT(axom::fmt::format("Saved quadrature point mesh to '{}'.", filename));
      return;
    }
#endif

#if defined(AXOM_USE_CONDUIT)
    if(m_bp_state != nullptr)
    {
      constexpr const char* quadName = "quadrature_points";
      const conduit::Node& bpMesh = m_bp_state->m_internal_node;

      if(!bpMesh.has_path(axom::fmt::format("coordsets/{}", quadName)) ||
         !bpMesh.has_path(axom::fmt::format("topologies/{}", quadName)))
      {
        SLIC_WARNING("No Blueprint quadrature point mesh is available to save.");
        return;
      }

      n_mesh["coordsets"][quadName].update(bpMesh.fetch_existing(axom::fmt::format("coordsets/{}", quadName)));
      n_mesh["topologies"][quadName].update(bpMesh.fetch_existing(axom::fmt::format("topologies/{}", quadName)));

      if(bpMesh.has_path("fields"))
      {
        const conduit::Node& fields = bpMesh.fetch_existing("fields");
        for(conduit::index_t i = 0; i < fields.number_of_children(); ++i)
        {
          const conduit::Node& field = fields.child(i);
          if(field.has_path("topology") && field.fetch_existing("topology").as_string() == quadName)
          {
            n_mesh["fields"][field.name()].update(field);
          }
        }
      }

  #ifdef CONDUIT_RELAY_MPI_ENABLED
      conduit::relay::mpi::io::blueprint::save_mesh(n_mesh, filename, "hdf5", MPI_COMM_WORLD);
  #else
      conduit::relay::io::blueprint::save_mesh(n_mesh, filename, "hdf5");
  #endif
      SLIC_INFO_ROOT(axom::fmt::format("Saved quadrature point mesh to '{}'.", filename));
      return;
    }
#endif

    SLIC_WARNING("No mesh state is available for quadrature-point export.");
#else
    AXOM_UNUSED_VAR(filename);
    SLIC_WARNING("Quadrature-point export requires Conduit Relay HDF5 support.");
#endif
  }

private:
  std::unique_ptr<shaping::MFEMState> createMFEMState() override
  {
    return std::make_unique<shaping::SamplingMFEMState>();
  }

  void initializeSamplingMFEMState()
  {
    // Shaper constructs its MFEM state in the base constructor, so upgrade it
    // here rather than relying on virtual dispatch during base construction.
    auto samplingState = std::make_unique<shaping::SamplingMFEMState>();
    if(m_mfem_state != nullptr)
    {
      samplingState->m_dc = m_mfem_state->m_dc;
    }
    m_mfem_state = std::move(samplingState);
  }

  shaping::SamplingMFEMState& samplingMFEMState()
  {
    SLIC_ASSERT(m_mfem_state != nullptr);
    return static_cast<shaping::SamplingMFEMState&>(*m_mfem_state);
  }

  const shaping::SamplingMFEMState& samplingMFEMState() const
  {
    SLIC_ASSERT(m_mfem_state != nullptr);
    return static_cast<const shaping::SamplingMFEMState&>(*m_mfem_state);
  }

  shaping::QFunctionCollection& shapeQFuncs() { return samplingMFEMState().m_inoutShapeQFuncs; }
  const shaping::QFunctionCollection& shapeQFuncs() const
  {
    return samplingMFEMState().m_inoutShapeQFuncs;
  }

  shaping::QFunctionCollection& materialQFuncs()
  {
    return samplingMFEMState().m_inoutMaterialQFuncs;
  }
  const shaping::QFunctionCollection& materialQFuncs() const
  {
    return samplingMFEMState().m_inoutMaterialQFuncs;
  }

  shaping::DenseTensorCollection& tensors() { return samplingMFEMState().m_inoutTensors; }
  const shaping::DenseTensorCollection& tensors() const
  {
    return samplingMFEMState().m_inoutTensors;
  }

  shaping::MFEMArrayCollection& arrays() { return samplingMFEMState().m_inoutArrays; }
  const shaping::MFEMArrayCollection& arrays() const
  {
    return samplingMFEMState().m_inoutArrays;
  }

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
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      applyReplacementRulesImpl(samplingMFEMState(), shape);
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT)
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

    auto& mfemState = samplingMFEMState();
    auto* mesh = mfemState.m_dc->GetMesh();
    ensureSamplingPositions(mfemState);
    auto* positionsQSpace = mfemState.m_inoutShapeQFuncs.Get("positions")->GetSpace();

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
      materialQFuncs().Register(matName, matQFunc, true);
    }
  }

  void adjustVolumeFractions() override
  {
    AXOM_ANNOTATE_SCOPE("adjustVolumeFractions");

    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    for(const auto& materialName : m_knownMaterials)
    {
      const auto matName = axom::fmt::format("mat_inout_{}", materialName);
      SLIC_INFO_ROOT(
        axom::fmt::format("Generating volume fraction fields for '{}' material", matName));

      switch(m_vfSampling)
      {
      case shaping::VolFracSampling::SAMPLE_AT_QPTS:
        this->computeVolumeFractionsForMaterial(matName);
        break;
      case shaping::VolFracSampling::SAMPLE_AT_DOFS:
        break;
      }
    }
  }

  /// Prints out the names of the registered fields related to shapes and materials
  /// This function is intended to help with debugging
  void printRegisteredFieldNames(const std::string& initialMessage)
  {
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      shaping::printRegisteredFieldNames(samplingMFEMState(),
                                         m_knownMaterials,
                                         m_vfSampling,
                                         initialMessage);
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bp_state != nullptr)
    {
      shaping::printRegisteredFieldNames(*m_bp_state,
                                         m_knownMaterials,
                                         m_vfSampling,
                                         initialMessage);
      return;
    }
#endif
    SLIC_INFO_ROOT(axom::fmt::format("SamplingShaper {} has no registered fields.",
                                     initialMessage));
  }

private:
  void ensureSamplingPositions(shaping::SamplingMFEMState& mfemState)
  {
    shaping::generateSamplingPositions(mfemState, m_sampleResolution, m_quadratureType);
  }

#if defined(AXOM_USE_CONDUIT)
  void ensureSamplingPositions(shaping::BlueprintState& bpState)
  {
    shaping::generateSamplingPositions(bpState, m_sampleResolution, m_quadratureType);
  }
#endif

  static int meshDimension(const shaping::SamplingMFEMState& mfemState)
  {
    return mfemState.m_dc->GetMesh()->Dimension();
  }

#if defined(AXOM_USE_CONDUIT)
  int meshDimension(const shaping::BlueprintState& bpState) const
  {
    AXOM_UNUSED_VAR(bpState);
    return getBlueprintMeshDimension();
  }
#endif

  // Handles 2D or 3D shaping for compatible samplers, based on the template and associated parameter
  template <typename MeshState, typename SamplerType>
  void runShapeQueryImplSampler(SamplerType* sampler, MeshState& meshState)
  {
    // Sample the InOut field at the mesh quadrature points
    if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
    {
      ensureSamplingPositions(meshState);
    }

    const int meshDim = meshDimension(meshState);
    switch(m_vfSampling)
    {
    case shaping::VolFracSampling::SAMPLE_AT_QPTS:
      switch(SamplerType::DIM)
      {
      case 2:
        if(meshDim == 2)
        {
          sampler->template sampleInOutField<2, 2>(meshState,
                                                   m_sampleResolution,
                                                   static_cast<int>(m_quadratureType),
                                                   m_projector22);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 2>(meshState,
                                                   m_sampleResolution,
                                                   static_cast<int>(m_quadratureType),
                                                   m_projector32);
        }
        break;
      case 3:
        if(meshDim == 2)
        {
          sampler->template sampleInOutField<2, 3>(meshState,
                                                   m_sampleResolution,
                                                   static_cast<int>(m_quadratureType),
                                                   m_projector23);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 3>(meshState,
                                                   m_sampleResolution,
                                                   static_cast<int>(m_quadratureType),
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
      runShapeQueryImplSampler(sampler, samplingMFEMState());
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bp_state != nullptr)
    {
      runShapeQueryImplSampler(sampler, *m_bp_state);
      return;
    }
#endif
    SLIC_ERROR("No mesh state is available for SamplingShaper.");
  }

  // Handles 2D or 3D shaping for InOutSampler, based on the template and associated parameter
  template <int DIM>
  void runShapeQueryImpl(shaping::WindingNumberSampler<DIM>* sampler)
  {
 #if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      runShapeQueryImplSampler(sampler, samplingMFEMState());
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bp_state != nullptr)
    {
      runShapeQueryImplSampler(sampler, *m_bp_state);
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
      const int meshDim = meshDimension(meshState);
      if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
      {
        ensureSamplingPositions(meshState);
      }

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
          sampler->template sampleInOutField<2, 3>(meshState,
                                                   m_sampleResolution,
                                                   static_cast<int>(m_quadratureType),
                                                   m_projector23);
        }
        else if(meshDim == 3)
        {
          sampler->template sampleInOutField<3, 3>(meshState,
                                                   m_sampleResolution,
                                                   static_cast<int>(m_quadratureType),
                                                   m_projector33);
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
      runImpl(samplingMFEMState());
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bp_state != nullptr)
    {
      runImpl(*m_bp_state);
      return;
    }
#endif
    SLIC_ERROR("No mesh state is available for SamplingShaper.");
  }

  template <typename MeshState>
  void applyReplacementRulesImpl(MeshState& meshState, const klee::Shape& shape)
  {
    const auto& shapeName = shape.getName();
    const auto& thisMatName = shape.getMaterial();

    SLIC_INFO_ROOT(
      axom::fmt::format("{:-^80}",
                        axom::fmt::format("Applying replacement rules for shape '{}'", shapeName)));

    auto* shapeFunc = shape.getGeometry().hasGeometry()
      ? meshState.getShapeFunction(axom::fmt::format("inout_{}", shapeName))
      : meshState.getMaterialFunction(axom::fmt::format("mat_inout_{}", thisMatName));

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
        meshState.getMaterialFunction(axom::fmt::format("mat_inout_{}", otherMatName));
      SLIC_ERROR_IF(otherMatFunc == nullptr,
                    axom::fmt::format("Missing inout samples for material '{}' while applying "
                                      "replacement rules for shape '{}'.",
                                      otherMatName,
                                      shapeName));

      quest::shaping::replaceMaterial(shapeFuncCopy, otherMatFunc, shouldReplace);
    }

    const std::string materialFunctionName = axom::fmt::format("mat_inout_{}", thisMatName);
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

    delete shapeFuncCopy;
    shapeFuncCopy = nullptr;

    m_knownMaterials.insert(thisMatName);
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
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      shaping::computeVolumeFractionsForMaterial(
        samplingMFEMState(),
        matField,
        m_volfracOrder,
        m_sampleResolution,
        m_quadratureType);
      return;
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(m_bp_state != nullptr)
    {
      shaping::computeVolumeFractionsForMaterial(*m_bp_state, matField);
      return;
    }
#endif
    SLIC_ERROR("No mesh state is available for SamplingShaper.");
  }

  bool usesAnisotropicCustomTensorQuadrature(const mfem::Mesh& mesh) const
  {
    if(m_quadratureType == axom::numerics::QuadratureType::Invalid)
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
  int m_sampleResolution[3] = {5, 5, 5};
  int m_volfracOrder {2};
  SamplingMethod m_samplingMethod {SamplingMethod::InOut};
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_SAMPLING_SHAPER__HPP_
