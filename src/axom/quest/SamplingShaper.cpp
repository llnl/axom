// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/quest/SamplingShaper.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"

namespace axom
{
namespace quest
{

  void SamplingShaper::setQuadratureType(axom::numerics::QuadratureType qtype)
  {
    if(m_bp_state != nullptr)
    {
      // For Blueprint, we rely on Axom quadrature types and not all are implementd yet.
      if(axom::numerics::is_valid_quadrature_type(static_cast<int>(qtype)))
      {
std::cout << "Setting m_quadratureType = " << static_cast<int>(qtype) << std::endl;
        m_quadratureType = qtype;
      }
      else
      {
        SLIC_ERROR(axom::fmt::format("Invalid quadrature type value {}", static_cast<int>(qtype)));
      }
    }
  }

  void SamplingShaper::setSamplingResolution(int sampleRes)
  {
    SLIC_ERROR_IF(sampleRes < 1, "Invalid sample resolution");
    m_samplingResolution.clear();
    const auto dim = meshDimension();
    for(int d = 0; d < dim; d++)
    {
      m_samplingResolution.push_back(sampleRes);
    }
  }

  void SamplingShaper::setSamplingResolution(axom::ArrayView<int> sampleRes)
  {
    const auto dim = meshDimension();
    SLIC_ERROR_IF(static_cast<axom::IndexType>(dim) != sampleRes.size(),
                  "Number of sample resolutions does not match mesh dimension.");
    m_samplingResolution.clear();
    for(int d = 0; d < dim; d++)
    {
      SLIC_ERROR_IF(sampleRes[d] < 1, "Invalid sample resolution");
      m_samplingResolution.push_back(sampleRes[d]);
    }
  }

  void SamplingShaper::initializeSamplingResolution()
  {
    // Initialize the default number of samples based on the mesh dimension.
    const int dim = meshDimension();
    for(int d = 0; d < dim; d++)
    {
      m_samplingResolution.push_back(5);
    }
  }

  bool SamplingShaper::verifyInputMeshImpl(std::string& whyBad) const
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

#if defined(AXOM_USE_CONDUIT)
  void SamplingShaper::saveBlueprintFile(const conduit::Node &n_mesh, const std::string &filename) const
  {                                             
  #ifdef CONDUIT_RELAY_MPI_ENABLED
      conduit::relay::mpi::io::blueprint::save_mesh(n_mesh, filename, outputProtocol(), m_comm);
  #else
      conduit::relay::io::blueprint::save_mesh(n_mesh, filename, outputProtocol());
  #endif
  }
#endif

  void SamplingShaper::saveQuadraturePoints(const std::string& filename) const
  {
#if defined(AXOM_USE_CONDUIT)
    conduit::Node n_mesh;

    // Save the quadrature points from MFEM as a Blueprint file.
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      auto* positions = shapeQFuncs().Get("positions");
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

      saveBlueprintFile(n_mesh, filename);
      SLIC_INFO_ROOT(axom::fmt::format("Saved quadrature point mesh to '{}'.", filename));
      return;
    }
#endif

    // Save the Blueprint quadrature point mesh as a Blueprint file.
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

      saveBlueprintFile(n_mesh, filename);
      SLIC_INFO_ROOT(axom::fmt::format("Saved quadrature point mesh to '{}'.", filename));
      return;
    }
    SLIC_WARNING("No mesh state is available for quadrature-point export.");
#else
    AXOM_UNUSED_VAR(filename);
    SLIC_WARNING("Quadrature-point export requires Conduit Relay HDF5 support.");
#endif
  }

  void SamplingShaper::loadShape(const klee::Shape& shape)
  {
#if defined(AXOM_USE_MFEM)
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
      return;
    }
#else
    SLIC_ERROR_IF(useWindingNumberSampler(shape),
                  "SamplingShaper winding-number sampling for MFEM shapes requires MFEM support.");
#endif

    Shaper::loadShape(shape);
  }

void SamplingShaper::prepareShapeQuery(klee::Dimensions shapeDimension, const klee::Shape& shape)
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

#if defined(AXOM_USE_MFEM)
  /**
   * \brief Import an initial set of material volume fractions before shaping
   *
   * \param [in] initialGridFuncions The input data as a map from material names to grid functions
   * 
   * The imported grid functions are interpolated at quadrature points and registered
   * with the supplied names as material-based quadrature fields
   */
  void SamplingShaper::importInitialVolumeFractions(const std::map<std::string, mfem::GridFunction*>& initialGridFunctions)
  {
    internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                      : slic::message::Warning);

    auto& mfemState = samplingMFEMState();
    auto* mesh = mfemState.m_dc->GetMesh();
    // Sample the InOut field at the mesh quadrature points
    if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
    {
      shaping::generateSamplingPositions(mfemState, m_samplingResolution.view(), m_quadratureType);
    }
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

      if(shaping::usesAnisotropicCustomTensorQuadrature(*mesh, m_samplingResolution, m_quadratureType))
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
#endif

  void SamplingShaper::printRegisteredFieldNames(const std::string& initialMessage)
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

  void SamplingShaper::saveResults(bool extra)
  {
    Shaper::saveResults(extra);
    if(extra)
    {
      saveQuadraturePoints("shaping_quadrature");
    }
  }

  void SamplingShaper::computeVolumeFractionsForMaterial(const std::string& matField)
  {
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state != nullptr)
    {
      // NOTE: We pass the m_samplingResolution and m_quadratureType values to this
      //       version of the function so we can detect whether we have anisotropic
      //       sampling, which is handled differently.
      shaping::computeVolumeFractionsForMaterial(
        samplingMFEMState(),
        matField,
        m_volfracOrder,
        m_samplingResolution,
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

  void SamplingShaper::adjustVolumeFractions()
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

  int SamplingShaper::meshDimension() const
  {
    const int InvalidDimension = -1;
    int dim = InvalidDimension;
#if defined(AXOM_USE_MFEM)
    if(m_mfem_state)
    {
      dim = m_mfem_state->meshDimension();
    }
#endif
#if defined(AXOM_USE_CONDUIT)
    if(dim == InvalidDimension && m_bp_state)
    {
      dim = m_bp_state->meshDimension();
    }
#endif
    return dim;
  }

} // end namespace quest
} // end namespace axom
