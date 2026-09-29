// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/quest/SamplingShaper.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"
#if defined(AXOM_USE_CONDUIT)
  #include "axom/quest/detail/shaping/shaping_helpers_blueprint.hpp"
  #include "conduit/conduit_relay_io.hpp"
  #ifdef CONDUIT_RELAY_IO_HDF5_ENABLED
    #ifdef CONDUIT_RELAY_MPI_ENABLED
      #include "conduit/conduit_relay_mpi_io_blueprint.hpp"
    #else
      #include "conduit/conduit_relay_io_blueprint.hpp"
    #endif
  #endif
#endif

namespace axom
{
namespace quest
{

void SamplingShaper::setQuadratureType(axom::numerics::QuadratureType qtype)
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

void SamplingShaper::setVolumeFractionOrder(int volfracOrder)
{
#if defined(AXOM_USE_CONDUIT)
  if(m_bp_state != nullptr)
  {
    SLIC_INFO("setVolumeFractionOrder is ignored for Blueprint meshes.");
    return;
  }
#endif
  m_volfracOrder = axom::utilities::max(1, volfracOrder);
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
void SamplingShaper::saveBlueprintFile(const conduit::Node& n_mesh, const std::string& filename) const
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
    auto* positions = m_mfem_state->shapeQFuncs().Get("positions");
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
    const conduit::Node& bpMesh = m_bp_state->getBlueprintMeshNode();

    if(!bpMesh.has_path(axom::fmt::format("coordsets/{}", shaping::QUADRATURE_COORDSET_NAME)) ||
       !bpMesh.has_path(axom::fmt::format("topologies/{}", shaping::QUADRATURE_TOPOLOGY_NAME)))
    {
      SLIC_WARNING("No Blueprint quadrature point mesh is available to save.");
      return;
    }

    n_mesh["coordsets"][shaping::QUADRATURE_COORDSET_NAME].update(
      bpMesh.fetch_existing(axom::fmt::format("coordsets/{}", shaping::QUADRATURE_COORDSET_NAME)));
    n_mesh["topologies"][shaping::QUADRATURE_TOPOLOGY_NAME].update(
      bpMesh.fetch_existing(axom::fmt::format("topologies/{}", shaping::QUADRATURE_TOPOLOGY_NAME)));

    if(bpMesh.has_path("fields"))
    {
      for(const auto& fieldName : m_bp_state->fieldNames())
      {
        const conduit::Node& field = m_bp_state->getField(fieldName);
        if(field.has_path("topology") &&
           field.fetch_existing("topology").as_string() == shaping::QUADRATURE_TOPOLOGY_NAME)
        {
          n_mesh["fields"][fieldName].update(field);
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

    SLIC_ERROR_IF(
      rc != quest::MFEMReader::READ_SUCCESS,
      axom::fmt::format("Failed to read MFEM shape '{}' from file '{}'.", shape.getName(), shapePath));
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
        sampler->initSpatialIndex(this->m_vertexWeldThreshold,
                                  m_inoutOctreeVtkOutputEnabled,
                                  m_inoutOctreeVtkOutputDirectory);
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
      mint::write_vtk(m_surfaceMesh.get(), axom::fmt::format("melded_shape_mesh_{}.vtk", shapeName));
    }
    else if(!m_contours.empty())
    {
      SLIC_INFO(axom::fmt::format("Contours contain {} curved polygons.", m_contours.size()));
    }
  }
}

#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
void SamplingShaper::importInitialVolumeFractions(
  const std::map<std::string, conduit::Node*>& initialVolumeFractions)
{
  internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                    : slic::message::Warning);
  SLIC_ERROR_IF(m_bp_state == nullptr, "This method requires Blueprint inputs.");
  // Generate the quadrature points.
  if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
  {
    shaping::generateSamplingPositions(*m_bp_state, m_samplingResolution.view(), m_quadratureType);
  }
  shaping::importInitialVolumeFractions(*m_bp_state, initialVolumeFractions);
}
#endif

#if defined(AXOM_USE_MFEM)
void SamplingShaper::importInitialVolumeFractions(
  const std::map<std::string, mfem::GridFunction*>& initialGridFunctions)
{
  internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                    : slic::message::Warning);

  SLIC_ERROR_IF(m_mfem_state == nullptr, "This method requires MFEM inputs.");

  auto& mfemState = *m_mfem_state;
  auto* mesh = mfemState.m_dc->GetMesh();
  // Generate the quadrature points.
  if(m_vfSampling == shaping::VolFracSampling::SAMPLE_AT_QPTS)
  {
    shaping::generateSamplingPositions(mfemState, m_samplingResolution.view(), m_quadratureType);
  }
  const bool anisotropic =
    shaping::usesAnisotropicCustomTensorQuadrature(*mesh, m_samplingResolution, m_quadratureType);
  shaping::importInitialVolumeFractions(mfemState, initialGridFunctions, anisotropic);
}
#endif

void SamplingShaper::printRegisteredFieldNames(const std::string& initialMessage)
{
#if defined(AXOM_USE_MFEM)
  if(m_mfem_state != nullptr)
  {
    shaping::printRegisteredFieldNames(*m_mfem_state, m_knownMaterials, m_vfSampling, initialMessage);
    return;
  }
#endif
#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
  if(m_bp_state != nullptr)
  {
    shaping::printRegisteredFieldNames(*m_bp_state, m_knownMaterials, m_vfSampling, initialMessage);
    return;
  }
#endif
  SLIC_INFO_ROOT(axom::fmt::format("SamplingShaper {} has no registered fields.", initialMessage));
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
    shaping::computeVolumeFractionsForMaterial(*m_mfem_state,
                                               matField,
                                               m_volfracOrder,
                                               m_samplingResolution,
                                               m_quadratureType);
    return;
  }
#endif
#if defined(AXOM_USE_CONDUIT) && defined(AXOM_USE_BUMP)
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
    const auto matName = shaping::materialInOutFieldName(materialName);
    SLIC_INFO_ROOT(axom::fmt::format("Generating volume fraction fields for '{}' material", matName));

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

}  // end namespace quest
}  // end namespace axom
