// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Shaper.hpp"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal/operators/split.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/Shaper.hpp"
#include "axom/quest/DiscreteShape.hpp"
#include "axom/quest/util/mesh_helpers.hpp"
#include "conduit_blueprint_mesh.hpp"

#include "axom/fmt.hpp"

#include "conduit/conduit_relay_io.hpp"
#ifdef CONDUIT_RELAY_IO_HDF5_ENABLED
  #ifdef CONDUIT_RELAY_MPI_ENABLED
    #include "conduit/conduit_relay_mpi_io_blueprint.hpp"
  #else
    #include "conduit/conduit_relay_io_blueprint.hpp"
  #endif
#endif

namespace axom
{
namespace quest
{

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
namespace
{
bool mpiIsActive()
{
  int initialized = 0;
  MPI_Initialized(&initialized);
  if(!initialized)
  {
    return false;
  }

  int finalized = 0;
  MPI_Finalized(&finalized);
  return finalized == 0;
}
}  // namespace
#endif

// These were needed for linking - but why? They are constexpr.
constexpr int Shaper::DEFAULT_SAMPLES_PER_KNOT_SPAN;
constexpr double Shaper::MINIMUM_PERCENT_ERROR;
constexpr double Shaper::MAXIMUM_PERCENT_ERROR;
constexpr double Shaper::DEFAULT_VERTEX_WELD_THRESHOLD;

#if defined(AXOM_USE_MFEM)
Shaper::Shaper(RuntimePolicy execPolicy,
               int allocatorId,
               const klee::ShapeSet& shapeSet,
               sidre::MFEMSidreDataCollection* dc)
  : m_execPolicy(execPolicy)
  , m_allocatorId(allocatorId != axom::INVALID_ALLOCATOR_ID
                    ? allocatorId
                    : axom::policyToDefaultAllocatorID(execPolicy))
  , m_shapeSet(shapeSet)
  , m_mfem_state()
  #if defined(AXOM_USE_CONDUIT)
  , m_bp_state()
  #endif
{
  m_mfem_state = createMFEMState();
  m_mfem_state->m_dc = dc;

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  m_comm = m_mfem_state->m_dc->GetComm();
  #endif
  m_cellCount = m_mfem_state->m_dc->GetMesh()->GetNE();

  setFilePath(shapeSet.getPath());
}
#endif

#if defined(AXOM_USE_CONDUIT)
Shaper::Shaper(RuntimePolicy execPolicy,
               int allocatorId,
               const klee::ShapeSet& shapeSet,
               sidre::Group* bpGrp,
               const std::string& topo)
  : m_execPolicy(execPolicy)
  , m_allocatorId(allocatorId != axom::INVALID_ALLOCATOR_ID
                    ? allocatorId
                    : axom::policyToDefaultAllocatorID(execPolicy))
  , m_shapeSet(shapeSet)
  #if defined(AXOM_USE_MFEM)
  , m_mfem_state()
  #endif
  , m_bp_state()
  #if defined(AXOM_USE_MPI)
  , m_comm(MPI_COMM_WORLD)
  #endif
{
  m_bp_state = createBlueprintState();
  bpGrp->setDefaultArrayAllocator(m_allocatorId);
  m_bp_state->initialize(bpGrp, m_allocatorId, resolveBlueprintTopologyName(bpGrp, topo));

  SLIC_ASSERT(m_bp_state->isSidreBacked());

  refreshBlueprintMeshState();

  setFilePath(shapeSet.getPath());
}
#endif

#if defined(AXOM_USE_CONDUIT)
Shaper::Shaper(RuntimePolicy execPolicy,
               int allocatorId,
               const klee::ShapeSet& shapeSet,
               conduit::Node& bpNode,
               const std::string& topo)
  : m_execPolicy(execPolicy)
  , m_allocatorId(allocatorId != axom::INVALID_ALLOCATOR_ID
                    ? allocatorId
                    : axom::policyToDefaultAllocatorID(execPolicy))
  , m_shapeSet(shapeSet)
  #if defined(AXOM_USE_MFEM)
  , m_mfem_state()
  #endif
  , m_bp_state()
  #if defined(AXOM_USE_MPI)
  , m_comm(MPI_COMM_WORLD)
  #endif
{
  AXOM_ANNOTATE_SCOPE("Shaper::Shaper_Node");

  m_bp_state = createBlueprintState();
  m_bp_state->initialize(&bpNode, m_allocatorId, resolveBlueprintTopologyName(bpNode, topo));

  refreshBlueprintMeshState();

  setFilePath(shapeSet.getPath());
}
#endif

Shaper::~Shaper() { }

void Shaper::setFilePath(const std::string& filePath)
{
  if(filePath.empty())
  {
    m_prefixPath.clear();
  }
  else
  {
    m_prefixPath = utilities::filesystem::getParentPath(filePath);
  }
}

void Shaper::setSamplesPerKnotSpan(int nSamples)
{
  using axom::utilities::clampLower;
  SLIC_WARNING_ROOT_IF(
    nSamples < 1,
    axom::fmt::format("Samples per knot span must be at least 1. Provided value was {}", nSamples));

  m_samplesPerKnotSpan = clampLower(nSamples, 1);
}

void Shaper::setVertexWeldThreshold(double threshold)
{
  SLIC_WARNING_ROOT_IF(
    threshold <= 0.,
    axom::fmt::format("Vertex weld threshold should be positive Provided value was {}", threshold));

  m_vertexWeldThreshold = threshold;
}

void Shaper::setPercentError(double percent)
{
  using axom::utilities::clampVal;
  SLIC_WARNING_ROOT_IF(percent <= MINIMUM_PERCENT_ERROR,
                       axom::fmt::format("Percent error must be greater than {}. Provided value "
                                         "was {}. Dynamic refinement will not be used.",
                                         MINIMUM_PERCENT_ERROR,
                                         percent));
  SLIC_WARNING_ROOT_IF(
    percent > MAXIMUM_PERCENT_ERROR,
    axom::fmt::format("Percent error must be less than {}. Provided value was {}",
                      MAXIMUM_PERCENT_ERROR,
                      percent));
  if(percent <= MINIMUM_PERCENT_ERROR)
  {
    m_refinementType = DiscreteShape::RefinementUniformSegments;
  }
  m_percentError = clampVal(percent, MINIMUM_PERCENT_ERROR, MAXIMUM_PERCENT_ERROR);
}

void Shaper::setRefinementType(Shaper::RefinementType t) { m_refinementType = t; }

bool Shaper::isValidFormat(const std::string& format) const
{
  static const char* formats[] = {
#if defined(AXOM_USE_MFEM)
    "mfem",
#endif
    "stl",
    "proe",
    "c2c",
    "blueprint-tets",
    "tet3D",
    "hex3D",
    "plane3D",
    "sphere3D",
    "sor3D",
    "none"};
  constexpr auto numFormats = sizeof(formats) / sizeof(const char*);
  const auto formats_end = formats + numFormats;
  return std::find(formats, formats + numFormats, format) != formats_end;
}

void Shaper::loadShape(const klee::Shape& shape)
{
  AXOM_ANNOTATE_SCOPE("loadShape");

  // Do not save the revolved volume in the default shaper.
  double revolved = 0.;
  loadShapeInternal(shape, m_percentError, revolved);
}

void Shaper::loadShapeInternal(const klee::Shape& shape, double percentError, double& revolvedVolume)
{
  AXOM_ANNOTATE_SCOPE("Shaper::loadShapeInternal");
  internal::ScopedLogLevelChanger logLevelChanger(this->isVerbose() ? slic::message::Debug
                                                                    : slic::message::Warning);

  SLIC_INFO_ROOT(
    axom::fmt::format("{:-^80}", axom::fmt::format(" Loading shape '{}' ", shape.getName())));

  SLIC_ERROR_ROOT_IF(
    !this->isValidFormat(this->shapeFormat(shape)),
    axom::fmt::format("Shape has unsupported format: '{}", this->shapeFormat(shape)));

  // Code for discretizing shapes has been factored into DiscreteShape class.
  DiscreteShape discreteShape(shape, m_dataStore.getRoot(), m_prefixPath);
  discreteShape.setSamplesPerKnotSpan(m_samplesPerKnotSpan);
  discreteShape.setVertexWeldThreshold(m_vertexWeldThreshold);
  discreteShape.setRefinementType(m_refinementType);
  if(percentError > 0)
  {
    discreteShape.setPercentError(percentError);
  }
  m_surfaceMesh = discreteShape.createMeshRepresentation();
  revolvedVolume = discreteShape.getRevolvedVolume();
}

bool Shaper::verifyInputMesh(std::string& whyBad) const { return verifyInputMeshImpl(whyBad); }

#if defined(AXOM_USE_CONDUIT)
std::string Shaper::resolveBlueprintTopologyName(const sidre::Group* bpMesh,
                                                 const std::string& topo) const
{
  SLIC_ASSERT(bpMesh != nullptr);
  auto* topologiesGrp = bpMesh->getGroup("topologies");
  SLIC_ERROR_IF(topologiesGrp == nullptr, "Blueprint mesh is missing a 'topologies' group.");

  const std::string topologyName = topo.empty() ? topologiesGrp->getGroupName(0) : topo;
  SLIC_ERROR_IF(topologyName == sidre::InvalidName,
                "Blueprint mesh does not contain any topology groups.");
  SLIC_ERROR_IF(!topologiesGrp->hasGroup(topologyName),
                axom::fmt::format("Blueprint mesh does not contain topology '{}'.", topologyName));

  return topologyName;
}

std::string Shaper::resolveBlueprintTopologyName(const conduit::Node& bpMesh,
                                                 const std::string& topo) const
{
  SLIC_ERROR_IF(!bpMesh.has_path("topologies"), "Blueprint mesh is missing a 'topologies' node.");

  const conduit::Node& topologies = bpMesh.fetch_existing("topologies");
  SLIC_ERROR_IF(topologies.number_of_children() == 0,
                "Blueprint mesh does not contain any topology nodes.");

  const std::string topologyName = topo.empty() ? topologies.child(0).name() : topo;
  SLIC_ERROR_IF(!topologies.has_child(topologyName),
                axom::fmt::format("Blueprint mesh does not contain topology '{}'.", topologyName));

  return topologyName;
}

void Shaper::refreshBlueprintMeshState()
{
  SLIC_ASSERT(m_bp_state != nullptr);
  m_bp_state->refreshBlueprintMeshNode();
  m_cellCount = conduit::blueprint::mesh::topology::length(getBlueprintTopologyNode());
}

const conduit::Node& Shaper::getBlueprintTopologyNode() const
{
  SLIC_ASSERT(m_bp_state != nullptr);
  return m_bp_state->getBlueprintTopologyNode();
}

const conduit::Node& Shaper::getBlueprintCoordsetNode() const
{
  SLIC_ASSERT(m_bp_state != nullptr);
  const std::string coordsetName = getBlueprintTopologyNode().fetch_existing("coordset").as_string();
  return m_bp_state->getBlueprintCoordsetNode(coordsetName);
}

std::string Shaper::getBlueprintCellShape() const
{
  return shaping::getBlueprintCellShape(getBlueprintTopologyNode());
}

int Shaper::getBlueprintMeshDimension() const
{
  const std::string shapeType = getBlueprintCellShape();
  if(shapeType == "quad")
  {
    return 2;
  }
  if(shapeType == "hex")
  {
    return 3;
  }

  SLIC_ERROR(axom::fmt::format("Unsupported Blueprint cell shape '{}'.", shapeType));
  return -1;
}

bool Shaper::verifyBlueprintMeshIsStructuredOrUnstructuredQuadHex(std::string& whyBad) const
{
  bool rval = true;

  if(m_bp_state != nullptr)
  {
    conduit::Node info;
    // Conduit's verify should work even if m_internal_node has array data on
    // devices. because the verification doesn't dereference array data.
    // If this changes in the future, more care must be taken.
    rval = conduit::blueprint::mesh::verify(m_bp_state->getBlueprintMeshNode(), info);
    if(rval)
    {
      const std::string topoType = getBlueprintTopologyNode().fetch_existing("type").as_string();
      rval = topoType == "unstructured" || topoType == "structured";
      info[0].set_string("Topology is not structured or unstructured.");
    }
    if(rval)
    {
      const std::string elemShape = getBlueprintCellShape();
      rval = (elemShape == "hex") || (elemShape == "quad");
      info[0].set_string("Topology elements are not hex or quad.");
    }
    if(rval)
    {
      const std::string coordsetType = getBlueprintCoordsetNode().fetch_existing("type").as_string();
      rval = coordsetType == "explicit";
      info[0].set_string("Coordset is not explicit.");
    }
    whyBad = info.to_summary_string();
  }

  return rval;
}

void Shaper::ensureBlueprintMeshIsUnstructured()
{
  if(m_bp_state == nullptr)
  {
    return;
  }

  const conduit::Node& topoNode = getBlueprintTopologyNode();
  const std::string topoType = topoNode.fetch_existing("type").as_string();
  if(topoType != "structured")
  {
    return;
  }

  if(!m_bp_state->isSidreBacked())
  {
    m_bp_state->ensureUnstructured(m_execPolicy);
    return;
  }

  AXOM_ANNOTATE_SCOPE("Shaper::convertStructured");
  m_bp_state->ensureUnstructured(m_execPolicy);
  m_cellCount = conduit::blueprint::mesh::topology::length(getBlueprintTopologyNode());
}
#endif

std::string Shaper::outputProtocol() const
{
#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
  return "hdf5";
#else
  return "yaml";
#endif
}

#if defined(AXOM_USE_MFEM)
bool Shaper::verifyMFEMInputMesh(std::string& whyBad) const
{
  AXOM_UNUSED_VAR(whyBad);

  if(getDC() != nullptr)
  {
    // No specific requirements for MFEM mesh.
  }

  return true;
}
#endif

void Shaper::saveResults(bool AXOM_UNUSED_PARAM(extra))
{
#ifdef MFEM_USE_MPI
  // If the target mesh was MFEM, save it.
  if(getDC() != nullptr)
  {
    getDC()->Save();
  }
#endif
#if defined(AXOM_USE_CONDUIT)
  // If the target mesh was Blueprint, save it.
  if(m_bp_state != nullptr)
  {
    const std::string filename("shaping");
  #if defined(CONDUIT_RELAY_MPI_ENABLED)
    conduit::relay::mpi::io::blueprint::save_mesh(m_bp_state->getBlueprintMeshNode(),
                                                  filename,
                                                  outputProtocol(),
                                                  m_comm);
  #else
    conduit::relay::io::blueprint::save_mesh(m_bp_state->getBlueprintMeshNode(),
                                             filename,
                                             outputProtocol());
  #endif
  }
#endif
}
// ----------------------------------------------------------------------------

int Shaper::getRank() const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(!mpiIsActive())
  {
    return 0;
  }
  int rank = -1;
  MPI_Comm_rank(m_comm, &rank);
  return rank;
#endif
  return 0;
}

double Shaper::allReduceSum(double val) const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(!mpiIsActive())
  {
    return val;
  }
  double global;
  MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_SUM, m_comm);
  return global;
#else
  return val;
#endif
}

double Shaper::allReduceMin(double val) const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(!mpiIsActive())
  {
    return val;
  }
  double global;
  MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_MIN, m_comm);
  return global;
#else
  return val;
#endif
}

double Shaper::allReduceMax(double val) const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(!mpiIsActive())
  {
    return val;
  }
  double global;
  MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_MAX, m_comm);
  return global;
#else
  return val;
#endif
}

}  // end namespace quest
}  // end namespace axom
