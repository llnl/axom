// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
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

namespace axom
{
namespace quest
{

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
  , m_dc(dc)
  #if defined(AXOM_USE_CONDUIT)
  , m_bpGrp(nullptr)
  , m_bpTopo()
  , m_bpNodeExt(nullptr)
  , m_bpNodeInt()
  #endif
{
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  m_comm = m_dc->GetComm();
  #endif
  m_cellCount = m_dc->GetMesh()->GetNE();

  setFilePath(shapeSet.getPath());
}
#endif

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
#if defined(AXOM_USE_CONDUIT)
  , m_bpGrp(bpGrp)
  , m_bpTopo(topo.empty() ? bpGrp->getGroup("topologies")->getGroupName(0) : topo)
  , m_bpNodeExt(nullptr)
  , m_bpNodeInt()
#endif
#if defined(AXOM_USE_MPI)
  , m_comm(MPI_COMM_WORLD)
#endif
{
  SLIC_ASSERT(m_bpTopo != sidre::InvalidName);

  // This may take too long if there are repeated construction.
  m_bpGrp->createNativeLayout(m_bpNodeInt);

  m_cellCount = conduit::blueprint::mesh::topology::length(
    m_bpNodeInt.fetch_existing("topologies").fetch_existing(m_bpTopo));

  setFilePath(shapeSet.getPath());
}

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
#if defined(AXOM_USE_CONDUIT)
  , m_bpGrp(nullptr)
  , m_bpTopo(topo.empty() ? bpNode.fetch_existing("topologies").child(0).name() : topo)
  , m_bpNodeExt(&bpNode)
  , m_bpNodeInt()
#endif
#if defined(AXOM_USE_MPI)
  , m_comm(MPI_COMM_WORLD)
#endif
{
  AXOM_ANNOTATE_SCOPE("Shaper::Shaper_Node");
  m_bpGrp = m_dataStore.getRoot()->createGroup("internalGrp");
  m_bpGrp->setDefaultArrayAllocator(m_allocatorId);
  m_bpGrp->importConduitTreeExternal(bpNode);

  // We want unstructured topo but can accomodate structured.
  const std::string topoType =
    bpNode.fetch_existing("topologies").fetch_existing(m_bpTopo).fetch_existing("type").as_string();

  if(topoType == "structured")
  {
    AXOM_ANNOTATE_SCOPE("Shaper::convertStructured");
    const std::string shapeType = bpNode.fetch_existing("topologies/mesh/elements/shape").as_string();

    if(shapeType == "hex")
    {
      axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_3d(m_bpGrp,
                                                                                  m_bpTopo,
                                                                                  m_execPolicy);
    }
    else if(shapeType == "quad")
    {
      axom::quest::util::convert_blueprint_structured_explicit_to_unstructured_2d(m_bpGrp,
                                                                                  m_bpTopo,
                                                                                  m_execPolicy);
    }
    else
    {
      SLIC_ERROR("Axom Internal error: Unhandled shape type.");
    }
  }

  m_bpGrp->createNativeLayout(m_bpNodeInt);

  m_cellCount = conduit::blueprint::mesh::topology::length(
    bpNode.fetch_existing("topologies").fetch_existing(m_bpTopo));

  setFilePath(shapeSet.getPath());
}

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

bool Shaper::verifyInputMesh(std::string& whyBad) const
{
  bool rval = true;

#if defined(AXOM_USE_CONDUIT)
  if(m_bpGrp != nullptr)
  {
    conduit::Node info;
    // Conduit's verify should work even if m_bpNodeInt has array data on
    // devices. because the verification doesn't dereference array data.
    // If this changes in the future, more care must be taken.
    rval = conduit::blueprint::mesh::verify(m_bpNodeInt, info);
    if(rval)
    {
      std::string topoType = m_bpNodeInt.fetch("topologies")[m_bpTopo]["type"].as_string();
      rval = topoType == "unstructured";
      info[0].set_string("Topology is not unstructured.");
    }
    if(rval)
    {
      std::string elemShape =
        m_bpNodeInt.fetch("topologies")[m_bpTopo]["elements"]["shape"].as_string();
      rval = (elemShape == "hex") || (elemShape == "quad");
      info[0].set_string("Topology elements are not hex or quad.");
    }
    whyBad = info.to_summary_string();
  }
#endif

#if defined(AXOM_USE_MFEM)
  if(m_dc != nullptr)
  {
    // No specific requirements for MFEM mesh.
  }
#endif

  return rval;
}

// ----------------------------------------------------------------------------

int Shaper::getRank() const
{
#if defined(AXOM_USE_MPI)
  int rank = -1;
  MPI_Comm_rank(m_comm, &rank);
  return rank;
#endif
  return 0;
}

double Shaper::allReduceSum(double val) const
{
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
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
  double global;
  MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_MAX, m_comm);
  return global;
#else
  return val;
#endif
}

}  // end namespace quest
}  // end namespace axom
