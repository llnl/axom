// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file shaping_helpers.hpp
 *
 * \brief Free-standing helper functions in support of shaping query
 */

#ifndef AXOM_QUEST_SHAPING_HELPERS__HPP_
#define AXOM_QUEST_SHAPING_HELPERS__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"

#if defined(AXOM_USE_MFEM)
  #include "mfem.hpp"
  #include "mfem/linalg/dtensor.hpp"
#endif
#if defined(AXOM_USE_CONDUIT)
  #include "conduit_node.hpp"
  #include "axom/bump/utilities/conduit_memory.hpp"
  #include "axom/bump/views/dispatch_coordset.hpp"
#endif

#include <set>
#include <vector>

namespace axom
{

template <typename Signature, size_t MaxSize = 16>
class function;

/**
 * \brief Basic implementation of a host/device compatible analogue to std::function
 *  
 * \tparam R The return type of the callable object
 * \tparam Args The parameter types of the callable object
 * \tparam MaxSize The maximum size of the callable (including its captured variables)
 * 
 * \note We will extend this and move it to the core component
 */
template <typename R, typename... Args, size_t MaxSize>
class function<R(Args...), MaxSize>
{
private:
  using Storage = typename std::aligned_storage<MaxSize>::type;

public:
  AXOM_HOST_DEVICE function() : invoke(nullptr) { }

  /**
   * \brief Constructs a function object from a callable object
   *
   * \tparam Callable The type of the callable object
   * \param callable The callable object to store and invoke
   *
   * This constructor stores the callable object in the internal storage
   * and sets up the invoke function pointer to call the stored object.
   * The callable object must be trivially copyable and its size must not
   * exceed the maximum storage size.
   */
  template <typename Callable>
  AXOM_HOST_DEVICE function(Callable callable)
  {
    static_assert(sizeof(Callable) <= MaxSize, "Callable object too large!");
    static_assert(std::is_trivially_copyable<Callable>::value,
                  "Callable must be trivially copyable!");
    //SLIC_WARNING("sizeof(Callable): " << sizeof(Callable));

    invoke = [](const void* storage, Args... args) -> R {
      return (*reinterpret_cast<const Callable*>(storage))(std::forward<Args>(args)...);
    };
    new(&storage) Callable(std::move(callable));
  }

  /**
   * \brief invoke the stored callable object
   *
   * \param args The arguments to be forwarded to the callable object
   * 
   * \return The result of invoking the callable object with the provided arguments.
   *         If the callable object is not set (i.e., `invoke` is null), a default-constructed
   *         value of type R is returned.
   */
  AXOM_HOST_DEVICE R operator()(Args... args) const
  {
    if(!invoke)
    {
      return R();
    }
    return invoke(&storage, std::forward<Args>(args)...);
  }

  /**
   * \brief Explicit conversion operator to check the validity of the object
   * 
   * \return True if `invoke` is not null, false otherwise
   */
  AXOM_HOST_DEVICE explicit operator bool() const { return invoke != nullptr; }

private:
  Storage storage;

  R (*invoke)(const void*, Args... args) = nullptr;
};

template <typename Lambda>
auto make_host_device_function(Lambda&& lambda)
{
  using Signature = decltype(&Lambda::operator());
  return function<Signature>(std::forward<Lambda>(lambda));
}
namespace quest
{

// clang-format off
using seq_exec = axom::SEQ_EXEC;

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  using omp_exec = axom::OMP_EXEC;
#else
  using omp_exec = seq_exec;
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  constexpr int CUDA_BLOCK_SIZE = 256;
  using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
#else
  using cuda_exec = seq_exec;
#endif

#if defined(AXOM_USE_HIP) && defined(AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  constexpr int HIP_BLOCK_SIZE = 64;
  using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
#else
  using hip_exec = seq_exec;
#endif
// clang-format on

namespace shaping
{

/// Alias to function pointer that projects a \a FromDim dimensional input point to
/// a \a ToDim dimensional query point when sampling the InOut field
template <int FromDim, int ToDim>
using PointProjector =
  axom::function<primal::Point<double, ToDim>(const primal::Point<double, FromDim>&)>;

enum class VolFracSampling : int
{
  SAMPLE_AT_DOFS,
  SAMPLE_AT_QPTS
};

#if defined(AXOM_USE_MFEM)

/*!
 * \brief Converts an Axom quadrature family to the corresponding MFEM
 *        `Quadrature1D` value.
 *
 * \note All `axom::numerics::QuadratureType` enumerators currently map 1:1 to
 *       MFEM names, even when Axom core numerics does not yet implement the
 *       corresponding rule family.
 */
int to_mfem_quadrature_type(axom::numerics::QuadratureType quadratureType);

using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;
using MFEMArrayCollection = mfem::NamedFieldsMap<mfem::Array<int>>;

/*!
 * \brief Contains the mesh and state used for shaping.
 */
struct MFEMState
{
  virtual ~MFEMState() = default;

  // For mesh represented as MFEMSidreDataCollection
  sidre::MFEMSidreDataCollection* m_dc {nullptr};
};

/*!
 * \brief An MFEMState subclass that contains additional data for sampling.
 */
struct SamplingMFEMState : public MFEMState
{
  ~SamplingMFEMState() override
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

  int meshDimension() const
  {
    return m_dc->GetMesh()->Dimension();
  }

  mfem::QuadratureFunction* getShapeFunction(const std::string& name)
  {
    return m_inoutShapeQFuncs.Get(name);
  }

  const mfem::QuadratureFunction* getShapeFunction(const std::string& name) const
  {
    return m_inoutShapeQFuncs.Get(name);
  }

  void deleteShapeFunction(const std::string& AXOM_UNUSED_PARAM(name))
  {
    // TODO: remove the function from m_inoutShapeQFuncs if it exists.
  }

  mfem::QuadratureFunction* getMaterialFunction(const std::string& name)
  {
    return m_inoutMaterialQFuncs.Get(name);
  }

  const mfem::QuadratureFunction* getMaterialFunction(const std::string& name) const
  {
    return m_inoutMaterialQFuncs.Get(name);
  }

  mfem::QuadratureFunction* createMaterialFunction(const std::string& name)
  {
    auto* positions = m_inoutShapeQFuncs.Get("positions");
    SLIC_ERROR_IF(positions == nullptr,
                  std::string("Cannot create material function '") + name +
                    "' without positions.");

    auto* qfunc = new mfem::QuadratureFunction(positions->GetSpace(), 1);
    qfunc->HostWrite();
    *qfunc = 0.;
    m_inoutMaterialQFuncs.Register(name, qfunc, true);
    return qfunc;
  }

  QFunctionCollection m_inoutShapeQFuncs;
  QFunctionCollection m_inoutMaterialQFuncs;
  DenseTensorCollection m_inoutTensors;
  MFEMArrayCollection m_inoutArrays;
};

/**
 * \brief Prints the registered sampling-related field names for an MFEM-backed
 *        sampling state.
 */
void printRegisteredFieldNames(const SamplingMFEMState& mfemState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling vfSampling,
                               const std::string& initialMessage);

/**
 * \brief Utility function to either return a grid function from the DataCollection \a dc, 
 * or to allocate the grud function through the dc, ensuring the memory doesn't leak
 * 
 * \return A pointer to the (allocated) grid function. nullptr if it cannot be allocated
 */
mfem::GridFunction* getOrAllocateL2GridFunction(mfem::DataCollection* dc,
                                                const std::string& gf_name,
                                                int order,
                                                int dim,
                                                const int basis);

/**
 * Utility function to zero out inout quadrature points for a material replaced by a shape
 *
 * Each location in space can only be covered by one material.
 * When \a shouldReplace is true, we clear all values in \a materialQFunc 
 * that are set in \a shapeQFunc. When it is false, we do the opposite.
 *
 * \param shapeQFunc The inout quadrature function for the shape samples
 * \param materialQFunc The inout quadrature function for the material samples
 * \param shapeReplacesMaterial Flag for whether the shape replaces the material 
 *   or whether the material remains and we should zero out the shape sample (when false)
 */
void replaceMaterial(mfem::QuadratureFunction* shapeQFunc,
                     mfem::QuadratureFunction* materialQFunc,
                     bool shouldReplace);

/**
 * \brief Utility function to copy inout quadrature point values from \a shapeQFunc to \a materialQFunc
 *
 * \param shapeQFunc The inout samples for the current shape
 * \param materialQFunc The inout samples for the material we're writing into
 * \param reuseExisting When a value is not set in \a shapeQFunc, should we retain existing values 
 * from \a materialQFunc or overwrite them based on \a shapeQFunc. The default is to retain values
 */
void copyShapeIntoMaterial(const mfem::QuadratureFunction* shapeQFunc,
                           mfem::QuadratureFunction* materialQFunc,
                           bool reuseExisting = true);

mfem::QuadratureFunction* cloneInOutFunction(const mfem::QuadratureFunction* qfunc);

/**
 * \brief Generates a "position" quadrature function corresponding to the mesh positions and
 *        store it in \a inoutQFuncs.
 *
 * \param mesh The mesh
 * \param inoutQFuncs A collection of quadrature functions where the new "position" function will be added.
 * \param sampleResolution The sample resolution in each logical dimension.
 * \param quadratureType An int corresponding to mfem::Quadrature1D enum values. If
 *                       Invalid is used then the default quadrature is constructed.
 *                       Otherwise, custom quadrature is constructed using the supplied
 *                       quadratureType -- the same type per dimension but the sampling
 *                       can vary.
 */
void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                int sampleResolution[3],
                                axom::numerics::QuadratureType quadratureType);

/**
 * \brief Generates a "position" quadrature function for the supplied MFEM state.
 */
void generateSamplingPositions(SamplingMFEMState& mfemState,
                               int sampleResolution[3],
                               axom::numerics::QuadratureType quadratureType);

void computeVolumeFractionsForMaterial(SamplingMFEMState& mfemState,
                                       const std::string& matField,
                                       int volfracOrder,
                                       int sampleResolution[3],
                                       axom::numerics::QuadratureType quadratureType);
/*!
 * \brief Identity transform for volume fractions from inout samples
 *
 * Copies \a inout samples from the quadrature function directly into volume fraction DOFs.
 * \param dc The data collection to which we will add the volume fractions
 * \param inout The inout samples
 * \param name The name of the generated volume fraction function
 * \note Assumes that the inout samples are co-located with the grid function DOFs.
 */
void computeVolumeFractionsIdentity(mfem::DataCollection* dc,
                                    mfem::QuadratureFunction* inout,
                                    const std::string& name);

/*!
  * \brief Samples the inout field over the indexed geometry, possibly using a
  * callback function to project the input points (from the computational mesh)
  * to query points on the spatial index
  *
  * \tparam FromDim The dimension of points from the input mesh
  * \tparam ToDim The dimension of points on the indexed shape
  * \tparam InsideFunc A function that takes a point and returns a bool indicating whether the
  *                    point is inside or outside of relevant shapes.
  *
  * \param [in] shapeName The name of the shape used in making data array names.
  * \param [in] mfemState The MFEM state containing the mesh and associated query points
  * \param [inout] inoutQFuncs A collection of quadrature functions for the shape and material
  * inout samples
  * \param [in] sampleRes The sampling resolution in each logical direction.
  * For custom quadrature families, these values specify the per-direction
  * sample counts directly, which in turn determine the quadrature rule used
  * in each logical direction.
  * \param [in] quadratureType The quadrature type to use to construct the sample point locations.
  * \param [in] checkInside The function that determines whether a point is inside.
  * \param [in] projector A callback function to apply to points from the input mesh
  * before querying them on the spatial index
  *
  * \note A projector callback must be supplied when \a FromDim is not equal
  *       to \a ToDim.
  */
template <int FromDim, int ToDim, typename InsideFunc>
void sampleInOutField(const std::string shapeName,
                      shaping::SamplingMFEMState& mfemState,
                      int AXOM_UNUSED_PARAM(sampleRes)[3],
                      int AXOM_UNUSED_PARAM(quadratureType),
                      InsideFunc&& checkInside,
                      PointProjector<FromDim, ToDim> projector = {})
{
  using FromPoint = primal::Point<double, FromDim>;
  using ToPoint = primal::Point<double, ToDim>;
  AXOM_ANNOTATE_SCOPE("sampleInOutField");

  SLIC_ERROR_IF(FromDim != ToDim && !projector,
                "A projector callback function is required when FromDim != ToDim");

  auto* mesh = mfemState.m_dc->GetMesh();
  SLIC_ASSERT(mesh != nullptr);
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  auto& inoutQFuncs = mfemState.m_inoutShapeQFuncs;
  SLIC_ASSERT(inoutQFuncs.Has("positions"));

  // Access the positions QFunc and associated QuadratureSpace
  mfem::QuadratureFunction* pos_coef = inoutQFuncs.Get("positions");
  auto* sp = pos_coef->GetSpace();
  const int nq = sp->GetIntRule(0).GetNPoints();
  const int numQueryPoints = sp->GetSize();
  SLIC_ASSERT(numQueryPoints == NE * nq);

  const auto pos = mfem::Reshape(pos_coef->HostRead(), dim, nq, NE);

  // Sample the in/out field at each point
  // store in QField which we register with the QFunc collection
  const std::string inoutName = axom::fmt::format("inout_{}", shapeName);
  const int vdim = 1;
  auto* inout = new mfem::QuadratureFunction(sp, vdim);
  inoutQFuncs.Register(inoutName, inout, true);
  auto inout_vals = mfem::Reshape(inout->HostWrite(), nq, NE);

  axom::utilities::Timer timer(true);
  if(projector)
  {
    for(int i = 0; i < NE; ++i)
    {
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt = projector(FromPoint(&pos(0, p, i), dim));
        inout_vals(p, i) = checkInside(pt) ? 1. : 0.;
      }
    }
  }
  else
  {
    for(int i = 0; i < NE; ++i)
    {
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt(&pos(0, p, i), dim);
        inout_vals(p, i) = checkInside(pt) ? 1. : 0.;
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
    static_cast<int>(numQueryPoints / timer.elapsed())));
}

/*!
  * \brief Samples the inout field over the indexed geometry, possibly using a
  * callback function to project the input points (from the computational mesh)
  * to query points on the spatial index
  *
  * \tparam FromDim The dimension of points from the input mesh
  * \tparam ToDim The dimension of points on the indexed shape
  * \tparam InsideFunc A function that takes a point and returns a bool indicating whether the
  *                    point is inside or outside of relevant shapes.
  *
  * \param [in] shapeName The name of the shape used in making data array names.
  * \param [in] dc The data collection containing the mesh and associated query points
  * \param [in] outputOrder The order of the output inout field
  * \param [in] checkInside The function that determines whether a point is inside.
  * \param [in] projector A callback function to apply to points from the input mesh
  *             before querying them on the spatial index
  *
  * \note A projector callback must be supplied when \a FromDim is not equal
  *       to \a ToDim.
  */
template <int FromDim, int ToDim, typename InsideFunc>
void computeVolumeFractionsBaseline(const std::string& shapeName,
                                    shaping::SamplingMFEMState& mfemState,
                                    int outputOrder,
                                    InsideFunc&& checkInside,
                                    PointProjector<FromDim, ToDim> projector = {})
{
  using FromPoint = primal::Point<double, FromDim>;
  using ToPoint = primal::Point<double, ToDim>;
  AXOM_ANNOTATE_SCOPE("computeVolumeFractionsBaseline");

  // Step 1 -- generate a QField w/ the spatial coordinates
  mfem::DataCollection* dc = mfemState.m_dc;
  mfem::Mesh* mesh = dc->GetMesh();
  const int NE = mesh->GetNE();
  const int dim = mesh->Dimension();

  if(NE < 1)
  {
    SLIC_WARNING("Mesh has no elements!");
    return;
  }

  const auto volFracName = axom::fmt::format("vol_frac_{}", shapeName);
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
  if(projector)
  {
    for(int i = 0; i < NE; ++i)
    {
      const mfem::DenseMatrix& m = pos_coef(i);
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt = projector(FromPoint(m.GetColumn(p), dim));
        res(p) = checkInside(pt) ? 1. : 0.;
      }

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, res);
    }
  }
  else
  {
    for(int i = 0; i < NE; ++i)
    {
      const mfem::DenseMatrix& m = pos_coef(i);
      for(int p = 0; p < nq; ++p)
      {
        const ToPoint pt(m.GetColumn(p), dim);
        res(p) = checkInside(pt) ? 1. : 0.;
      }

      fes->GetElementDofs(i, dofs);
      volFrac->SetSubVector(dofs, res);
    }
  }
}
#endif

#if defined(AXOM_USE_CONDUIT)
//------------------------------------------------------------------------------
/**
 * \brief Returns the element shape for a supported Blueprint topology node.
 *
 * Structured topologies may omit `elements/shape`, in which case the shape is
 * inferred from `elements/dims`.
 */
std::string getBlueprintCellShape(const conduit::Node& topoNode);

/*!
 * \brief A Blueprint-based state class used for shaping.
 */
struct BlueprintState
{
  virtual ~BlueprintState() = default;

  //! @brief Version of the mesh for computations.
  axom::sidre::Group* m_group_ptr {nullptr};
  int m_allocator_id {axom::getDefaultAllocatorID()};
  std::string m_topology_name;
  //! @brief Mesh in an external Node, when provided as a Node.
  conduit::Node* m_external_node_ptr {nullptr};
  //! @brief Internal Node representation used for blueprint operations.
  conduit::Node m_internal_node;

  int meshDimension() const
  {
    const std::string shapeType = shaping::getBlueprintCellShape(getBlueprintTopologyNode());

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

  const conduit::Node& getBlueprintTopologyNode() const
  {
    return m_internal_node.fetch_existing("topologies").fetch_existing(m_topology_name);
  }

  conduit::Node* getShapeFunction(const std::string& name)
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

  const conduit::Node* getShapeFunction(const std::string& name) const
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

  void deleteShapeFunction(const std::string& name)
  {
    // This method lets us delete the shape functions as we go
    if(m_internal_node.has_path("fields"))
    {
      conduit::Node &n_fields = m_internal_node["fields"];
      if(n_fields.has_path(name))
      {
        n_fields.remove(name);
      }
    }
  }

  conduit::Node* getMaterialFunction(const std::string& name)
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

  const conduit::Node* getMaterialFunction(const std::string& name) const
  {
    return m_internal_node.has_path("fields/" + name) ? &m_internal_node["fields/" + name]
                                                      : nullptr;
  }

  conduit::Node* createMaterialFunction(const std::string& name)
  {
    constexpr const char* quadratureTopologyName = "quadrature_points";
    SLIC_ERROR_IF(!m_internal_node.has_path("coordsets/quadrature_points/values"),
                  std::string("Cannot create material function '") + name +
                    "' without quadrature points.");

    conduit::Node& fieldNode = m_internal_node["fields/" + name];
    fieldNode.reset();
    fieldNode["association"] = "element";
    fieldNode["topology"] = quadratureTopologyName;

    const auto conduitAllocatorId =
      axom::sidre::ConduitMemory::axomAllocIdToConduit(m_allocator_id);
    conduit::Node& valuesNode = fieldNode["values"];
    valuesNode.set_allocator(conduitAllocatorId);

    const conduit::Node& values =
      m_internal_node["coordsets/quadrature_points"].fetch_existing("values");
    const auto numValues = values.child(0).dtype().number_of_elements();
    valuesNode.set(conduit::DataType::float64(numValues));

    auto fieldValues = axom::bump::utilities::make_array_view<double>(valuesNode);
    for(axom::IndexType i = 0; i < fieldValues.size(); ++i)
    {
      fieldValues[i] = 0.;
    }

    return &fieldNode;
  }
};

/**
 * \brief Prints the registered sampling-related field names for a Blueprint-backed
 *        sampling state.
 */
void printRegisteredFieldNames(const BlueprintState& bpState,
                               const std::set<std::string>& knownMaterials,
                               VolFracSampling vfSampling,
                               const std::string& initialMessage);

void replaceMaterial(conduit::Node* shapeNode, conduit::Node* materialNode, bool shouldReplace);

void copyShapeIntoMaterial(const conduit::Node* shapeNode,
                           conduit::Node* materialNode,
                           bool reuseExisting = true);

conduit::Node* cloneInOutFunction(const conduit::Node* node);

/**
 * \brief Generates a derived Blueprint quadrature point mesh within the
 *        supplied Blueprint mesh node.
 *
 * \param bpMeshNode The Blueprint mesh node to augment.
 * \param topologyName The source topology name to sample.
 * \param allocatorID Allocator id used for generated storage.
 * \param sampleResolution The sample resolution in each logical dimension.
 * \param quadratureType An int corresponding to `mfem::Quadrature1D` when MFEM
 *        is enabled, or to `axom::numerics::QuadratureType` otherwise.
 */
void generateQuadraturePointMesh(conduit::Node& bpMeshNode,
                                 const std::string& topologyName,
                                 int allocatorID,
                                 int sampleResolution[3],
                                 axom::numerics::QuadratureType quadratureType);

/**
 * \brief Generates a derived Blueprint quadrature point mesh for the supplied
 *        Blueprint state.
 */
void generateSamplingPositions(BlueprintState& bpState,
                               int sampleResolution[3],
                               axom::numerics::QuadratureType quadratureType);

void computeVolumeFractionsForMaterial(BlueprintState& bpState, const std::string& matField);

/*!
  * \brief Samples the inout field over the indexed geometry, possibly using a
  * callback function to project the input points (from the computational mesh)
  * to query points on the spatial index
  *
  * \tparam FromDim The dimension of points from the input mesh
  * \tparam ToDim The dimension of points on the indexed shape
  * \tparam InsideFunc A function that takes a point and returns a bool indicating whether the
  *                    point is inside or outside of relevant shapes.
  *
  * \param [in] shapeName The name of the shape used in making data array names.
  * \param [in] bpState The Blueprint state containing the mesh and associated query points
  * \param [in] sampleRes The sampling resolution in each logical direction.
  * For custom quadrature families, these values specify the per-direction
  * sample counts directly, which in turn determine the quadrature rule used
  * in each logical direction.
  * \param [in] quadratureType The quadrature type to use to construct the sample point locations.
  * \param [in] checkInside The function that determines whether a point is inside.
  * \param [in] projector A callback function to apply to points from the input mesh
  * before querying them on the spatial index
  *
  * \note A projector callback must be supplied when \a FromDim is not equal
  *       to \a ToDim.
  */
template <int FromDim, int ToDim, typename InsideFunc>
void sampleInOutField(const std::string& shapeName,
                      shaping::BlueprintState& bpState,
                      int AXOM_UNUSED_PARAM(sampleRes)[3],
                      int AXOM_UNUSED_PARAM(quadratureType),
                      InsideFunc&& checkInside,
                      PointProjector<FromDim, ToDim> projector = {})
{
  using FromPoint = primal::Point<double, FromDim>;
  using ToPoint = primal::Point<double, ToDim>;
  AXOM_ANNOTATE_SCOPE("sampleInOutField");

  SLIC_ERROR_IF(FromDim != ToDim && !projector,
                "A projector callback function is required when FromDim != ToDim");

  constexpr const char* quadratureCoordsetName = "quadrature_points";
  constexpr const char* quadratureTopologyName = "quadrature_points";
  const std::string inoutName = axom::fmt::format("inout_{}", shapeName);

  conduit::Node& bpMeshNode = bpState.m_internal_node;
  SLIC_ERROR_IF(!bpMeshNode.has_path("coordsets/quadrature_points"),
                "Missing Blueprint quadrature coordset. Generate sampling positions first.");
  SLIC_ERROR_IF(!bpMeshNode.has_path("topologies/quadrature_points"),
                "Missing Blueprint quadrature topology. Generate sampling positions first.");

  conduit::Node& inoutNode = bpMeshNode["fields/" + inoutName];
  inoutNode.reset();
  inoutNode["association"] = "element";
  inoutNode["topology"] = quadratureTopologyName;

  namespace utils = axom::bump::utilities;
  const auto conduitAllocatorId =
    axom::sidre::ConduitMemory::axomAllocIdToConduit(bpState.m_allocator_id);
  conduit::Node& valuesNode = inoutNode["values"];
  valuesNode.set_allocator(conduitAllocatorId);

  axom::utilities::Timer timer(true);
  axom::IndexType numQueryPoints = 0;
  axom::bump::views::dispatch_explicit_coordset(
    bpMeshNode["coordsets/" + std::string(quadratureCoordsetName)], [&](auto coordsetView) {
      using CoordsetView = typename std::decay<decltype(coordsetView)>::type;

      SLIC_ERROR_IF(CoordsetView::dimension() != FromDim,
                    axom::fmt::format("Expected {}D quadrature point coordset, got {}D.",
                                      FromDim,
                                      CoordsetView::dimension()));

      numQueryPoints = coordsetView.size();
      valuesNode.set(conduit::DataType::float64(numQueryPoints));
      auto inoutValues = utils::make_array_view<double>(valuesNode);

      for(axom::IndexType i = 0; i < numQueryPoints; ++i)
      {
        // Make a FromPoint from the coordsetView. The coordsetView might have
        // float or double, depending on the Blueprint data.
        FromPoint fromPt;
        const auto coordsetPoint = coordsetView[i];
        for(int d = 0; d < FromDim; ++d)
        {
          fromPt[d] = coordsetPoint[d];
        }

        // Sample at the query point.
        const ToPoint queryPt = projector ? projector(fromPt) : ToPoint(fromPt.data());
        inoutValues[i] = checkInside(queryPt) ? 1. : 0.;
      }
    });
  timer.stop();

  SLIC_INFO_ROOT(axom::fmt::format(
    axom::utilities::locale(),
    "\t Sampling inout field '{}' took {:.3Lf} seconds (@ {:L} queries per second)",
    inoutName,
    timer.elapsed(),
    static_cast<int>(numQueryPoints / timer.elapsed())));
}
#endif

/** 
 * Implements flux-corrected transport (FCT) to correct the solution obtained
 * when converting from inout samples (ones and zeros) to a grid function 
 * on the degrees of freedom such that the volume fractions are doubles
 * between 0 and 1 ( \a y_min and \a y_max )
 */
void FCT_correct(const double* M,
                 const int s,
                 const double* m,
                 const double y_min,  // 0
                 const double y_max,  // 1
                 double* xy,
                 double* fct_mat);  // scratch buffer

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPING_HELPERS__HPP_
