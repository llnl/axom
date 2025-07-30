// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_MESH_VIEW_UTIL_H_
#define QUEST_MESH_VIEW_UTIL_H_

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT

  #include "axom/core.hpp"
  #include "axom/core/NumericLimits.hpp"
  #include "axom/fmt.hpp"
  #include "axom/slic.hpp"
  #include "conduit_blueprint.hpp"
  #include "conduit_blueprint_mcarray.hpp"
  #ifdef AXOM_USE_MPI
    #include "conduit_blueprint_mpi.hpp"
    #include "conduit_relay_mpi_io_blueprint.hpp"
  #endif

  #include <memory>
  #include <cstdlib>
  #include <cmath>
  #include <vector>

namespace axom
{
namespace quest
{
namespace internal
{
template <typename T, int DIM>
inline axom::StackArray<T, DIM> makeStackArray(T v = axom::numeric_limits<T>::max())
{
  axom::StackArray<T, DIM> rval;
  for(int d = 0; d < DIM; ++d)
  {
    rval[d] = v;
  }
  return rval;
}
template <typename T, int DIM, typename U>
inline axom::StackArray<T, DIM> makeStackArray(const U* v)
{
  axom::StackArray<T, DIM> rval;
  for(int d = 0; d < DIM; ++d)
  {
    rval[d] = v[d];
  }
  return rval;
}
template <typename T, int DIM>
static inline T product(const axom::StackArray<T, DIM>& a)
{
  T prod = a[0];
  for(int d = 1; d < DIM; ++d)
  {
    prod *= a[d];
  }
  return prod;
}

///@{
//!@name Conversions between shape specifications and strides-and-offsets.

/*!
 * @brief Convert shape specifications to blueprint-style offsets and strides.
 *
 * @tparam IType Index type
 * @tparam DIM Spatial dimension
 *
 * @param realShape [i]
 * @param loPads [i] Ghost padding amount on low side.
 * @param hiPads [i] Ghost padding amount ont high side.
 * @param strideOrder [i] Fastest-to-slowest advancing index directions.
 * @param minStride [i] Stride of fastest advancing index direction.
 * @param offsets [o] Blueprint-style index offsets.
 * @param strides [o] Blueprint-style strides.
 * @param valuesCount [o] Number of values in ghost-padded data.
 */
template <typename IType, int DIM>
static void shapesToStridesAndOffsets(const axom::StackArray<IType, DIM>& realShape,
                                      const axom::StackArray<IType, DIM>& loPads,
                                      const axom::StackArray<IType, DIM>& hiPads,
                                      const axom::StackArray<IType, DIM>& strideOrder,
                                      IndexType minStride,
                                      axom::StackArray<IType, DIM>& offsets,
                                      axom::StackArray<IType, DIM>& strides,
                                      IndexType& valuesCount)
{
  axom::StackArray<IType, DIM> paddedShape;
  for(int d = 0; d < DIM; ++d)
  {
    offsets[d] = loPads[d];
    paddedShape[d] = realShape[d] + loPads[d] + hiPads[d];
  }

  strides[strideOrder[0]] = minStride;
  for(int nd = 1; nd < DIM; ++nd)
  {
    const int& curDir = strideOrder[nd];
    const int& prevDir = strideOrder[nd - 1];
    strides[curDir] = strides[prevDir] * paddedShape[prevDir];
  }
  auto slowestDir = strideOrder[DIM - 1];
  valuesCount = strides[slowestDir] * paddedShape[slowestDir];
}

/*!
 * @brief Convert blueprint-style offsets and strides to shape specifications.
 *
 * @tparam IType Index type
 * @tparam DIM Spatial dimension
 *
 * @param realShape [i]
 * @param offsets [i] Blueprint-style index offsets.
 * @param strides [i] Blueprint-style strides.
 * @param valuesCount [i] Number of values in ghost-padded data.
 * @param paddedShape [o] \a realShape + \a loPads + \a hiPads
 * @param loPads [o] Ghost padding amount on low side.
 * @param hiPads [o] Ghost padding amount ont high side.
 * @param minStride [i] Stride of fastest advancing index direction.
 * @param strideOrder [i] Fastest-to-slowest advancing index directions.
 */
template <typename IType, int DIM>
static void stridesAndOffsetsToShapes(const axom::StackArray<IType, DIM>& realShape,
                                      const axom::StackArray<IType, DIM>& offsets,
                                      const axom::StackArray<IType, DIM>& strides,
                                      const IndexType& valuesCount,
                                      axom::StackArray<IType, DIM>& paddedShape,
                                      axom::StackArray<IType, DIM>& loPads,
                                      axom::StackArray<IType, DIM>& hiPads,
                                      axom::StackArray<IType, DIM>& strideOrder)
{
  // Sort directions from fastest to slowest.
  for(int d = 0; d < DIM; ++d)
  {
    strideOrder[d] = d;
  }
  for(int s = 0; s < DIM; ++s)
  {
    for(int d = s; d < DIM; ++d)
    {
      if(strides[strideOrder[d]] < strides[strideOrder[s]])
      {
        std::swap(strideOrder[s], strideOrder[d]);
      }
    }
  }

  for(int nd = 0; nd < DIM - 1; ++nd)
  {
    const int& curDir = strideOrder[nd];
    const int& nextDir = strideOrder[nd + 1];
    paddedShape[curDir] = strides[nextDir] / strides[curDir];
  }
  paddedShape[strideOrder[DIM - 1]] = valuesCount / strides[strideOrder[DIM - 1]];

  for(int d = 0; d < DIM; ++d)
  {
    loPads[d] = offsets[d];
    hiPads[d] = paddedShape[d] - realShape[d] - loPads[d];
  }
}
///@}

}  // namespace internal

/**
 * \brief Utility for high-level access into a blueprint mesh,
 * for structured mesh with explicit coordinates.
 *
 * Note: This class was written for a specific use and supports
 * only structured domains with explicit node values.
 *
 * Blueprint mesh data is sufficient but sparse, leaving users to
 * compute some intermediate data to get to high-level data of
 * interest, such as views into array data.  This class encapsulates
 * those common functions.
 *
 * Views are single-domain-specific.  They don't apply to multi-domain
 * meshes.  They are also topology specific, with the topology name
 * given to the constructor.  They are valid only while their domain
 * exists with no change to its data layout.
 *
 * This class recognizes potential ghost (a.k.a. phony, image) data
 * layers around the domain.  Some methods and paramenters names refer
 * to the data with ghosts, while others refer to the data without
 * ghosts (a.k.a. real data).
 *
 * TODO: Figure out if there's a better place for this utility.
 * It's only in axom/quest because the initial need was there.
 */
template <int DIM, axom::MemorySpace MemSpace = MemorySpace::Dynamic>
class MeshViewUtil
{
public:
  using MdIndices = axom::StackArray<axom::IndexType, DIM>;
  using CoordsViewsType = axom::StackArray<axom::ArrayView<double, DIM, MemSpace>, DIM>;
  using ConstCoordsViewsType = axom::StackArray<axom::ArrayView<const double, DIM, MemSpace>, DIM>;

  ///@{
  /// @name Setting up

  //! @brief Default constructor doesn't set up the view utility.
  MeshViewUtil()
    : m_dom(nullptr)
    , m_topology(nullptr)
    , m_coordset(nullptr)
    , m_cdom(nullptr)
    , m_ctopology(nullptr)
    , m_ccoordset(nullptr)
  { }

  /*!
   * @brief Construct view of a non-const domain.
   *
   * @param [in] bpDomain Blueprint single domain.
   * @param [in] topologyName Name of topology in the domain.
   *
   * If \a topologyName is omitted, use the first topology.
   *
   * @pre The topology dimension must match DIM.
  */
  MeshViewUtil(conduit::Node& bpDomain, const std::string& topologyName = "")
    : m_dom(&bpDomain)
    , m_topologyName(topologyName)
  {
    if(m_topologyName.empty())
    {
      m_topologyName = bpDomain.fetch_existing("topologies")[0].name();
    }

    const std::string topoPath = "topologies/" + m_topologyName;
    m_topology = &m_dom->fetch_existing(topoPath);

    const std::string coordsetPath =
      "coordsets/" + m_dom->fetch_existing(topoPath + "/coordset").as_string();
    m_coordset = &m_dom->fetch_existing(coordsetPath);

    m_cdom = m_dom;
    m_ctopology = m_topology;
    m_ccoordset = m_coordset;

    if(!isValid())
    {
      SLIC_ERROR(
        "Invalid domain in MeshViewUtil.  Domain must be blueprint-valid,"
        " structured, with exlicit coordinates and dimension equal to "
        " the template DIM parameter.");
    }

    computeCoordsDataLayout();
  }

  /*!
   * @brief Construct view of a const domain.
   *
   * @param [in] bpDomain const Blueprint single domain.
   * @param [in] topologyName Name op topology in the domain.
   *
   * If \a topologyName is omitted, use the first topology.
   * @pre The topology dimension must match DIM.
  */
  MeshViewUtil(const conduit::Node& bpDomain, const std::string& topologyName)
    : m_cdom(&bpDomain)
    , m_topologyName(topologyName)
  {
    const std::string topoPath = "topologies/" + m_topologyName;
    SLIC_ASSERT(m_cdom->has_path(topoPath));
    m_ctopology = &m_cdom->fetch_existing(topoPath);

    const std::string coordsetPath =
      "coordsets/" + m_cdom->fetch_existing(topoPath + "/coordset").as_string();
    m_ccoordset = &m_cdom->fetch_existing(coordsetPath);

    if(!isValid())
    {
      SLIC_ERROR(
        "Invalid domain in MeshViewUtil.  Domain must be blueprint-valid,"
        " structured, with exlicit coordinates. and dimension equal to "
        " the template DIM parameter.");
    }

    computeCoordsDataLayout();
  }

  /*!
  * @brief Whether the domain is a valid one for this utility.
  * In addition to checking for blueprint validity (unless \a checkBlueprint is false), 
  * this method checks against its own requirements.
  */
  bool isValid(bool checkBlueprint = false, bool printFailureCause = false) const
  {
    bool valA = true;
    if(checkBlueprint)
    {
      conduit::Node info;
  #ifdef AXOM_USE_MPI
      valA = conduit::blueprint::mpi::verify("mesh", *m_cdom, info, MPI_COMM_WORLD);
  #else
      valA = conduit::blueprint::verify("mesh", *m_cdom, info);
  #endif
      if(printFailureCause && !valA)
      {
        info.print();
      }
    }
    bool valB = m_ctopology->fetch_existing("type").as_string() == "structured";
    if(printFailureCause && !valB)
    {
      SLIC_INFO("MeshViewUtil domain is not structured");
    }

    bool valC = m_ccoordset->fetch_existing("type").as_string() == "explicit";
    if(printFailureCause && !valC)
    {
      SLIC_INFO("MeshViewUtil domain coords is not explicit");
    }

    bool valD = conduit::blueprint::mesh::coordset::dims(*m_ccoordset) == DIM;
    if(printFailureCause && !valD)
    {
      SLIC_INFO("MeshViewUtil domain has wrong dimension");
    }

    return valA && valB && valC && valD;
  }
  ///@}

  ///@{
  //!@name Blueprint data organization
  /*!
   * @brief Get the named topology.
   * 
   * @return the Conduit node at "topologies/<topologyName>".
  */
  const conduit::Node& getTopology() { return *m_topology; }
  const conduit::Node& getTopology() const { return *m_ctopology; }

  //! @brief Get the coordinates conduit::Node for the named topology.
  conduit::Node& getCoordSet() { return *m_coordset; }
  const conduit::Node& getCoordSet() const { return *m_ccoordset; }

  //! @brief Get the spatial dimension of the named topology.
  conduit::index_t getTopologyDim() const { return DIM; }
  ///@}

  ///@{
  //! @name General sizes and shapes of domain

  //! @brief Get the number of cells in each direction of the domain.
  MdIndices getCellShape() const { return m_cellShape; }

  //! @brief Get the number of nodes in each direction of the domain.
  MdIndices getNodeShape() const { return m_nodeShape; }

  //! @brief Return number of (real) cells.
  axom::IndexType getCellCount() const { return axom::quest::internal::product(m_cellShape); }

  //! @brief Return number of (real) nodes.
  axom::IndexType getNodeCount() const { return axom::quest::internal::product(m_nodeShape); }

  //! @brief Return the real (ghost-free) shape of mesh data.
  const MdIndices& getRealShape(const std::string& association)
  {
    if(association == "vertex")
    {
      return m_nodeShape;
    }
    else if(association == "element")
    {
      return m_cellShape;
    }

    SLIC_ERROR(
      axom::fmt::format("MeshVieuUtil only supports element and vertex data "
                        "association for now, not '{}'.",
                        association));

    return m_cellShape;
  }
  ///@}

  ///@{
  //! @name Coordinates and field data sizes and shapes.

  /// @brief Return the array strides for ghost-free nodal coordinates.
  MdIndices getCoordsStrides() const { return m_coordsStrides; }

  /// @brief Return the array index offsets for ghost-free nodal coordinates.
  MdIndices getCoordsOffsets() const { return m_coordsOffsets; }

  //! @brief Return number of points, excluding ghosts.
  axom::IndexType getCoordsCount() const { return getNodeCount(); }

  //! @brief Return number of points, including ghosts.
  axom::IndexType getCoordsCountWithGhosts() const
  {
    const conduit::Node& valuesNode = m_ccoordset->fetch_existing("values");
    axom::IndexType rval = valuesNode[0].dtype().number_of_elements();
    return rval;
  }
  ///@}

  ///@{
  //! @name Coords and data views

  //! @brief Return the views of the DIM coordinates component data.
  CoordsViewsType getCoordsViews(bool withGhosts = false)
  {
    conduit::Node& valuesNode = getCoordSet().fetch_existing("values");
    axom::StackArray<axom::ArrayView<double, DIM, MemSpace>, DIM> rval;
    for(int d = 0; d < DIM; ++d)
    {
      auto* dataPtr = valuesNode[d].as_double_ptr();
      rval[d] = axom::ArrayView<double, DIM, MemSpace>(dataPtr, m_coordsPaddedShape, m_coordsStrides);
    }

    if(withGhosts == false)
    {
      // Compute a view without ghosts.
      MdIndices offsets =
        conduitIndicesToStackArray(m_ctopology->fetch_existing("elements/dims"), "offsets", 0);
      MdIndices counts = m_nodeShape;
      for(int d = 0; d < DIM; ++d)
      {
        auto rval1 = rval[d];
        rval[d] = rval1.subspan(offsets, counts);
      }
    }

    return rval;
  }

  //! @brief Return the views of the DIM coordinates component data.
  ConstCoordsViewsType getConstCoordsViews(bool withGhosts = false) const
  {
    const conduit::Node& valuesNode = getCoordSet().fetch_existing("values");
    axom::StackArray<axom::ArrayView<const double, DIM, MemSpace>, DIM> rval;
    for(int d = 0; d < DIM; ++d)
    {
      auto* dataPtr = valuesNode[d].as_double_ptr();
      rval[d] =
        axom::ArrayView<const double, DIM, MemSpace>(dataPtr, m_coordsPaddedShape, m_coordsStrides);
    }

    if(withGhosts == false)
    {
      MdIndices offsets =
        conduitIndicesToStackArray(m_ctopology->fetch_existing("elements/dims"), "offsets", 0);
      MdIndices counts = m_cellShape;
      for(int d = 0; d < DIM; ++d)
      {
        ++counts[d];
      }
      for(int d = 0; d < DIM; ++d)
      {
        auto rval1 = rval[d];
        rval[d] = rval1.subspan(offsets, counts);
      }
    }

    return rval;
  }

  /*!
   * @brief Return view to a scalar field variable.
   * 
   * @warning view returned has an allocator id determined by
   * \a MemSpace, regardless of the memory type.
   * 
   * @warning Assuming, without checking, that the field contains
   * data of type \a T.  User is responsible for using the correct type.
   */
  template <typename T>
  axom::ArrayView<T, DIM, MemSpace> getFieldView(const std::string& fieldName, bool withGhosts = false)
  {
    if(!m_dom)
    {
      SLIC_ERROR(
        "Cannot use getFieldView from MeshViewUtil initialized with a const "
        "conduit::Node.  Use getConstFieldView.");
    }
    conduit::Node& fieldNode = m_dom->fetch_existing("fields/" + fieldName);
    const std::string association = fieldNode.fetch_existing("association").as_string();
    conduit::Node& valuesNode = fieldNode.fetch_existing("values");
    const axom::IndexType valuesCount = valuesNode.dtype().number_of_elements();

    SLIC_ASSERT_MSG(association == "vertex" || association == "element",
                    "MeshViewUtil only supports vertex and element-based fields right now.");
    const bool onVertex = association == "vertex";

    const auto& realShape = onVertex ? m_nodeShape : m_cellShape;

    MdIndices strides = conduitIndicesToStackArray(fieldNode, "strides");
    if(!fieldNode.has_child("strides"))
    {
      axom::IndexType tmpStride = 1;
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = tmpStride;
        tmpStride *= m_cellShape[d] + onVertex;
      }
    }

    MdIndices offsets = conduitIndicesToStackArray(fieldNode, "offsets", 0);

    MdIndices loPads, hiPads, paddedShape, strideOrder;
    axom::quest::internal::stridesAndOffsetsToShapes(realShape,
                                                     offsets,
                                                     strides,
                                                     valuesCount,
                                                     paddedShape,
                                                     loPads,
                                                     hiPads,
                                                     strideOrder);

    T* dataPtr = static_cast<T*>(valuesNode.data_ptr());
    axom::ArrayView<T, DIM, MemSpace> rval(dataPtr, paddedShape, strides);

    if(withGhosts == false)
    {
      auto rval1 = rval;
      rval = rval1.subspan(offsets, realShape);
    }

    return rval;
  }

  /*!
   * @brief Return view to a scalar field variable.
   * 
   * @warning The view returned has an allocator id determined by
   * \a MemSpace, regardless of the memory type.
   * 
   * @warning Assuming, without checking, that the field contains
   * data of type \a T.  User is responsible for using the correct type.
   */
  template <typename T>
  axom::ArrayView<const T, DIM, MemSpace> getConstFieldView(const std::string& fieldName,
                                                            bool withGhosts = false)
  {
    const conduit::Node& fieldNode = m_cdom->fetch_existing("fields/" + fieldName);
    const std::string association = fieldNode.fetch_existing("association").as_string();
    const conduit::Node& valuesNode = fieldNode.fetch_existing("values");
    const axom::IndexType valuesCount = valuesNode.dtype().number_of_elements();

    SLIC_ASSERT_MSG(association == "vertex" || association == "element",
                    "MeshViewUtil only supports vertex and element-based fields right now.");
    const bool onVertex = association == "vertex";

    const auto& realShape = onVertex ? m_nodeShape : m_cellShape;

    MdIndices strides = conduitIndicesToStackArray(fieldNode, "strides");
    if(!fieldNode.has_child("strides"))
    {
      axom::IndexType tmpStride = 1;
      for(int d = 0; d < DIM; ++d)
      {
        strides[d] = tmpStride;
        tmpStride *= m_cellShape[d] + onVertex;
      }
    }

    MdIndices offsets = conduitIndicesToStackArray(fieldNode, "offsets", 0);

    MdIndices loPads, hiPads, paddedShape, strideOrder;
    axom::quest::internal::stridesAndOffsetsToShapes(realShape,
                                                     offsets,
                                                     strides,
                                                     valuesCount,
                                                     paddedShape,
                                                     loPads,
                                                     hiPads,
                                                     strideOrder);

    const T* dataPtr = static_cast<const T*>(valuesNode.data_ptr());
    axom::ArrayView<const T, DIM, MemSpace> rval(dataPtr, paddedShape, strides);

    if(withGhosts == false)
    {
      auto rval1 = rval;
      rval = rval1.subspan(offsets, realShape);
    }

    return rval;
  }
  ///@}

  ///@{
  //!@name Creating data

  /*!
   * @brief Create a new scalar nodal data field with specified
   * strides and offsets.
   *
   * @param [in] fieldName
   * @param [in] association "vertex" or "element"
   * @param [in] dtype Conduit data type to put in the field.  Must be at least
   *             big enough for the strides and offsets specified.
   * @param [in] strides Data strides.  Set to zero for no ghosts and default strides.
   * @param [in] offsets Data index offsets.  Set to zero for no ghosts.
   *
   * Field data allocation is done by Conduit, so the data lives in
   * host memory.  Conduit currently doesn't provide a means to allocate
   * the array in device memory.
   *
   * Creating a field with given strides and offsets may only be useful
   * for matching the strides and offsets of existing data.  It's more
   * natural to create the field based on ghost layer thickness and
   * index advancement order (row-major, column-major or some other).
   * Use createFieldPadded() for this.
   */
  void createField(const std::string& fieldName,
                   const std::string& association,
                   const conduit::DataType& dtype,
                   const MdIndices& strides,
                   const MdIndices& offsets)
  {
    SLIC_ERROR_IF(m_dom == nullptr,
                  axom::fmt::format("Cannot create field {}."
                                    "  MeshViewUtil was not constructed with a non-const domain.",
                                    fieldName));
    SLIC_ERROR_IF(m_dom->has_path("fields/" + fieldName),
                  axom::fmt::format("Cannot create field {}.  It already exists.", fieldName));

    SLIC_ERROR_IF(
      association != "vertex" && association != "element",
      axom::fmt::format("MeshViewUtil doesn't support association '{}' yet.", association));

    const auto& realShape = getRealShape(association);
    MdIndices loPads, hiPads, paddedShape, strideOrder;
    axom::quest::internal::stridesAndOffsetsToShapes(realShape,
                                                     offsets,
                                                     strides,
                                                     dtype.number_of_elements(),
                                                     paddedShape,
                                                     loPads,
                                                     hiPads,
                                                     strideOrder);
    for(int d = 0; d < DIM; ++d)
    {
      SLIC_ERROR_IF(offsets[d] < 0 || offsets[d] > paddedShape[d] - 1,
                    axom::fmt::format("Bad offsets {} for paddedShape {}", offsets, paddedShape));
    }

    conduit::Node& fieldNode = m_dom->fetch("fields/" + fieldName);
    fieldNode["association"] = association;
    fieldNode["topology"] = m_topologyName;

    constexpr bool isInt32 = std::is_same<axom::IndexType, std::int32_t>::value;
    const conduit::DataType conduitDtype =
      isInt32 ? conduit::DataType::int32(DIM) : conduit::DataType::int64(DIM);
    // Make temporary non-const copies for the "set" methods.
    auto tmpStrides = strides, tmpOffsets = offsets;
    fieldNode["strides"].set(conduitDtype, &tmpStrides[0]);
    fieldNode["offsets"].set(conduitDtype, &tmpOffsets[0]);

    axom::IndexType slowDir = strideOrder[DIM - 1];
    auto extras =
      dtype.number_of_elements() - strides[slowDir] * (offsets[slowDir] + realShape[slowDir]);
    if(extras < 0)
    {
      SLIC_ERROR(
        axom::fmt::format("Insufficient space allocated for data, given the "
                          "offsets and strides."));
    }

    fieldNode["values"].set(dtype);
  }

  /*!
   * @brief Create a new scalar nodal data field with specified ghost paddings.
   * 
   * @param [in] fieldName
   * @param [in] association "vertex" or "element"
   * @param [in] dtype Conduit data type to put in the field.
   *   If dtype has too few elements, the minimum sufficient
   *   size will be allocated.
   * @param loPads [i] Ghost padding amount on low side.
   * @param hiPads [i] Ghost padding amount ont high side.
   * @param strideOrder [i] Fastest-to-slowest advancing index directions.
   * 
   * Field data allocation is done by Conduit, so the data lives in
   * host memory.  Conduit currently doesn't provide a means to allocate
   * the array in device memory.
   */
  void createField(const std::string& fieldName,
                   const std::string& association,
                   const conduit::DataType& dtype,
                   const MdIndices& loPads,
                   const MdIndices& hiPads,
                   const MdIndices& strideOrder)
  {
    SLIC_ERROR_IF(m_dom == nullptr,
                  axom::fmt::format("Cannot create field {}."
                                    "  MeshViewUtil was not constructed with a non-const domain.",
                                    fieldName));
    SLIC_ERROR_IF(m_dom->has_path("fields/" + fieldName),
                  axom::fmt::format("Cannot create field {}.  It already exists.", fieldName));

    SLIC_ERROR_IF(association != "vertex" && association != "element",
                  axom::fmt::format("Not yet supporting association '{}'.", association));

    axom::StackArray<axom::IndexType, DIM> offsets;
    axom::StackArray<axom::IndexType, DIM> strides;
    axom::IndexType valuesCount;
    const auto& realShape = getRealShape(association);
    axom::quest::internal::shapesToStridesAndOffsets(realShape,
                                                     loPads,
                                                     hiPads,
                                                     strideOrder,
                                                     1,
                                                     offsets,
                                                     strides,
                                                     valuesCount);

    conduit::Node& fieldNode = m_dom->fetch("fields/" + fieldName);
    fieldNode["association"] = association;
    fieldNode["topology"] = m_topologyName;

    constexpr bool isInt32 = std::is_same<axom::IndexType, std::int32_t>::value;
    const conduit::DataType conduitDtype =
      isInt32 ? conduit::DataType::int32(DIM) : conduit::DataType::int64(DIM);
    // Make temporary non-const copies for the "set" methods.
    auto tmpStrides = strides, tmpOffsets = offsets;
    fieldNode["strides"].set(conduitDtype, &tmpStrides[0]);
    fieldNode["offsets"].set(conduitDtype, &tmpOffsets[0]);

    conduit::DataType bumpedDtype = dtype;
    if(bumpedDtype.number_of_elements() < valuesCount)
    {
      bumpedDtype.set_number_of_elements(valuesCount);
    }
    fieldNode["values"].set(bumpedDtype);
  }
  ///@}

private:
  conduit::Node* m_dom = nullptr;
  conduit::Node* m_topology = nullptr;
  conduit::Node* m_coordset = nullptr;
  const conduit::Node* m_cdom = nullptr;
  const conduit::Node* m_ctopology = nullptr;
  const conduit::Node* m_ccoordset = nullptr;
  std::string m_topologyName;

  MdIndices m_cellShape;
  MdIndices m_nodeShape;

  MdIndices m_coordsStrides;
  MdIndices m_coordsOffsets;
  MdIndices m_coordsPaddedShape;

  MdIndices conduitIndicesToStackArray(
    const conduit::Node& node,
    const std::string& path,
    axom::IndexType defaultVal = axom::numeric_limits<axom::IndexType>::max()) const
  {
    if(node.has_path(path))
    {
      const auto& child = node.fetch_existing(path);
      if(child.dtype().is_int32())
      {
        const auto* ptr = node.fetch_existing(path).as_int32_ptr();
        return internal::makeStackArray<axom::IndexType, DIM>(ptr);
      }
      else if(child.dtype().is_int64())
      {
        const auto* ptr = node.fetch_existing(path).as_int64_ptr();
        return internal::makeStackArray<axom::IndexType, DIM>(ptr);
      }
      else
      {
        SLIC_ERROR("MeshViewUtil internal error: Unanticipated type.");
      }
    }
    return internal::makeStackArray<axom::IndexType, DIM>(defaultVal);
  }

  void computeCoordsDataLayout()
  {
    const conduit::Node& topologyDims = m_ctopology->fetch_existing("elements/dims");

    const std::string coordsetName = m_ctopology->fetch_existing("coordset").as_string();
    const conduit::Node& coordsValues =
      m_cdom->fetch_existing(axom::fmt::format("coordsets/{}/values", coordsetName));

    const axom::IndexType coordValuesCount = coordsValues[0].dtype().number_of_elements();

    const bool coordsAreInterleaved = conduit::blueprint::mcarray::is_interleaved(coordsValues);

    for(int i = 0; i < DIM; ++i)
    {
      m_cellShape[i] = topologyDims[i].to_value();
      m_nodeShape[i] = m_cellShape[i] + 1;
    }

    m_coordsOffsets = conduitIndicesToStackArray(topologyDims, "offsets", 0);
    m_coordsStrides = conduitIndicesToStackArray(topologyDims, "strides", -1);

    if(!topologyDims.has_child("strides"))
    {
      // Compute strides manually, assuming (as Conduit does) that
      // direction 0 is fastest.
      axom::IndexType tmpStride = coordsAreInterleaved ? DIM : 1;
      for(int d = 0; d < DIM; ++d)
      {
        m_coordsStrides[d] = tmpStride;
        tmpStride *= 1 + m_cellShape[d];
      }
    }

    // Compute the padded shape of coords data.
    MdIndices loPads;
    MdIndices hiPads;
    MdIndices strideOrder;
    axom::quest::internal::stridesAndOffsetsToShapes(m_nodeShape,
                                                     m_coordsOffsets,
                                                     m_coordsStrides,
                                                     coordValuesCount,
                                                     m_coordsPaddedShape,
                                                     loPads,
                                                     hiPads,
                                                     strideOrder);
  }

  bool isEqual(const MdIndices& a, const MdIndices& b)
  {
    for(int d = 0; d < DIM; ++d)
    {
      if(a[d] != b[d])
      {
        return false;
      }
    }
    return true;
  }
};

}  // end namespace quest
}  // end namespace axom

#endif  //  AXOM_USE_CONDUIT
#endif  //  QUEST_MESH_VIEW_UTIL_H_
