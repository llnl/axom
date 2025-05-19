// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_BLUEPRINT_UTILITIES_HPP_
#define AXOM_MIR_BLUEPRINT_UTILITIES_HPP_

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/mir/views/NodeArrayView.hpp"
#include "axom/export/mir.h"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <utility>
#include <string>
#include <vector>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
//------------------------------------------------------------------------------
/*!
 * \brief This class provides a couple of type traits that let us map C++ types
 *        to types / values useful in Conduit.
 */
template <typename T>
struct cpp2conduit
{ };

template <>
struct cpp2conduit<conduit::int8>
{
  using type = conduit::int8;
  static constexpr conduit::index_t id = conduit::DataType::INT8_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::int16>
{
  using type = conduit::int16;
  static constexpr conduit::index_t id = conduit::DataType::INT16_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::int32>
{
  using type = conduit::int32;
  static constexpr conduit::index_t id = conduit::DataType::INT32_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::int64>
{
  using type = conduit::int64;
  static constexpr conduit::index_t id = conduit::DataType::INT64_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint8>
{
  using type = conduit::uint8;
  static constexpr conduit::index_t id = conduit::DataType::UINT8_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint16>
{
  using type = conduit::uint16;
  static constexpr conduit::index_t id = conduit::DataType::UINT16_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint32>
{
  using type = conduit::uint32;
  static constexpr conduit::index_t id = conduit::DataType::UINT32_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::uint64>
{
  using type = conduit::uint64;
  static constexpr conduit::index_t id = conduit::DataType::UINT64_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::float32>
{
  using type = conduit::float32;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT32_ID;
  AXOM_MIR_EXPORT static const char *name;
};

template <>
struct cpp2conduit<conduit::float64>
{
  using type = conduit::float64;
  static constexpr conduit::index_t id = conduit::DataType::FLOAT64_ID;
  AXOM_MIR_EXPORT static const char *name;
};

//------------------------------------------------------------------------------

/*!
 * \brief Make an axom::ArrayView from a Conduit node.
 *
 * \tparam T The type for the array view elements.
 *
 * \param n The conduit node for which we want an array view.
 *
 * \return An axom::ArrayView that wraps the data in the Conduit node.
 */
/// @{
template <typename T>
inline axom::ArrayView<T> make_array_view(conduit::Node &n)
{
  SLIC_ASSERT_MSG(cpp2conduit<T>::id == n.dtype().id(),
                  axom::fmt::format("Cannot create ArrayView<{}> for Conduit {} data.",
                                    cpp2conduit<T>::name,
                                    n.dtype().name()));
  return axom::ArrayView<T>(static_cast<T *>(n.data_ptr()), n.dtype().number_of_elements());
}

template <typename T>
inline axom::ArrayView<T> make_array_view(const conduit::Node &n)
{
  SLIC_ASSERT_MSG(cpp2conduit<T>::id == n.dtype().id(),
                  axom::fmt::format("Cannot create ArrayView<{}> for Conduit {} data.",
                                    cpp2conduit<T>::name,
                                    n.dtype().name()));
  return axom::ArrayView<T>(static_cast<T *>(const_cast<void *>(n.data_ptr())),
                            n.dtype().number_of_elements());
}
/// @}

//------------------------------------------------------------------------------
/*!
 * \brief This class registers a Conduit allocator that can make Conduit allocate
 *        through Axom's allocate/deallocate functions using a specific allocator.
 *        This permits Conduit to allocate through Axom's UMPIRE logic.
 *
 * \tparam ExecSpace The execution space.
 */
template <typename ExecSpace>
class ConduitAllocateThroughAxom
{
public:
  /*!
   * \brief Get the Conduit allocator ID for this ExecSpace.
   *
   * \return The Conduit allocator ID for this ExecSpace.
   */
  static conduit::index_t getConduitAllocatorID()
  {
    constexpr conduit::index_t NoAllocator = -1;
    static conduit::index_t conduitAllocatorID = NoAllocator;
    if(conduitAllocatorID == NoAllocator)
    {
      conduitAllocatorID = conduit::utils::register_allocator(internal_allocate, internal_free);
    }
    return conduitAllocatorID;
  }

private:
  /*!
   * \brief A function we register with Conduit to allocate memory.
   *
   * \param items The number of items to allocate.
   * \param item_size The size of each item in bytes.
   *
   * \brief A block of newly allocated memory large enough for the requested items.
   */
  static void *internal_allocate(size_t items, size_t item_size)
  {
    const auto axomAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
    void *ptr = static_cast<void *>(axom::allocate<std::uint8_t>(items * item_size, axomAllocatorID));
    //std::cout << axom::execution_space<ExecSpace>::name()
    //  << ": Allocated for Conduit via axom: items=" << items
    //  << ", item_size=" << item_size << ", ptr=" << ptr << std::endl;
    return ptr;
  }

  /*!
   * \brief A deallocation function we register with Conduit.
   */
  static void internal_free(void *ptr)
  {
    //std::cout << axom::execution_space<ExecSpace>::name()
    //  << ": Dellocating for Conduit via axom: ptr=" << ptr << std::endl;
    axom::deallocate(ptr);
  }
};

//------------------------------------------------------------------------------
/*!
 * \brief Copies a Conduit tree in the \a src node to a new Conduit \a dest node,
 *        making sure to allocate array data in the appropriate memory space for
 *        the execution space.
 *
 * \tparam ExecSpace The destination execution space (e.g. axom::SEQ_EXEC).
 *
 * \param dest The conduit node that will receive the copied data.
 * \param src The source data to be copied.
 */
template <typename ExecSpace>
void copy(conduit::Node &dest, const conduit::Node &src)
{
  ConduitAllocateThroughAxom<ExecSpace> c2a;
  dest.reset();
  if(src.number_of_children() > 0)
  {
    for(conduit::index_t i = 0; i < src.number_of_children(); i++)
    {
      copy<ExecSpace>(dest[src[i].name()], src[i]);
    }
  }
  else
  {
    const int allocatorID = axom::getAllocatorIDFromPointer(src.data_ptr());
    bool deviceAllocated = isDeviceAllocator(allocatorID);
    if(deviceAllocated || (!src.dtype().is_string() && src.dtype().number_of_elements() > 1))
    {
      // Allocate the node's memory in the right place.
      dest.reset();
      dest.set_allocator(c2a.getConduitAllocatorID());
      dest.set(conduit::DataType(src.dtype().id(), src.dtype().number_of_elements()));

      // Copy the data to the destination node. Axom uses Umpire to manage that.
      if(src.is_compact())
        axom::copy(dest.data_ptr(), src.data_ptr(), src.dtype().bytes_compact());
      else
      {
        // NOTE: This assumes that src is on the host.
        conduit::Node tmp;
        src.compact_to(tmp);
        axom::copy(dest.data_ptr(), tmp.data_ptr(), tmp.dtype().bytes_compact());
      }
    }
    else
    {
      // The data fits in the node or is a string. It's on the host.
      dest.set(src);
    }
  }
}

//------------------------------------------------------------------------------
/*!
 * \brief Fill an array with int values from a Conduit node.
 *
 * \tparam ArrayType The array type being filled.
 *
 * \param n The node that contains the data.
 * \param key The name of the node that contains the data in \a n.
 * \param[out] arr The array being filled.
 * \param moveToHost Sometimes data are on device and need to be moved to host first.
 */
template <typename ArrayType>
bool fillFromNode(const conduit::Node &n, const std::string &key, ArrayType &arr, bool moveToHost = false)
{
  bool found = false;
  if((found = n.has_path(key)) == true)
  {
    if(moveToHost)
    {
      // Make sure data are on host.
      conduit::Node hostNode;
      copy<axom::SEQ_EXEC>(hostNode, n.fetch_existing(key));

      const auto acc = hostNode.as_int_accessor();
      for(int i = 0; i < arr.size(); i++)
      {
        arr[i] = acc[i];
      }
    }
    else
    {
      const auto acc = n.fetch_existing(key).as_int_accessor();
      for(int i = 0; i < arr.size(); i++)
      {
        arr[i] = acc[i];
      }
    }
  }
  return found;
}
//------------------------------------------------------------------------------

/*!
 * \brief Return the names of the axes for a coordset.
 *
 * \param n_input A Conduit node containing a coordset.
 *
 * \return A vector containing the names of the coordset's axes.
 */
std::vector<std::string> coordsetAxes(const conduit::Node &n_input);

//------------------------------------------------------------------------------
/*!
 * \brief Returns the input index (no changes).
 */
struct DirectIndexing
{
  /*!
   * \brief Return the input index (no changes).
   * \param index The input index.
   * \return The input index.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const { return index; }
};

//------------------------------------------------------------------------------
/*!
 * \brief Help turn slice data zone indices into strided structured element field indices.
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSElementFieldIndexing
{
  /*!
   * \brief Update the indexing offsets/strides from a Conduit node.
   * \param field The Conduit node for a field.
   *
   * \note Executes on the host.
   */
  void update(const conduit::Node &field)
  {
    fillFromNode(field, "offsets", m_indexing.m_offsets, true);
    fillFromNode(field, "strides", m_indexing.m_strides, true);
  }

  /*!
   * \brief Transforms the index from local to global through an indexing object.
   * \param index The local index
   * \return The global index for the field.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const
  {
    return m_indexing.LocalToGlobal(index);
  }

  Indexing m_indexing {};
};

//------------------------------------------------------------------------------
/*!
 * \brief Help turn blend group node indices (global) into vertex field indices.
 * \tparam Indexing A StridedStructuredIndexing of some dimension.
 */
template <typename Indexing>
struct SSVertexFieldIndexing
{
  /*!
   * \brief Update the indexing offsets/strides from a Conduit node.
   * \param field The Conduit node for a field.
   *
   * \note Executes on the host.
   */
  void update(const conduit::Node &field)
  {
    fillFromNode(field, "offsets", m_fieldIndexing.m_offsets, true);
    fillFromNode(field, "strides", m_fieldIndexing.m_strides, true);
  }

  /*!
   * \brief Transforms the index from local to global through an indexing object.
   * \param index The global index
   * \return The global index for the field.
   */
  AXOM_HOST_DEVICE
  inline axom::IndexType operator[](axom::IndexType index) const
  {
    // Make the global index into a global logical in the topo.
    const auto topoGlobalLogical = m_topoIndexing.GlobalToGlobal(index);
    // Make the global logical into a local logical in the topo.
    const auto topoLocalLogical = m_topoIndexing.GlobalToLocal(topoGlobalLogical);
    // Make the global logical index in the field.
    const auto fieldGlobalLogical = m_fieldIndexing.LocalToGlobal(topoLocalLogical);
    // Make the global index in the field.
    const auto fieldGlobalIndex = m_fieldIndexing.GlobalToGlobal(fieldGlobalLogical);
    return fieldGlobalIndex;
  }

  Indexing m_topoIndexing {};
  Indexing m_fieldIndexing {};
};

//------------------------------------------------------------------------------
/*!
 * \brief Get the min/max values for the data in a Conduit node or ArrayView.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam ReturnType The data type of the returned min/max values.
 */
template <typename ExecSpace, typename ReturnType>
struct MinMax
{
  /*!
   * \brief Get the min/max values for the data in a Conduit node.
   *
   * \param[in] n The Conduit node whose data we're checking.
   *
   * \return A pair containing the min,max values in the node.
   */
  static std::pair<ReturnType, ReturnType> execute(const conduit::Node &n)
  {
    SLIC_ASSERT(n.dtype().number_of_elements() > 0);
    std::pair<ReturnType, ReturnType> retval;

    axom::mir::views::Node_to_ArrayView(n, [&](auto nview) { retval = execute(nview); });
    return retval;
  }

  /*!
   * \brief Get the min/max values for the data in an ArrayView.
   *
   * \param[in] n The Conduit node whose data we're checking.
   *
   * \return A pair containing the min,max values in the node.
   */
  template <typename T>
  static std::pair<ReturnType, ReturnType> execute(const axom::ArrayView<T> nview)
  {
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

    RAJA::ReduceMin<reduce_policy, T> vmin(axom::numeric_limits<T>::max());
    RAJA::ReduceMax<reduce_policy, T> vmax(axom::numeric_limits<T>::min());

    axom::for_all<ExecSpace>(
      nview.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        vmin.min(nview[index]);
        vmax.max(nview[index]);
      });

    return std::pair<ReturnType, ReturnType> {static_cast<ReturnType>(vmin.get()),
                                              static_cast<ReturnType>(vmax.get())};
  }
};

//------------------------------------------------------------------------------
/*!
 * \brief Base template for computing a shape's area or volume.
 */
template <int NDIMS>
struct ComputeShapeAmount
{ };

/*!
 * \brief 2D specialization for shapes to compute area.
 */
template <>
struct ComputeShapeAmount<2>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.area();
  }
};

/*!
 * \brief 3D specialization for shapes to compute volume.
 */
template <>
struct ComputeShapeAmount<3>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.volume();
  }
};

//------------------------------------------------------------------------------
/*!
 * \brief Save a Blueprint mesh to a legacy ASCII VTK file.
 *
 * \param node The node that contains the mesh data.
 * \param path The file path to save.
 *
 * \note This function currently handles only unstructured topos with explicit coordsets.
 */
void save_vtk(const conduit::Node &node, const std::string &path);

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
