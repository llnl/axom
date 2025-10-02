// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_UTILITIES_UNIQUE_HPP_
#define AXOM_BUMP_UTILITIES_UNIQUE_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

// Uncomment this to print debug data.
// #define AXOM_DEBUG_UNIQUE

namespace axom
{
namespace bump
{

#if defined(AXOM_DEBUG_UNIQUE)
namespace detail
{
/*!
 * \brief Print a container.
 *
 * \param name The name to print for the container.
 * \param container The container whose data will be printed.
 */
template <typename ContainerType>
void printContainer(const std::string &name, const ContainerType &container)
{
  std::cout << name << "=[";
  for(axom::IndexType i = 0; i < container.size(); i++)
  {
    if(i > 0)
    {
      std::cout << ", ";
    }
    std::cout << container[i];
  }
  std::cout << "]\n";
}

/*!
 * \brief Print a map.
 *
 * \param name The name to print for the container.
 * \param container The container whose data will be printed.
 */
template <typename MapType>
void printMap(const std::string &name, const MapType &container, bool printKey)
{
  std::cout << name << "=[";
  for(auto it = container.begin(); it != container.end(); it++)
  {
    if(it != container.begin())
    {
      std::cout << ", ";
    }
    if(printKey)
      std::cout << it->first;
    else
      std::cout << it->second;
  }
  std::cout << "]\n";
}

}  // end namespace detail
#endif

/*!
 * \brief Makes a unique array of values from values that could contain multiple
 *        instances of a key.
 *
 * \tparam ExecSpace The execution space.
 * \tparam KeyType   The data type for the keys.
 *
 */
template <typename ExecSpace, typename KeyType>
struct Unique
{
  /*!
   * \brief This function makes a unique array of values from an input list of keys. The
   *        output set of unique keys is sorted.
   *
   * \param[in] keys_orig_view The input view that contains the input keys to be made unique.
   * \param[out] skeys     A sorted unique array of keys produced from \a keys_orig_view.
   *                       If there were duplicates in \a keys_orig_view then the size
   *                       will be smaller since duplicates will have been removed.
   * \param[out] sindices  An array of indices that indicate where in the original
   *                       view the keys came from (the index of the key that was used
   *                       the unique key). This array contains the same number of elements
   *                       as \a skeys. This array is useful when the calling algorithm
   *                       needs to select data out of other arrays associated with the
   *                       keys.
   *
   * \note key_orig_view is passed by value so it does not require a local copy to capture it.
   */
  static void execute(const axom::ArrayView<KeyType> keys_orig_view,
                      axom::Array<KeyType> &skeys,
                      axom::Array<axom::IndexType> &sindices)
  {
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Make a copy of the keys and make original indices.
    const auto n = keys_orig_view.size();
    axom::Array<KeyType> keys(n, n, allocatorID);
    axom::Array<axom::IndexType> indices(n, n, allocatorID);
    auto keys_view = keys.view();
    auto indices_view = indices.view();
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) {
        keys_view[i] = keys_orig_view[i];
        indices_view[i] = i;
      });
#if defined(AXOM_DEBUG_UNIQUE)
    // Input values
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      detail::printContainer("keys", keys_view);
      detail::printContainer("indices", indices_view);
    }
#endif
    // Sort the keys, indices in place.
    axom::stable_sort_pairs<ExecSpace>(keys_view, indices_view);

    // Make a mask array for where differences occur.
    axom::Array<axom::IndexType> mask(n, n, allocatorID);
    auto mask_view = mask.view();
    axom::ReduceSum<ExecSpace, axom::IndexType> mask_sum(0);
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) {
        const axom::IndexType m = (i >= 1) ? ((keys_view[i] != keys_view[i - 1]) ? 1 : 0) : 1;
        mask_view[i] = m;
        mask_sum += m;
      });

    // Do a scan on the mask array to build an offset array.
    axom::Array<axom::IndexType> offsets(n, n, allocatorID);
    auto offsets_view = offsets.view();
    axom::exclusive_scan<ExecSpace>(mask_view, offsets_view);

#if defined(AXOM_DEBUG_UNIQUE)
    // Post-sorting values.
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      detail::printContainer("sorted_keys", keys_view);
      detail::printContainer("sorted_indices", indices_view);
      detail::printContainer("mask", mask_view);
      detail::printContainer("offsets", offsets_view);
    }
#endif

    // Allocate the output arrays.
    const axom::IndexType newsize = mask_sum.get();
    skeys = axom::Array<KeyType>(newsize, newsize, allocatorID);
    sindices = axom::Array<axom::IndexType>(newsize, newsize, allocatorID);

    // Iterate over the mask/offsets to store values at the right
    // offset in the new array.
    auto skeys_view = skeys.view();
    auto sindices_view = sindices.view();
    axom::for_all<ExecSpace>(
      n,
      AXOM_LAMBDA(axom::IndexType i) {
        if(mask_view[i])
        {
          skeys_view[offsets_view[i]] = keys_view[i];
          sindices_view[offsets_view[i]] = indices_view[i];
        }
      });

#if defined(AXOM_DEBUG_UNIQUE)
    // Output values
    if(!axom::execution_space<ExecSpace>::onDevice())
    {
      detail::printContainer("skeys", skeys_view);
      detail::printContainer("sindices", sindices_view);
    }
#endif
  }
};

//------------------------------------------------------------------------------
/// Partial specialization for SEQ_EXEC.
template <typename KeyType>
struct Unique<axom::SEQ_EXEC, KeyType>
{
  /*!
   * \brief This function makes a unique array of values from an input list of keys.
   *
   * \param[in] keys_orig_view The input view that contains the input keys to be made unique.
   * \param[out] skeys     A sorted unique array of keys produced from keys_orig_view.
   * \param[out] sindices  An array of indices that indicate where in the original view the keys came from.
   *
   */
  static void execute(const axom::ArrayView<KeyType> &keys_orig_view,
                      axom::Array<KeyType> &skeys,
                      axom::Array<axom::IndexType> &sindices)
  {
    // Make unique values and store the indices.
    std::unordered_map<KeyType, axom::IndexType> unique_map;
    const axom::IndexType n = keys_orig_view.size();
    for(axom::IndexType index = 0; index < n; ++index)
    {
      const auto k = keys_orig_view[index];
      if(unique_map.find(k) == unique_map.end())
      {
        unique_map[k] = index;
      }
    }

    // Store the unordered_map into a vector.
    std::vector<std::pair<KeyType, axom::IndexType>> unique_vector(unique_map.begin(),
                                                                   unique_map.end());

    // Sort the vector by the keys.
    std::sort(unique_vector.begin(),
              unique_vector.end(),
              [](const std::pair<KeyType, axom::IndexType> &a,
                 const std::pair<KeyType, axom::IndexType> &b) { return a.first < b.first; });

    // Allocate the output arrays and populate them
    const axom::IndexType newsize = unique_vector.size();
    const int allocatorID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    skeys = axom::Array<KeyType>(newsize, newsize, allocatorID);
    sindices = axom::Array<axom::IndexType>(newsize, newsize, allocatorID);

    for(axom::IndexType index = 0; index < newsize; ++index)
    {
      skeys[index] = unique_vector[index].first;
      sindices[index] = unique_vector[index].second;
    }
  }
};

}  // end namespace bump
}  // end namespace axom

#endif
