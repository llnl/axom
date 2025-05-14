// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_UTILITIES_SORTING_HPP
#define AXOM_CORE_UTILITIES_SORTING_HPP
#include <axom/core/utilities/Utilities.hpp>
#include <axom/core/NumericLimits.hpp>

namespace axom
{
namespace utilities
{
namespace detail
{

/*!
 * \brief Computes stack size for qsort.
 * \param N The largest number of items to support.
 * \return A number of elements for an array-based stack.
 */
template <typename T>
AXOM_HOST_DEVICE constexpr static T stack_size(T N)
{
  T v = N;
  T i = 0;
  while(v > 0)
  {
    i++;
    v /= 2;
  }
  return i * 2;
}

/*!
 * \brief Swap 2 elements if b < a.
 * param a The first value.
 * param b The second value.
 */
template <typename T>
AXOM_HOST_DEVICE inline void ifswap(T &a, T &b)
{
  if(a > b)
  {
    T tmp = a;
    a = b;
    b = tmp;
  }
}

/*!
 * \brief This method is called when there are no more values to swap.
 */
AXOM_HOST_DEVICE
inline void internal_swap(axom::IndexType, axom::IndexType) { }

/*!
 * \brief Swap one set of values, move to next ones.
 *
 * \param idx1 The first array index.
 * \param idx2 The second array index.
 * \param values The current array values being swapped.
 * \param args A parameter pack of additional arrays that need to be swapped.
 */
template <typename T, typename... Args>
AXOM_HOST_DEVICE inline void internal_swap(axom::IndexType idx1,
                                           axom::IndexType idx2,
                                           T values,
                                           Args... args)
{
  axom::utilities::swap(values[idx1], values[idx2]);
  internal_swap(idx1, idx2, args...);
}

template <typename T, typename Predicate, typename... Args>
AXOM_HOST_DEVICE inline void ifswap_multiple(Predicate &&predicate,
                                             axom::IndexType idx1,
                                             axom::IndexType idx2,
                                             T values,
                                             Args... args)
{
  if(predicate(values[idx2], values[idx1]))
  {
    internal_swap(idx1, idx2, values, args...);
  }
}

/*!
 * \brief Create a partition for qsort.
 *
 * \param[inout] values The array to be sorted.
 * \param low The lower bound of the input partition.
 * \param high The upper bound of the input partition.
 *
 * \return A new pivot.
 */
template <typename T, typename Predicate, typename... Args>
AXOM_HOST_DEVICE static axom::IndexType partitionMultiple(Predicate &&predicate,
                                                          axom::IndexType low,
                                                          axom::IndexType high,
                                                          T values,
                                                          Args... args)
{
  const auto pivot = values[high];
  axom::IndexType i = low - 1;
  for(axom::IndexType j = low; j < high; j++)
  {
    if(predicate(values[j], pivot) || values[j] == pivot)
    {
      i++;
      internal_swap(i, j, values, args...);
    }
  }
  internal_swap(i + 1, high, values, args...);
  return i + 1;
}

/*!
 * \brief Sort the input array using qsort.
 *
 * \param[inout] values The array to be sorted.
 * \param n The number of values in the array.
 */
template <typename T, typename Predicate, typename... Args>
AXOM_HOST_DEVICE static void qsortMultiple(Predicate &&predicate,
                                           axom::IndexType n,
                                           T values,
                                           Args... args)
{
  if(n <= 1) return;
  axom::IndexType stack[stack_size(axom::numeric_limits<axom::IndexType>::max())][2];
  axom::IndexType stack_count = 1;
  stack[0][0] = 0;
  stack[0][1] = n - 1;

  while(stack_count > 0)
  {
    stack_count--;
    const auto low = stack[stack_count][0];
    const auto high = stack[stack_count][1];

    if(low < high)
    {
      const auto pivot_index = partitionMultiple(predicate, low, high, values, args...);

      stack[stack_count][0] = low;
      stack[stack_count][1] = pivot_index - 1;
      stack_count++;

      stack[stack_count][0] = pivot_index + 1;
      stack[stack_count][1] = high;
      stack_count++;
    }
  }
}

/*!
 * \brief Perform actual sort on multiple data arrays.
 *
 * \param n The length of the data arrays.
 * \param values The array being sorted.
 * \param args The rest of the arrays to be sorted the same way.
 */
template <typename T, typename Predicate, typename... Args>
AXOM_HOST_DEVICE static void insertionSortMultiple(Predicate &&predicate,
                                                   axom::IndexType n,
                                                   T values,
                                                   Args... args)
{
  for(axom::IndexType i = 1; i < n; i++)
  {
    axom::IndexType j = i;
    // Keep swapping elements until we're not out-of-order.
    while(j > 0 && predicate(values[j], values[j - 1]))
    {
      internal_swap(j, j - 1, values, args...);
      j--;
    }
  }
}

// This struct is used in shifting template arguments.
struct Delimiter
{ };

/*!
 * \brief Terminal case for argument reordering where we know all of the
 *        array arguments as well as the array length argument.
 *
 * \param predicate The predicate to use for comparing elements,
 * \param first The first data array used for sorting.
 * \param d Delimiter
 * \param n The length of the arrays.
 * \param args A parameter pack containing all other arrays to be sorted.
 */
template <typename T, typename Predicate, typename... Args>
AXOM_HOST_DEVICE inline static void sort_multiple_internal(Predicate &&predicate,
                                                           axom::IndexType n,
                                                           T first,
                                                           Delimiter,
                                                           Args... args)
{
  constexpr int SortingCutoffSize = 11;
  if(n < SortingCutoffSize)
  {
    insertionSortMultiple(predicate, n, first, args...);
  }
  else
  {
    qsortMultiple(predicate, n, first, args...);
  }
}

/*!
 * \brief Shift args to end until the array length is at the desired position.
 */
template <typename T, typename Predicate, typename... Args>
AXOM_HOST_DEVICE inline static void sort_multiple_internal(Predicate &&predicate,
                                                           T first,
                                                           Delimiter d,
                                                           axom::IndexType n,
                                                           Args... args)
{
  sort_multiple_internal(predicate, n, first, d, args...);
}

/*!
 * \brief Shift args to end until the array length is at the desired position.
 */
template <typename T, typename Predicate, typename S, typename... Args>
AXOM_HOST_DEVICE inline static void sort_multiple_internal(Predicate &&predicate,
                                                           T first,
                                                           Delimiter d,
                                                           S second,
                                                           Args... args)
{
  sort_multiple_internal(predicate, first, d, args..., second);
}

/*!
 * \brief Predicate that implements less than.
 */
template <typename T>
struct less_than
{
  AXOM_HOST_DEVICE inline bool operator()(const T &a, const T &b) const { return a < b; }
};

/*!
 * \brief Predicate that implements greater than.
 */
template <typename T>
struct greater_than
{
  AXOM_HOST_DEVICE inline bool operator()(const T &a, const T &b) const { return a > b; }
};

}  // end namespace detail

/*!
 * \brief Sort multiple arrays, using the first array as the key for sorting.
 *        All other arrays are sorted the same.
 *
 * \tparam T The type of the first array.
 * \tparam Args Parameter pack containing multiple array types and an int.
 *
 * \param first The first data array.
 * \param args A parameter pack of all other args.
 *
 * \note This function is not part of the Sorting struct because CUDA did not
 *       like variadic templates on member functions.
 */
template <typename T, typename... Args>
AXOM_HOST_DEVICE inline static void sort_multiple(T first, Args... args)
{
  using ElementType = typename std::remove_pointer<T>::type;
  detail::sort_multiple_internal(detail::less_than<ElementType> {},
                                 first,
                                 detail::Delimiter {},
                                 args...);
}

/*!
 * \brief Sort multiple arrays in reverse order, using the first array as the
 *        key for sorting. All other arrays are sorted the same.
 *
 * \tparam T The type of the first array.
 * \tparam Args Parameter pack containing multiple array types and an int.
 *
 * \param first The first data array.
 * \param args A parameter pack of all other args.
 *
 * \note This function is not part of the Sorting struct because CUDA did not
 *       like variadic templates on member functions.
 */
template <typename T, typename... Args>
AXOM_HOST_DEVICE inline static void reverse_sort_multiple(T first, Args... args)
{
  using ElementType = typename std::remove_pointer<T>::type;
  detail::sort_multiple_internal(detail::greater_than<ElementType> {},
                                 first,
                                 detail::Delimiter {},
                                 args...);
}

/*!
 * \brief This is a template suitable for sorting small arrays on device.
 *
 * \tparam T The data type being sorted.
 * \tparam N The biggest array size that will be used.
 *
 * \note For very short arrays, a simpler sort is faster. As the array size
 *       increases, the algorithm switches to qsort. Also, this is designed as
 *       a template class so it can be specialized.
 */
template <typename T, int N = axom::numeric_limits<int>::max()>
struct Sorting
{
  static constexpr int MAX_SIZE = N;
  static constexpr int SORT_SIZE_CUTOFF = 11;

  /*!
   * \brief Sort an array of values in place.
   *
   * \param[inout] values The array of values to sort.
   * \param n The number of values in the array. 
   */
  AXOM_HOST_DEVICE
  inline static void sort(T *values, int n)
  {
    if(n < SORT_SIZE_CUTOFF)
    {
      insertionSort(values, n);
    }
    else
    {
      qsort(values, n);
    }
  }

  //----------------------------------------------------------------------------
  // Internal methods (sort)
  //----------------------------------------------------------------------------

  /*!
   * \brief Sort the input array using qsort.
   *
   * \param[inout] values The array to be sorted.
   * \param n The number of values in the array.
   */
  AXOM_HOST_DEVICE
  static void qsort(T *values, int n)
  {
    if(n <= 1)
    {
      return;
    }
    int stack[detail::stack_size(N)][2];
    int stack_count = 1;
    stack[0][0] = 0;
    stack[0][1] = n - 1;

    while(stack_count > 0)
    {
      stack_count--;
      int low = stack[stack_count][0];
      int high = stack[stack_count][1];

      if(low < high)
      {
        const int pivot_index = partition(values, low, high);

        stack[stack_count][0] = low;
        stack[stack_count][1] = pivot_index - 1;
        stack_count++;

        stack[stack_count][0] = pivot_index + 1;
        stack[stack_count][1] = high;
        stack_count++;
      }
    }
  }

  /*!
   * \brief Create a partition for qsort.
   *
   * \param[inout] values The array to be sorted.
   * \param low The lower bound of the input partition.
   * \param high The upper bound of the input partition.
   *
   * \return A new pivot.
   */
  AXOM_HOST_DEVICE
  static int partition(T *values, int low, int high)
  {
    // Median-of-three pivot selection
    const int mid = low + (high - low) / 2;

    // Find the median of values[low], values[mid], values[high]
    detail::ifswap(values[low], values[mid]);
    detail::ifswap(values[low], values[high]);
    detail::ifswap(values[mid], values[high]);

    // Use the median as the pivot
    axom::utilities::swap(values[mid], values[high]);

    const T pivot = values[high];

    int i = low - 1;
    for(int j = low; j < high; j++)
    {
      if(values[j] <= pivot)
      {
        i++;
        axom::utilities::swap(values[i], values[j]);
      }
    }
    axom::utilities::swap(values[i + 1], values[high]);
    return i + 1;
  }

  /*!
   * \brief Sort the input array using insertion sort.
   *
   * \param[inout] values The array to be sorted.
   * \param n The number of values in the array.
   */
  AXOM_HOST_DEVICE
  inline static void insertionSort(T *values, int n)
  {
    for(int i = 1; i < n; i++)
    {
      int j = i;
      // Keep swapping elements until we're not out-of-order.
      while(j > 0 && (values[j] < values[j - 1]))
      {
        axom::utilities::swap(values[j], values[j - 1]);
        j--;
      }
    }
  }
};

/*!
 * \brief Template specialization for sorting arrays with 3 elements.
 * \note Specializing resulted in a small speedup over general sorting.
 */
template <typename T>
struct Sorting<T, 3>
{
  /*!
   * \brief Sort the input array.
   *
   * \param[inout] values The array to be sorted.
   */
  AXOM_HOST_DEVICE
  inline static void sort(T *values, int AXOM_UNUSED_PARAM(n))
  {
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[1], values[2]);
    detail::ifswap(values[0], values[1]);
  }
};

/*!
 * \brief Template specialization for sorting arrays with 4 elements.
 * \note Specializing resulted in a small speedup over general sorting.
 */
template <typename T>
struct Sorting<T, 4>
{
  /*!
   * \brief Sort the input array.
   *
   * \param[inout] values The array to be sorted.
   */
  AXOM_HOST_DEVICE
  inline static void sort(T *values, int AXOM_UNUSED_PARAM(n))
  {
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[2], values[3]);
    detail::ifswap(values[1], values[2]);
    detail::ifswap(values[0], values[1]);
    detail::ifswap(values[2], values[3]);
    detail::ifswap(values[1], values[2]);
  }
};

}  // end namespace utilities
}  // end namespace axom

#endif
