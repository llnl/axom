// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file SetBuilders.hpp
 *
 * \brief Free-function "make" helpers that construct SLAM sets from a buffer or
 *  a range while deducing the full policy stack.
 *
 * SLAM's sets are configured by a long list of orthogonal policy template
 * parameters. Spelling the full stack at every construction site is verbose:
 *
 * \code
 *   using Set = slam::ArrayViewIndirectionSet<int, double>;
 *   Set s(Set::SetBuilder().size(v.size()).data(v));
 * \endcode
 *
 * The helpers here collapse that to a single call that deduces the element type
 * from the buffer:
 *
 * \code
 *   auto s = slam::make_array_view_set(v);   // -> ArrayViewIndirectionSet<.., double>
 * \endcode
 *
 * \note On CTAD vs. helpers. A class-template-argument deduction guide cannot
 *  recover a set's policy stack from a SetBuilder argument: a guide parameter of
 *  the form `typename OrderedSet<P,...>::SetBuilder` is a non-deduced context
 *  (the template arguments appear only as a nested-name-specifier),
 *  so `OrderedSet s(builder)` can never deduce. Deduction only works from a
 *  directly-named argument type such as `axom::ArrayView<T>`. These free
 *  functions are therefore the portable way to get stack-deducing construction in C++17.
 */

#ifndef SLAM_SET_BUILDERS_H_
#define SLAM_SET_BUILDERS_H_

#include "axom/core/Array.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/IndirectionSet.hpp"

#include <cstddef>
#include <vector>

namespace axom::slam
{
namespace detail
{
template <typename T>
struct type_identity
{
  using type = T;
};

template <typename T>
using type_identity_t = typename type_identity<T>::type;
}  // namespace detail

/// \name Set construction helpers
/// \brief Construct a SLAM set while deducing its policy stack from the buffer
///  or range. \a PosType defaults to slam's default position type and may be
///  supplied explicitly as the leading template argument.
/// \{

/*!
 * \brief Make a contiguous range set \f$[0, size)\f$.
 * \param size the number of elements
 * \return a RangeSet<PosType, ElemType>
 */
template <typename PosType = DefaultPositionType, typename ElemType = DefaultElementType>
RangeSet<PosType, ElemType> make_range_set(detail::type_identity_t<PosType> size)
{
  return RangeSet<PosType, ElemType>(size);
}

/*!
 * \brief Make a contiguous range set \f$[lower, upper)\f$.
 * \param lower the first element of the range
 * \param upper one past the last element of the range
 * \return a RangeSet<PosType, ElemType>
 */
template <typename PosType = DefaultPositionType, typename ElemType = DefaultElementType>
RangeSet<PosType, ElemType> make_range_set(detail::type_identity_t<PosType> lower,
                                           detail::type_identity_t<PosType> upper)
{
  return RangeSet<PosType, ElemType>(lower, upper);
}

/*!
 * \brief Make a set whose elements indirect through an axom::ArrayView.
 *
 * The element type is deduced from \a view; the set is device-capable
 * (ArrayView indirection is host-device). The set's size matches the view.
 *
 * \param view the backing array view
 * \return an ArrayViewIndirectionSet<PosType, T>
 */
template <typename PosType = DefaultPositionType, typename T>
ArrayViewIndirectionSet<PosType, T> make_array_view_set(axom::ArrayView<T> view)
{
  using SetType = ArrayViewIndirectionSet<PosType, T>;
  return SetType(typename SetType::SetBuilder().size(static_cast<PosType>(view.size())).data(view));
}

/*!
 * \brief Make a set whose elements indirect through an std::vector.
 *
 * The element type is deduced from \a vec. The set's size matches the vector.
 * \param vec the backing vector (must outlive the set)
 * \return a VectorIndirectionSet<PosType, T>
 */
template <typename PosType = DefaultPositionType, typename T>
VectorIndirectionSet<PosType, T> make_vector_set(std::vector<T>& vec)
{
  using SetType = VectorIndirectionSet<PosType, T>;
  return SetType(typename SetType::SetBuilder().size(static_cast<PosType>(vec.size())).data(&vec));
}

/*!
 * \brief Make a set whose elements indirect through an axom::Array.
 *
 * The element type is deduced from \a arr; the set is device-capable
 * (axom::Array indirection is host-device). The set's size matches the array.
 *
 * \param arr the backing array (must outlive the set)
 * \return an ArrayIndirectionSet<PosType, T>
 */
template <typename PosType = DefaultPositionType, typename T, int DIM, MemorySpace SPACE, typename StoragePolicy>
ArrayIndirectionSet<PosType, T> make_array_set(axom::Array<T, DIM, SPACE, StoragePolicy>& arr)
{
  using SetType = ArrayIndirectionSet<PosType, T>;
  return SetType(typename SetType::SetBuilder().size(static_cast<PosType>(arr.size())).data(&arr));
}

/*!
 * \brief Make a set whose elements indirect through a C array.
 *
 * \param data pointer to the backing buffer (must outlive the set)
 * \param size the number of elements
 * \return a CArrayIndirectionSet<PosType, T>
 */
template <typename PosType = DefaultPositionType, typename T>
CArrayIndirectionSet<PosType, T> make_c_array_set(T* data, detail::type_identity_t<PosType> size)
{
  using SetType = CArrayIndirectionSet<PosType, T>;
  return SetType(typename SetType::SetBuilder().size(size).data(data));
}

/// \}

}  // end namespace axom::slam

#endif  // SLAM_SET_BUILDERS_H_
