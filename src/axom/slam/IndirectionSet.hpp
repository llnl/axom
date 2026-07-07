// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file IndirectionSet.hpp
 *
 * \brief Defines alias templates for OrderedSets with indirection.
 *
 * Use \c ArrayIndirectionSet for a set indexed through an \c axom::Array it manages,
 * and \c ArrayViewIndirectionSet for a set indexed through an \c axom::ArrayView of a buffer managed elsewhere. 
 * \c CArrayIndirectionSet and \c VectorIndirectionSet index raw-pointer and \c std::vector storage,
 * for interoperation and as reference examples for custom indirection policies.
 */

#include <cstddef>
#include <vector>

#include "axom/slam/OrderedSet.hpp"

namespace axom
{
namespace slam
{
/**
 * \brief Alias template for an OrderedSet with indirection over a C array.
 *
 * \note Indexes a raw pointer, for interoperation with C-style array storage.
 *  For an \c axom::ArrayView of a buffer managed elsewhere, use \c ArrayViewIndirectionSet.
 *
 * \tparam PosType The position type for indexing into the set
 * \tparam ElemType The type for the set's elements
 * \sa OrderedSet
 */
template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
using CArrayIndirectionSet = OrderedSet<PosType,
                                        ElemType,
                                        policies::RuntimeSize<PosType>,
                                        policies::ZeroOffset<PosType>,
                                        policies::StrideOne<PosType>,
                                        policies::CArrayIndirection<PosType, ElemType>>;

/**
 * \brief Alias template for an OrderedSet with indirection over a std::vector.
 *
 * \note Indexes a (host-only) \c std::vector, for interoperation with existing \c std::vector storage. 
 *  For \c axom::Array / \c axom::ArrayView storage, use \c ArrayIndirectionSet or \c ArrayViewIndirectionSet.
 *
 * \tparam PosType The position type for indexing into the set
 * \tparam ElemType The type for the set's elements
 * \sa OrderedSet
 */
template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
using VectorIndirectionSet = OrderedSet<PosType,
                                        ElemType,
                                        policies::RuntimeSize<PosType>,
                                        policies::ZeroOffset<PosType>,
                                        policies::StrideOne<PosType>,
                                        policies::STLVectorIndirection<PosType, ElemType>>;

/**
 * \brief Alias template for an OrderedSet with indirection over an axom::Array.
 *
 * This is the canonical Axom alias for binding a set to an \c axom::Array buffer.
 * The array object is managed outside the set and must outlive it.
 *
 * \tparam PosType The position type for indexing into the set
 * \tparam ElemType The type for the set's elements
 * \sa OrderedSet
 */
template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
using ArrayIndirectionSet = OrderedSet<PosType,
                                       ElemType,
                                       policies::RuntimeSize<PosType>,
                                       policies::ZeroOffset<PosType>,
                                       policies::StrideOne<PosType>,
                                       policies::ArrayIndirection<PosType, ElemType>>;

/**
 * \brief Alias template for an OrderedSet with indirection over an axom::ArrayView.
 *
 * A set indexed through an \c axom::ArrayView of a buffer managed elsewhere.
 * The backing allocation must outlive the set. Because \c axom::ArrayView is
 * trivially copyable, it can be captured by value into device kernels.
 *
 * \tparam PosType The position type for indexing into the set
 * \tparam ElemType The type for the set's elements
 * \sa OrderedSet
 */
template <typename PosType = slam::DefaultPositionType, typename ElemType = slam::DefaultElementType>
using ArrayViewIndirectionSet = OrderedSet<PosType,
                                           ElemType,
                                           policies::RuntimeSize<PosType>,
                                           policies::ZeroOffset<PosType>,
                                           policies::StrideOne<PosType>,
                                           policies::ArrayViewIndirection<PosType, ElemType>>;

}  // end namespace slam
}  // end namespace axom
