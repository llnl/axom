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
 * New Axom code should prefer \c ArrayIndirectionSet when binding an
 * \c axom::Array buffer and \c ArrayViewIndirectionSet when viewing storage
 * owned elsewhere. \c CArrayIndirectionSet and \c VectorIndirectionSet are
 * retained for raw-pointer/std::vector interop and as reference examples for
 * custom indirection policies.
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
 * \note This is a low-level compatibility/reference alias. Prefer
 *  \c ArrayViewIndirectionSet for new non-owning buffers when possible.
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
 * \note This is a host-only compatibility/reference alias. Prefer
 *  \c ArrayIndirectionSet or \c ArrayViewIndirectionSet for new Axom code.
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
 * This is the canonical Axom alias for binding a set to an \c axom::Array
 * buffer. The array object is managed outside the set and must outlive it.
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
 * This is the canonical Axom alias for a non-owning, view-backed indirection
 * set. The view's backing allocation must outlive the set. Because
 * \c axom::ArrayView is trivially copyable, this alias is the preferred shape
 * for sets captured by value into device kernels.
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
