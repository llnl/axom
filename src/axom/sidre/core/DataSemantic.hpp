// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file DataSemantic.hpp
 *
 * \brief   Header file containing definition of Semantic type.
 *
 ******************************************************************************
 */

#ifndef SIDRE_DATA_SEMANTIC_HPP_
#define SIDRE_DATA_SEMANTIC_HPP_

namespace axom
{
namespace sidre
{
  /*!
    The term semantic is borrowed from "copy semantics."  Here, it affects
    how the data stored for usage as well as how it may be copied.

    Note: The term "semantic" may not be the best.  What are alternatives?

    VALUEs are always copied by value.  REFERENCEs may be copied by reference
    or explicitly deep-copied like value-copying.

    In usage, it's helpful (but not required) to have separate
    allocators for VALUEs and REFERENCEs.  There are different default
    allocators for the hierarchy.  @see Group::setDefaultArrayAllocator()
    and @see Group::setDefaultTupleAllocator().

    The allocators can be overridden individually for each node in the
    hierarchy.
  */
  enum DataSemantic
  {
    UNKNOWN,   // Semantic is unknown for empty Views
    VALUE,     // For scalars, tuples and strings
    REFERENCE, // For BUFFER and EXTERNAL states.
  };

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_DATA_SEMANTIC_HPP_ */
