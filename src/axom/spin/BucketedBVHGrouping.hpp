// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file BucketedBVHGrouping.hpp
 * \brief Defines leaf-grouping policies for a bucketed BVH.
 */

namespace axom::spin
{
/*! \brief Selects how spatially sorted primitives are grouped into BVH leaves. */
enum class BucketedBVHGrouping
{
  /// Keep each bucket within one shared spatial-key-prefix range.
  SharedPrefix,

  /// Form consecutive fixed-size buckets in sorted code order.
  Sequential
};
}  // namespace axom::spin
