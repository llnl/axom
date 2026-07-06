// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file MIRMeshTypes.hpp
 * 
 * \brief Contains the specifications for types aliases used throughout the MIR component.
 */

#include "axom/core.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"

namespace axom
{
namespace mir
{
enum Shape
{
  Triangle,
  Quad,
  Tetrahedron,
  Triangular_Prism,
  Pyramid,
  Hexahedron
};

using Point2 = primal::Point<double, 3>;

// SET TYPE ALIASES
using PosType = slam::DefaultPositionType;
using ElemType = slam::DefaultElementType;

using ArrayIndir = slam::policies::CArrayIndirection<PosType, ElemType>;

using VertSet = slam::PositionSet<PosType, ElemType>;
using ElemSet = slam::PositionSet<PosType, ElemType>;

// RELATION TYPE ALIASES
using VarCard = slam::policies::VariableCardinality<PosType, ArrayIndir>;

// Note: This is the actual relation type, which takes in a bunch of policies and data entries to relate to each other.
// Note: It is the relation of the elements to the vertices.

using ElemToVertRelation =
  slam::StaticRelation<PosType, ElemType, VarCard, ArrayIndir, ElemSet, VertSet>;
using VertToElemRelation =
  slam::StaticRelation<PosType, ElemType, VarCard, ArrayIndir, VertSet, ElemSet>;

// MAP TYPE ALIASES
using BaseSet = slam::Set<PosType, ElemType>;
using ScalarMap = slam::ArrayMap<BaseSet, axom::float64>;
using PointMap = slam::ArrayMap<BaseSet, Point2>;
using IntMap = slam::ArrayMap<BaseSet, int>;
}  // namespace mir
}  // namespace axom
