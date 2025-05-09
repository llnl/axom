// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MIRMeshTypes.hpp
 * 
 * \brief Contains the specifications for types aliases used throughout the MIR component.
 */

#ifndef __MIR_MESH_TYPES_H__
#define __MIR_MESH_TYPES_H__

#include "axom/core.hpp"  // for axom macros
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
using ScalarMap = slam::Map<axom::float64, BaseSet>;
using PointMap = slam::Map<Point2, BaseSet>;
using IntMap = slam::Map<int, BaseSet>;
}  // namespace mir
}  // namespace axom
#endif
