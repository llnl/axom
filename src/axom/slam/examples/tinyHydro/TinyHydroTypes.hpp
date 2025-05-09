// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Part class holds the material data for a single material.

#ifndef __TINY_HYDRO_TYPES_H__
#define __TINY_HYDRO_TYPES_H__


#include "VectorXY.hpp"

#include "axom/fmt.hpp"

#include "axom/slam.hpp"

namespace slam = axom::slam;

namespace tinyHydro {

  using PositionType = int;
  using IndexType = PositionType;
  using ElementType = PositionType;

  using SetBase = slam::Set<PositionType, ElementType>;

  using IndexField = slam::Map<int>;

  using ScalarField = slam::Map<double,SetBase>;
  using NodalScalarField = ScalarField;
  using ZonalScalarField = ScalarField;

  using VectorField = slam::Map<VectorXY,SetBase>;
  using NodalVectorField = VectorField;
  using ZonalVectorField = VectorField;
  using FaceVectorField = VectorField;

  using BoundaryEdgeVectorField = VectorField;

  using ZoneSet = slam::PositionSet<PositionType, ElementType>;
  using NodeSet = slam::PositionSet<PositionType, ElementType>;
  using FaceSet = slam::PositionSet<PositionType, ElementType>;
  using CornerSet = slam::PositionSet<PositionType, ElementType>;

  using ZoneSubset = slam::VectorIndirectionSet<PositionType, ElementType>;
  using NodeSubset = slam::VectorIndirectionSet<PositionType, ElementType>;

  enum
  {
    NODES_PER_ZONE        = 4,
    FACES_PER_ZONE        = 4,
    BD_BOTTOM             = 0,         // lower boundary nodes
    BD_RIGHT              = 1,         // right boundary nodes
    BD_TOP                = 2,         // top boundary nodes
    BD_LEFT               = 3,         // left boundary nodes
    NUM_DOMAIN_BOUNDARIES = 4
  };

  using STLIndirection = slam::policies::STLVectorIndirection<PositionType, PositionType>;

  using ZNStride = slam::policies::CompileTimeStride<PositionType, NODES_PER_ZONE>;
  using ZNCard = slam::policies::ConstantCardinality<PositionType, ZNStride>;
  using ZoneToNodeRelation = slam::StaticRelation<PositionType, ElementType,
                                                  ZNCard, STLIndirection, ZoneSet, NodeSet>;
  using ZNodeSet = ZoneToNodeRelation::RelationSubset;

  using ZFStride = slam::policies::CompileTimeStride<PositionType, FACES_PER_ZONE>;
  using ZFCard = slam::policies::ConstantCardinality<PositionType, ZFStride>;
  using ZoneToFaceRelation = slam::StaticRelation<PositionType, ElementType,
                                                  ZFCard, STLIndirection, ZoneSet, FaceSet>;
  using ZFaceSet = ZoneToFaceRelation::RelationSubset;

  using NUM_BD_SZ = slam::policies::CompileTimeSize<ZoneSet::PositionType, NUM_DOMAIN_BOUNDARIES>;
  using BoundaryEdgeSet = slam::OrderedSet< PositionType, ElementType, NUM_BD_SZ>;

  using IndexMap = slam::Map<IndexType>;

  using IndexRegistry = slam::FieldRegistry<SetBase, ZoneSet::PositionType>;
  using IndexBuffer = IndexRegistry::BufferType;

  struct DataRegistry
  {
    static IndexRegistry setRegistry;
  };

} // end namespace tinyHydro

#endif
