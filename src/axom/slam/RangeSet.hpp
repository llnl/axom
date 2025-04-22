// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file RangeSet.hpp
 *
 * \brief Basic API for an ordered set of entities in a simulation
 *
 */

#ifndef SLAM_RANGE_SET_H_
#define SLAM_RANGE_SET_H_

#include "axom/slam/OrderedSet.hpp"

namespace axom
{
namespace slam
{
/**
 * \class GenericRangeSet
 * \brief Models a set whose elements belong to a contiguous range
 *  \f$ \in [lowerIndex,upperIndex) \f$
 *
 * \note The \a ElementType here needs to be computable as offsets
 *  (of \a PositionType) from the lowerIndex.
 *  Examples include: signed and unsigned integral types.
 *  This version of a range set still allows you to have
 *  different policies on striding, indirection and subsetting
 *  \sa OrderedSet, PositionSet, RangeSet
 */
template <typename P = slam::DefaultPositionType,
          typename E = slam::DefaultElementType,
          typename OffsetPolicy = policies::RuntimeOffset<P>,
          typename StridingPolicy = policies::StrideOne<P>,
          typename IndirectionPolicy = policies::NoIndirection<P, E>,
          typename SubsettingPolicy = policies::NoSubset,
          typename InterfacePolicy = policies::VirtualInterface>
class GenericRangeSet
  : public OrderedSet<P, E, policies::RuntimeSize<P>, OffsetPolicy, StridingPolicy, IndirectionPolicy, SubsettingPolicy, InterfacePolicy>
{
public:
  using PositionType = P;
  using ElementType = E;

private:
  using OrderedSetType =
    OrderedSet<P, E, policies::RuntimeSize<P>, OffsetPolicy, StridingPolicy, IndirectionPolicy, SubsettingPolicy, InterfacePolicy>;

  template <typename OtherIndirectionPolicy, typename OtherInterfaceType>
  using ConvertibleRangeSet =
    GenericRangeSet<P, E, OffsetPolicy, StridingPolicy, OtherIndirectionPolicy, SubsettingPolicy, OtherInterfaceType>;

public:
  using ConcreteSet = ConvertibleRangeSet<IndirectionPolicy, policies::ConcreteInterface>;

  using VirtualSet = ConvertibleRangeSet<IndirectionPolicy, policies::VirtualInterface>;

  using OtherIfaceSet =
    std::conditional_t<std::is_same<InterfacePolicy, policies::VirtualInterface>::value, ConcreteSet, VirtualSet>;

  GenericRangeSet(const OtherIfaceSet& otherSet) : OrderedSetType(otherSet) { }

public:
  GenericRangeSet(PositionType size = OrderedSetType::SizePolicyType::DEFAULT_VALUE)
    : OrderedSetType(size,
                     OrderedSetType::OffsetPolicyType::DEFAULT_VALUE,
                     OrderedSetType::StridePolicyType::DEFAULT_VALUE)
  { }

  AXOM_HOST_DEVICE GenericRangeSet(const typename OrderedSetType::SetBuilder& builder)
    : OrderedSetType(builder)
  { }

  GenericRangeSet(PositionType lowerIndex, PositionType upperIndex)
    : OrderedSetType(upperIndex - lowerIndex,
                     lowerIndex,
                     OrderedSetType::StridePolicyType::DEFAULT_VALUE)
  { }
};

/**
 * \class PositionSet
 * \brief Alias template for an OrderedSet whose elements belong
 * to a contiguous range \f$ \in [0,size) \f$
 *
 * \tparam P The PositionType
 * \tparam E The ElementType
 * \sa OrderedSet
 */
template <typename P = slam::DefaultPositionType, typename E = slam::DefaultElementType>
using PositionSet = GenericRangeSet<P, E, policies::ZeroOffset<P>>;

/**
 * \class RangeSet
 * \brief An alias of GenericRangeSet with stride 1 and no indirection
 *
 * \tparam P The PositionType
 * \tparam E The ElementType
 *  \sa GenericRangeSet, OrderedSet, PositionSet
 *
 * \note The \a ElementType needs to be computable as offsets (of
 * \a PositionType) from the lowerIndex.
 * Examples include: signed and unsigned integral types
 */
template <typename P = slam::DefaultPositionType, typename E = slam::DefaultElementType>
using RangeSet = GenericRangeSet<P, E>;

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_RANGE_SET_H_
