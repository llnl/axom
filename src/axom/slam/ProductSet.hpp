// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file ProductSet.hpp
 *
 * \brief Basic API for a SLAM Cartesian product set
 *
 */

#ifndef SLAM_PRODUCT_SET_H_
#define SLAM_PRODUCT_SET_H_

#include "axom/core/IteratorBase.hpp"
#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/RangeSet.hpp"

#include "axom/slam/policies/BivariateSetInterfacePolicies.hpp"

#include <cassert>

namespace axom
{
namespace slam
{
/**
 * \class ProductSet
 *
 * \brief Models a set whose element is the Cartesian product of two sets. The
 *        number of elements in this set is the product of the sizes of the two
 *        input sets.
 *
 *        Users should refer to the BivariateSet documentation for descriptions
 *        of the different indexing names (SparseIndex, DenseIndex, FlatIndex).
 *
 * \see   BivariateSet
 */
template <typename SetType1 = slam::Set<>,
          typename SetType2 = slam::Set<>,
          typename InterfaceType = policies::VirtualInterface>
class ProductSet final : public policies::BivariateSetInterface<InterfaceType, SetType1, SetType2>
{
private:
  using BaseType = policies::BivariateSetInterface<InterfaceType, SetType1, SetType2>;

public:
  using RangeSetType = typename BaseType::RangeSetType;
  using FirstSetType = SetType1;
  using SecondSetType = SetType2;
  using PositionType = typename BaseType::PositionType;
  using ElementType = typename BaseType::ElementType;
  using ProductSetType = ProductSet<SetType1, SetType2>;

  using IteratorType = BivariateSetIterator<ProductSet>;

private:
  template <typename Dummy, typename SetType>
  struct RowSet
  {
    using Type = SetType;
    RowSet(PositionType secondSetSize) : m_data(secondSetSize)
    {
      //fill in the row data now for getElements(i) function,
      //since every row is the same, a call to getElements() returns the same set.
      //
      // HACK -- this should actually be returning a PositionSet since it always
      //         goes from 0 to secondSetSize()
      // This requires a change to the return type of BivariateSet::getElements()
      std::iota(m_data.begin(), m_data.end(), 0);
      m_set = typename SetType::SetBuilder().size(secondSetSize).offset(0).data(m_data.view());
    }

    Type get(PositionType) const { return m_set; }

    axom::Array<PositionType> m_data;
    SetType m_set;
  };

  // This dummy parameter works around a C++ bug where partial specializations,
  // but not explicit instantiations, are allowed in class scope.
  // This is fixed in C++17: https://wg21.cmeerw.net/cwg/issue727
  template <typename Dummy>
  struct RowSet<Dummy, void>
  {
    using Type = PositionSet<PositionType, ElementType>;

    RowSet(PositionType) { }

    Type get(PositionType secondSetSize) const { return Type(secondSetSize); }
  };

public:
  using ConcreteSet = ProductSet<SetType1, SetType2, policies::ConcreteInterface>;
  using VirtualSet = ProductSet<SetType1, SetType2, policies::VirtualInterface>;

  using OtherSet =
    std::conditional_t<std::is_same<InterfaceType, policies::VirtualInterface>::value, ConcreteSet, VirtualSet>;

  ProductSet(const OtherSet& other)
    : BaseType(other.getFirstSet(), other.getSecondSet())
    , m_rowSet(this->secondSetSize())
  { }

public:
  using SubsetType = typename RowSet<void, typename BaseType::SubsetType>::Type;

  /** \brief Default constructor */
  ProductSet() : m_rowSet(0) { }

  /**
   * \brief Constructor taking in pointers of two Sets.
   *
   * \param set1  Pointer to the first Set.
   * \param set2  Pointer to the second Set.
   */

  ProductSet(const FirstSetType* set1, const SecondSetType* set2)
    : BaseType(set1, set2)
    , m_rowSet(this->secondSetSize())
  { }

  /**
   * \brief Return the element SparseIndex. Since ProductSet is the full
   *        Cartesian product of the two sets, the SparseIndex is the same as
   *        its DenseIndex, which is the same as the \a pos2 parameter.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The element's SparseIndex, which is the same as pos2
   */
  PositionType findElementIndex(PositionType pos1, PositionType pos2) const
  {
    verifyPositionImpl(pos1, pos2);
    return pos2;
  }

  /**
   * \brief Returns an element's FlatIndex given its DenseIndex. Since
   *        ProductSet is the full Cartesian product of the two sets, an
   *        element's FlatIndex is equal to `pos1*secondSetSize()+pos2`.
   *
   * \param pos1  The first set position.
   * \param pos2  The second set position.
   *
   * \return  The element's FlatIndex.
   */
  AXOM_HOST_DEVICE PositionType findElementFlatIndex(PositionType pos1, PositionType pos2) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(pos1, pos2);
#endif
    PositionType size2 = this->secondSetSize();
    PositionType pos = size2 * pos1 + pos2;

    return pos;
  }

  /**
   * \brief Returns the FlatIndex of the first element in the specified row.
   *        This is equal to `pos1*secondSetSize()`.
   *
   * \param pos1  The first set position that specifies the row.
   */
  PositionType findElementFlatIndex(PositionType pos1) const
  {
    return findElementFlatIndex(pos1, 0);
  }

  /**
   * \brief Given the flat index, return the associated to-set index in the
   *        relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos2  The to-set index.
   */
  AXOM_HOST_DEVICE PositionType flatToSecondIndex(PositionType flatIndex) const
  {
    if(flatIndex < 0 || flatIndex > size())
    {
      SLIC_ASSERT("Flat index out of bounds of the relation set.");
    }
    return flatIndex % this->secondSetSize();
  }

  /**
   * \brief Given the flat index, return the associated from-set index in the
   *        relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos1  The from-set index.
   */
  AXOM_HOST_DEVICE PositionType flatToFirstIndex(PositionType flatIndex) const
  {
    if(flatIndex < 0 || flatIndex > size())
    {
      SLIC_ASSERT("Flat index out of bounds of the relation set.");
    }
    return flatIndex / this->secondSetSize();
  }

  /**
   * \brief Return all elements from the second set associated with position
   *        \a pos1 in the first set.
   *
   * \param pos1   The first set position that specifies the row.
   *
   * \return  An OrderedSet of the elements in the row.
   */
  SubsetType getElements(PositionType AXOM_DEBUG_PARAM(pos1)) const
  {
    SLIC_ASSERT(pos1 >= 0 && pos1 < this->firstSetSize());

    return m_rowSet.get(this->secondSetSize());
  }

  AXOM_HOST_DEVICE ElementType at(PositionType pos) const { return pos % this->secondSetSize(); }

  AXOM_HOST_DEVICE PositionType size() const
  {
    return this->firstSetSize() * this->secondSetSize();
  }

  PositionType size(PositionType) const { return this->secondSetSize(); }

  /*!
   * \brief Return an iterator to the first pair of set elements in the
   *  relation.
   */
  IteratorType begin() const { return IteratorType(this, 0); }

  /*!
   * \brief Return an iterator to one past the last pair of set elements in the
   *  relation.
   */
  IteratorType end() const { return IteratorType(this, size()); }

  AXOM_HOST_DEVICE RangeSetType elementRangeSet(PositionType pos1) const
  {
    const auto sz = this->secondSetSize();
    return typename RangeSetType::SetBuilder().size(sz).offset(sz * pos1);
  }

  bool isValidIndex(PositionType s1, PositionType s2) const
  {
    PositionType size1 = this->firstSetSize();
    PositionType size2 = this->secondSetSize();
    return s1 >= 0 && s1 < size1 && s2 >= 0 && s2 < size2;
  }

  bool isValid(bool verboseOutput = false) const { return BaseType::isValid(verboseOutput); }

private:
  /** \brief verify the FlatIndex \a pos is within the valid range. */
  void verifyPosition(PositionType pos) const
  {  //from RangeSet, overloading to avoid warning in compiler
    verifyPositionImpl(pos);
  }

  /** \brief implementation for verifyPosition */
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos)) const
  {  //from RangeSet, overloading to avoid warning in compiler
    SLIC_ASSERT_MSG(pos >= 0 && pos < size(),
                    "SLAM::ProductSet -- requested out-of-range element at position "
                      << pos << ", but set only has " << size() << " elements.");
  }

  /** \brief verify the SparseIndex (which is the same as its DenseIndex) is
   *         within the valid range. */
  void verifyPosition(PositionType pos1, PositionType pos2) const
  {
    verifyPositionImpl(pos1, pos2);
  }

  /** \brief verify the SparseIndex (which is the same as its DenseIndex) is
   *         within the valid range. */
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos1),
                                 PositionType AXOM_DEBUG_PARAM(pos2)) const
  {
    SLIC_ASSERT_MSG(isValidIndex(pos1, pos2),
                    "SLAM::ProductSet -- requested out-of-range element at position ("
                      << pos1 << "," << pos2 << "), but set only has " << this->firstSetSize()
                      << "x" << this->secondSetSize() << " elements.");
  }

private:
  RowSet<void, typename BaseType::SubsetType> m_rowSet;
};

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_PRODUCT_SET_H
