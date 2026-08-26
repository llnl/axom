// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file ProductSet.hpp
 *
 * \brief Basic API for a SLAM Cartesian product set
 */

#include "axom/core/IteratorBase.hpp"
#include "axom/slam/BivariateSet.hpp"
#include "axom/slam/RangeSet.hpp"

#include "axom/slam/policies/BivariateSetInterfacePolicies.hpp"

#include <cassert>
#include <optional>

namespace axom::slam
{
/**
 * \class ProductSet
 *
 * \brief Models a set whose element is the Cartesian product of two sets. 
 *        The number of elements in this set is the product of the sizes of the two input sets.
 *
 *        Users should refer to the BivariateSet documentation for descriptions
 *        of the different indexing names (SparseIndex, DenseIndex, FlatIndex).
 *
 * \see   BivariateSet
 */
template <typename SetType1 = slam::Set<>,
          typename SetType2 = slam::Set<>,
          typename InterfaceType = policies::VirtualInterface,
          typename FlatPosition = detail::default_flat_position_t<typename SetType1::PositionType,
                                                                  typename SetType2::PositionType>>
class ProductSet final
  : public policies::BivariateSetInterface<InterfaceType, SetType1, SetType2, FlatPosition>
{
private:
  using BaseType = policies::BivariateSetInterface<InterfaceType, SetType1, SetType2, FlatPosition>;

public:
  using RangeSetType = typename BaseType::RangeSetType;
  using FirstSetType = SetType1;
  using SecondSetType = SetType2;
  using FirstPositionType = typename BaseType::FirstPositionType;
  using SecondPositionType = typename BaseType::SecondPositionType;
  using PositionType = typename BaseType::PositionType;
  using ElementType = typename BaseType::ElementType;
  using ProductSetType = ProductSet<SetType1, SetType2, InterfaceType, FlatPosition>;

  using IteratorType = BivariateSetIterator<ProductSet>;

private:
  template <typename Dummy, typename SetType>
  struct RowSet
  {
    using Type = SetType;
    RowSet(SecondPositionType secondSetSize) : m_data(static_cast<axom::IndexType>(secondSetSize))
    {
      //fill in the row data now for getElements(i) function,
      //since every row is the same, a call to getElements() returns the same set.
      //
      // HACK -- this should actually be returning a PositionSet since it always
      //         goes from 0 to secondSetSize()
      // This requires a change to the return type of BivariateSet::getElements()
      std::iota(m_data.begin(), m_data.end(), 0);
      m_set = typename SetType::SetBuilder()
                .size(static_cast<PositionType>(secondSetSize))
                .offset(0)
                .data(m_data.view());
    }

    Type get(SecondPositionType) const { return m_set; }

    axom::Array<SecondPositionType> m_data;
    SetType m_set;
  };

  // This dummy parameter works around a C++ bug where partial specializations,
  // but not explicit instantiations, are allowed in class scope.
  // This is fixed in C++17: https://wg21.cmeerw.net/cwg/issue727
  template <typename Dummy>
  struct RowSet<Dummy, void>
  {
    using Type = PositionSet<PositionType, SecondPositionType>;

    RowSet(SecondPositionType) { }

    Type get(SecondPositionType secondSetSize) const
    {
      return Type(static_cast<PositionType>(secondSetSize));
    }
  };

public:
  using ConcreteSet = ProductSet<SetType1, SetType2, policies::ConcreteInterface, FlatPosition>;
  using VirtualSet = ProductSet<SetType1, SetType2, policies::VirtualInterface, FlatPosition>;

  using OtherSet =
    std::conditional_t<std::is_same<InterfaceType, policies::VirtualInterface>::value, ConcreteSet, VirtualSet>;

  ProductSet(const OtherSet& other)
    : BaseType(other.getFirstSet(), other.getSecondSet())
    , m_rowSet(this->secondSetSize())
  { }

public:
  using SubsetType = typename RowSet<void, typename BaseType::SubsetType>::Type;

  /// \brief Default constructor
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
  PositionType findElementIndex(FirstPositionType pos1, SecondPositionType pos2) const
  {
    verifyPositionImpl(pos1, pos2);
    return static_cast<PositionType>(pos2);
  }

  /*!
   * \brief Optional-returning wrapper for `findElementIndex`.
   *
   * \return An engaged `std::optional` with the SparseIndex if valid, else empty.
   */
  [[nodiscard]] std::optional<PositionType> findElementIndexOptional(FirstPositionType pos1,
                                                                     SecondPositionType pos2) const
  {
    const auto idx = findElementIndex(pos1, pos2);
    return idx != BaseType::INVALID_POS ? std::optional<PositionType>(idx)
                                        : std::optional<PositionType> {};
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
  AXOM_HOST_DEVICE PositionType findElementFlatIndex(FirstPositionType pos1,
                                                     SecondPositionType pos2) const
  {
#ifndef AXOM_DEVICE_CODE
    verifyPositionImpl(pos1, pos2);
#endif
    const PositionType size2 = static_cast<PositionType>(this->secondSetSize());
    const PositionType pos =
      size2 * static_cast<PositionType>(pos1) + static_cast<PositionType>(pos2);

    return pos;
  }

  /*!
   * \brief Optional-returning wrapper for `findElementFlatIndex(pos1, pos2)`.
   *
   * \return An engaged `std::optional` with the FlatIndex if valid, else empty.
   */
  [[nodiscard]] AXOM_HOST_DEVICE std::optional<PositionType> findElementFlatIndexOptional(
    FirstPositionType pos1,
    SecondPositionType pos2) const
  {
    const auto idx = findElementFlatIndex(pos1, pos2);
    return idx != BaseType::INVALID_POS ? std::optional<PositionType>(idx)
                                        : std::optional<PositionType> {};
  }

  /**
   * \brief Returns the FlatIndex of the first element in the specified row.
   *        This is equal to `pos1*secondSetSize()`.
   *
   * \param pos1  The first set position that specifies the row.
   */
  PositionType findElementFlatIndex(FirstPositionType pos1) const
  {
    SLIC_ASSERT_MSG(pos1 >= 0 && pos1 < this->firstSetSize(),
                    "SLAM::ProductSet -- requested out-of-range first-set position "
                      << pos1 << ", but set only has " << this->firstSetSize() << " rows.");

    if(this->secondSetSize() == 0)
    {
      return BaseType::INVALID_POS;
    }

    return findElementFlatIndex(pos1, 0);
  }

  /*!
   * \brief Optional-returning wrapper for `findElementFlatIndex(pos1)`.
   *
   * \return An engaged `std::optional` with the FlatIndex if valid, else empty.
   */
  [[nodiscard]] std::optional<PositionType> findElementFlatIndexOptional(FirstPositionType pos1) const
  {
    const auto idx = findElementFlatIndex(pos1);
    return idx != BaseType::INVALID_POS ? std::optional<PositionType>(idx)
                                        : std::optional<PositionType> {};
  }

  /**
   * \brief Given the flat index, return the associated to-set index in the relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos2  The to-set index.
   */
  AXOM_HOST_DEVICE SecondPositionType flatToSecondIndex(PositionType flatIndex) const
  {
#ifndef AXOM_DEVICE_CODE
    SLIC_ASSERT_MSG(flatIndex >= 0 && flatIndex < size(),
                    "SLAM::ProductSet -- requested out-of-range flat index "
                      << flatIndex << "; set has " << size() << " elements.");
#endif
    return static_cast<SecondPositionType>(flatIndex %
                                           static_cast<PositionType>(this->secondSetSize()));
  }

  /**
   * \brief Given the flat index, return the associated from-set index in the relation pair.
   *
   * \param flatIndex The FlatIndex of the from-set/to-set pair.
   *
   * \return pos1  The from-set index.
   */
  AXOM_HOST_DEVICE FirstPositionType flatToFirstIndex(PositionType flatIndex) const
  {
#ifndef AXOM_DEVICE_CODE
    SLIC_ASSERT_MSG(flatIndex >= 0 && flatIndex < size(),
                    "SLAM::ProductSet -- requested out-of-range flat index "
                      << flatIndex << "; set has " << size() << " elements.");
#endif
    return static_cast<FirstPositionType>(flatIndex /
                                          static_cast<PositionType>(this->secondSetSize()));
  }

  /**
   * \brief Return all elements from the second set associated with position
   *        \a pos1 in the first set.
   *
   * \param pos1   The first set position that specifies the row.
   *
   * \return  An OrderedSet of the elements in the row.
   */
  SubsetType getElements(FirstPositionType AXOM_DEBUG_PARAM(pos1)) const
  {
    SLIC_ASSERT(pos1 >= 0 && pos1 < this->firstSetSize());

    return m_rowSet.get(this->secondSetSize());
  }

  [[nodiscard]] AXOM_HOST_DEVICE ElementType at(PositionType pos) const
  {
    return {flatToFirstIndex(pos), flatToSecondIndex(pos)};
  }

  AXOM_HOST_DEVICE PositionType size() const
  {
    return static_cast<PositionType>(this->firstSetSize()) *
      static_cast<PositionType>(this->secondSetSize());
  }

  PositionType size(FirstPositionType) const
  {
    return static_cast<PositionType>(this->secondSetSize());
  }

  /// \brief Return an iterator to the first pair of set elements in the relation.
  IteratorType begin() const { return IteratorType(this, 0); }

  /// \brief Return an iterator to one past the last pair of set elements in the relation.
  IteratorType end() const { return IteratorType(this, size()); }

  AXOM_HOST_DEVICE RangeSetType elementRangeSet(FirstPositionType pos1) const
  {
    const auto sz = static_cast<PositionType>(this->secondSetSize());
    return typename RangeSetType::SetBuilder().size(sz).offset(sz * static_cast<PositionType>(pos1));
  }

  [[nodiscard]] bool isValidIndex(FirstPositionType s1, SecondPositionType s2) const
  {
    FirstPositionType size1 = this->firstSetSize();
    SecondPositionType size2 = this->secondSetSize();
    return s1 >= 0 && s1 < size1 && s2 >= 0 && s2 < size2;
  }

  [[nodiscard]] bool isValid(bool verboseOutput = false) const
  {
    return BaseType::isValid(verboseOutput);
  }

private:
  /// \brief verify the FlatIndex \a pos is within the valid range.
  void verifyPosition(PositionType pos) const
  {  //from RangeSet, overloading to avoid warning in compiler
    verifyPositionImpl(pos);
  }

  /// \brief implementation for verifyPosition
  inline void verifyPositionImpl(PositionType AXOM_DEBUG_PARAM(pos)) const
  {  //from RangeSet, overloading to avoid warning in compiler
    SLIC_ASSERT_MSG(pos >= 0 && pos < size(),
                    "SLAM::ProductSet -- requested out-of-range element at position "
                      << pos << ", but set only has " << size() << " elements.");
  }

  /// \brief verify the SparseIndex (which is the same as its DenseIndex) is within the valid range.
  void verifyPosition(FirstPositionType pos1, SecondPositionType pos2) const
  {
    verifyPositionImpl(pos1, pos2);
  }

  /// \brief verify the SparseIndex (which is the same as its DenseIndex) is within the valid range.
  inline void verifyPositionImpl(FirstPositionType AXOM_DEBUG_PARAM(pos1),
                                 SecondPositionType AXOM_DEBUG_PARAM(pos2)) const
  {
    SLIC_ASSERT_MSG(isValidIndex(pos1, pos2),
                    "SLAM::ProductSet -- requested out-of-range element at position ("
                      << pos1 << "," << pos2 << "), but set only has " << this->firstSetSize()
                      << "x" << this->secondSetSize() << " elements.");
  }

private:
  RowSet<void, typename BaseType::SubsetType> m_rowSet;
};

}  // end namespace axom::slam
