// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_DYNAMIC_MAP_HPP_
#define SLAM_DYNAMIC_MAP_HPP_

#include <vector>
#include <sstream>

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/slam/DynamicSet.hpp"
#include "axom/slam/Map.hpp"

namespace axom
{
namespace slam
{
/**
 * \class DynamicMap
 * \brief A slam map class that supports adding and removing entries.
 *
 * \detail An entry in the map is considered valid if
 * its corresponding set's entry is valid
 */
template <typename SetType, typename DataType>
class DynamicMap
{
public:
  using SetPosition = typename SetType::PositionType;
  using SetElement = typename SetType::ElementType;

  using OrderedMap = std::vector<DataType>;

public:
  /** \brief Default constructor   */
  DynamicMap() : m_set(nullptr) { }

  /**
   * \brief Constructor from a set pointer
   *
   * \param theSet A pointer to the map's set

   * The map will be allocated with theSet->size() entries.
   * There is no guarantee that the values will be initialized
   */
  DynamicMap(SetType* theSet) : m_set(theSet)
  {
    if(m_set != nullptr)
    {
      m_data.resize(m_set->size());
    }
  }

  /**
   * \brief Constructor from a set pointer
   *
   * \param theSet A pointer to the map's set
   * \param defaultValue The value that each entry in the map will
   * be initialized
   *
   * The map will be allocated with \a theSet->size() entries.
   * Each entry will have value \a defaultValue
   */
  DynamicMap(SetType* theSet, DataType defaultValue) : m_set(theSet)
  {
    if(m_set != nullptr)
    {
      m_data.resize(m_set->size(), defaultValue);
    }
  }

  ~DynamicMap() { }

public:
  /** \brief Returns a pointer to the map's underlying set */
  const SetType* set() const { return m_set; }

  /// \name DynamicMap individual access functions
  /// @{
  ///

  /** \brief Return the value at set index \a setIndex   */
  const DataType& operator[](SetPosition setIndex) const
  {
    verifyPosition(setIndex);
    return m_data[setIndex];
  }

  /// @}

  /** \brief Access to underlying data */
  OrderedMap& data() { return m_data; }

  /** \brief Const access to underlying data */
  const OrderedMap& data() const { return m_data; }

  /// \name DynamicMap cardinality functions
  /// @{

  /** \brief Returns the size of map's set */
  SetPosition size() const { return m_data.size(); }

  /**
   * \brief Return the number of valid entries
   *
   * An entry at a given index is considered valid if corresponding
   * set element is valid.
   */
  SetPosition numberOfValidEntries() const
  {
    return (m_set != nullptr) ? m_set->numberOfValidEntries() : 0;
  }

  /// @}

  /// \name DynamicMap validity check functions
  /// @{

  bool isValidEntry(SetPosition pos) const
  {
    return (m_set != nullptr) ? m_set->isValidEntry(pos) : false;
  }

  /** \brief Predicate to check if this DynamicMap instance is valid */
  bool isValid(bool verboseOutput = false) const;

  /// @}

private:
  /** \brief Debug check that the index is not out of range    */
  inline void verifyPosition(SetPosition AXOM_DEBUG_PARAM(setIndex)) const
  {
    SLIC_ASSERT_MSG(
      setIndex >= 0 && setIndex < (int)m_data.size(),
      "Attempted to access entry " << setIndex << " but map's set has size " << m_data.size());

    SLIC_ASSERT_MSG(isValidEntry(setIndex), "Attempted to access invalid set entry " << setIndex);
  }

public:
  /// \name Functions that modify the map's cardinality
  /// @{

  /**
   * \brief Get the value at position \a position
   *
   * \note Increases the map size if position is out of range
   */
  DataType& operator[](SetPosition position)
  {
    if(size() < position + 1)
    {
      resize(position + 1);
    }

    return m_data[position];
  }

  /**
   * \brief Insert \a value into \a position in the map
   *
   * \note Increases the map size if position is out of range
   */
  void insert(SetPosition position, DataType value) { operator[](position) = value; }

  /**
   * \brief Resizes the map to have at least \a s positions
   * \param s The minimum necessary capacity for resizing
   * \pre s >= 0
   */
  void resize(SetPosition s)
  {
    // Note (KW): For this to be a valid DynamicMap operation,
    // we would need to also increase the size of the map's set!

    // Note (KW): Do we want this to shrink the size of the map
    // when s < size() ?

    SLIC_ASSERT_MSG(s >= 0, "Attempted to resize vector with a negative size " << s);

    m_data.resize(s);
  }

  /// @}

private:
  SetType* m_set;
  OrderedMap m_data;
};

template <typename SetType, typename DataType>
bool DynamicMap<SetType, DataType>::isValid(bool verboseOutput) const
{
  bool bValid = true;

  std::stringstream errStr;

  if(m_set == nullptr)
  {
    if(!m_data.empty())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set was never provided,"
               << " but its associated data is not empty"
               << " , data has size " << m_data.size();
      }
      bValid = false;
    }
  }
  else
  {
    // Check the data array and set data have equal size
    if(static_cast<SetPosition>(m_data.size()) != m_set->size())
    {
      if(verboseOutput)
      {
        errStr << "\n\t* the underlying set and its associated mapped data"
               << " have different sizes, underlying set has size " << m_set->size()
               << " , data has size " << m_data.size();
        ;
      }

      bValid = false;
    }
  }

  if(verboseOutput)
  {
    if(bValid)
    {
      SLIC_DEBUG("Map was valid.");
    }
    else
    {
      SLIC_DEBUG("Map was not valid. " << errStr.str());
    }
  }

  return bValid;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_DYNAMIC_MAP_HPP_
