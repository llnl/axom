// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SLAM_FIELD_REGISTRY_H_
#define SLAM_FIELD_REGISTRY_H_

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/Map.hpp"

#include <functional>
#include <map>
#include <optional>
#include <string>
#include <string_view>

namespace axom::slam
{
/**
 * \brief Simple container for fields of type DataType w/ minimal error checking
 *
 * \note We are using concrete instances for int and double in the code below.
 *       This should eventually be replaced with the sidre datastore.
 *
 * \note FieldRegistry is a host-only facility: it stores std::map tables
 *       keyed by std::string and its find APIs return std::optional. The
 *       three lookup tables use transparent comparison (std::less<>) so callers may
 *       query with a std::string_view without constructing a temporary std::string.
 */
template <typename SetType, typename TheDataType>
class FieldRegistry
{
public:
  using PositionType = typename SetType::PositionType;
  using DataType = TheDataType;
  using KeyType = std::string;
  using MapType = slam::Map<DataType, SetType>;
  using BufferType = typename MapType::OrderedMap;

  // Transparent comparator (std::less<>) enables heterogeneous lookup: a
  // std::string_view (or const char*) can be used as a key without allocating.
  using DataVecMap = std::map<KeyType, MapType, std::less<>>;
  using DataBufferMap = std::map<KeyType, BufferType, std::less<>>;
  using DataAttrMap = std::map<KeyType, DataType, std::less<>>;

public:
  [[nodiscard]] bool hasField(std::string_view key) const
  {
    return m_maps.find(key) != m_maps.end();
  }

  MapType& addField(KeyType key, const SetType* theSet)
  {
    return m_maps[std::move(key)] = MapType(theSet);
  }

  MapType& addNamelessField(const SetType* theSet)
  {
    static int cnt = 0;
    return m_maps[axom::fmt::format("__field_{}", cnt++)] = MapType(theSet);
  }

  MapType& getField(std::string_view key)
  {
    verifyFieldsKey(key);
    return m_maps.find(key)->second;
  }

  const MapType& getField(std::string_view key) const
  {
    verifyFieldsKey(key);
    return m_maps.find(key)->second;
  }

  /*!
   * \brief Find a field by name without inserting or asserting.
   * \return an optional referencing the field if present, else empty.
   */
  [[nodiscard]] std::optional<std::reference_wrapper<MapType>> findField(std::string_view key)
  {
    auto it = m_maps.find(key);
    return it != m_maps.end() ? std::optional<std::reference_wrapper<MapType>>(it->second)
                              : std::nullopt;
  }

  [[nodiscard]] std::optional<std::reference_wrapper<const MapType>> findField(std::string_view key) const
  {
    auto it = m_maps.find(key);
    return it != m_maps.end() ? std::optional<std::reference_wrapper<const MapType>>(it->second)
                              : std::nullopt;
  }

  [[nodiscard]] bool hasBuffer(std::string_view key) const
  {
    return m_buff.find(key) != m_buff.end();
  }

  BufferType& addBuffer(KeyType key, int size = 0)
  {
    return m_buff[std::move(key)] = BufferType(size);
  }

  BufferType& addNamelessBuffer(int size = 0)
  {
    static int cnt = 0;
    return m_buff[axom::fmt::format("__buffer_{}", cnt++)] = BufferType(size);
  }

  BufferType& getBuffer(std::string_view key)
  {
    verifyBufferKey(key);
    return m_buff.find(key)->second;
  }
  const BufferType& getBuffer(std::string_view key) const
  {
    verifyBufferKey(key);
    return m_buff.find(key)->second;
  }

  /*!
   * \brief Find a buffer by name without inserting or asserting.
   * \return an optional referencing the buffer if present, else empty.
   */
  [[nodiscard]] std::optional<std::reference_wrapper<BufferType>> findBuffer(std::string_view key)
  {
    auto it = m_buff.find(key);
    return it != m_buff.end() ? std::optional<std::reference_wrapper<BufferType>>(it->second)
                              : std::nullopt;
  }

  [[nodiscard]] bool hasScalar(std::string_view key) const
  {
    return m_scal.find(key) != m_scal.end();
  }

  DataType& addScalar(KeyType key, DataType val) { return m_scal[std::move(key)] = val; }

  DataType& getScalar(std::string_view key)
  {
    verifyScalarKey(key);
    return m_scal.find(key)->second;
  }

  const DataType& getScalar(std::string_view key) const
  {
    verifyScalarKey(key);
    return m_scal.find(key)->second;
  }

  /*!
   * \brief Find a scalar by name without inserting or asserting.
   * \return an engaged optional with the scalar's value if present, else empty.
   */
  [[nodiscard]] std::optional<DataType> findScalar(std::string_view key) const
  {
    auto it = m_scal.find(key);
    return it != m_scal.end() ? std::optional<DataType>(it->second) : std::nullopt;
  }

private:
  inline void verifyFieldsKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
    SLIC_ASSERT_MSG(hasField(key), "Didn't find field named " << key);
  }

  inline void verifyBufferKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
    SLIC_ASSERT_MSG(hasBuffer(key), "Didn't find buffer named " << key);
  }

  inline void verifyScalarKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
    SLIC_ASSERT_MSG(hasScalar(key), "Didn't find scalar named " << key);
  }

private:
  DataVecMap m_maps;
  DataBufferMap m_buff;
  DataAttrMap m_scal;
};

}  // end namespace axom::slam

#endif  // SLAM_FIELD_REGISTRY_H_
