// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include "axom/slam/Utilities.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/Aliases.hpp"
#include "axom/slam/Map.hpp"
#include "axom/slam/MapBuilders.hpp"

#include <functional>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

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
 *       lookup tables use transparent comparison (std::less<>) so callers may
 *       query with a std::string_view without constructing a temporary std::string.
 *
 * \note FieldRegistry supports two field storage modes under a single set of field keys:
 *       - fields that manage their own buffer, stored as a `slam::Map` with the
 *         default `axom::Array` indirection
 *       - fields that refer to an externally-managed buffer, stored as a `slam::Map`
 *         with an `axom::ArrayView` indirection
 *         (useful when generating Slam objects from externally-managed arrays)
 *
 * \note Owning buffers and owning fields now use `axom::Array` rather than `std::vector`.
 *       Prefer `auto`, `FieldRegistry::BufferType`, and `FieldRegistry::MapType` 
 *       at registry boundaries. Code that needs a non-owning view
 *       should call `buffer.view()` or register storage with `addFieldView()` 
 *       or `addBufferView()` instead of assuming a `std::vector` interface.
 */
template <typename SetType, typename TheDataType>
class FieldRegistry
{
public:
  /// \brief The set's position type.
  using PositionType = typename SetType::PositionType;

  /// \brief Data type for scalars, buffers, and map values stored in this registry.
  using DataType = TheDataType;

  /// \brief Key type used for registry entries.
  using KeyType = std::string;

  /*!
   * \brief Field type stored by this registry, a `slam::Map` from \a SetType to \a DataType.
   *
   * Uses Slam's default indirection (`axom::Array`), so registry-managed
   * fields hold their own buffer.  See \ref FieldRegistry::BufferType.
   */
  using MapType = slam::Map<DataType, SetType>;

  /// \brief Buffer type backing `MapType` (uses the default Map buffer type: `axom::Array<DataType>`
  using BufferType = typename MapType::OrderedMap;

  /*!
   * \brief View-backed field type stored by this registry.
   *
   * This is the return type of `slam::make_map(set, axom::ArrayView<DataType>{...})`
   * and uses the canonical non-owning `policies::ArrayViewIndirection` policy.
   * The view's backing storage must outlive the registered field.
   */
  using ViewMapType = decltype(slam::make_map(std::declval<const SetType*>(),
                                              std::declval<axom::ArrayView<DataType>>()));

  /// \brief View-backed buffer type stored by this registry.
  using ViewBufferType = axom::ArrayView<DataType>;

  /// \brief Storage for one named field, either owning or view-backed.
  using FieldStorageType = std::variant<MapType, ViewMapType>;

  // Transparent comparator (std::less<>) enables heterogeneous lookup: a
  // std::string_view (or const char*) can be used as a key without allocating.
  using DataVecMap = std::map<KeyType, FieldStorageType, std::less<>>;
  using DataBufferMap = std::map<KeyType, BufferType, std::less<>>;
  using DataViewBufferMap = std::map<KeyType, ViewBufferType, std::less<>>;
  using DataAttrMap = std::map<KeyType, DataType, std::less<>>;

public:
  /// \name Fields
  /// @{

  /*!
   * \brief Returns true if a field with the given key exists.
   *
   * This function performs heterogeneous lookup, so callers may pass a
   * `std::string_view` (or `const char*`) without allocating.
   */
  [[nodiscard]] bool hasField(std::string_view key) const
  {
    return m_maps.find(key) != m_maps.end();
  }

  /*!
   * \brief Adds (or replaces) an owning field with the given key.
   *
   * \param key The field name.
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \return A mutable reference to the stored field.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  MapType& addField(std::string_view key, const SetType* theSet)
  {
    auto [it, _] = m_maps.insert_or_assign(KeyType(key), FieldStorageType(MapType(theSet)));
    return std::get<MapType>(it->second);
  }

  /*!
   * \brief Adds (or replaces) a view-backed field from an `axom::ArrayView`.
   *
   * \param key The field name.
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \param data View of externally-owned storage (must outlive the map).
   * \return A mutable reference to the stored field.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  ViewMapType& addField(std::string_view key, const SetType* theSet, axom::ArrayView<DataType> data)
  {
    auto [it, _] =
      m_maps.insert_or_assign(KeyType(key), FieldStorageType(slam::make_map(theSet, data)));
    return std::get<ViewMapType>(it->second);
  }

  /*!
   * \brief Adds (or replaces) a view-backed field from a raw pointer buffer.
   *
   * \param key The field name.
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \param data Pointer to externally-owned storage (must outlive the map).
   *             The helper wraps this as an `axom::ArrayView<DataType>` sized to `theSet->size()`.
   * \return A mutable reference to the stored field.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  ViewMapType& addField(std::string_view key, const SetType* theSet, DataType* data)
  {
    auto [it, _] =
      m_maps.insert_or_assign(KeyType(key), FieldStorageType(slam::make_map(theSet, data)));
    return std::get<ViewMapType>(it->second);
  }

  /*!
   * \brief Adds a new owning field using an auto-generated unique key.
   *
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \return A mutable reference to the stored field.
   *
   * \note The generated key is unique within this translation unit and template instantiation, 
   *       but this API is not intended to be used concurrently from multiple threads.
   */
  MapType& addNamelessField(const SetType* theSet)
  {
    return addField(makeNamelessFieldKey(), theSet);
  }

  /*!
   * \brief Adds a new view-backed field using an auto-generated unique key.
   *
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \param data View of externally-owned storage (must outlive the map).
   * \return A mutable reference to the stored field.
   */
  ViewMapType& addNamelessField(const SetType* theSet, axom::ArrayView<DataType> data)
  {
    return addField(makeNamelessFieldKey(), theSet, data);
  }

  /*!
   * \brief Returns a mutable reference to the owning field for \a key.
   *
   * \param key Field name.
   * \return A mutable reference to the stored field.
   *
   * \pre `findField(key)` has a value (the key exists and is owning, not view-backed).
   * \note In debug builds, this asserts if the key is missing or is view-backed.
   * \note For view-backed fields, use `getFieldView()`. Requesting an owning field for a
   *  view-backed key is a precondition violation: it asserts in debug builds, and in release
   *  builds the underlying `std::get` throws `std::bad_variant_access`.
   */
  MapType& getField(std::string_view key)
  {
    verifyFieldsKey(key);
    return std::get<MapType>(m_maps.find(key)->second);
  }

  /*!
   * \brief Returns a const reference to the owning field for \a key.
   *
   * \param key Field name.
   * \return A const reference to the stored field.
   *
   * \pre `findField(key)` has a value (the key exists and is owning, not view-backed).
   * \note In debug builds, this asserts if the key is missing or is view-backed.
   * \note For view-backed fields, use `getFieldView()`. Requesting an owning field for a
   *  view-backed key is a precondition violation: it asserts in debug builds, and in release
   *  builds the underlying `std::get` throws `std::bad_variant_access`.
   */
  const MapType& getField(std::string_view key) const
  {
    verifyFieldsKey(key);
    return std::get<MapType>(m_maps.find(key)->second);
  }

  /*!
   * \brief Finds an owning field by name without inserting or asserting.
   *
   * \param key Field name.
   * \return An optional referencing the field if present, else empty.
   *
   * \note This function never inserts a new entry.
   */
  [[nodiscard]] std::optional<std::reference_wrapper<MapType>> findField(std::string_view key)
  {
    auto it = m_maps.find(key);
    if(it == m_maps.end())
    {
      return std::nullopt;
    }

    auto* field = std::get_if<MapType>(&it->second);
    return field != nullptr ? std::optional<std::reference_wrapper<MapType>>(*field) : std::nullopt;
  }

  /*!
   * \brief Finds an owning field by name (const overload).
   *
   * \param key Field name.
   * \return An optional referencing the field if present, else empty.
   *
   * \note This function never inserts a new entry.
   */
  [[nodiscard]] std::optional<std::reference_wrapper<const MapType>> findField(std::string_view key) const
  {
    auto it = m_maps.find(key);
    if(it == m_maps.end())
    {
      return std::nullopt;
    }

    const auto* field = std::get_if<MapType>(&it->second);
    return field != nullptr ? std::optional<std::reference_wrapper<const MapType>>(*field)
                            : std::nullopt;
  }

  /// @}

  /// \name View-backed Fields
  /// @{

  /// \brief Returns true if a view-backed field with the given key exists.
  [[nodiscard]] bool hasFieldView(std::string_view key) const
  {
    auto it = m_maps.find(key);
    return it != m_maps.end() && std::holds_alternative<ViewMapType>(it->second);
  }

  /*!
   * \brief Adds (or replaces) a view-backed field from an `axom::ArrayView`.
   *
   * \param key The field name.
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \param data View of externally-owned storage (must outlive the map).
   * \return A mutable reference to the stored field.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  ViewMapType& addFieldView(std::string_view key, const SetType* theSet, axom::ArrayView<DataType> data)
  {
    return addField(key, theSet, data);
  }

  /*!
   * \brief Adds (or replaces) a view-backed field from a raw pointer buffer.
   *
   * \param key The field name.
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \param data Pointer to externally-owned storage (must outlive the map).
   *             The helper wraps this as an `axom::ArrayView<DataType>` sized to `theSet->size()`.
   * \return A mutable reference to the stored field.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  ViewMapType& addFieldView(std::string_view key, const SetType* theSet, DataType* data)
  {
    return addField(key, theSet, data);
  }

  /*!
   * \brief Adds a new view-backed field using an auto-generated unique key.
   *
   * \param theSet Pointer to the set associated with the field (must outlive the map).
   * \param data View of externally-owned storage (must outlive the map).
   * \return A mutable reference to the stored field.
   */
  ViewMapType& addNamelessFieldView(const SetType* theSet, axom::ArrayView<DataType> data)
  {
    return addNamelessField(theSet, data);
  }

  /*!
   * \brief Returns a mutable reference to the view-backed field for \a key.
   *
   * \pre `hasFieldView(key)` is true (the key exists and is view-backed, not owning).
   * \note For owning fields, use `getField()`. Requesting a view-backed field for an owning
   *  key is a precondition violation: it asserts in debug builds, and in release builds the
   *  underlying `std::get` throws `std::bad_variant_access`.
   */
  ViewMapType& getFieldView(std::string_view key)
  {
    verifyFieldViewKey(key);
    return std::get<ViewMapType>(m_maps.find(key)->second);
  }

  /*!
   * \brief Returns a const reference to the view-backed field for \a key.
   *
   * \pre `hasFieldView(key)` is true (the key exists and is view-backed, not owning).
   * \note For owning fields, use `getField()`. Requesting a view-backed field for an owning
   *  key is a precondition violation: it asserts in debug builds, and in release builds the
   *  underlying `std::get` throws `std::bad_variant_access`.
   */
  const ViewMapType& getFieldView(std::string_view key) const
  {
    verifyFieldViewKey(key);
    return std::get<ViewMapType>(m_maps.find(key)->second);
  }

  /*!
   * \brief Finds a view-backed field by name without inserting or asserting.
   *
   * \return An optional referencing the field if present, else empty.
   */
  [[nodiscard]] std::optional<ViewMapType> findFieldView(std::string_view key)
  {
    auto it = m_maps.find(key);
    if(it == m_maps.end())
    {
      return std::nullopt;
    }

    auto* field = std::get_if<ViewMapType>(&it->second);
    return field != nullptr ? std::optional<ViewMapType>(*field) : std::nullopt;
  }

  /*!
   * \brief Finds a view-backed field by name (const overload).
   *
   * \return An optional referencing the field if present, else empty.
   */
  [[nodiscard]] std::optional<ViewMapType> findFieldView(std::string_view key) const
  {
    auto it = m_maps.find(key);
    if(it == m_maps.end())
    {
      return std::nullopt;
    }

    const auto* field = std::get_if<ViewMapType>(&it->second);
    return field != nullptr ? std::optional<ViewMapType>(*field) : std::nullopt;
  }

  /// @}

  /// \name Owning Buffers
  /// @{

  /// \brief Returns true if an owning buffer with the given key exists.
  [[nodiscard]] bool hasBuffer(std::string_view key) const
  {
    return m_buff.find(key) != m_buff.end();
  }

  /*!
   * \brief Adds (or replaces) an owning buffer with the given key.
   *
   * \param key Buffer name.
   * \param size Initial buffer size.
   * \return A mutable reference to the stored buffer.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  BufferType& addBuffer(std::string_view key, int size = 0)
  {
    auto [it, _] = m_buff.insert_or_assign(KeyType(key), BufferType(size));
    return it->second;
  }

  /*!
   * \brief Adds a new owning buffer using an auto-generated unique key.
   *
   * \param size Initial buffer size.
   * \return A mutable reference to the stored buffer.
   */
  BufferType& addNamelessBuffer(int size = 0)
  {
    static int cnt = 0;
    return addBuffer(axom::fmt::format("__buffer_{}", cnt++), size);
  }

  /*!
   * \brief Returns a mutable reference to the owning buffer for \a key.
   *
   * \pre `hasBuffer(key)` is true.
   */
  BufferType& getBuffer(std::string_view key)
  {
    verifyBufferKey(key);
    return m_buff.find(key)->second;
  }

  /*!
   * \brief Returns a const reference to the owning buffer for \a key.
   *
   * \pre `hasBuffer(key)` is true.
   */
  const BufferType& getBuffer(std::string_view key) const
  {
    verifyBufferKey(key);
    return m_buff.find(key)->second;
  }

  /*!
   * \brief Finds an owning buffer by name without inserting or asserting.
   *
   * \return An optional referencing the buffer if present, else empty.
   */
  [[nodiscard]] std::optional<std::reference_wrapper<BufferType>> findBuffer(std::string_view key)
  {
    auto it = m_buff.find(key);
    return it != m_buff.end() ? std::optional<std::reference_wrapper<BufferType>>(it->second)
                              : std::nullopt;
  }

  /*!
   * \brief Finds an owning buffer by name (const overload).
   *
   * \return An optional referencing the buffer if present, else empty.
   */
  [[nodiscard]] std::optional<std::reference_wrapper<const BufferType>> findBuffer(
    std::string_view key) const
  {
    auto it = m_buff.find(key);
    return it != m_buff.end() ? std::optional<std::reference_wrapper<const BufferType>>(it->second)
                              : std::nullopt;
  }

  /// @}

  /// \name View-backed Buffers
  /// @{

  /// \brief Returns true if a view-backed buffer with the given key exists.
  [[nodiscard]] bool hasBufferView(std::string_view key) const
  {
    return m_view_buff.find(key) != m_view_buff.end();
  }

  /*!
   * \brief Adds (or replaces) a view-backed buffer from an `axom::ArrayView`.
   *
   * \param key Buffer name.
   * \param view View of externally-owned storage (must outlive the registry entry).
   * \return A mutable reference to the stored view.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  ViewBufferType& addBufferView(std::string_view key, ViewBufferType view)
  {
    auto [it, _] = m_view_buff.insert_or_assign(KeyType(key), view);
    return it->second;
  }

  /*!
   * \brief Adds (or replaces) a view-backed buffer from a raw pointer and size.
   *
   * \param key Buffer name.
   * \param data Pointer to externally-owned storage (must outlive the registry entry).
   * \param size Number of elements.
   * \return A mutable reference to the stored view.
   */
  ViewBufferType& addBufferView(std::string_view key, DataType* data, PositionType size)
  {
    return addBufferView(key, axom::ArrayView<DataType>(data, static_cast<axom::IndexType>(size)));
  }

  /*!
   * \brief Adds a new view-backed buffer using an auto-generated unique key.
   *
   * \param view View of externally-owned storage (must outlive the registry entry).
   * \return A mutable reference to the stored view.
   */
  ViewBufferType& addNamelessBufferView(ViewBufferType view)
  {
    static int cnt = 0;
    return addBufferView(axom::fmt::format("__buffer_view_{}", cnt++), view);
  }

  /*!
   * \brief Returns a mutable reference to the view-backed buffer for \a key.
   *
   * \pre `hasBufferView(key)` is true.
   */
  ViewBufferType& getBufferView(std::string_view key)
  {
    verifyBufferViewKey(key);
    return m_view_buff.find(key)->second;
  }

  /*!
   * \brief Returns a const reference to the view-backed buffer for \a key.
   *
   * \pre `hasBufferView(key)` is true.
   */
  const ViewBufferType& getBufferView(std::string_view key) const
  {
    verifyBufferViewKey(key);
    return m_view_buff.find(key)->second;
  }

  /*!
   * \brief Finds a view-backed buffer by name without inserting or asserting.
   *
   * \return An optional referencing the buffer if present, else empty.
   */
  [[nodiscard]] std::optional<ViewBufferType> findBufferView(std::string_view key)
  {
    auto it = m_view_buff.find(key);
    return it != m_view_buff.end() ? std::optional<ViewBufferType>(it->second) : std::nullopt;
  }

  /*!
   * \brief Finds a view-backed buffer by name (const overload).
   *
   * \return An optional referencing the buffer if present, else empty.
   */
  [[nodiscard]] std::optional<ViewBufferType> findBufferView(std::string_view key) const
  {
    auto it = m_view_buff.find(key);
    return it != m_view_buff.end() ? std::optional<ViewBufferType>(it->second) : std::nullopt;
  }

  /// @}

  /// \name Scalars
  /// @{

  /// \brief Returns true if a scalar with the given key exists.
  [[nodiscard]] bool hasScalar(std::string_view key) const
  {
    return m_scal.find(key) != m_scal.end();
  }

  /*!
   * \brief Adds (or replaces) a scalar value with the given key.
   *
   * \param key Scalar name.
   * \param val Scalar value.
   * \return A mutable reference to the stored value.
   *
   * \note If an entry with the same key already exists, it is overwritten.
   */
  DataType& addScalar(std::string_view key, DataType val)
  {
    auto [it, _] = m_scal.insert_or_assign(KeyType(key), val);
    return it->second;
  }

  /*!
   * \brief Returns a mutable reference to the scalar for \a key.
   *
   * \pre `hasScalar(key)` is true.
   */
  DataType& getScalar(std::string_view key)
  {
    verifyScalarKey(key);
    return m_scal.find(key)->second;
  }

  /*!
   * \brief Returns a const reference to the scalar for \a key.
   *
   * \pre `hasScalar(key)` is true.
   */
  const DataType& getScalar(std::string_view key) const
  {
    verifyScalarKey(key);
    return m_scal.find(key)->second;
  }

  /*!
   * \brief Finds a scalar by name without inserting or asserting.
   *
   * \param key Scalar name.
   * \return An engaged optional with the scalar's value if present, else empty.
   */
  [[nodiscard]] std::optional<DataType> findScalar(std::string_view key) const
  {
    auto it = m_scal.find(key);
    return it != m_scal.end() ? std::optional<DataType>(it->second) : std::nullopt;
  }

  /// @}

private:
  static KeyType makeNamelessFieldKey()
  {
    static int cnt = 0;
    return axom::fmt::format("__field_{}", cnt++);
  }

  inline void verifyFieldsKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
#ifdef AXOM_DEBUG
    auto it = m_maps.find(key);
    SLIC_ASSERT_MSG(it != m_maps.end(), "Didn't find field named " << key);
    SLIC_ASSERT_MSG(std::holds_alternative<MapType>(it->second),
                    "Field named '" << key << "' is view-backed; use getFieldView()");
#endif
  }

  inline void verifyFieldViewKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
    SLIC_ASSERT_MSG(hasFieldView(key), "Didn't find view field named " << key);
  }

  inline void verifyBufferKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
    SLIC_ASSERT_MSG(hasBuffer(key), "Didn't find buffer named " << key);
  }

  inline void verifyBufferViewKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
    SLIC_ASSERT_MSG(hasBufferView(key), "Didn't find view buffer named " << key);
  }

  inline void verifyScalarKey(std::string_view AXOM_DEBUG_PARAM(key)) const
  {
    SLIC_ASSERT_MSG(hasScalar(key), "Didn't find scalar named " << key);
  }

private:
  DataVecMap m_maps;
  DataBufferMap m_buff;
  DataViewBufferMap m_view_buff;
  DataAttrMap m_scal;
};

}  // end namespace axom::slam
