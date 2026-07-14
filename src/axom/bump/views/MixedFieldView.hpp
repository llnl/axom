// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/core/ArrayView.hpp"
#include "axom/core/StaticArray.hpp"
#include "axom/bump/views/MaterialView.hpp"

namespace axom
{
namespace bump
{
namespace views
{

// Base template for MixedFieldTraits. These traits are used in conjunction with the
// MixedFieldView to access MatsetView-specific field data.
template <typename MatsetView, typename FieldT>
struct MixedFieldTraits
{ };

/*!
 * \brief Specialization for UnibufferMaterialView that can store field data organized
 *        in a compatible manner.
 */
template <typename IndexT, typename FloatT, typename FieldT, axom::IndexType MAXMATERIALS>
struct MixedFieldTraits<UnibufferMaterialView<IndexT, FloatT, MAXMATERIALS>, FieldT>
{
  static_assert(MAXMATERIALS > 0, "MAXMATERIALS must be greater than 0.");
  using MatsetView = UnibufferMaterialView<IndexT, FloatT, MAXMATERIALS>;
  using ValueView = axom::ArrayView<FieldT>;

  /*!
   * \brief Set the field values view.
   * \param values The field values.
   */
  void set(ValueView values) { m_values = values; }

  /*!
   * \brief Return the field value for the iterator index.
   *
   * \param index The IteratorIndex that provides the indices to use for looking up the value.
   *
   * \return The field value at the provided index.
   */
  AXOM_HOST_DEVICE FieldT value(const IteratorIndex &index) const
  {
    SLIC_ASSERT(axom::utilities::inBounds_0_N(index.m_localIndex, m_values.size()));
    return m_values[index.m_localIndex];
  }

  ValueView m_values;
};

/*!
 * \brief Specialization for ElementDominantMaterialView that can store field data organized
 *        in a compatible manner.
 */
template <typename IndexT, typename FloatT, typename FieldT, axom::IndexType MAXMATERIALS>
struct MixedFieldTraits<ElementDominantMaterialView<IndexT, FloatT, MAXMATERIALS>, FieldT>
{
  static_assert(MAXMATERIALS > 0, "MAXMATERIALS must be greater than 0.");
  using MatsetView = ElementDominantMaterialView<IndexT, FloatT, MAXMATERIALS>;
  using ValueView = axom::ArrayView<FieldT>;

  /*!
   * \brief Add the values to the list of buffers.
   * \param values The values to be added.
   */
  void add(ValueView values) { m_values.push_back(values); }

  /*!
   * \brief Return the field value for the iterator index.
   *
   * \param index The IteratorIndex that provides the indices to use for looking up the value.
   *
   * \return The field value at the provided index.
   */
  AXOM_HOST_DEVICE FieldT value(const IteratorIndex &index) const
  {
    SLIC_ASSERT(axom::utilities::inBounds_0_N(index.m_bufferIndex, m_values.size()));
    SLIC_ASSERT(
      axom::utilities::inBounds_0_N(index.m_localIndex, m_values[index.m_bufferIndex].size()));
    return m_values[index.m_bufferIndex][index.m_localIndex];
  }

  axom::StaticArray<ValueView, MAXMATERIALS> m_values;
};

/*!
 * \brief Specialization for MaterialDominantMaterialView that can store field data organized
 *        in a compatible manner.
 */
template <typename IndexT, typename FloatT, typename FieldT, axom::IndexType MAXMATERIALS>
struct MixedFieldTraits<MaterialDominantMaterialView<IndexT, FloatT, MAXMATERIALS>, FieldT>
{
  static_assert(MAXMATERIALS > 0, "MAXMATERIALS must be greater than 0.");
  using MatsetView = MaterialDominantMaterialView<IndexT, FloatT, MAXMATERIALS>;
  using ValueView = axom::ArrayView<FieldT>;

  /*!
   * \brief Add the values to the list of buffers.
   * \param values The values to be added.
   */
  void add(ValueView values) { m_values.push_back(values); }

  /*!
   * \brief Return the field value for the iterator index.
   *
   * \param index The IteratorIndex that provides the indices to use for looking up the value.
   *
   * \return The field value at the provided index.
   */
  AXOM_HOST_DEVICE FieldT value(const IteratorIndex &index) const
  {
    SLIC_ASSERT(axom::utilities::inBounds_0_N(index.m_bufferIndex, m_values.size()));
    SLIC_ASSERT(
      axom::utilities::inBounds_0_N(index.m_localIndex, m_values[index.m_bufferIndex].size()));
    return m_values[index.m_bufferIndex][index.m_localIndex];
  }

  axom::StaticArray<ValueView, MAXMATERIALS> m_values;
};

/*!
 * \param This view enables the user to traverse Blueprint mixed field data using
 *        data from a MatsetView const_iterator.
 */
template <typename MatsetView, typename FieldT>
class MixedFieldView
{
public:
  using Traits = MixedFieldTraits<MatsetView, FieldT>;

  Traits &traits() { return m_traits; }
  const Traits &traits() const { return m_traits; }

  /*!
   * \brief Given a MatsetView's const_iterator, use it to look up the typed field
   *        data in the field.
   */
  AXOM_HOST_DEVICE FieldT value(const typename MatsetView::const_iterator &it) const
  {
    return m_traits.value(it.index());
  }

private:
  Traits m_traits;
};

}  // end namespace views
}  // end namespace bump
}  // end namespace axom
