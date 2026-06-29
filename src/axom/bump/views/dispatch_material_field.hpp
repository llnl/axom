// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/bump/views/dispatch_material.hpp"
#include "axom/bump/views/MixedFieldView.hpp"

namespace axom
{
namespace bump
{
namespace views
{
namespace detail
{
inline void verifyMatchingMaterialOrder(const conduit::Node &mat_values,
                                        const conduit::Node &field_values)
{
  SLIC_ERROR_IF(
    mat_values.number_of_children() != field_values.number_of_children(),
    "The matset volume_fractions and field matset_values have different numbers of materials.");

  for(conduit::index_t i = 0; i < mat_values.number_of_children(); i++)
  {
    SLIC_ERROR_IF(
      mat_values[i].name() != field_values[i].name(),
      "The matset volume_fractions and field matset_values do not have matching material order.");
  }
}

/*!
 * \brief Dispatch a unibuffer matset_values field.
 *
 * \tparam MatsetView The type of matset view associated with the field. This has
 *                    implications for the field storage.
 * \tparam FuncType The function/lambda type to call on the MixedFieldView.
 */
template <typename MatsetView, typename FuncType>
bool dispatch_unibuffer_field(const conduit::Node &n_field, FuncType &&func)
{
  bool rv = false;
  detail::verifyMixedField(n_field);
  const conduit::Node &matset_values = n_field["matset_values"];
  SLIC_ERROR_IF(!matset_values.dtype().is_number(), "The matset_values must be a number.");
  // NOTE: For now support float, double types.
  axom::bump::views::floatNodeToArrayView(matset_values, [&](auto valuesView) {
    using FieldT = typename decltype(valuesView)::value_type;
    MixedFieldView<MatsetView, FieldT> mixedFieldView;
    mixedFieldView.traits().set(valuesView);
    func(mixedFieldView);
    rv = true;
  });
  return rv;
}

/*!
 * \brief Dispatch a multibuffer matset_values field.
 *
 * \tparam MatsetView The type of matset view associated with the field. This has
 *                    implications for the field storage.
 * \tparam FuncType The function/lambda type to call on the MixedFieldView.
 */
template <typename MatsetView, typename FuncType>
bool dispatch_multibuffer_field(const conduit::Node &n_field, FuncType &&func)
{
  bool rv = false;
  detail::verifyMixedField(n_field);
  const conduit::Node &matset_values = n_field["matset_values"];
  SLIC_ERROR_IF(matset_values.number_of_children() < 1, "Missing fields in matset_values.");
  // NOTE: For now support float, double types.
  axom::bump::views::floatNodeToArrayView(matset_values[0], [&](auto firstValuesView) {
    using FieldT = typename decltype(firstValuesView)::value_type;
    MixedFieldView<MatsetView, FieldT> mixedFieldView;
    mixedFieldView.traits().add(firstValuesView);
    for(conduit::index_t i = 1; i < matset_values.number_of_children(); i++)
    {
      mixedFieldView.traits().add(axom::bump::utilities::make_array_view<FieldT>(matset_values[i]));
    }
    func(mixedFieldView);
    rv = true;
  });
  return rv;
}

}  // end namespace detail

/*!
 * \brief Dispatch Conduit nodes containing a unibuffer matset and a values array
 *        to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset and mixed field view.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset and mixed field views.
 */
template <typename FuncType, axom::IndexType MAXMATERIALS = 20>
bool dispatch_material_unibuffer_field(const conduit::Node &matset,
                                       const conduit::Node &n_field,
                                       FuncType &&func)
{
  detail::verifyPositiveMaxMaterials<MAXMATERIALS>();
  verify(matset, "matset");
  detail::verifyMixedField(n_field);

  auto handleMatset = [&](auto matsetView) {
    using MatsetView = decltype(matsetView);
    detail::dispatch_unibuffer_field<MatsetView>(n_field, [&](auto mixedFieldView) {
      func(matsetView, mixedFieldView);
    });
  };

  return dispatch_material_unibuffer<decltype(handleMatset), MAXMATERIALS>(
    matset,
    std::forward<decltype(handleMatset)>(handleMatset));
}

/*!
 * \brief Dispatch Conduit nodes containing a element-dominant matset and related field
 *        to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset and mixed field view.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset and mixed field views.
 */
template <typename FuncType, axom::IndexType MAXMATERIALS = 20>
bool dispatch_material_element_dominant_field(const conduit::Node &matset,
                                              const conduit::Node &n_field,
                                              FuncType &&func)
{
  detail::verifyPositiveMaxMaterials<MAXMATERIALS>();
  verify(matset, "matset");
  detail::verifyMixedField(n_field);
  detail::verifyMatchingMaterialOrder(matset["volume_fractions"], n_field["matset_values"]);
  auto handleMatset = [&](auto matsetView) {
    using MatsetView = decltype(matsetView);
    detail::dispatch_multibuffer_field<MatsetView>(n_field, [&](auto mixedFieldView) {
      func(matsetView, mixedFieldView);
    });
  };

  return dispatch_material_element_dominant<decltype(handleMatset), MAXMATERIALS>(
    matset,
    std::forward<decltype(handleMatset)>(handleMatset));
}

/*!
 * \brief Dispatch Conduit nodes containing a material-dominant matset and related field
 *        to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset and mixed field view.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset and mixed field views.
 */
template <typename FuncType, axom::IndexType MAXMATERIALS = 20>
bool dispatch_material_material_dominant_field(const conduit::Node &matset,
                                               const conduit::Node &n_field,
                                               FuncType &&func)
{
  detail::verifyPositiveMaxMaterials<MAXMATERIALS>();
  verify(matset, "matset");
  detail::verifyMixedField(n_field);
  detail::verifyMatchingMaterialOrder(matset["volume_fractions"], n_field["matset_values"]);
  auto handleMatset = [&](auto matsetView) {
    using MatsetView = decltype(matsetView);
    detail::dispatch_multibuffer_field<MatsetView>(n_field, [&](auto mixedFieldView) {
      func(matsetView, mixedFieldView);
    });
  };

  return dispatch_material_material_dominant<decltype(handleMatset), MAXMATERIALS>(
    matset,
    std::forward<decltype(handleMatset)>(handleMatset));
}

/*!
 * \brief Dispatch Conduit nodes containing a matset and related field
 *        to a function as the appropriate type of matset view. The matset will be used
 *        to access the per-material field data.
 *
 * \tparam FuncType The function/lambda type that will take the matset and mixed field view.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset and mixed field views.
 */
template <typename FuncType, axom::IndexType MAXMATERIALS = 20>
bool dispatch_material_field(const conduit::Node &matset, const conduit::Node &n_field, FuncType &&func)
{
  detail::verifyPositiveMaxMaterials<MAXMATERIALS>();
  bool retval =
    dispatch_material_unibuffer_field<FuncType, MAXMATERIALS>(matset,
                                                              n_field,
                                                              std::forward<FuncType>(func));
  // Multibuffer
  if(!retval)
  {
    retval =
      dispatch_material_element_dominant_field<FuncType, MAXMATERIALS>(matset,
                                                                       n_field,
                                                                       std::forward<FuncType>(func));
  }
  if(!retval)
  {
    retval = dispatch_material_material_dominant_field<FuncType, MAXMATERIALS>(
      matset,
      n_field,
      std::forward<FuncType>(func));
  }
  // NOTE: Blueprint describes a unibuffer material-dominant matset type but does not technically implement it.
  //       https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#material-sets
  if(!retval && conduit::blueprint::mesh::matset::is_uni_buffer(matset) &&
     conduit::blueprint::mesh::matset::is_material_dominant(matset))
  {
    SLIC_ERROR("Unibuffer material dominant matsets are unsupported.");
  }
  return retval;
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom
