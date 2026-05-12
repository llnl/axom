// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_DISPATCH_MATERIAL_FIELD_HPP_
#define AXOM_BUMP_DISPATCH_MATERIAL_FIELD_HPP_
#include "axom/bump/views/dispatch_material.hpp"

namespace axom
{
namespace bump
{
namespace views
{

/*!
 * \brief Dispatch Conduit nodes containing a unibuffer matset and a values array
 *        to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset view.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_unibuffer_field(const conduit::Node &matset,
                                       const conduit::Node &n_field,
                                       FuncType &&func)
{
  verify(matset, "matset");
  detail::verifyMixedField(n_field);
  return detail::dispatch_material_unibuffer_with_values<FuncType, MAXMATERIALS>(
    matset,
    n_field["matset_values"],
    std::forward<FuncType>(func));
}

/*!
 * \brief Dispatch Conduit nodes containing a element-dominant matset and related field
 *        to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset view.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_element_dominant_field(const conduit::Node &matset,
                                              const conduit::Node n_field,
                                              FuncType &&func)
{
  verify(matset, "matset");
  detail::verifyMixedField(n_field);
  return detail::dispatch_material_element_dominant_with_values<FuncType, MAXMATERIALS>(
    matset,
    n_field["matset_values"],
    std::forward<FuncType>(func));
}

/*!
 * \brief Dispatch Conduit nodes containing a material-dominant matset and related field
 *        to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset view.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_material_dominant_field(const conduit::Node &matset,
                                               const conduit::Node n_field,
                                               FuncType &&func)
{
  verify(matset, "matset");
  detail::verifyMixedField(n_field);
  return detail::dispatch_material_material_dominant_with_values<FuncType, MAXMATERIALS>(
    matset,
    n_field["matset_values"],
    std::forward<FuncType>(func));
}

/*!
 * \brief Dispatch Conduit nodes containing a matset and related field
 *        to a function as the appropriate type of matset view. The matset will be used
 *        to access the per-material field data.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param n_field The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset view.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_field(const conduit::Node &matset, const conduit::Node &n_field, FuncType &&func)
{
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

#endif
