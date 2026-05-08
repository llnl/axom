// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_DISPATCH_MATERIAL_HPP_
#define AXOM_BUMP_DISPATCH_MATERIAL_HPP_

#include "axom/bump/views/MaterialView.hpp"
#include "axom/bump/views/NodeArrayView.hpp"
#include "axom/bump/views/dispatch_utilities.hpp"

#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace bump
{
namespace views
{
namespace detail
{

inline void verifyMixedField(const conduit::Node &n_field)
{
  SLIC_ERROR_IF(!n_field.has_path("matset_values"),
    "The mixed field does not contain matset_values");
}

/*!
 * \brief Dispatch Conduit nodes containing a unibuffer matset and a values array
 *        to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param values The node that contains the values to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_unibuffer_with_values(const conduit::Node &matset, const conduit::Node &values, FuncType &&func)
{
  bool retval = false;
  verify(matset, "matset");
  if(conduit::blueprint::mesh::matset::is_uni_buffer(matset))
  {
    indexNodeToArrayViewSame(
      matset["material_ids"],
      matset["sizes"],
      matset["offsets"],
      matset["indices"],
      [&](auto material_ids, auto sizes, auto offsets, auto indices) {
        floatNodeToArrayView(values, [&](auto typedValues) {
          using IndexType = typename decltype(material_ids)::value_type;
          using FloatType = typename decltype(typedValues)::value_type;

          UnibufferMaterialView<IndexType, FloatType, MAXMATERIALS> matsetView;
          matsetView.set(material_ids, typedValues, sizes, offsets, indices);
          func(matsetView);
        });
      });
    retval = true;
  }
  return retval;
}

/*!
 * \brief Use the material_map, if it exists, to get the material id of the named material.
 *
 * \param matset The node that contains the matset.
 * \param matName The name of the material.
 * \param defaultValue The default matno to use if the material_map is not present.
 *
 * \return The material id for the named material.
 */
template <typename IntElement>
IntElement getMaterialID(const conduit::Node &matset,
                         const std::string &matName,
                         IntElement defaultValue)
{
  IntElement matno = static_cast<IntElement>(defaultValue);
  if(matset.has_child("material_map"))
  {
    const conduit::Node &n_mm = matset["material_map"];
    if(n_mm.has_child(matName))
    {
      matno = static_cast<IntElement>(n_mm[matName].to_int());
    }
  }
  return matno;
}

/*!
 * \brief Dispatch a Conduit node containing a multibuffer matset to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param values_object The node that contains the values to use as volume fractions / field values and indices.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_multibuffer_with_values(const conduit::Node &matset, const conduit::Node &values_object, FuncType &&func)
{
  bool retval = false;
  verify(matset, "matset");
  if(conduit::blueprint::mesh::matset::is_multi_buffer(matset))
  {
    if(values_object.number_of_children() > 0)
    {
      const conduit::Node &n_firstValues = values_object[0].fetch_existing("values");
      const conduit::Node &n_firstIndices = values_object[0].fetch_existing("indices");
      indexNodeToArrayView(n_firstIndices, [&](auto firstIndices) {
        floatNodeToArrayView(n_firstValues, [&](auto firstValues) {
          using IntElement =
            typename std::remove_const<typename decltype(firstIndices)::value_type>::type;
          using FloatElement =
            typename std::remove_const<typename decltype(firstValues)::value_type>::type;
          using IntView = axom::ArrayView<IntElement>;
          using FloatView = axom::ArrayView<FloatElement>;

          MultiBufferMaterialView<IntElement, FloatElement, MAXMATERIALS> matsetView;

          for(conduit::index_t i = 0; i < values_object.number_of_children(); i++)
          {
            const conduit::Node &values = values_object[i].fetch_existing("values");
            const conduit::Node &indices = values_object[i].fetch_existing("indices");

            const IntElement *indices_ptr = indices.value();
            const FloatElement *values_ptr = values.value();

            IntView indices_view(const_cast<IntElement *>(indices_ptr),
                                 indices.dtype().number_of_elements());
            FloatView values_view(const_cast<FloatElement *>(values_ptr),
                                  values.dtype().number_of_elements());

            // Get the material number if we can.
            IntElement matno = getMaterialID<IntElement>(matset,
                                                         values_object[i].name(),
                                                         static_cast<IntElement>(i));

            matsetView.add(matno, indices_view, values_view);
          }

          func(matsetView);
        });
      });
      retval = true;
    }
  }
  return retval;
}

/*!
 * \brief Dispatch a Conduit node containing a element-dominant matset to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param values_object The node that contains the values object to be used as volume fractions / field.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_element_dominant_with_values(const conduit::Node &matset, const conduit::Node &values_object, FuncType &&func)
{
  bool retval = false;
  verify(matset, "matset");
  if(conduit::blueprint::mesh::matset::is_element_dominant(matset))
  {
    if(values_object.number_of_children() > 0)
    {
      const conduit::Node &n_firstValues = values_object[0];
      floatNodeToArrayView(n_firstValues, [&](auto firstValues) {
        using FloatElement =
          typename std::remove_const<typename decltype(firstValues)::value_type>::type;
        using FloatView = axom::ArrayView<FloatElement>;
        using IntElement = axom::IndexType;

        ElementDominantMaterialView<IntElement, FloatElement, MAXMATERIALS> matsetView;

        for(conduit::index_t i = 0; i < values_object.number_of_children(); i++)
        {
          const conduit::Node &values = values_object[i];
          const FloatElement *values_ptr = values.value();
          FloatView values_view(const_cast<FloatElement *>(values_ptr),
                                values.dtype().number_of_elements());

          // Get the material number if we can.
          IntElement matno =
            getMaterialID<IntElement>(matset, values_object[i].name(), static_cast<IntElement>(i));

          matsetView.add(matno, values_view);
        }

        func(matsetView);
      });
      retval = true;
    }
  }
  return retval;
}

/*!
 * \brief Dispatch a Conduit node containing a material-dominant matset to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param values_object The node that contains the values to use as volume fractions / field values and indices.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_material_dominant_with_values(const conduit::Node &matset, const conduit::Node &values_object, FuncType &&func)
{
  bool retval = false;
  verify(matset, "matset");
  if(conduit::blueprint::mesh::matset::is_material_dominant(matset))
  {
    const conduit::Node &element_ids = matset.fetch_existing("element_ids");
    if(values_object.number_of_children() > 0 &&
       values_object.number_of_children() == element_ids.number_of_children())
    {
      const conduit::Node &n_firstValues = values_object[0];
      const conduit::Node &n_firstIndices = element_ids[0];

      indexNodeToArrayView(n_firstIndices, [&](auto firstIndices) {
        floatNodeToArrayView(n_firstValues, [&](auto firstValues) {
          using FloatElement =
            typename std::remove_const<typename decltype(firstValues)::value_type>::type;
          using IntElement =
            typename std::remove_const<typename decltype(firstIndices)::value_type>::type;
          using FloatView = axom::ArrayView<FloatElement>;
          using IntView = axom::ArrayView<IntElement>;

          MaterialDominantMaterialView<IntElement, FloatElement, MAXMATERIALS> matsetView;

          for(conduit::index_t i = 0; i < values_object.number_of_children(); i++)
          {
            const conduit::Node &indices = element_ids[i];
            const conduit::Node &values = values_object[i];

            const IntElement *indices_ptr = indices.value();
            const FloatElement *values_ptr = values.value();

            IntView indices_view(const_cast<IntElement *>(indices_ptr),
                                 indices.dtype().number_of_elements());
            FloatView values_view(const_cast<FloatElement *>(values_ptr),
                                  values.dtype().number_of_elements());

            // Get the material number if we can.
            IntElement matno =
              getMaterialID<IntElement>(matset, values.name(), static_cast<IntElement>(i));

            matsetView.add(matno, indices_view, values_view);
          }

          func(matsetView);
        });
      });
      retval = true;
    }
  }
  return retval;
}

} // end namespace detail

/*!
 * \brief Make a unibuffer matset view from a Conduit node.
 */
template <typename IntType, typename FloatType, size_t MAXMATERIALS = 20>
struct make_unibuffer_matset
{
  using MatsetView = UnibufferMaterialView<IntType, FloatType, MAXMATERIALS>;

  /*!
   * \brief Wrap the Conduit node as a unibuffer matset view.
   *
   * \param n_matset The Conduit node that contains the matset.
   *
   * \return A UnibufferMaterialView.
   */
  static MatsetView view(const conduit::Node &n_matset)
  {
    namespace utils = axom::bump::utilities;
    verify(n_matset, "matset");
    MatsetView m;
    m.set(utils::make_array_view<IntType>(n_matset["material_ids"]),
          utils::make_array_view<FloatType>(n_matset["volume_fractions"]),
          utils::make_array_view<IntType>(n_matset["sizes"]),
          utils::make_array_view<IntType>(n_matset["offsets"]),
          utils::make_array_view<IntType>(n_matset["indices"]));
    return m;
  }

  /*!
   * \brief Wrap the Conduit matset and field nodes as a unibuffer matset view
   *        so we can traverse the field using material machinery. The field
   *        data will be accessible as the volume component in the matset view.
   *
   * \param n_matset The Conduit node that contains the matset.
   * \param n_field The Conduit node that contains the mixed field; its format
   *                must match that of the matset.
   *
   * \return A UnibufferMaterialView.
   */
  static MatsetView mixedFieldView(const conduit::Node &n_matset,
                                   const conduit::Node &n_field)
  {
    namespace utils = axom::bump::utilities;
    verify(n_matset, "matset");
    detail::verifyMixedField(n_field);
    // NOTE: further field length and type checking happens in the MatsetView.
    MatsetView m;
    m.set(utils::make_array_view<IntType>(n_matset["material_ids"]),
          utils::make_array_view<FloatType>(n_field["matset_values"]),
          utils::make_array_view<IntType>(n_matset["sizes"]),
          utils::make_array_view<IntType>(n_matset["offsets"]),
          utils::make_array_view<IntType>(n_matset["indices"]));
    return m;
  }
};

/*!
 * \brief Dispatch a Conduit node containing a unibuffer matset to a function as
 *        the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_unibuffer(const conduit::Node &matset, FuncType &&func)
{
  verify(matset, "matset");
  return detail::dispatch_material_unibuffer_with_values<FuncType, MAXMATERIALS>(matset,
           matset["volume_fractions"], std::forward<FuncType>(func));
}

/*!
 * \brief Dispatch a Conduit node containing a multibuffer matset to a function as
 *        the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param values_object The node that contains the values to use as volume fractions / field values and indices.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_multibuffer(const conduit::Node &matset, FuncType &&func)
{
  verify(matset, "matset");
  return detail::dispatch_material_multibuffer_with_values<FuncType, MAXMATERIALS>(matset, matset["volume_fractions"], std::forward<FuncType>(func));
}

/*!
 * \brief Dispatch a Conduit node containing a element-dominant matset to a function as
 *        the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_element_dominant(const conduit::Node &matset, FuncType &&func)
{
  verify(matset, "matset");
  return detail::dispatch_material_element_dominant_with_values(matset, matset["volume_fractions"], std::forward<FuncType>(func));
}

/*!
 * \brief Dispatch a Conduit node containing a material-dominant matset to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material_material_dominant(const conduit::Node &matset, FuncType &&func)
{
  verify(matset, "matset");
  return detail::dispatch_material_material_dominant_with_values<FuncType, MAXMATERIALS>(matset, matset["volume_fractions"], std::forward<FuncType>(func));
}

/*!
 * \brief Dispatch a Conduit node containing a matset to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param func   The function/lambda that will operate on the matset view.
 *
 * \return true if the dispatch worked, false otherwise.
 */
template <typename FuncType, size_t MAXMATERIALS = 20>
bool dispatch_material(const conduit::Node &matset, FuncType &&func)
{
  bool retval =
    dispatch_material_unibuffer<FuncType, MAXMATERIALS>(matset, std::forward<FuncType>(func));
  if(!retval)
  {
    retval =
      dispatch_material_multibuffer<FuncType, MAXMATERIALS>(matset, std::forward<FuncType>(func));
  }
  if(!retval)
  {
    retval = dispatch_material_element_dominant<FuncType, MAXMATERIALS>(matset, std::forward<FuncType>(func));
  }
  if(!retval)
  {
    retval = dispatch_material_material_dominant<FuncType, MAXMATERIALS>(matset, std::forward<FuncType>(func));
  }
  return retval;
}

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
