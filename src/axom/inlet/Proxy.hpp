// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 *******************************************************************************
 * \file Proxy.hpp
 *
 * \brief This file contains the class definition of Inlet's Proxy class.
 *******************************************************************************
 */

#include <type_traits>
#include <utility>

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Container.hpp"

namespace axom
{
namespace inlet
{

namespace detail
{
template <typename T, typename SFINAE = void>
struct has_ProxyFromInlet_specialization : std::false_type
{ };

template <typename T>
struct has_ProxyFromInlet_specialization<
  T,
  typename std::enable_if<
    std::is_same<T, decltype(std::declval<FromInlet<T>&>()(std::declval<const Proxy&>()))>::value>::type>
  : std::true_type
{ };

/*!
 *******************************************************************************
 * \brief Maps a Proxy::evaluate() argument to the type passed to an Inlet function
 *
 * Inlet converts numeric arguments to double and character strings to std::string.
 * This lets callers pass an int or a string literal to functions whose argument
 * types are FunctionTag::Double or FunctionTag::String.
 *******************************************************************************
 */
template <typename Arg>
using function_value_arg_t = std::conditional_t<
  std::is_arithmetic_v<std::decay_t<Arg>>,
  double,
  std::conditional_t<std::is_convertible_v<std::decay_t<Arg>, std::string>, std::string, std::decay_t<Arg>>>;
}  // namespace detail

/*!
 *******************************************************************************
 * \class Proxy
 *
 * \brief Provides a uniform interface for access and conversion to primitive
 * and user-defined types
 *
 * \see Inlet Field
 * \see Inlet Container
 *******************************************************************************
 */
class Proxy
{
public:
  Proxy() = default;
  /*!
   *******************************************************************************
   * \brief Constructs a proxy view onto a container
   *
   * \param [in] container The container to construct a proxy into
   *******************************************************************************
   */
  Proxy(const Container& container) : m_container(&container) { }

  /*!
   *******************************************************************************
   * \brief Constructs a proxy view onto a field
   *
   * \param [in] field The field to construct a proxy into
   *******************************************************************************
   */
  Proxy(const Field& field) : m_field(&field) { }

  /*!
   *******************************************************************************
   * \brief Constructs a proxy view onto a function
   * 
   * \param [in] func The function to construct a proxy into
   *******************************************************************************
   */
  Proxy(const Function& func) : m_func(&func) { }

  /*!
   *******************************************************************************
   * \brief Returns an object from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \note Implicit conversions are defined only for primitive types - to 
   * convert to a user-defined type, use an explicit conversion with static_cast
   * or T{}
   *******************************************************************************
   */
  template <typename T,
            typename SFINAE = typename std::enable_if<
              detail::is_inlet_primitive<T>::value || detail::is_inlet_primitive_array<T>::value ||
              detail::is_inlet_primitive_dict<T>::value || detail::is_std_function<T>::value ||
              detail::is_primitive_std_vector<T>::value>::type>
  operator T() const
  {
    return get<T>();
  }

  /*!
   *****************************************************************************
   * \brief Return whether a subobject with the given name is present in 
   * the container referred to by the calling proxy.
   *
   * \return Boolean value indicating whether this Container's subtree contains a
   * Field or Container with the given name.
   *****************************************************************************
   */
  bool contains(const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief Returns the type of the stored value
   * 
   * \return The type
   * \see InletType
   *****************************************************************************
   */
  InletType type() const;

  /*!
   *******************************************************************************
   * \brief Obtains a proxy view into the proxy for either a Field/Container subobject
   * 
   * Returns a reference via a lightweight proxy object to the element in the 
   * datastore at the index specified by the name.  This can be a field 
   * or a container.
   * 
   * \param [in] name The name of the subobject
   * \return A view onto the subobject
   *******************************************************************************
   */
  Proxy operator[](const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief Returns pointer to the Sidre Group for the underlying object
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for the underlying object.
   *****************************************************************************
   */
  const axom::sidre::Group* sidreGroup() const;

  /*!
   *****************************************************************************
   * \return The full name for the underlying object.
   *****************************************************************************
  */
  std::string name() const;

  /*!
   *******************************************************************************
   * \brief Returns a user-defined type from the proxy
   *
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \pre The Proxy must refer to a container object
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<!detail::is_inlet_primitive<T>::value && !detail::is_std_function<T>::value &&
                            detail::has_ProxyFromInlet_specialization<T>::value,
                          T>::type
  get() const
  {
    checkNotFunction();
    SLIC_ASSERT_MSG(m_container != nullptr,
                    "[Inlet] Tried to read a user-defined type from a Proxy "
                    "containing a single field or function");
    FromInlet<T> from_inlet;
    return from_inlet(*this);
  }

  /*!
   *******************************************************************************
   * \brief Returns a user-defined type from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \pre The Proxy must refer to a container object
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<!detail::is_inlet_primitive<T>::value && !detail::is_std_function<T>::value &&
                            !detail::has_ProxyFromInlet_specialization<T>::value,
                          T>::type
  get() const
  {
    checkNotFunction();
    SLIC_ASSERT_MSG(m_container != nullptr,
                    "[Inlet] Tried to read a user-defined type from a Proxy "
                    "containing a single field or function");
    return m_container->get<T>();
  }

  /*!
   *******************************************************************************
   * \brief Returns a function type from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \pre The Proxy must refer to a function object
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_std_function<T>::value, T>::type get() const
  {
    SLIC_ASSERT_MSG(m_func != nullptr,
                    "[Inlet] Tried to read a function from a Proxy "
                    "containing a field or container");
    return m_func->get<typename detail::std_function_signature<T>::type>();
  }

  /*!
   *******************************************************************************
   * \brief Calls the function
   * 
   * \param [in] args The parameter pack for the function's arguments
   * \tparam Ret The user-specified return type, needed to fully disambiguate the
   * function to call
   * \tparam Args The types of the user-specified arguments, deduced automatically
   * 
   * \return The function's result
   *******************************************************************************
   */
  template <typename Ret, typename... Args>
  Ret call(Args&&... args) const
  {
    SLIC_ASSERT_MSG(m_func != nullptr,
                    "[Inlet] Tried to call a Proxy "
                    "containing a field or container");
    return m_func->call<Ret>(std::forward<Args>(args)...);
  }

  /*!
   *******************************************************************************
   * \brief Evaluates an entry declared with orFunction().
   *
   * If the input supplied a function, this calls it with \a args.
   * If the input supplied a value, this returns the value and ignores \a args.
   *
   * \param [in] args Arguments for the function. Numbers are passed as double and
   * character strings as std::string, matching FunctionTag::Double and FunctionTag::String.
   * \tparam T The value type: a primitive type, or std::vector<double> for an array
   * \tparam Args The types of the arguments, deduced automatically
   *
   * \return The function's result or the value
   *
   * \pre The Proxy refers to a Field or an array, or to a function returning \a T
   * whose arguments match \a args. A function for an array of doubles returns
   * FunctionTag::Vector, whose components are returned as a std::vector<double>.
   *******************************************************************************
   */
  template <typename T, typename... Args>
  T evaluate(Args&&... args) const
  {
    constexpr bool is_array = std::is_same_v<T, std::vector<double>>;
    static_assert(detail::is_inlet_primitive<T>::value || is_array,
                  "evaluate() returns a primitive type (bool, int, double, std::string) "
                  "or std::vector<double>");
    if(m_func != nullptr)
    {
      if constexpr(is_array)
      {
        const auto result = m_func->call<FunctionType::Vector>(
          detail::function_value_arg_t<Args>(std::forward<Args>(args))...);
        return T(result.vec.data(), result.vec.data() + result.dim);
      }
      else
      {
        return m_func->call<T>(detail::function_value_arg_t<Args>(std::forward<Args>(args))...);
      }
    }
    return get<T>();
  }

  /*!
   *******************************************************************************
   * \brief Returns a primitive type from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \pre The Proxy must refer to a field object
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_inlet_primitive<T>::value, T>::type get() const
  {
    checkNotFunction();
    SLIC_ASSERT_MSG(m_field != nullptr,
                    "[Inlet] Tried to read a primitive type from a Proxy "
                    "containing a container or function");
    return m_field->get<T>();
  }

private:
  /*!
   *******************************************************************************
   * \brief Reports an error if this Proxy refers to a function, which can't be read as a value
   *******************************************************************************
   */
  void checkNotFunction() const
  {
    SLIC_ERROR_IF(m_func != nullptr,
                  fmt::format("[Inlet] '{0}' was supplied as a function. "
                              "Use evaluate() or check type() before reading it as a value",
                              detail::userFacingPath(m_func->name())));
  }

  const Container* const m_container = nullptr;
  const Field* const m_field = nullptr;
  const Function* const m_func = nullptr;
};

}  // end namespace inlet
}  // end namespace axom
