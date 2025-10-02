// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Function.hpp
 *
 * \brief This file contains the class definition of Inlet's Function class.
 *******************************************************************************
 */

#ifndef INLET_FUNCTION_HPP
#define INLET_FUNCTION_HPP

#include <memory>
#include <functional>
#include <tuple>
#include <type_traits>
#include <typeinfo>

#include "axom/fmt.hpp"

#include "axom/sidre.hpp"

#include "axom/inlet/Field.hpp"
#include "axom/inlet/InletVector.hpp"
#include "axom/inlet/Verifiable.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \brief The tags used to describe function signatures in the input file
 * 
 * \note Additions to this enumeration should be propagated to FunctionType, LuaReader
 * and func_signature_lists defined below (the mapping from enum to types)
 * 
 * \note Vector corresponds to a vector with max dimension of three, Double corresponds to
 * a floating-point scalar
 * 
 * \note A two-dimensional vector was intentionally excluded for simplicity as
 * a 3D vector can be used in its place (with the third component being empty/zero)
 *******************************************************************************
 */
enum class FunctionTag
{
  Vector,
  Double,
  Void,
  String
};

/*!
 *******************************************************************************
 * \brief The types used to describe function signatures in the input file
 * 
 * \note These are aliases intended to improve readability
 *******************************************************************************
 */
struct FunctionType
{
  using Vector = InletVector;
  using Double = double;
  using Void = void;
  using String = std::string;
};

namespace detail
{
/*!
 *******************************************************************************
 * \class inlet_function_arg_type
 *
 * \brief A type trait for modifying function argument types to enforce const
 * correctness and to avoid copies
 * \tparam Arg The argument type
 * 
 * Maps a type "Arg" to "const Arg&" if that type does not satisfy std::is_arithmetic,
 * which matches Inlet primitive types bool, int, and double, but not vectors or
 * std::strings - though only vectors and doubles are currently supported as arg
 * types
 *******************************************************************************
 */
template <typename Arg>
struct inlet_function_arg_type
{
  using type = typename std::conditional<
    std::is_arithmetic<Arg>::value,
    Arg,
    typename std::add_lvalue_reference<typename std::add_const<Arg>::type>::type>::type;
};

/*!
 *******************************************************************************
 * \class cleanup_function_signature
 *
 * \brief Takes a Ret(Args...) signature and applies cvref qualifiers to applicable
 * arguments
 * 
 * \tparam Ret The function's return type
 * \tparam Args... The function's arguments
 * 
 * Designed to be used with user-facing retrieval functions so the user is not 
 * required to specify the qualifiers in something like a get<std::function<Ret(Args...)>
 *******************************************************************************
 */
template <typename>
struct cleanup_function_signature;

template <typename Ret, typename... Args>
struct cleanup_function_signature<Ret(Args...)>
{
  using type = Ret(typename inlet_function_arg_type<Args>::type...);
};

static constexpr std::size_t MAX_NUM_ARGS = 2u;

// Represents a buffer of bytes aligned like an std::function.
struct FunctionBuffer
{
  static constexpr size_t Alignment = alignof(std::function<void()>);
  static constexpr size_t Size = sizeof(std::function<void()>);

  alignas(Alignment) std::uint8_t m_bytes[Size];
};

template <typename Func>
inline void destroy_func_inst(FunctionBuffer* function_storage)
{
  Func* function = reinterpret_cast<Func*>(function_storage);
  // Destroy underlying function object
  function->~Func();
  // Destroy allocated buffer
  delete function_storage;
}

template <>
inline void destroy_func_inst<void>(FunctionBuffer*)
{ }

}  // end namespace detail

/*!
 *******************************************************************************
 * \class FunctionWrapper
 *
 * \brief A sum type for callables with arbitrary signature
 * \tparam FunctionTags The types of supported functions
 * 
 * This provides an interface not templated on a specific function signature
 * for uniform retrieval through the Reader interface
 *******************************************************************************
 */
class FunctionWrapper
{
private:
  using StorageType = std::unique_ptr<detail::FunctionBuffer, void (*)(detail::FunctionBuffer*)>;

public:
  /*!
   *******************************************************************************
   * \brief Primary constructor, initializes a member of the "variant"
   * 
   * \param [in] func The function to initialize with
   * \tparam FuncType The function's signature
   * 
   * \note "Empty" functions are allowable and result in the object comparing
   * false when converted to bool
   *******************************************************************************
   */
  template <typename FuncType>
  FunctionWrapper(std::function<FuncType>&& func)
  {
    static_assert(alignof(std::function<FuncType>) == detail::FunctionBuffer::Alignment,
                  "std::function does not have the same alignment across all instances.");
    static_assert(sizeof(std::function<FuncType>) == detail::FunctionBuffer::Size,
                  "std::function does not have the same size across all instances.");

    m_function_valid = static_cast<bool>(func);

    // Use unique_ptr to hold the destructor for this function type
    m_func =
      StorageType {new detail::FunctionBuffer, &detail::destroy_func_inst<std::function<FuncType>>};
    // Construct function object in-place in byte storage
    new(m_func.get()) std::function<FuncType>(std::move(func));
    // Store the type information of the passed-in std::function
    m_func_type = typeid(std::function<FuncType>);
  }

  FunctionWrapper() = default;

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
    using FuncType = std::function<Ret(typename detail::inlet_function_arg_type<Args>::type...)>;
    SLIC_ERROR_IF(typeid(FuncType) != m_func_type.get(),
                  fmt::format("[Inlet] Attempted to call function '{0}' with incorrect type.\n"
                              " - Stored type: {1}\n"
                              " - Expected type: {2}\n",
                              m_name,
                              m_func_type.get().name(),
                              typeid(FuncType).name()));

    FuncType* ptr = reinterpret_cast<FuncType*>(m_func.get());
    SLIC_ERROR_IF(
      !m_function_valid || !(*ptr),
      fmt::format("[Inlet] No valid function '{0} assigned to function wrapper.", m_name));

    const auto& func = *ptr;
    return func(std::forward<typename detail::inlet_function_arg_type<Args>::type>(args)...);
  }

  template <typename FuncType>
  std::function<FuncType> get() const
  {
    using StoredFuncType = std::function<typename detail::cleanup_function_signature<FuncType>::type>;
    SLIC_ERROR_IF(typeid(StoredFuncType) != m_func_type.get(),
                  fmt::format("[Inlet] Attempted to get function '{0}' with incorrect type.\n"
                              " - Stored type: {1}\n"
                              " - Expected type: {2}\n",
                              m_name,
                              m_func_type.get().name(),
                              typeid(StoredFuncType).name()));

    StoredFuncType* ptr = reinterpret_cast<StoredFuncType*>(m_func.get());

    return *ptr;
  }

  /*!
   *******************************************************************************
   * \brief Checks whether the function exists
   *******************************************************************************
   */
  explicit operator bool() const { return m_function_valid; }

  /*!
   *******************************************************************************
   * \brief Sets the function's name
   * 
   * \note Needs to be separate from constructor to allow the compiler to deduce
   * template arguments correctly
   *******************************************************************************
   */
  void setName(std::string&& name) { m_name = std::move(name); }

private:
  StorageType m_func {nullptr, &detail::destroy_func_inst<void>};
  std::reference_wrapper<const std::type_info> m_func_type {typeid(void)};

  bool m_function_valid = false;
  std::string m_name;
};

using FunctionVariant = FunctionWrapper;

class Function : public Verifiable<Function>
{
public:
  Function(axom::sidre::Group* sidreGroup, axom::sidre::Group* root, FunctionVariant&& func)
    : m_sidreGroup(sidreGroup)
    , m_sidreRootGroup(root)
    , m_func(std::move(func))
  {
    m_func.setName(name());
  }

  /*!
   *****************************************************************************
   * \brief Returns a function of requested type
   * 
   * \return The value
   * \tparam T The type to retrieve
   *****************************************************************************
   */
  template <typename FuncType>
  std::function<FuncType> get() const
  {
    return m_func.get<FuncType>();
  }

  template <typename Ret, typename... Args>
  Ret call(Args&&... args) const
  {
    return m_func.call<Ret>(std::forward<Args>(args)...);
  }

  /*!
   *****************************************************************************
   * \brief Returns pointer to the Sidre Group class for this Function.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for this Function instance.
   *****************************************************************************
   */
  const axom::sidre::Group* sidreGroup() const { return m_sidreGroup; };

  /*!
   *******************************************************************************
   * \brief Checks whether the function exists
   *******************************************************************************
   */
  operator bool() const { return static_cast<bool>(m_func); }

  /*!
   *****************************************************************************
   * \return The full name for this Function.
   *****************************************************************************
  */
  std::string name() const;

  bool verify(std::vector<VerificationError>* errors = nullptr) const override;

  Function& required(bool isRequired = true) override;

  bool isRequired() const override;

  using Verifiable<Function>::registerVerifier;

  Function& registerVerifier(Verifier lambda) override;

private:
  // This function's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  Verifier m_verifier;
  FunctionVariant m_func;
};

}  // end namespace inlet
}  // end namespace axom

#endif  // INLET_FUNCTION_HPP
