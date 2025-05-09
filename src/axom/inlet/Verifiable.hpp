// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Verifiable.hpp
 *
 * \brief This file defines an interface for things that are verifiable
 *******************************************************************************
 */

#ifndef INLET_VERIFIABLE_HPP
#define INLET_VERIFIABLE_HPP

#include <functional>

#include "axom/inlet/inlet_utils.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class Verifiable
 *
 * \brief Interface for trivially verifiable objects - namely those that can
 * be marked as required or checked with a user-provided lambda
 * 
 * \tparam BaseType The "base" type of the object, used for the argument type
 * of a verifying predicate
 * 
 * In practice this interface is used for the Container and Function classes.
 * Currently the only supported means of verifying a composite type (container)
 * or function type are the methods exposed by this interface.
 *******************************************************************************
 */
template <typename BaseType>
class Verifiable
{
public:
  /**
   * A function which can verify the contents of the item being verifier.
   * It should report any errors via INLET_VERIFICATION_WARNING, passing
   * in the given array of errors.
   */
  using Verifier = std::function<bool(const BaseType&, std::vector<VerificationError>* errors)>;

  virtual ~Verifiable() = default;

  // Should not be reassignable
  Verifiable<BaseType>& operator=(const Verifiable<BaseType>&) = delete;
  /*!
   *****************************************************************************
   * \brief Set the required status of this object.
   *
   * Set whether this object is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether object is required
   *
   * \return Reference to calling object, for chaining
   *****************************************************************************
   */
  virtual Verifiable<BaseType>& required(bool isRequired = true) = 0;

  /*!
   *****************************************************************************
   * \brief Return the required status.
   *
   * Return that this object is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \return Boolean value of whether this object is required
   *****************************************************************************
   */
  virtual bool isRequired() const = 0;

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this object's contents
   * during the verification stage.
   * 
   * \param [in] verifier The function object.
   *****************************************************************************
  */
  Verifiable<BaseType>& registerVerifier(std::function<bool(const BaseType&)> verifier)
  {
    return registerVerifier(
      [verifier](const BaseType& item, std::vector<VerificationError>*) { return verifier(item); });
  };

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this object's contents
   * during the verification stage.
   *
   * \param [in] verifier The function which will verify the contents of
   * the container.
   *****************************************************************************
  */
  virtual Verifiable<BaseType>& registerVerifier(Verifier verifier) = 0;

  /*!
   *****************************************************************************
   * \brief Verifies the object to make sure it satisfies the imposed requirements
   * \param [in] errors An optional vector of errors to append to in the case
   * of verification failure
   * 
   * Ownership is not taken of @a errors, the raw pointer is only used for its
   * optional reference semantics, as opposed to something like
   * std::optional<std::reference_wrapper<T>>
   *****************************************************************************
  */
  virtual bool verify(std::vector<VerificationError>* errors = nullptr) const = 0;
};

/*!
 *******************************************************************************
 * \class AggregateVerifiable
 *
 * \brief Implementation of the Verifiable interface for aggregates of BaseTypes
 * 
 * \tparam BaseType The "base" type of the object, used for the argument type
 * of a verifying predicate
 * 
 * In practice this interface is used for the Container and Function classes.
 *******************************************************************************
 */
template <typename BaseType>
class AggregateVerifiable : public Verifiable<BaseType>
{
  using BaseVerifiable = Verifiable<BaseType>;

public:
  AggregateVerifiable(std::vector<std::reference_wrapper<BaseVerifiable>>&& verifiables)
    : m_verifiables(std::move(verifiables))
  { }

  AggregateVerifiable& required(bool isRequired = true) override
  {
    for(auto& verifiable : m_verifiables)
    {
      verifiable.get().required(isRequired);
    }
    return *this;
  }

  bool isRequired() const override
  {
    return std::any_of(m_verifiables.begin(),
                       m_verifiables.end(),
                       [](const BaseVerifiable& verifiable) { return verifiable.isRequired(); });
  }

  using Verifiable<BaseType>::registerVerifier;

  AggregateVerifiable& registerVerifier(typename Verifiable<BaseType>::Verifier lambda) override
  {
    for(auto& verifiable : m_verifiables)
    {
      verifiable.get().registerVerifier(lambda);
    }
    return *this;
  }

  bool verify(std::vector<VerificationError>* errors = nullptr) const override
  {
    return std::all_of(
      m_verifiables.begin(),
      m_verifiables.end(),
      [&errors](const BaseVerifiable& verifiable) { return verifiable.verify(errors); });
  }

private:
  std::vector<std::reference_wrapper<BaseVerifiable>> m_verifiables;
};

}  // namespace inlet
}  // namespace axom

#endif  // INLET_VERIFIABLE_HPP
