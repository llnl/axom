// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file VerifiableScalar.hpp
 *
 * \brief This file defines an interface for scalars that are verifiable
 *******************************************************************************
 */

#ifndef INLET_VERIFIABLE_SCALAR_HPP
#define INLET_VERIFIABLE_SCALAR_HPP

#include <memory>
#include <functional>
#include <vector>
#include <initializer_list>

#include "axom/inlet/inlet_utils.hpp"

namespace axom
{
namespace inlet
{
// Forward declarations
class Field;

/*!
 *******************************************************************************
 * \class VerifiableScalar
 *
 * \brief Basic interface for verifiable scalar values of Inlet primitive type,
 * namely int, double, bool, or std::string - implementations can use this
 * directly (inlet::Field) or forward to all elements of a collection (inlet::AggregrateField)
 * 
 * In practice this interface is used for the Field and AggregateField types.
 * In addition to the ability to mark things as required and use user-defined lambdas
 * (as in Verifiable), default values can be provided, as can ranges of valid
 * values and discrete sets of valid values.
 *******************************************************************************
 */
class VerifiableScalar
{
public:
  /**
   * A function which can verify the contents of the item being verified.
   * It should report any errors via the INLET_VERIFICATION_WARNING macro,
   * passing in the given array of errors.
   */
  using Verifier =
    std::function<bool(const axom::inlet::Field&, std::vector<VerificationError>* errors)>;

  virtual ~VerifiableScalar() = default;

  // Should not be reassignable
  VerifiableScalar& operator=(const VerifiableScalar&) = delete;

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
  virtual VerifiableScalar& required(bool isRequired = true) = 0;

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
   * \param [in] lambda The function object.
   *****************************************************************************
  */
  VerifiableScalar& registerVerifier(std::function<bool(const axom::inlet::Field&)> lambda);

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this object's contents
   * during the verification stage.
   *
   * \param [in] verifier The function object.
   *****************************************************************************
  */
  virtual VerifiableScalar& registerVerifier(Verifier verifier) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(const std::string& value) = 0;
  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(const char* value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default boolean value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(bool value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default integer value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(int value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default double value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(double value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the range of this object.
   *
   * Set the continuous range for the object in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& range(double startVal, double endVal) = 0;

  /*!
   *****************************************************************************
   * \brief Set the range of this object.
   *
   * Set the continuous range for the object in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& range(int startVal, int endVal) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An vector containing the set of allowed integer values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::vector<int>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An vector containing the set of allowed double values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::vector<double>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set A vector containing the set of allowed string values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::vector<std::string>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An initializer list containing the set of allowed C-string 
   * values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::initializer_list<const char*>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An initializer list containing the valid integer values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::initializer_list<int>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values.
   *
   * \param [in] set An initializer list containing the valid double values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::initializer_list<double>& set) = 0;

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

}  // namespace inlet
}  // namespace axom

#endif  // INLET_VERIFIABLE_SCALAR_HPP
