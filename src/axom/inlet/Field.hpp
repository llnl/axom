// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 *******************************************************************************
 * \file Field.hpp
 *
 * \brief This file contains the class definition of Inlet's Field and AggregateField classes.
 *******************************************************************************
 */

#include "axom/sidre.hpp"
#include "axom/inlet/VariantKey.hpp"
#include "axom/inlet/VerifiableScalar.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <functional>
#include <vector>

namespace axom
{
namespace inlet
{
// Forward declarations
class Function;
class Reader;
enum class FunctionTag;

/*!
 *******************************************************************************
 * \class Field
 *
 * \brief Provides functions to help define how individual field variables in an
 * input file are expected to behave.  It also holds the Sidre Group to 
 * the individual field.
 *
 * \see Inlet Container
 *******************************************************************************
 */
class Field : public VerifiableScalar
{
public:
  /*!
   *****************************************************************************
   * \brief Constructor for the Field class.
   *
   * This class provides functions to define the behavior of the Field
   * data already read and stored in the given Sidre Group.
   *
   * \param [in] sidreGroup Pointer to the already created Sidre Group.
   * \param [in] root Pointer to the sidreRootGroup containing this Field
   * \param [in] type FieldType specifying the data type of this Field instance.
   * Default is FieldType::UNSPECIFIED.
   * \param [in] docEnabled Boolean indicating whether or not documentation
   * generation is enabled for Input file this Field instance belongs to.
   * \param [in] reader The Reader for the input file, used by orFunction().
   * Fields without a reader, e.g. those reconstructed from a datastore,
   * cannot read a function alternative.
   * \param [in] lookupPath The path of this Field within the input file
   *****************************************************************************
   */
  Field(axom::sidre::Group* sidreGroup,
        axom::sidre::Group* root,
        axom::sidre::DataTypeId type = axom::sidre::DataTypeId::NO_TYPE_ID,
        bool docEnabled = true,
        Reader* reader = nullptr,
        const std::string& lookupPath = "");

  ~Field() override;

  /*!
   *****************************************************************************
   * \brief Returns pointer to the Sidre Group class for this Field.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for this Field class.
   *
   * \return Pointer to the Sidre Group class for this Field
   *****************************************************************************
   */
  const axom::sidre::Group* sidreGroup() const { return m_sidreGroup; };

  Field& required(bool isRequired = true) override;

  bool isRequired() const override;

  Field& defaultValue(const std::string& value) override;

  Field& defaultValue(const char* value) override;

  Field& defaultValue(bool value) override;

  Field& defaultValue(int value) override;

  Field& defaultValue(double value) override;

  Field& range(double startVal, double endVal) override;

  Field& range(int startVal, int endVal) override;

  Field& validValues(const std::vector<int>& set) override;

  Field& validValues(const std::vector<double>& set) override;

  Field& validValues(const std::vector<std::string>& set) override;

  Field& validValues(const std::initializer_list<const char*>& set) override;

  Field& validValues(const std::initializer_list<int>& set) override;

  Field& validValues(const std::initializer_list<double>& set) override;

  /*!
   *****************************************************************************
   * \brief Allows this Field to accept a value or a function.
   *
   * The Field's type sets the function's return type.
   * A double Field accepts a function returning FunctionTag::Double.
   * A string Field accepts a function returning FunctionTag::String.
   *
   * A supplied function makes the Field present and user-provided, and satisfies required().
   * In that case, the Field doesn't have a concrete value or a default.
   * Use Proxy::evaluate() to read either form, or Proxy::type()
   * to determine which form the input supplied.
   *
   * \param [in] arg_types The argument types of the function
   *
   * \return Reference to this Field, for chaining
   *
   * \pre The Field holds a double or a string
   * \note Range and valid-value constraints apply only to a concrete value.
   *****************************************************************************
   */
  Field& orFunction(const std::vector<FunctionTag>& arg_types = {}) override;

  using VerifiableScalar::registerVerifier;

  Field& registerVerifier(Verifier lambda) override;

  bool verify(std::vector<VerificationError>* errors = nullptr) const override;

  /*!
   *****************************************************************************
   * \return The full name for this Field.
   *****************************************************************************
  */
  std::string name() const;

  /*!
   *****************************************************************************
   * \brief Returns a value of primitive type
   * 
   * \return The value
   * \tparam T The type to retrieve
   *****************************************************************************
   */
  template <typename T>
  T get() const;

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
   *****************************************************************************
   * \brief Returns whether a value for the Field exists, i.e., if a value 
   * was provided in the input file or if a default was provided
   *
   * A function supplied through orFunction() also counts.
   *****************************************************************************
   */
  bool exists() const;

  /*!
   *****************************************************************************
   * \brief Returns whether a value was provided in the input file
   *
   * A function supplied through orFunction() also counts.
   *****************************************************************************
   */
  bool isUserProvided() const;

  /*!
   *****************************************************************************
   * \brief Returns whether the input supplied a function for this Field.
   *
   * \see orFunction()
   *****************************************************************************
   */
  bool isFunction() const;

  /*!
   *****************************************************************************
   * \brief Returns the function the input supplied for this Field.
   *
   * \pre isFunction() is true
   *****************************************************************************
   */
  const Function& getFunction() const;

private:
  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set A vector containing the set of allowed scalar values
   *****************************************************************************
  */
  template <typename T>
  void setScalarValidValues(std::vector<T> set);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *****************************************************************************
  */
  template <typename T>
  void setRange(T startVal, T endVal);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default value
   *****************************************************************************
  */
  template <typename T>
  void setDefaultValue(T value);

  /*!
   *****************************************************************************
   * \brief Checks the validity of a field value
   *
   * \param [in] view The view to verify
   * \param [out] errors The list of errors to record to, when non-null
   *
   * \return Whether the value satisfied all constraints
   *****************************************************************************
  */
  bool verifyValue(const axom::sidre::View& view,
                   std::vector<VerificationError>* errors = nullptr) const;

  /*!
   *****************************************************************************
   * \brief Checks if the given value is within the range.
   * 
   * \param [in] view The view containing the value that will be checked.
   * \param [out] errors The list of errors to record to, when non-null
   * 
   * \return true if the given value was within its respective range, else false.
   * \pre T must define bool operator<=(T, T)
   *****************************************************************************
   */
  template <typename T>
  bool checkRange(const axom::sidre::View& view,
                  std::vector<VerificationError>* errors = nullptr) const;

  /*!
   *****************************************************************************
   * \brief Checks if the given value is found in the list of valid values.
   * 
   * \param [in] view The view containing the value that will be checked.
   * \param [out] errors The list of errors to record to, when non-null
   * 
   * \return true if the given target was found in its respective valid values, 
   *  else false.
   *****************************************************************************
   */
  template <typename T>
  bool searchValidValues(const axom::sidre::View& view,
                         std::vector<VerificationError>* errors = nullptr) const;

  /*!
   *****************************************************************************
   * \brief Checks the existence and type of the value for the field
   *
   * \param [in] expected The expected type for the value
   *
   * \return Non-owning pointer to the Sidre view containing the value
   * \note Treats a nonexistent value or type mismatch as an error and will
   * emit a SLIC_ERROR accordingly
   *****************************************************************************
  */
  const axom::sidre::View* checkExistenceAndType(const axom::sidre::DataTypeId expected) const;

  // This Field's sidre group
  axom::sidre::Group* m_sidreGroup {nullptr};
  axom::sidre::Group* m_sidreRootGroup {nullptr};
  axom::sidre::DataTypeId m_type {axom::sidre::DataTypeId::NO_TYPE_ID};
  bool m_docEnabled;
  Verifier m_verifier;
  // Used by orFunction() to read a function alternative
  Reader* m_reader {nullptr};
  std::string m_lookupPath;
  // The function supplied in place of a value, if any
  std::unique_ptr<Function> m_function;
};

// Prototypes for template specializations
template <>
bool Field::get<bool>() const;

template <>
int Field::get<int>() const;

template <>
double Field::get<double>() const;

template <>
std::string Field::get<std::string>() const;

template <>
bool Field::searchValidValues<std::string>(const axom::sidre::View& view,
                                           std::vector<VerificationError>* errors) const;

/*!
   *****************************************************************************
   * \brief A wrapper class that enables constraints on groups of Fields
   *****************************************************************************
  */
class AggregateField : public VerifiableScalar
{
public:
  AggregateField(std::vector<std::reference_wrapper<VerifiableScalar>>&& fields)
    : m_fields(std::move(fields))
  { }

  virtual ~AggregateField() = default;

  bool verify(std::vector<VerificationError>* errors = nullptr) const override;

  AggregateField& required(bool isRequired) override;

  bool isRequired() const override;

  AggregateField& defaultValue(const std::string& value) override;

  AggregateField& defaultValue(const char* value) override;

  AggregateField& defaultValue(bool value) override;

  AggregateField& defaultValue(int value) override;

  AggregateField& defaultValue(double value) override;

  AggregateField& range(double startVal, double endVal) override;

  AggregateField& range(int startVal, int endVal) override;

  AggregateField& validValues(const std::vector<int>& set) override;

  AggregateField& validValues(const std::vector<double>& set) override;

  AggregateField& validValues(const std::vector<std::string>& set) override;

  AggregateField& validValues(const std::initializer_list<const char*>& set) override;

  AggregateField& validValues(const std::initializer_list<int>& set) override;

  AggregateField& validValues(const std::initializer_list<double>& set) override;

  AggregateField& orFunction(const std::vector<FunctionTag>& arg_types = {}) override;

  using VerifiableScalar::registerVerifier;

  AggregateField& registerVerifier(Verifier lambda) override;

private:
  std::vector<std::reference_wrapper<VerifiableScalar>> m_fields;
};

}  // end namespace inlet
}  // end namespace axom
