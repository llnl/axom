// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Field.hpp
 *
 * \brief This file contains the class definition of Inlet's Field and AggregateField classes.
 *******************************************************************************
 */

#ifndef INLET_FIELD_HPP
#define INLET_FIELD_HPP

#include "axom/sidre.hpp"
#include "axom/inlet/VariantKey.hpp"
#include "axom/inlet/VerifiableScalar.hpp"

#include <memory>
#include <type_traits>
#include <functional>

namespace axom
{
namespace inlet
{
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
   *****************************************************************************
   */
  Field(axom::sidre::Group* sidreGroup,
        axom::sidre::Group* root,
        axom::sidre::DataTypeId type = axom::sidre::DataTypeId::NO_TYPE_ID,
        bool docEnabled = true)
    : m_sidreGroup(sidreGroup)
    , m_sidreRootGroup(root)
    , m_type(type)
    , m_docEnabled(docEnabled)
  { }

  virtual ~Field() = default;

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
   *****************************************************************************
   */
  bool exists() const;

  /*!
   *****************************************************************************
   * \brief Returns whether a value was provided in the input file
   *****************************************************************************
   */
  bool isUserProvided() const;

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
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  axom::sidre::DataTypeId m_type = axom::sidre::DataTypeId::NO_TYPE_ID;
  bool m_docEnabled;
  Verifier m_verifier;
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

  using VerifiableScalar::registerVerifier;

  AggregateField& registerVerifier(Verifier lambda) override;

private:
  std::vector<std::reference_wrapper<VerifiableScalar>> m_fields;
};

}  // end namespace inlet
}  // end namespace axom

#endif
