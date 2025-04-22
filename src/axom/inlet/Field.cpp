// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Field.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "axom/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
Field& Field::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "[Inlet] Field specific Sidre Datastore Group not set");
  setFlag(*m_sidreGroup, *m_sidreRootGroup, detail::REQUIRED_FLAG, isRequired);
  return *this;
}

bool Field::isRequired() const
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "[Inlet] Field specific Sidre Datastore Group not set");
  return checkFlag(*m_sidreGroup, *m_sidreRootGroup, detail::REQUIRED_FLAG);
}

template <typename T>
void Field::setDefaultValue(T value)
{
  if(m_sidreGroup->hasView("defaultValue"))
  {
    const std::string msg = fmt::format(
      "[Inlet] Field has already defined "
      "default value: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    if(m_docEnabled)
    {
      m_sidreGroup->createViewScalar("defaultValue", value);
    }
    if(!m_sidreGroup->hasView("value"))
    {
      m_sidreGroup->createViewScalar("value", value);
    }
  }
}

template <>
void Field::setDefaultValue<std::string>(std::string value)
{
  if(m_sidreGroup->hasView("defaultValue"))
  {
    const std::string msg = fmt::format(
      "[Inlet] Field has already defined "
      "default value: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    if(m_docEnabled)
    {
      m_sidreGroup->createViewString("defaultValue", value);
    }
    if(!m_sidreGroup->hasView("value"))
    {
      m_sidreGroup->createViewString("value", value);
    }
  }
}

Field& Field::defaultValue(const char* value)
{
  std::string str = "";
  if(value)
  {
    str = value;
  }
  setDefaultValue(std::string(value));
  return *this;
}

Field& Field::defaultValue(const std::string& value)
{
  if(m_type != axom::sidre::DataTypeId::CHAR8_STR_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match STRING");
    setWarningFlag(m_sidreRootGroup);
  }
  setDefaultValue(value);
  return *this;
}

Field& Field::defaultValue(bool value)
{
  if(m_type != axom::sidre::DataTypeId::INT8_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match BOOL");
    setWarningFlag(m_sidreRootGroup);
  }
  setDefaultValue((std::int8_t)value);
  return *this;
}

Field& Field::defaultValue(int value)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::INT_ID:
    setDefaultValue(value);
    break;
  case axom::sidre::DataTypeId::DOUBLE_ID:
    setDefaultValue((double)value);
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match INT/DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return *this;
}

Field& Field::defaultValue(double value)
{
  if(m_type != axom::sidre::DataTypeId::DOUBLE_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
  }
  setDefaultValue(value);
  return *this;
}

template <typename T>
void Field::setRange(T startVal, T endVal)
{
  if(m_sidreGroup->hasView("range"))
  {
    const std::string msg =
      fmt::format("[Inlet] Inlet Field has already defined range: {0}", m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else if(m_sidreGroup->hasView("validValues") || m_sidreGroup->hasView("validStringValues"))
  {
    const std::string msg = fmt::format(
      "[Inlet] Cannot set range for Inlet Field "
      "after setting valid values: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    auto* view = m_sidreGroup->createViewAndAllocate("range", m_type, 2);
    T* pair = view->getArray();
    pair[0] = startVal;
    pair[1] = endVal;
  }
}

Field& Field::range(int startVal, int endVal)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::INT_ID:
    this->setRange(startVal, endVal);
    break;
  case axom::sidre::DataTypeId::DOUBLE_ID:
    this->setRange(static_cast<double>(startVal), static_cast<double>(endVal));
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match INT or DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return *this;
}

Field& Field::range(double startVal, double endVal)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::DOUBLE_ID:
    this->setRange(startVal, endVal);
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return *this;
}

template <>
bool Field::get<bool>() const
{
  const auto valueView = checkExistenceAndType(axom::sidre::INT8_ID);
  // There is no boolean type in conduit/sidre so we use std::int8_t
  const std::int8_t intValue = valueView->getScalar();
  if(intValue < 0 || intValue > 1)
  {
    const std::string msg = fmt::format(
      "[Inlet] Invalid integer value stored in "
      " boolean value named {0}",
      name());
    SLIC_ERROR(msg);
  }

  return static_cast<bool>(intValue);
}

template <>
double Field::get<double>() const
{
  const auto valueView = checkExistenceAndType(axom::sidre::DOUBLE_ID);
  return valueView->getScalar();
}

template <>
int Field::get<int>() const
{
  const auto valueView = checkExistenceAndType(axom::sidre::INT_ID);
  return valueView->getScalar();
}

template <>
std::string Field::get<std::string>() const
{
  const auto valueView = checkExistenceAndType(axom::sidre::CHAR8_STR_ID);

  const char* valueStr = valueView->getString();
  return (valueStr == nullptr) ? "" : valueStr;
}

const axom::sidre::View* Field::checkExistenceAndType(const axom::sidre::DataTypeId expected) const
{
  const axom::sidre::View* valueView = m_sidreGroup->getView("value");

  if(valueView == nullptr)
  {
    SLIC_ERROR("[Inlet] Field does not contain a value");
  }

  if(valueView->getTypeID() != expected)
  {
    const std::string msg = fmt::format(
      "[Inlet] Field with name '{0}' was expected to be of type '{1}'"
      " but was actually of type '{2}'",
      name(),
      conduit::DataType::id_to_name(expected),
      conduit::DataType::id_to_name(valueView->getTypeID()));
    SLIC_ERROR(msg);
  }

  return valueView;
}

InletType Field::type() const
{
  axom::sidre::View* valueView = m_sidreGroup->getView("value");
  if(valueView == nullptr)
  {
    return InletType::Nothing;
  }

  switch(valueView->getTypeID())
  {
  case axom::sidre::NO_TYPE_ID:
    return InletType::Nothing;
  case axom::sidre::INT8_ID:
    return InletType::Bool;
  case axom::sidre::INT_ID:
    return InletType::Integer;
  case axom::sidre::CHAR8_STR_ID:
    return InletType::String;
  case axom::sidre::DOUBLE_ID:
    return InletType::Double;
  default:
    std::string msg =
      fmt::format("Type ID {0} for field not recognized, returning InletType::Nothing",
                  valueView->getTypeID());
    SLIC_WARNING(msg);
    return InletType::Nothing;
  }
}

bool Field::exists() const { return m_sidreGroup->hasView("value"); }

bool Field::isUserProvided() const
{
  if(m_sidreGroup->hasView("retrieval_status"))
  {
    auto status = static_cast<ReaderResult>(
      static_cast<int>(m_sidreGroup->getView("retrieval_status")->getData()));
    // Even if it wasn't of the right type, we still say that it was user-provided
    return (status != ReaderResult::NotFound);
  }
  return exists();
}

template <typename T>
void Field::setScalarValidValues(std::vector<T> set)
{
  if(m_sidreGroup->hasView("validValues") || m_sidreGroup->hasView("validStringValues"))
  {
    std::string msg = fmt::format(
      "[Inlet] Inlet Field has already "
      "defined valid values: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else if(m_sidreGroup->hasView("range"))
  {
    std::string msg = fmt::format(
      "[Inlet] Cannot set valid values "
      "after defining range: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    auto view = m_sidreGroup->createViewAndAllocate("validValues", m_type, set.size());
    view->getBuffer()->copyBytesIntoBuffer(&set[0], set.size() * sizeof(T));
  }
}

Field& Field::validValues(const std::vector<int>& set)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::INT_ID:
    setScalarValidValues(set);
    break;
  case axom::sidre::DataTypeId::DOUBLE_ID:
    setScalarValidValues(std::vector<double>(set.begin(), set.end()));
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match INT OR DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return *this;
}

Field& Field::validValues(const std::vector<double>& set)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::DOUBLE_ID:
    setScalarValidValues(set);
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return *this;
}

Field& Field::validValues(const std::vector<std::string>& set)
{
  if(m_type != axom::sidre::DataTypeId::CHAR8_STR_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match STRING");
    setWarningFlag(m_sidreRootGroup);
  }
  if(m_sidreGroup->hasView("validValues") || m_sidreGroup->hasView("validStringValues"))
  {
    std::string msg = fmt::format(
      "[Inlet] Inlet Field has already defined "
      "valid values: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else if(m_sidreGroup->hasView("range"))
  {
    std::string msg = fmt::format(
      "[Inlet] Cannot set valid values after defining "
      "range: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    auto group = m_sidreGroup->createGroup("validStringValues", true);
    for(std::string str : set)
    {
      group->createViewString("", str);
    }
  }
  return *this;
}

Field& Field::validValues(const std::initializer_list<const char*>& set)
{
  return validValues(std::vector<std::string>(set.begin(), set.end()));
}

Field& Field::validValues(const std::initializer_list<int>& set)
{
  return validValues(std::vector<int>(set));
}

Field& Field::validValues(const std::initializer_list<double>& set)
{
  return validValues(std::vector<double>(set));
}

Field& Field::registerVerifier(Verifier lambda)
{
  SLIC_WARNING_IF(m_verifier,
                  fmt::format("[Inlet] Verifier for Field "
                              "already set: {0}",
                              m_sidreGroup->getPathName()));
  m_verifier = lambda;
  return *this;
}

bool Field::verify(std::vector<VerificationError>* errors) const
{
  // If this field was required, make sure something was defined in it
  if(!verifyRequired(*m_sidreGroup, m_sidreGroup->hasView("value"), "Field", errors))
  {
    return false;
  }
  // Verify the provided value
  if(m_sidreGroup->hasView("value") && !verifyValue(*m_sidreGroup->getView("value"), errors))
  {
    return false;
  }

  // Verify the default value
  if(m_sidreGroup->hasView("defaultValue") &&
     !verifyValue(*m_sidreGroup->getView("defaultValue"), errors))
  {
    return false;
  }

  // Lambda verification step
  if(m_verifier && !m_verifier(*this, errors))
  {
    const std::string msg =
      fmt::format("[Inlet] Field failed lambda verification: {0}", m_sidreGroup->getPathName());
    INLET_VERIFICATION_WARNING(m_sidreGroup->getPathName(), msg, errors);
    return false;
  }
  return true;
}

bool Field::verifyValue(const axom::sidre::View& view, std::vector<VerificationError>* errors) const
{
  const auto type = view.getTypeID();
  if(m_sidreGroup->hasView("validValues"))
  {
    if(type == axom::sidre::INT_ID)
    {
      return searchValidValues<int>(view, errors);
    }
    else
    {
      return searchValidValues<double>(view, errors);
    }
  }
  else if(m_sidreGroup->hasView("range"))
  {
    if(type == axom::sidre::INT_ID)
    {
      return checkRange<int>(view, errors);
    }
    else
    {
      return checkRange<double>(view, errors);
    }
  }
  else if(m_sidreGroup->hasGroup("validStringValues"))
  {
    return searchValidValues<std::string>(view, errors);
  }
  return true;
}

template <typename T>
bool Field::checkRange(const axom::sidre::View& view, std::vector<VerificationError>* errors) const
{
  const T val = view.getScalar();
  const T* range = m_sidreGroup->getView("range")->getArray();
  const bool in_range = range[0] <= val && val <= range[1];
  if(!in_range)
  {
    std::string msg = fmt::format(
      "[Inlet] Given value '{0}' for '{1}' did not meet range "
      " constraints: [{2}, {3}]",
      val,
      view.getPath(),
      range[0],
      range[1]);
    INLET_VERIFICATION_WARNING(view.getPath(), msg, errors);
  }
  return in_range;
}

template <typename T>
bool Field::searchValidValues(const axom::sidre::View& view,
                              std::vector<VerificationError>* errors) const
{
  const T target = view.getScalar();
  const auto valid_vals = m_sidreGroup->getView("validValues");
  const T* valuesArray = valid_vals->getArray();
  const size_t size = valid_vals->getBuffer()->getNumElements();
  const T* result = std::find(valuesArray, valuesArray + size, target);
  const bool is_valid = result != valuesArray + size;
  if(!is_valid)
  {
    std::string msg = fmt::format(
      "[Inlet] Given value '{0}' for '{1}' did not meet valid value(s) "
      " constraints: {2}",
      target,
      view.getPath(),
      fmt::join(valuesArray, valuesArray + size, ","));
    INLET_VERIFICATION_WARNING(view.getPath(), msg, errors);
  }
  return is_valid;
}

template <>
bool Field::searchValidValues<std::string>(const axom::sidre::View& view,
                                           std::vector<VerificationError>* errors) const
{
  const auto string_group = m_sidreGroup->getGroup("validStringValues");
  const std::string value = view.getString();
  bool is_valid = false;
  std::vector<std::string> valid_values;
  for(auto& view : string_group->views())
  {
    // Store the valid values so we can print them in the error message
    valid_values.push_back(view.getString());
    if(valid_values.back() == value)
    {
      is_valid = true;
    }
  }
  if(!is_valid)
  {
    std::string msg = fmt::format(
      "[Inlet] Given value '{0}' for '{1}' did not meet valid value(s) "
      " constraints: {2}",
      value,
      view.getPath(),
      fmt::join(valid_values, ","));
    INLET_VERIFICATION_WARNING(view.getPath(), msg, errors);
  }
  return is_valid;
}

std::string Field::name() const
{
  return utilities::string::removePrefix(m_sidreRootGroup->getPathName(),
                                         m_sidreGroup->getPathName());
}

bool AggregateField::verify(std::vector<VerificationError>* errors) const
{
  return std::all_of(m_fields.begin(), m_fields.end(), [&errors](const VerifiableScalar& field) {
    return field.verify(errors);
  });
}

AggregateField& AggregateField::required(bool isRequired)
{
  for(auto& field : m_fields)
  {
    field.get().required(isRequired);
  }
  return *this;
}

bool AggregateField::isRequired() const
{
  return std::any_of(m_fields.begin(), m_fields.end(), [](const VerifiableScalar& field) {
    return field.isRequired();
  });
}

AggregateField& AggregateField::defaultValue(const std::string& value)
{
  for(auto& field : m_fields)
  {
    field.get().defaultValue(value);
  }
  return *this;
}

AggregateField& AggregateField::defaultValue(const char* value)
{
  for(auto& field : m_fields)
  {
    field.get().defaultValue(value);
  }
  return *this;
}

AggregateField& AggregateField::defaultValue(bool value)
{
  for(auto& field : m_fields)
  {
    field.get().defaultValue(value);
  }
  return *this;
}

AggregateField& AggregateField::defaultValue(int value)
{
  for(auto& field : m_fields)
  {
    field.get().defaultValue(value);
  }
  return *this;
}

AggregateField& AggregateField::defaultValue(double value)
{
  for(auto& field : m_fields)
  {
    field.get().defaultValue(value);
  }
  return *this;
}

AggregateField& AggregateField::range(double startVal, double endVal)
{
  for(auto& field : m_fields)
  {
    field.get().range(startVal, endVal);
  }
  return *this;
}

AggregateField& AggregateField::range(int startVal, int endVal)
{
  for(auto& field : m_fields)
  {
    field.get().range(startVal, endVal);
  }
  return *this;
}

AggregateField& AggregateField::validValues(const std::vector<int>& set)
{
  for(auto& field : m_fields)
  {
    field.get().validValues(set);
  }
  return *this;
}

AggregateField& AggregateField::validValues(const std::vector<double>& set)
{
  for(auto& field : m_fields)
  {
    field.get().validValues(set);
  }
  return *this;
}

AggregateField& AggregateField::validValues(const std::vector<std::string>& set)
{
  for(auto& field : m_fields)
  {
    field.get().validValues(set);
  }
  return *this;
}

AggregateField& AggregateField::validValues(const std::initializer_list<const char*>& set)
{
  for(auto& field : m_fields)
  {
    field.get().validValues(set);
  }
  return *this;
}

AggregateField& AggregateField::validValues(const std::initializer_list<int>& set)
{
  for(auto& field : m_fields)
  {
    field.get().validValues(set);
  }
  return *this;
}

AggregateField& AggregateField::validValues(const std::initializer_list<double>& set)
{
  for(auto& field : m_fields)
  {
    field.get().validValues(set);
  }
  return *this;
}

AggregateField& AggregateField::registerVerifier(Verifier lambda)
{
  for(auto& field : m_fields)
  {
    field.get().registerVerifier(lambda);
  }
  return *this;
}

}  // end namespace inlet
}  // end namespace axom
