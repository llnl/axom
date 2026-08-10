// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "GeometryOperatorsIO.hpp"
#include "IOUtil.hpp"

#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"

#include "axom/fmt.hpp"

#include <functional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace axom
{
namespace klee
{
namespace internal
{
namespace
{
using OpPtr = CompositeOperator::OpPtr;
using OperatorParser =
  std::function<OpPtr(const SingleOperatorData &,
                      const TransformableGeometryProperties &,
                      const std::string &)>;
using internal::toDoubleVector;
using primal::Point3D;
using primal::Vector3D;
using FieldSet = std::unordered_set<std::string>;

std::string childName(const inlet::Container& container, const std::string& name)
{
  std::string result = axom::utilities::string::removePrefix(container.name(), name);
  if(axom::utilities::string::startsWith(result, '/'))
  {
    result = result.substr(1);
  }
  return result;
}

bool hasCallback(const inlet::Container &container, char const *fieldName)
{
  return container.containsFunctionValueAlternative(fieldName);
}

bool containsFieldOrCallback(const inlet::Container &container, char const *fieldName)
{
  return container.contains(fieldName) || hasCallback(container, fieldName);
}

Path fieldPath(const inlet::Container &container, char const *fieldName)
{
  return Path::join({Path {container.name()}, Path {std::string {fieldName}}});
}

std::string callbackContext(const inlet::Container &container,
                            char const *fieldName,
                            const std::string &ownerLabel)
{
  Path path {container.name()};
  std::string operatorIndex = path.baseName();
  if(operatorIndex == "slice")
  {
    operatorIndex = path.parent().baseName();
  }

  const auto operatorLabel = operatorIndex.empty() ? std::string {"operator at "} + container.name()
                                                   : std::string {"operator "} + operatorIndex;
  if(ownerLabel.empty())
  {
    return axom::fmt::format("Error evaluating callback for '{}' in {}", fieldName, operatorLabel);
  }
  return axom::fmt::format(
    "Error evaluating callback for '{}' in {} {}",
    fieldName,
    ownerLabel,
    operatorLabel);
}

[[noreturn]] void throwCallbackAwareValidationError(const inlet::Container &container,
                                                    char const *fieldName,
                                                    const std::string &ownerLabel,
                                                    const Path &fallbackPath,
                                                    const std::string &message)
{
  if(hasCallback(container, fieldName))
  {
    throw KleeError(
      {fieldPath(container, fieldName),
       axom::fmt::format("{}: {}", callbackContext(container, fieldName, ownerLabel), message)});
  }

  throw KleeError({fallbackPath, message});
}

template <typename Result, typename Func>
Result wrapCallbackErrors(const inlet::Container &container,
                          char const *fieldName,
                          const std::string &ownerLabel,
                          Func &&func)
{
  // Convert generic Inlet/Lua callback failures into Klee diagnostics at the
  // boundary where the shape, operator, and field context are all available.
  try
  {
    return func();
  }
  catch(const KleeError &)
  {
    throw;
  }
  catch(const std::exception &ex)
  {
    throw KleeError(
      {fieldPath(container, fieldName),
       axom::fmt::format("{}: {}", callbackContext(container, fieldName, ownerLabel), ex.what())});
  }
}

double getScalar(const inlet::Container &container,
                 char const *fieldName,
                 const std::string &ownerLabel)
{
  if(hasCallback(container, fieldName))
  {
    return wrapCallbackErrors<double>(container, fieldName, ownerLabel, [&]() {
      return container.getFunctionValueAlternative(fieldName).call<inlet::FunctionType::Double>();
    });
  }
  return container[fieldName].get<double>();
}

std::string getString(const inlet::Container &container,
                      char const *fieldName,
                      const std::string &ownerLabel)
{
  if(hasCallback(container, fieldName))
  {
    return wrapCallbackErrors<std::string>(container, fieldName, ownerLabel, [&]() {
      return container.getFunctionValueAlternative(fieldName).call<inlet::FunctionType::String>();
    });
  }
  return container[fieldName].get<std::string>();
}

std::vector<double> callbackVectorToDoubleVector(const inlet::FunctionType::Vector &value)
{
  std::vector<double> result;
  result.reserve(value.dim);
  for(int i = 0; i < value.dim; ++i)
  {
    result.push_back(value.vec[i]);
  }
  return result;
}

std::vector<double> getDoubleVector(const inlet::Container &container,
                                    char const *fieldName,
                                    Dimensions expectedDims,
                                    const std::string &ownerLabel)
{
  if(hasCallback(container, fieldName))
  {
    auto values = wrapCallbackErrors<std::vector<double>>(container, fieldName, ownerLabel, [&]() {
      return callbackVectorToDoubleVector(
        container.getFunctionValueAlternative(fieldName).call<inlet::FunctionType::Vector>());
    });
    auto actualSize = values.size();
    auto expectedSize = static_cast<std::size_t>(expectedDims);
    if(actualSize != expectedSize)
    {
      throw KleeError({fieldPath(container, fieldName),
                       fmt::format("{}: Wrong size for {}. Expected {}. Got {}.",
                                   callbackContext(container, fieldName, ownerLabel),
                                   fieldName,
                                   expectedSize,
                                   actualSize)});
    }
    return values;
  }
  return toDoubleVector(container[fieldName], expectedDims, fieldName);
}

template <typename T>
T toArrayLike(const inlet::Container &parent,
              char const *fieldName,
              Dimensions expectedDims,
              const std::string &ownerLabel)
{
  auto values = getDoubleVector(parent, fieldName, expectedDims, ownerLabel);
  return T {values.data(), static_cast<int>(expectedDims)};
}

template <typename T>
T toArrayLike(const inlet::Container &parent,
              char const *fieldName,
              Dimensions expectedDims,
              const T &defaultValue,
              const std::string &ownerLabel)
{
  if(containsFieldOrCallback(parent, fieldName))
  {
    return toArrayLike<T>(parent, fieldName, expectedDims, ownerLabel);
  }
  return defaultValue;
}

Point3D getPoint(const inlet::Container &parent,
                 char const *fieldName,
                 Dimensions expectedDims,
                 const std::string &ownerLabel)
{
  return toArrayLike<Point3D>(parent, fieldName, expectedDims, ownerLabel);
}

Point3D getPoint(const inlet::Container &parent,
                 char const *fieldName,
                 Dimensions expectedDims,
                 const Point3D &defaultValue,
                 const std::string &ownerLabel)
{
  return toArrayLike(parent, fieldName, expectedDims, defaultValue, ownerLabel);
}

Vector3D getVector(const inlet::Container &parent,
                   char const *fieldName,
                   Dimensions expectedDims,
                   const std::string &ownerLabel)
{
  return toArrayLike<Vector3D>(parent, fieldName, expectedDims, ownerLabel);
}

Vector3D getVector(const inlet::Container &parent,
                   char const *fieldName,
                   Dimensions expectedDims,
                   const Vector3D &defaultValue,
                   const std::string &ownerLabel)
{
  return toArrayLike(parent, fieldName, expectedDims, defaultValue, ownerLabel);
}

/**
 * Get the names of all the children in the given container.
 *
 * @param container the Container whose children to get
 * @return the names of all the children
 */
std::unordered_set<std::string> getChildNames(const inlet::Container& container)
{
  std::unordered_set<std::string> allChildren;

  std::vector<std::string> unexpectedNames = container.unexpectedNames();
  allChildren.insert(unexpectedNames.begin(), unexpectedNames.end());

  for(auto& child : container.getChildContainers())
  {
    if(child.second->exists())
    {
      allChildren.insert(childName(container, child.first));
    }
  }

  for(auto& child : container.getChildFields())
  {
    if(child.second->exists())
    {
      allChildren.insert(childName(container, child.first));
    }
  }

  for(const auto &name : container.getFunctionValueAlternativeNames())
  {
    allChildren.insert(name);
  }

  return allChildren;
}

/**
 * Verify that a Container has the correct fields.
 *
 * While Inlet can do a lot of the verification, there are some situations
 * where we need some extra manual checks. For example, operators can look like this:
 *
 * \code{.yaml}
 *
 * operators:
 *   - translate: [10, 20, 30]
 *   - rotate: 90
 *     center: [1, 2, 3]
 *     axis: [4, 5, 6]
 *
 * \endcode
 *
 * In the above, "translate", "rotate", "center", and "axis" are all valid
 * entries, but not in arbitrary combinations. You can't specify both
 * "translate" and "axis", for example, or "translate" and "rotate" within
 * the same entry.
 *
 * This function can be used to handle cases like the above.
 *
 * \param containerToTest the Container to test
 * \param name the name of the container. This must be one of its fields.
 * \param additionalRequiredFields any additional required fields
 * \param optionalFields any additional optional fields
 * \throws KleeError if a required field is missing or an unexpected field is present
 */
void verifyObjectFields(const inlet::Container& containerToTest,
                        const std::string& name,
                        const FieldSet& additionalRequiredFields,
                        const FieldSet& optionalFields)
{
  std::unordered_set<std::string> requiredFields {additionalRequiredFields};
  requiredFields.insert(name);

  for(auto& requiredField : requiredFields)
  {
    if(!containsFieldOrCallback(containerToTest, requiredField.c_str()))
    {
      throw KleeError(
        {containerToTest.name(),
         axom::fmt::format("Missing required parameter '{}' for operator '{}'", requiredField, name)});
    }
  }

  for(auto& child : getChildNames(containerToTest))
  {
    if(requiredFields.find(child) != requiredFields.end())
    {
      continue;
    }
    if(optionalFields.find(child) != optionalFields.end())
    {
      continue;
    }

    throw KleeError({containerToTest.name(),
                     axom::fmt::format("Unexpected parameter '{}' for operator '{}'", name, child)});
  }
}

/**
 * Parse a "translate" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 * \throws KleeError if the operator fields or vector dimensions are invalid
 */
OpPtr parseTranslate(const SingleOperatorData &data,
                     const TransformableGeometryProperties &startProperties,
                     const std::string &ownerLabel)
{
  const auto &opContainer = *data.m_container;
  verifyObjectFields(opContainer, "translate", FieldSet {}, FieldSet {});
  return std::make_shared<Translation>(
    getVector(opContainer, "translate", startProperties.dimensions, ownerLabel),
    startProperties);
}

/**
 * Parse a "rotate" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 * \throws KleeError if the rotation is invalid for the start dimensions or operator fields
 */
OpPtr parseRotate(const SingleOperatorData &data,
                  const TransformableGeometryProperties &startProperties,
                  const std::string &ownerLabel)
{
  const auto &opContainer = *data.m_container;
  switch(startProperties.dimensions)
  {
  case Dimensions::Two:
  {
    verifyObjectFields(opContainer, "rotate", FieldSet {}, {"center"});
    auto angle = getScalar(opContainer, "rotate", ownerLabel);
    auto center =
      getPoint(opContainer, "center", Dimensions::Two, Point3D {0, 0, 0}, ownerLabel);
    Vector3D axis {0, 0, 1};
    return std::make_shared<Rotation>(angle, center, axis, startProperties);
  }
  break;
  case Dimensions::Three:
  {
    verifyObjectFields(opContainer, "rotate", {"axis"}, {"center"});
    auto angle = getScalar(opContainer, "rotate", ownerLabel);
    auto center =
      getPoint(opContainer, "center", Dimensions::Three, Point3D {0, 0, 0}, ownerLabel);
    auto axis = getVector(opContainer, "axis", Dimensions::Three, ownerLabel);
    if(axis.is_zero())
    {
      auto message = std::string {"The 'axis' vector must not be a zero vector"};
      if(hasCallback(opContainer, "axis"))
      {
        message = axom::fmt::format("{}: {}",
                                    callbackContext(opContainer, "axis", ownerLabel),
                                    message);
      }
      throw KleeError({fieldPath(opContainer, "axis"), message});
    }
    return std::make_shared<Rotation>(angle, center, axis, startProperties);
  }
  break;
  default:
  case Dimensions::Unspecified:
    throw KleeError({opContainer.name(), "Rotations can only be applied to 2D or 3D shapes"});
  }
}

/**
 * Make a slice, ensuring all the values are valid.
 *
 * \param origin the origin of the coordinate system
 * \param normal a vector normal to the plane
 * \param up a vector which defines the positive Y direction
 * \param startProperties the properties before the slice
 * \param sliceContainer the Inlet container describing the slice
 * \param ownerLabel a description of the owning shape or named operator
 * \return the created operator
 * \throws KleeError if any value is invalid
 */
OpPtr makeCheckedSlice(Point3D origin,
                       Vector3D normal,
                       Vector3D up,
                       const TransformableGeometryProperties &startProperties,
                       const inlet::Container &sliceContainer,
                       const std::string &ownerLabel)
{
  if(normal.is_zero())
  {
    throwCallbackAwareValidationError(sliceContainer,
                                      "normal",
                                      ownerLabel,
                                      Path {sliceContainer.name()},
                                      "The 'normal' vector must not be a zero vector");
  }
  if(!utilities::isNearlyEqual(normal.dot(up), 0.))
  {
    const std::string message = "The 'normal' and 'up' vectors must be perpendicular";
    if(hasCallback(sliceContainer, "up"))
    {
      throwCallbackAwareValidationError(
        sliceContainer,
        "up",
        ownerLabel,
        Path {sliceContainer.name()},
        message);
    }
    if(hasCallback(sliceContainer, "normal"))
    {
      throwCallbackAwareValidationError(
        sliceContainer,
        "normal",
        ownerLabel,
        Path {sliceContainer.name()},
        message);
    }
    throw KleeError({sliceContainer.name(), message});
  }
  return std::make_shared<SliceOperator>(origin, normal, up, startProperties);
}

/**
 * Get the origin to use for a perpendicular slice.
 *
 * \param sliceContainer the Container describing the slice
 * \param planeName the name of the plane ("x", "y", or "z")
 * \param defaultNormal the default normal vector
 * \return the point to use as the origin
 * \throws KleeError if the specified origin is not on the slice plane
 */
primal::Point3D getPerpendicularSliceOrigin(const inlet::Container &sliceContainer,
                                            char const *planeName,
                                            const primal::Vector3D &defaultNormal,
                                            const std::string &ownerLabel)
{
  double axisIntercept = getScalar(sliceContainer, planeName, ownerLabel);

  primal::Point3D defaultOrigin;
  int nonZeroIndex = -1;
  for(int i = 0; i < 3; ++i)
  {
    defaultOrigin[i] = axisIntercept * defaultNormal[i];
    if(!utilities::isNearlyEqual(defaultNormal[i], 0.0))
    {
      nonZeroIndex = i;
    }
  }

  if(!containsFieldOrCallback(sliceContainer, "origin"))
  {
    return defaultOrigin;
  }

  primal::Point3D givenOrigin = getPoint(sliceContainer, "origin", Dimensions::Three, ownerLabel);
  if(givenOrigin[nonZeroIndex] != axisIntercept)
  {
    throwCallbackAwareValidationError(sliceContainer,
                                      "origin",
                                      ownerLabel,
                                      Path {sliceContainer["origin"].name()},
                                      "The origin must be on the slice plane");
  }
  return givenOrigin;
}

/**
 * Get the normal vector to use for a perpendicular slice.
 *
 * \param sliceContainer the Container describing the slice
 * \param defaultNormal the default normal vector
 * \return the vector to use as the normal
 * \throws KleeError if the specified normal is not parallel to the slice plane normal
 */
primal::Vector3D getPerpendicularSliceNormal(const inlet::Container &sliceContainer,
                                             const primal::Vector3D &defaultNormal,
                                             const std::string &ownerLabel)
{
  if(!containsFieldOrCallback(sliceContainer, "normal"))
  {
    return defaultNormal;
  }

  primal::Vector3D givenNormal = getVector(sliceContainer, "normal", Dimensions::Three, ownerLabel);
  auto cross = primal::Vector3D::cross_product(givenNormal, defaultNormal);
  bool parallel = cross.is_zero();
  if(!parallel)
  {
    throwCallbackAwareValidationError(sliceContainer,
                                      "normal",
                                      ownerLabel,
                                      Path {sliceContainer["normal"].name()},
                                      "Invalid normal");
  }
  return givenNormal;
}

/**
 * Read a perpendicular slice.
 *
 * \param sliceContainer the Container describing the slice
 * \param planeName the name of the plane ("x", "y", or "z")
 * \param defaultNormal the default normal vector for the type of plane being parsed
 * \param defaultUp the default up vector for the plane being parsed
 * \param startProperties the properties prior to this operator
 * \return the parsed plane
 * \throws KleeError if the slice fields or values are invalid
 */
OpPtr readPerpendicularSlice(const inlet::Container &sliceContainer,
                             char const *planeName,
                             Vector3D const &defaultNormal,
                             Vector3D const &defaultUp,
                             const TransformableGeometryProperties &startProperties,
                             const std::string &ownerLabel)
{
  verifyObjectFields(sliceContainer, planeName, FieldSet {}, {"origin", "normal", "up"});
  const primal::Vector3D defaultNormalVec {defaultNormal.data()};

  auto origin = getPerpendicularSliceOrigin(sliceContainer, planeName, defaultNormalVec, ownerLabel);
  auto normal = getPerpendicularSliceNormal(sliceContainer, defaultNormalVec, ownerLabel);
  auto up = getVector(sliceContainer, "up", Dimensions::Three, defaultUp, ownerLabel);

  return makeCheckedSlice(origin, normal, up, startProperties, sliceContainer, ownerLabel);
}

/**
 * Parse a "slice" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 * \throws KleeError if the slice fields or values are invalid
 */
OpPtr parseSlice(const SingleOperatorData &data,
                 const TransformableGeometryProperties &startProperties,
                 const std::string &ownerLabel)
{
  const auto &opContainer = *data.m_container;
  if(startProperties.dimensions != Dimensions::Three)
  {
    throw KleeError({opContainer.name(), "Cannot do a slice from 2D"});
  }
  verifyObjectFields(opContainer, "slice", FieldSet {}, FieldSet {});
  auto &sliceContainer = *opContainer.getChildContainers().at(opContainer.name() + "/slice").get();
  if(containsFieldOrCallback(sliceContainer, "x"))
  {
    return readPerpendicularSlice(sliceContainer,
                                  "x",
                                  {1, 0, 0},
                                  {0, 0, 1},
                                  startProperties,
                                  ownerLabel);
  }
  else if(containsFieldOrCallback(sliceContainer, "y"))
  {
    return readPerpendicularSlice(sliceContainer,
                                  "y",
                                  {0, 1, 0},
                                  {1, 0, 0},
                                  startProperties,
                                  ownerLabel);
  }
  else if(containsFieldOrCallback(sliceContainer, "z"))
  {
    return readPerpendicularSlice(sliceContainer,
                                  "z",
                                  {0, 0, 1},
                                  {0, 1, 0},
                                  startProperties,
                                  ownerLabel);
  }

  verifyObjectFields(sliceContainer, "origin", {"normal", "up"}, FieldSet {});

  auto origin = getPoint(sliceContainer, "origin", Dimensions::Three, ownerLabel);
  auto normal = getVector(sliceContainer, "normal", Dimensions::Three, ownerLabel);
  auto up = getVector(sliceContainer, "up", Dimensions::Three, ownerLabel);
  return makeCheckedSlice(origin, normal, up, startProperties, sliceContainer, ownerLabel);
}

/**
 * Parse a "scale" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 * \throws KleeError if the scale fields or vector dimensions are invalid
 */
OpPtr parseScale(const SingleOperatorData &data,
                 const TransformableGeometryProperties &startProperties,
                 const std::string &ownerLabel)
{
  const auto &opContainer = *data.m_container;
  verifyObjectFields(opContainer, "scale", FieldSet {}, FieldSet {"center"});
  auto factors = hasCallback(opContainer, "scale")
    ? wrapCallbackErrors<std::vector<double>>(
        opContainer,
        "scale",
        ownerLabel,
        [&]() {
          return callbackVectorToDoubleVector(
            opContainer.getFunctionValueAlternative("scale").call<inlet::FunctionType::Vector>());
        })
    : opContainer["scale"].get<std::vector<double>>();

  const bool isUniform = factors.size() == 1;
  if(!isUniform && hasCallback(opContainer, "scale"))
  {
    auto actualSize = factors.size();
    auto expectedSize = static_cast<std::size_t>(startProperties.dimensions);
    if(actualSize != expectedSize)
    {
      throw KleeError({fieldPath(opContainer, "scale"),
                        fmt::format("{}: Wrong size for scale. Expected {}. Got {}.",
                                    callbackContext(opContainer, "scale", ownerLabel),
                                    expectedSize,
                                    actualSize)});
    }
  }
  else if(!isUniform)
  {
    factors = toDoubleVector(opContainer["scale"], startProperties.dimensions, "scale");
  }
  if(!isUniform && startProperties.dimensions == Dimensions::Two)
  {
    factors.emplace_back(1.0);
  }

  Point3D center {0., 0., 0.};
  if(containsFieldOrCallback(opContainer, "center"))
  {
    center = getPoint(opContainer,
                      "center",
                      startProperties.dimensions,
                      Point3D {0, 0, 0},
                      ownerLabel);
  }

  if(isUniform)
  {
    return std::make_shared<Scale>(
      factors[0],
      factors[0],
      factors[0],
      center,
      startProperties);
  }

  return std::make_shared<Scale>(factors[0], factors[1], factors[2], center, startProperties);
}

/**
 * Parse a "convert_units_to" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 * \throws KleeError if the unit string or operator fields are invalid
 */
OpPtr parseConvertUnits(const SingleOperatorData &data,
                        const TransformableGeometryProperties &startProperties,
                        const std::string &ownerLabel)
{
  const auto &opContainer = *data.m_container;
  verifyObjectFields(opContainer, "convert_units_to", FieldSet {}, FieldSet {});
  const auto unitName = getString(opContainer, "convert_units_to", ownerLabel);
  const auto path = fieldPath(opContainer, "convert_units_to");
  LengthUnit endUnits;
  try
  {
    endUnits = internal::parseLengthUnits(unitName, static_cast<std::string>(path));
  }
  catch(const KleeError &err)
  {
    if(!hasCallback(opContainer, "convert_units_to"))
    {
      throw;
    }
    throw KleeError(
      {path,
       axom::fmt::format(
         "{}: {}",
         callbackContext(opContainer, "convert_units_to", ownerLabel),
         err.what())});
  }
  return std::make_shared<UnitConverter>(endUnits, startProperties);
}

/**
 * Parse an operator specified via the "ref" command.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties before the "ref" command
 * \param namedOperators a map of named operators from which to get referenced operators
 * \return the created operator
 * \throws KleeError if the reference is missing or the operator fields are invalid
 */
OpPtr parseRef(const SingleOperatorData &data,
               const TransformableGeometryProperties &startProperties,
               const NamedOperatorMap &namedOperators,
               const std::string &ownerLabel)
{
  const auto &opContainer = *data.m_container;
  verifyObjectFields(opContainer, "ref", FieldSet {}, FieldSet {});
  const auto operatorName = getString(opContainer, "ref", ownerLabel);
  auto opIter = namedOperators.find(operatorName);
  if(opIter == namedOperators.end())
  {
    std::string message = "No operator named '";
    message += operatorName;
    message += '\'';
    if(hasCallback(opContainer, "ref"))
    {
      message = axom::fmt::format("{}: {}",
                                  callbackContext(opContainer, "ref", ownerLabel),
                                  message);
    }
    throw KleeError({fieldPath(opContainer, "ref"), message});
  }
  auto referencedOperator = opIter->second;
  bool startUnitsMatch = startProperties.units == referencedOperator->getStartProperties().units;
  bool endUnitsMatch = startProperties.units == referencedOperator->getEndProperties().units;

  if(startUnitsMatch && endUnitsMatch)
  {
    return referencedOperator;
  }

  auto compositeWithConversions = std::make_shared<CompositeOperator>(startProperties);
  if(!startUnitsMatch)
  {
    compositeWithConversions->addOperator(
      std::make_shared<UnitConverter>(referencedOperator->getStartProperties().units,
                                      startProperties));
  }
  compositeWithConversions->addOperator(referencedOperator);
  if(!endUnitsMatch)
  {
    compositeWithConversions->addOperator(
      std::make_shared<UnitConverter>(startProperties.units, referencedOperator->getEndProperties()));
  }
  return compositeWithConversions;
}

/**
 * Convert a single operator.
 *
 * \param data the data from which to convert the operator
 * \param startProperties the properties before the operator
 * \param namedOperators a map of named operators from which to get referenced operators
 * \return the created operator
 * \throws KleeError if the operator type or fields are invalid
 */
OpPtr convertOperator(SingleOperatorData const& data,
                      TransformableGeometryProperties startProperties,
                      const NamedOperatorMap &namedOperators,
                      const std::string &ownerLabel)
{
  std::unordered_map<std::string, OperatorParser> parsers {
    {"translate", parseTranslate},
    {"rotate", parseRotate},
    {"slice", parseSlice},
    {"scale", parseScale},
    {"convert_units_to", parseConvertUnits},
    {"ref",
     [&namedOperators](const SingleOperatorData &opData,
                       const TransformableGeometryProperties &startProperties,
                       const std::string &ownerLabel) {
       return parseRef(opData, startProperties, namedOperators, ownerLabel);
     }},
  };

  for(auto& entry : parsers)
  {
    if(containsFieldOrCallback(*data.m_container, entry.first.c_str()))
    {
      return entry.second(data, startProperties, ownerLabel);
    }
  }

  auto childNames = getChildNames(*data.m_container);
  std::string message = axom::fmt::format("Invalid transformation at {}", data.m_container->name());
  if(!childNames.empty())
  {
    message += ". Found parameters:";
    for(const auto &name : childNames)
    {
      message += " ";
      message += name;
    }
  }
  throw KleeError({data.m_container->name(), message});
}

}  // namespace

GeometryOperatorData::GeometryOperatorData(const Path& path)
  : m_path {path}
  , m_singleOperatorData {}
{ }

GeometryOperatorData::GeometryOperatorData(const Path& path,
                                           std::vector<SingleOperatorData>&& singleOperatorData)
  : m_path {path}
  , m_singleOperatorData {std::move(singleOperatorData)}
{ }

inlet::Container &GeometryOperatorData::defineSchema(inlet::Container &parent,
                                                     const std::string &fieldName,
                                                     const std::string &description,
                                                     bool enableLuaCallbacks)
{
  auto &opContainer = parent.addStructArray(fieldName, description).strict();
  auto &slice = opContainer.addStruct("slice");

  if(enableLuaCallbacks)
  {
    const auto addCallbackAlternative =
      [](inlet::Container &container, const char *fieldName, inlet::FunctionTag returnType) {
        container.addFunctionAsValueAlternative(
          fieldName,
          returnType,
          {});
      };

    addCallbackAlternative(opContainer, "translate", inlet::FunctionTag::Vector);
    addCallbackAlternative(opContainer, "rotate", inlet::FunctionTag::Double);
    addCallbackAlternative(opContainer, "center", inlet::FunctionTag::Vector);
    addCallbackAlternative(opContainer, "axis", inlet::FunctionTag::Vector);
    addCallbackAlternative(opContainer, "scale", inlet::FunctionTag::Vector);
    addCallbackAlternative(opContainer, "convert_units_to", inlet::FunctionTag::String);
    addCallbackAlternative(opContainer, "ref", inlet::FunctionTag::String);

    addCallbackAlternative(slice, "x", inlet::FunctionTag::Double);
    addCallbackAlternative(slice, "y", inlet::FunctionTag::Double);
    addCallbackAlternative(slice, "z", inlet::FunctionTag::Double);
    addCallbackAlternative(slice, "origin", inlet::FunctionTag::Vector);
    addCallbackAlternative(slice, "normal", inlet::FunctionTag::Vector);
    addCallbackAlternative(slice, "up", inlet::FunctionTag::Vector);
  }

  opContainer.addDoubleArray("translate");

  opContainer.addDouble("rotate");
  opContainer.addDoubleArray("center");
  opContainer.addDoubleArray("axis");

  opContainer.addDoubleArray("scale");

  opContainer.addString("convert_units_to");

  slice.addDouble("x");
  slice.addDouble("y");
  slice.addDouble("z");
  slice.addDoubleArray("origin");
  slice.addDoubleArray("normal");
  slice.addDoubleArray("up");

  opContainer.addString("ref");
  return opContainer;
}

std::shared_ptr<GeometryOperator> GeometryOperatorData::makeOperator(
  const TransformableGeometryProperties &startProperties,
  const NamedOperatorMap &namedOperators,
  const std::string &ownerLabel) const
{
  if(m_singleOperatorData.empty())
  {
    return nullptr;
  }
  if(startProperties.units == LengthUnit::unspecified)
  {
    throw KleeError({m_singleOperatorData[0].m_container->name(),
                     "Cannot specify operators without specifying units"});
  }
  auto composite = std::make_shared<CompositeOperator>(startProperties);
  for(auto& data : m_singleOperatorData)
  {
    composite->addOperator(
      convertOperator(data, composite->getEndProperties(), namedOperators, ownerLabel));
  }
  return composite;
}

void NamedOperatorData::defineSchema(inlet::Container &container, bool enableLuaCallbacks)
{
  container.addString("name").required();
  defineDimensionsField(container, "start_dimensions", "The initial dimensions of the operator");
  defineUnitsSchema(container,
                    "The units (both start and end) of the operator",
                    "The start units of the operator",
                    "The end units of the operator");
  GeometryOperatorData::defineSchema(container,
                                     "value",
                                     "The operation to apply",
                                     enableLuaCallbacks);  //.required();
}

NamedOperatorMapData::NamedOperatorMapData(std::vector<NamedOperatorData>&& operatorData)
  : m_operatorData {operatorData}
{ }

void NamedOperatorMapData::defineSchema(inlet::Container &parent,
                                        const std::string &name,
                                        bool enableLuaCallbacks)
{
  auto &container = parent.addStructArray(name);
  NamedOperatorData::defineSchema(container, enableLuaCallbacks);
}

NamedOperatorMap NamedOperatorMapData::makeNamedOperatorMap(Dimensions fileDimensions) const
{
  NamedOperatorMap namedOperators;

  for(auto& opData : m_operatorData)
  {
    Dimensions dimensions = fileDimensions;
    if(opData.startDimsSet)
    {
      dimensions = opData.startDims;
    }

    TransformableGeometryProperties startProperties {
      dimensions,
      opData.startUnits,
    };
    const auto ownerLabel = axom::fmt::format("named operator '{}'", opData.name);
    auto op = opData.value.makeOperator(startProperties, namedOperators, ownerLabel);

    if(op->getEndProperties().units != opData.endUnits)
    {
      throw KleeError({opData.value.getPath(), "Specified end units did not match actual units"});
    }
    namedOperators.insert({opData.name, op});
  }
  return namedOperators;
}

}  // namespace internal
}  // namespace klee
}  // namespace axom

template <>
struct FromInlet<axom::klee::internal::SingleOperatorData>
{
  axom::klee::internal::SingleOperatorData operator()(const axom::inlet::Container& base)
  {
    return axom::klee::internal::SingleOperatorData {&base};
  }
};

axom::klee::internal::GeometryOperatorData
FromInlet<axom::klee::internal::GeometryOperatorData>::operator()(const axom::inlet::Container& base)
{
  std::vector<axom::klee::internal::SingleOperatorData> v =
    base.get<std::vector<axom::klee::internal::SingleOperatorData>>();
  return axom::klee::internal::GeometryOperatorData {base.name(), std::move(v)};
}

axom::klee::internal::NamedOperatorData FromInlet<axom::klee::internal::NamedOperatorData>::operator()(
  const axom::inlet::Container& base)
{
  axom::klee::internal::NamedOperatorData data;
  std::tie(data.startUnits, data.endUnits) = axom::klee::internal::getStartAndEndUnits(base);
  data.name = base["name"];
  data.value = base["value"].get<axom::klee::internal::GeometryOperatorData>();
  if(base.contains("start_dimensions"))
  {
    data.startDimsSet = true;
    data.startDims = axom::klee::internal::toDimensions(base["start_dimensions"]);
  }
  else
  {
    data.startDimsSet = false;
  }
  return data;
}

axom::klee::internal::NamedOperatorMapData
FromInlet<axom::klee::internal::NamedOperatorMapData>::operator()(const axom::inlet::Container& base)
{
  return axom::klee::internal::NamedOperatorMapData {
    base.get<std::vector<axom::klee::internal::NamedOperatorData>>()};
}
