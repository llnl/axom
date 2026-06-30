// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "GeometryOperatorsIO.hpp"
#include "IOUtil.hpp"

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"

#include "axom/fmt.hpp"

#include <functional>
#include <stdexcept>
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
  std::function<OpPtr(const SingleOperatorData&, const TransformableGeometryProperties&)>;
using internal::toDoubleVector;
using primal::Point3D;
using primal::Vector3D;
using FieldSet = std::unordered_set<std::string>;

// Callback schema entries are internal aliases: the public input path remains
// the ordinary operator field name, while Inlet stores the function separately.
constexpr char const* LUA_CALLBACK_SUFFIX = "__klee_lua_callback";

std::string callbackName(char const* fieldName)
{
  return std::string(fieldName) + LUA_CALLBACK_SUFFIX;
}

std::string publicNameForCallback(const std::string& childName)
{
  const std::string suffix = LUA_CALLBACK_SUFFIX;
  if(childName.size() > suffix.size() &&
     childName.compare(childName.size() - suffix.size(), suffix.size(), suffix) == 0)
  {
    return childName.substr(0, childName.size() - suffix.size());
  }
  return childName;
}

bool hasCallback(const inlet::Container& container, char const* fieldName)
{
  const auto name = callbackName(fieldName);
  return container.contains(name);
}

bool containsFieldOrCallback(const inlet::Container& container, char const* fieldName)
{
  return container.contains(fieldName) || hasCallback(container, fieldName);
}

Path fieldPath(const inlet::Container& container, char const* fieldName)
{
  return Path::join({Path {container.name()}, Path {std::string {fieldName}}});
}

std::string callbackContext(const inlet::Container& container,
                            char const* fieldName,
                            const std::string& shapeName)
{
  Path path {container.name()};
  std::string operatorIndex = path.baseName();
  if(operatorIndex == "slice")
  {
    operatorIndex = path.parent().baseName();
  }

  const auto operatorLabel = operatorIndex.empty() ? std::string {"operator at "} + container.name()
                                                   : std::string {"operator "} + operatorIndex;
  if(shapeName.empty())
  {
    return axom::fmt::format("Error evaluating callback for '{}' in {}", fieldName, operatorLabel);
  }
  return axom::fmt::format("Error evaluating callback for '{}' in shape '{}' {}",
                           fieldName,
                           shapeName,
                           operatorLabel);
}

template <typename Result, typename Func>
Result wrapCallbackErrors(const inlet::Container& container,
                          char const* fieldName,
                          const std::string& shapeName,
                          Func&& func)
{
  // Convert generic Inlet/Lua callback failures into Klee diagnostics at the
  // boundary where the shape, operator, and field context are all available.
  try
  {
    return func();
  }
  catch(const KleeError&)
  {
    throw;
  }
  catch(const std::exception& ex)
  {
    throw KleeError(
      {fieldPath(container, fieldName),
       axom::fmt::format("{}: {}", callbackContext(container, fieldName, shapeName), ex.what())});
  }
}

double getScalar(const inlet::Container& container, char const* fieldName, const std::string& shapeName)
{
  if(hasCallback(container, fieldName))
  {
    return wrapCallbackErrors<double>(container, fieldName, shapeName, [&]() {
      return container[callbackName(fieldName)].call<inlet::FunctionType::Double>();
    });
  }
  return container[fieldName].get<double>();
}

std::vector<double> callbackVectorToDoubleVector(const inlet::FunctionType::Vector& value)
{
  std::vector<double> result;
  result.reserve(value.dim);
  for(int i = 0; i < value.dim; ++i)
  {
    result.push_back(value.vec[i]);
  }
  return result;
}

Point3D inletVectorToPoint(const inlet::FunctionType::Vector& value, Dimensions expectedDims)
{
  const auto expectedSize = static_cast<int>(expectedDims);
  if(value.dim != expectedSize)
  {
    throw std::runtime_error(
      fmt::format("Wrong size for transform. Expected {}. Got {}.", expectedSize, value.dim));
  }
  return Point3D {value.vec[0], value.vec[1], expectedDims == Dimensions::Three ? value.vec[2] : 0.};
}

std::vector<double> getDoubleVector(const inlet::Container& container,
                                    char const* fieldName,
                                    Dimensions expectedDims,
                                    const std::string& shapeName)
{
  if(hasCallback(container, fieldName))
  {
    auto values = wrapCallbackErrors<std::vector<double>>(container, fieldName, shapeName, [&]() {
      return callbackVectorToDoubleVector(
        container[callbackName(fieldName)].call<inlet::FunctionType::Vector>());
    });
    auto actualSize = values.size();
    auto expectedSize = static_cast<std::size_t>(expectedDims);
    if(actualSize != expectedSize)
    {
      throw KleeError({fieldPath(container, fieldName),
                       fmt::format("{}: Wrong size for {}. Expected {}. Got {}.",
                                   callbackContext(container, fieldName, shapeName),
                                   fieldName,
                                   expectedSize,
                                   actualSize)});
    }
    return values;
  }
  return toDoubleVector(container[fieldName], expectedDims, fieldName);
}

template <typename T>
T toArrayLike(const inlet::Container& parent,
              char const* fieldName,
              Dimensions expectedDims,
              const std::string& shapeName)
{
  auto values = getDoubleVector(parent, fieldName, expectedDims, shapeName);
  return T {values.data(), static_cast<int>(expectedDims)};
}

template <typename T>
T toArrayLike(const inlet::Container& parent,
              char const* fieldName,
              Dimensions expectedDims,
              const T& defaultValue,
              const std::string& shapeName)
{
  if(containsFieldOrCallback(parent, fieldName))
  {
    return toArrayLike<T>(parent, fieldName, expectedDims, shapeName);
  }
  return defaultValue;
}

Point3D getPoint(const inlet::Container& parent,
                 char const* fieldName,
                 Dimensions expectedDims,
                 const std::string& shapeName)
{
  return toArrayLike<Point3D>(parent, fieldName, expectedDims, shapeName);
}

Point3D getPoint(const inlet::Container& parent,
                 char const* fieldName,
                 Dimensions expectedDims,
                 const Point3D& defaultValue,
                 const std::string& shapeName)
{
  return toArrayLike(parent, fieldName, expectedDims, defaultValue, shapeName);
}

Vector3D getVector(const inlet::Container& parent,
                   char const* fieldName,
                   Dimensions expectedDims,
                   const std::string& shapeName)
{
  return toArrayLike<Vector3D>(parent, fieldName, expectedDims, shapeName);
}

Vector3D getVector(const inlet::Container& parent,
                   char const* fieldName,
                   Dimensions expectedDims,
                   const Vector3D& defaultValue,
                   const std::string& shapeName)
{
  return toArrayLike(parent, fieldName, expectedDims, defaultValue, shapeName);
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

  // Add 1 for the "/" separator
  auto prefixLength = container.name().size() + 1;

  for(auto& child : container.getChildContainers())
  {
    if(child.second->exists())
    {
      allChildren.insert(child.first.substr(prefixLength));
    }
  }

  for(auto& child : container.getChildFields())
  {
    if(child.second->exists())
    {
      allChildren.insert(child.first.substr(prefixLength));
    }
  }

  for(auto& child : container.getChildFunctions())
  {
    if(*child.second)
    {
      allChildren.insert(publicNameForCallback(child.first.substr(prefixLength)));
    }
  }

  return allChildren;
}

/**
 * Verify that a Container has the correct fields.
 *
 * While Inlet can do a lot of the verification, there are some situations
 * where we need some extra manual checks. For example, operators can look
 * like this:
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
 */
OpPtr parseTranslate(const SingleOperatorData& data,
                     const TransformableGeometryProperties& startProperties)
{
  const auto& opContainer = *data.m_container;
  verifyObjectFields(opContainer, "translate", FieldSet {}, FieldSet {});
  return std::make_shared<Translation>(
    getVector(opContainer, "translate", startProperties.dimensions, data.m_shapeName),
    startProperties);
}

/**
 * Parse a "rotate" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseRotate(const SingleOperatorData& data,
                  const TransformableGeometryProperties& startProperties)
{
  const auto& opContainer = *data.m_container;
  switch(startProperties.dimensions)
  {
  case Dimensions::Two:
  {
    verifyObjectFields(opContainer, "rotate", FieldSet {}, {"center"});
    Vector3D axis {0, 0, 1};
    return std::make_shared<Rotation>(
      getScalar(opContainer, "rotate", data.m_shapeName),
      getPoint(opContainer, "center", Dimensions::Two, Point3D {0, 0, 0}, data.m_shapeName),
      axis,
      startProperties);
  }
  break;
  case Dimensions::Three:
  {
    verifyObjectFields(opContainer, "rotate", {"axis"}, {"center"});
    return std::make_shared<Rotation>(
      getScalar(opContainer, "rotate", data.m_shapeName),
      getPoint(opContainer, "center", Dimensions::Three, Point3D {0, 0, 0}, data.m_shapeName),
      getVector(opContainer, "axis", Dimensions::Three, data.m_shapeName),
      startProperties);
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
 * \param path the path where the slice is specified, for error reporting
 * \return the created operator
 * \throws KleeError if any value is invalid
 */
OpPtr makeCheckedSlice(Point3D origin,
                       Vector3D normal,
                       Vector3D up,
                       const TransformableGeometryProperties& startProperties,
                       const Path& path)
{
  if(normal.is_zero())
  {
    throw KleeError({path, "The 'normal' vector must not be a zero vector"});
  }
  if(!utilities::isNearlyEqual(normal.dot(up), 0.))
  {
    throw KleeError({path, "The 'normal' and 'up' vectors must be perpendicular"});
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
 */
primal::Point3D getPerpendicularSliceOrigin(const inlet::Container& sliceContainer,
                                            char const* planeName,
                                            const primal::Vector3D& defaultNormal,
                                            const std::string& shapeName)
{
  double axisIntercept = getScalar(sliceContainer, planeName, shapeName);

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

  primal::Point3D givenOrigin = getPoint(sliceContainer, "origin", Dimensions::Three, shapeName);
  if(givenOrigin[nonZeroIndex] != axisIntercept)
  {
    throw KleeError({sliceContainer["origin"].name(), "The origin must be on the slice plane"});
  }
  return givenOrigin;
}

/**
 * Get the normal vector to use for a perpendicular slice.
 *
 * \param sliceContainer the Container describing the slice
 * \param defaultNormal the default normal vector
 * \return the vector to use as the normal
 */
primal::Vector3D getPerpendicularSliceNormal(const inlet::Container& sliceContainer,
                                             const primal::Vector3D& defaultNormal,
                                             const std::string& shapeName)
{
  if(!containsFieldOrCallback(sliceContainer, "normal"))
  {
    return defaultNormal;
  }

  primal::Vector3D givenNormal = getVector(sliceContainer, "normal", Dimensions::Three, shapeName);
  auto cross = primal::Vector3D::cross_product(givenNormal, defaultNormal);
  bool parallel = cross.is_zero();
  if(!parallel)
  {
    throw KleeError({sliceContainer["normal"].name(), "Invalid normal"});
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
 */
OpPtr readPerpendicularSlice(const inlet::Container& sliceContainer,
                             char const* planeName,
                             Vector3D const& defaultNormal,
                             Vector3D const& defaultUp,
                             const TransformableGeometryProperties& startProperties,
                             const std::string& shapeName)
{
  verifyObjectFields(sliceContainer, planeName, FieldSet {}, {"origin", "normal", "up"});
  const primal::Vector3D defaultNormalVec {defaultNormal.data()};

  auto origin = getPerpendicularSliceOrigin(sliceContainer, planeName, defaultNormalVec, shapeName);
  auto normal = getPerpendicularSliceNormal(sliceContainer, defaultNormalVec, shapeName);
  auto up = getVector(sliceContainer, "up", Dimensions::Three, defaultUp, shapeName);

  return makeCheckedSlice(origin, normal, up, startProperties, sliceContainer.name());
}

/**
 * Parse a "slice" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseSlice(const SingleOperatorData& data,
                 const TransformableGeometryProperties& startProperties)
{
  const auto& opContainer = *data.m_container;
  if(startProperties.dimensions != Dimensions::Three)
  {
    throw KleeError({opContainer.name(), "Cannot do a slice from 2D"});
  }
  verifyObjectFields(opContainer, "slice", FieldSet {}, FieldSet {});
  auto& sliceContainer = *opContainer.getChildContainers().at(opContainer.name() + "/slice").get();
  if(containsFieldOrCallback(sliceContainer, "x"))
  {
    return readPerpendicularSlice(sliceContainer,
                                  "x",
                                  {1, 0, 0},
                                  {0, 0, 1},
                                  startProperties,
                                  data.m_shapeName);
  }
  else if(containsFieldOrCallback(sliceContainer, "y"))
  {
    return readPerpendicularSlice(sliceContainer,
                                  "y",
                                  {0, 1, 0},
                                  {1, 0, 0},
                                  startProperties,
                                  data.m_shapeName);
  }
  else if(containsFieldOrCallback(sliceContainer, "z"))
  {
    return readPerpendicularSlice(sliceContainer,
                                  "z",
                                  {0, 0, 1},
                                  {0, 1, 0},
                                  startProperties,
                                  data.m_shapeName);
  }

  verifyObjectFields(sliceContainer, "origin", {"normal", "up"}, FieldSet {});

  return makeCheckedSlice(getPoint(sliceContainer, "origin", Dimensions::Three, data.m_shapeName),
                          getVector(sliceContainer, "normal", Dimensions::Three, data.m_shapeName),
                          getVector(sliceContainer, "up", Dimensions::Three, data.m_shapeName),
                          startProperties,
                          sliceContainer.name());
}

/**
 * Parse a "scale" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseScale(const SingleOperatorData& data,
                 const TransformableGeometryProperties& startProperties)
{
  const auto& opContainer = *data.m_container;
  verifyObjectFields(opContainer, "scale", FieldSet {}, FieldSet {"center"});
  auto factors = hasCallback(opContainer, "scale")
    ? wrapCallbackErrors<std::vector<double>>(
        opContainer,
        "scale",
        data.m_shapeName,
        [&]() {
          return callbackVectorToDoubleVector(
            opContainer[callbackName("scale")].call<inlet::FunctionType::Vector>());
        })
    : opContainer["scale"].get<std::vector<double>>();
  if(factors.size() == 1)
  {
    return std::make_shared<Scale>(factors[0], factors[0], factors[0], startProperties);
  }
  if(hasCallback(opContainer, "scale"))
  {
    auto actualSize = factors.size();
    auto expectedSize = static_cast<std::size_t>(startProperties.dimensions);
    if(actualSize != expectedSize)
    {
      throw KleeError({fieldPath(opContainer, "scale"),
                       fmt::format("{}: Wrong size for scale. Expected {}. Got {}.",
                                   callbackContext(opContainer, "scale", data.m_shapeName),
                                   expectedSize,
                                   actualSize)});
    }
  }
  else
  {
    factors = toDoubleVector(opContainer["scale"], startProperties.dimensions, "scale");
  }
  if(startProperties.dimensions == Dimensions::Two)
  {
    factors.emplace_back(1.0);
  }
  Point3D center {0., 0., 0.};
  if(containsFieldOrCallback(opContainer, "center"))
  {
    center =
      getPoint(opContainer, "center", startProperties.dimensions, Point3D {0, 0, 0}, data.m_shapeName);
  }

  return std::make_shared<Scale>(factors[0], factors[1], factors[2], center, startProperties);
}

/**
 * Parse a point transform.
 *
 * \param data the operator data to parse
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseTransform(const SingleOperatorData& data,
                     const TransformableGeometryProperties& startProperties)
{
  const auto& opContainer = *data.m_container;
  verifyObjectFields(opContainer, "transform", FieldSet {}, FieldSet {});
  auto transform = opContainer[callbackName("transform")]
                     .get<std::function<inlet::FunctionType::Vector(inlet::FunctionType::Vector)>>();
  const auto expectedDims = startProperties.dimensions;
  auto path = static_cast<std::string>(fieldPath(opContainer, "transform"));
  auto context = callbackContext(opContainer, "transform", data.m_shapeName);

  return std::make_shared<PointTransform>(
    [transform, expectedDims](const Point3D& point) {
      inlet::FunctionType::Vector input {point.data(), static_cast<int>(expectedDims)};
      return inletVectorToPoint(transform(input), expectedDims);
    },
    startProperties,
    path,
    context);
}

/**
 * Parse a "convert_units_to" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseConvertUnits(const SingleOperatorData& data,
                        const TransformableGeometryProperties& startProperties)
{
  const auto& opContainer = *data.m_container;
  verifyObjectFields(opContainer, "convert_units_to", FieldSet {}, FieldSet {});
  auto endUnits = parseLengthUnits(opContainer["convert_units_to"]);
  return std::make_shared<UnitConverter>(endUnits, startProperties);
}

/**
 * Parse an operator specified via the "ref" command.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties before the "ref" command
 * \param namedOperators a map of named operators from which to get referenced operators
 * \return the created operator
 */
OpPtr parseRef(const SingleOperatorData& data,
               const TransformableGeometryProperties& startProperties,
               const NamedOperatorMap& namedOperators)
{
  const auto& opContainer = *data.m_container;
  verifyObjectFields(opContainer, "ref", FieldSet {}, FieldSet {});
  std::string const& operatorName = opContainer["ref"];
  auto opIter = namedOperators.find(operatorName);
  if(opIter == namedOperators.end())
  {
    std::string message = "No operator named '";
    message += operatorName;
    message += '\'';
    throw KleeError({opContainer["ref"].name(), message});
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
 */
OpPtr convertOperator(SingleOperatorData const& data,
                      TransformableGeometryProperties startProperties,
                      const NamedOperatorMap& namedOperators)
{
  std::unordered_map<std::string, OperatorParser> parsers {
    {"translate", parseTranslate},
    {"rotate", parseRotate},
    {"slice", parseSlice},
    {"scale", parseScale},
    {"transform", parseTransform},
    {"convert_units_to", parseConvertUnits},
    {"ref",
     [&namedOperators](const SingleOperatorData& opData,
                       const TransformableGeometryProperties& startProperties) {
       return parseRef(opData, startProperties, namedOperators);
     }},
  };

  for(auto& entry : parsers)
  {
    if(containsFieldOrCallback(*data.m_container, entry.first.c_str()))
    {
      return entry.second(data, startProperties);
    }
  }

  auto childNames = getChildNames(*data.m_container);
  std::string message = axom::fmt::format("Invalid transformation at {}", data.m_container->name());
  if(!childNames.empty())
  {
    message += ". Found parameters:";
    for(const auto& name : childNames)
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
  , m_singleOperatorData {singleOperatorData}
{ }

void GeometryOperatorData::setShapeName(std::string shapeName)
{
  m_shapeName = std::move(shapeName);
  for(auto& data : m_singleOperatorData)
  {
    data.m_shapeName = m_shapeName;
  }
}

inlet::Container& GeometryOperatorData::defineSchema(inlet::Container& parent,
                                                     const std::string& fieldName,
                                                     const std::string& description,
                                                     bool enableLuaCallbacks)
{
  auto& opContainer = parent.addStructArray(fieldName, description).strict();

  opContainer.addDoubleArray("translate");

  opContainer.addDouble("rotate");
  opContainer.addDoubleArray("center");
  opContainer.addDoubleArray("axis");

  opContainer.addDoubleArray("scale");

  opContainer.addString("convert_units_to");

  auto& slice = opContainer.addStruct("slice");
  slice.addDouble("x");
  slice.addDouble("y");
  slice.addDouble("z");
  slice.addDoubleArray("origin");
  slice.addDoubleArray("normal");
  slice.addDoubleArray("up");

  opContainer.addString("ref");
  if(enableLuaCallbacks)
  {
    // These Lua-only function alternatives read from the public field paths via
    // pathOverride, leaving YAML and concrete Lua field parsing unchanged.
    opContainer.addFunction(callbackName("translate"),
                            inlet::FunctionTag::Vector,
                            {},
                            "",
                            "translate");
    opContainer.addFunction(callbackName("rotate"), inlet::FunctionTag::Double, {}, "", "rotate");
    opContainer.addFunction(callbackName("center"), inlet::FunctionTag::Vector, {}, "", "center");
    opContainer.addFunction(callbackName("axis"), inlet::FunctionTag::Vector, {}, "", "axis");
    opContainer.addFunction(callbackName("scale"), inlet::FunctionTag::Vector, {}, "", "scale");
    opContainer.addFunction(callbackName("transform"),
                            inlet::FunctionTag::Vector,
                            {inlet::FunctionTag::Vector},
                            "",
                            "transform");

    slice.addFunction(callbackName("x"), inlet::FunctionTag::Double, {}, "", "x");
    slice.addFunction(callbackName("y"), inlet::FunctionTag::Double, {}, "", "y");
    slice.addFunction(callbackName("z"), inlet::FunctionTag::Double, {}, "", "z");
    slice.addFunction(callbackName("origin"), inlet::FunctionTag::Vector, {}, "", "origin");
    slice.addFunction(callbackName("normal"), inlet::FunctionTag::Vector, {}, "", "normal");
    slice.addFunction(callbackName("up"), inlet::FunctionTag::Vector, {}, "", "up");
  }
  return opContainer;
}

std::shared_ptr<GeometryOperator> GeometryOperatorData::makeOperator(
  const TransformableGeometryProperties& startProperties,
  const NamedOperatorMap& namedOperators) const
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
    composite->addOperator(convertOperator(data, composite->getEndProperties(), namedOperators));
  }
  return composite;
}

void NamedOperatorData::defineSchema(inlet::Container& container, bool enableLuaCallbacks)
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

void NamedOperatorMapData::defineSchema(inlet::Container& parent,
                                        const std::string& name,
                                        bool enableLuaCallbacks)
{
  auto& container = parent.addStructArray(name);
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
    auto op = opData.value.makeOperator(startProperties, namedOperators);

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
    return axom::klee::internal::SingleOperatorData {&base, ""};
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
