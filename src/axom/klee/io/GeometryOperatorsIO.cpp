// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "GeometryOperatorsIO.hpp"
#include "IOUtil.hpp"

#include "axom/klee/Geometry.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"

#include "axom/fmt.hpp"

#include <functional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

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
  std::function<OpPtr(const inlet::Container &, const TransformableGeometryProperties &)>;
using internal::toDoubleVector;
using primal::Point3D;
using primal::Vector3D;
using FieldSet = std::unordered_set<std::string>;

/**
 * Get the names of all the children in the given container.
 *
 * @param container the Container whose children to get
 * @return the names of all the children
 */
std::unordered_set<std::string> getChildNames(const inlet::Container &container)
{
  std::unordered_set<std::string> allChildren;

  std::vector<std::string> unexpectedNames = container.unexpectedNames();
  allChildren.insert(unexpectedNames.begin(), unexpectedNames.end());

  // Add 1 for the "/" separator
  auto prefixLength = container.name().size() + 1;

  for(auto &child : container.getChildContainers())
  {
    if(child.second->exists())
    {
      allChildren.insert(child.first.substr(prefixLength));
    }
  }

  for(auto &child : container.getChildFields())
  {
    if(child.second->exists())
    {
      allChildren.insert(child.first.substr(prefixLength));
    }
  }

  for(auto &child : container.getChildFunctions())
  {
    if(*child.second)
    {
      allChildren.insert(child.first.substr(prefixLength));
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
void verifyObjectFields(const inlet::Container &containerToTest,
                        const std::string &name,
                        const FieldSet &additionalRequiredFields,
                        const FieldSet &optionalFields)
{
  std::unordered_set<std::string> requiredFields {additionalRequiredFields};
  requiredFields.insert(name);

  for(auto &requiredField : requiredFields)
  {
    if(!containerToTest.contains(requiredField))
    {
      throw KleeError(
        {containerToTest.name(),
         axom::fmt::format("Missing required parameter '{}' for operator '{}'", requiredField, name)});
    }
  }

  for(auto &child : getChildNames(containerToTest))
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
OpPtr parseTranslate(const inlet::Container &opContainer,
                     const TransformableGeometryProperties &startProperties)
{
  verifyObjectFields(opContainer, "translate", FieldSet {}, FieldSet {});
  return std::make_shared<Translation>(toVector(opContainer, "translate", startProperties.dimensions),
                                       startProperties);
}

/**
 * Parse a "rotate" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseRotate(const inlet::Container &opContainer,
                  const TransformableGeometryProperties &startProperties)
{
  switch(startProperties.dimensions)
  {
  case Dimensions::Two:
  {
    verifyObjectFields(opContainer, "rotate", FieldSet {}, {"center"});
    Vector3D axis {0, 0, 1};
    return std::make_shared<Rotation>(
      opContainer["rotate"].get<double>(),
      toPoint(opContainer, "center", Dimensions::Two, Point3D {0, 0, 0}),
      axis,
      startProperties);
  }
  break;
  case Dimensions::Three:
  {
    verifyObjectFields(opContainer, "rotate", {"axis"}, {"center"});
    return std::make_shared<Rotation>(
      opContainer["rotate"].get<double>(),
      toPoint(opContainer, "center", Dimensions::Three, Point3D {0, 0, 0}),
      toVector(opContainer, "axis", Dimensions::Three),
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
                       const TransformableGeometryProperties &startProperties,
                       const Path &path)
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
primal::Point3D getPerpendicularSliceOrigin(const inlet::Container &sliceContainer,
                                            char const *planeName,
                                            const primal::Vector3D &defaultNormal)
{
  double axisIntercept = sliceContainer[planeName];

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

  if(!sliceContainer.contains("origin"))
  {
    return defaultOrigin;
  }

  primal::Point3D givenOrigin = toPoint(sliceContainer, "origin", Dimensions::Three);
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
primal::Vector3D getPerpendicularSliceNormal(const inlet::Container &sliceContainer,
                                             const primal::Vector3D &defaultNormal)
{
  if(!sliceContainer.contains("normal"))
  {
    return defaultNormal;
  }

  primal::Vector3D givenNormal = toVector(sliceContainer, "normal", Dimensions::Three);
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
OpPtr readPerpendicularSlice(const inlet::Container &sliceContainer,
                             char const *planeName,
                             Vector3D const &defaultNormal,
                             Vector3D const &defaultUp,
                             const TransformableGeometryProperties &startProperties)
{
  verifyObjectFields(sliceContainer, planeName, FieldSet {}, {"origin", "normal", "up"});
  const primal::Vector3D defaultNormalVec {defaultNormal.data()};

  auto origin = getPerpendicularSliceOrigin(sliceContainer, planeName, defaultNormalVec);
  auto normal = getPerpendicularSliceNormal(sliceContainer, defaultNormalVec);
  auto up = toVector(sliceContainer, "up", Dimensions::Three, defaultUp);

  return makeCheckedSlice(origin, normal, up, startProperties, sliceContainer.name());
}

/**
 * Parse a "slice" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseSlice(const inlet::Container &opContainer,
                 const TransformableGeometryProperties &startProperties)
{
  if(startProperties.dimensions != Dimensions::Three)
  {
    throw KleeError({opContainer.name(), "Cannot do a slice from 2D"});
  }
  verifyObjectFields(opContainer, "slice", FieldSet {}, FieldSet {});
  auto &sliceContainer = *opContainer.getChildContainers().at(opContainer.name() + "/slice").get();
  if(sliceContainer.contains("x"))
  {
    return readPerpendicularSlice(sliceContainer, "x", {1, 0, 0}, {0, 0, 1}, startProperties);
  }
  else if(sliceContainer.contains("y"))
  {
    return readPerpendicularSlice(sliceContainer, "y", {0, 1, 0}, {1, 0, 0}, startProperties);
  }
  else if(sliceContainer.contains("z"))
  {
    return readPerpendicularSlice(sliceContainer, "z", {0, 0, 1}, {0, 1, 0}, startProperties);
  }

  verifyObjectFields(sliceContainer, "origin", {"normal", "up"}, FieldSet {});

  return makeCheckedSlice(toPoint(sliceContainer, "origin", Dimensions::Three),
                          toVector(sliceContainer, "normal", Dimensions::Three),
                          toVector(sliceContainer, "up", Dimensions::Three),
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
OpPtr parseScale(const inlet::Container &opContainer,
                 const TransformableGeometryProperties &startProperties)
{
  verifyObjectFields(opContainer, "scale", FieldSet {}, FieldSet {});
  auto factors = opContainer["scale"].get<std::vector<double>>();
  if(factors.size() == 1)
  {
    return std::make_shared<Scale>(factors[0], factors[0], factors[0], startProperties);
  }
  factors = toDoubleVector(opContainer["scale"], startProperties.dimensions, "scale");
  if(startProperties.dimensions == Dimensions::Two)
  {
    factors.emplace_back(1.0);
  }
  return std::make_shared<Scale>(factors[0], factors[1], factors[2], startProperties);
}

/**
 * Parse a "convert_units_to" operator.
 *
 * \param opContainer the Container from which to read the operator
 * \param startProperties the properties prior to this operator
 * \return the created operator
 */
OpPtr parseConvertUnits(const inlet::Container &opContainer,
                        const TransformableGeometryProperties &startProperties)
{
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
OpPtr parseRef(const inlet::Container &opContainer,
               const TransformableGeometryProperties &startProperties,
               const NamedOperatorMap &namedOperators)
{
  verifyObjectFields(opContainer, "ref", FieldSet {}, FieldSet {});
  std::string const &operatorName = opContainer["ref"];
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
OpPtr convertOperator(SingleOperatorData const &data,
                      TransformableGeometryProperties startProperties,
                      const NamedOperatorMap &namedOperators)
{
  std::unordered_map<std::string, OperatorParser> parsers {
    {"translate", parseTranslate},
    {"rotate", parseRotate},
    {"slice", parseSlice},
    {"scale", parseScale},
    {"convert_units_to", parseConvertUnits},
    {"ref",
     [&namedOperators](const inlet::Container &opNode,
                       const TransformableGeometryProperties &startProperties) {
       return parseRef(opNode, startProperties, namedOperators);
     }},
  };

  for(auto &entry : parsers)
  {
    if(data.m_container->contains(entry.first))
    {
      return entry.second(*data.m_container, startProperties);
    }
  }

  throw KleeError({data.m_container->name(), "Invalid transformation"});
}

}  // namespace

GeometryOperatorData::GeometryOperatorData(const Path &path)
  : m_path {path}
  , m_singleOperatorData {}
{ }

GeometryOperatorData::GeometryOperatorData(const Path &path,
                                           std::vector<SingleOperatorData> &&singleOperatorData)
  : m_path {path}
  , m_singleOperatorData {singleOperatorData}
{ }

inlet::Container &GeometryOperatorData::defineSchema(inlet::Container &parent,
                                                     const std::string &fieldName,
                                                     const std::string &description)
{
  auto &opContainer = parent.addStructArray(fieldName, description).strict();

  opContainer.addDoubleArray("translate");

  opContainer.addDouble("rotate");
  opContainer.addDoubleArray("center");
  opContainer.addDoubleArray("axis");

  opContainer.addDoubleArray("scale");

  opContainer.addString("convert_units_to");

  auto &slice = opContainer.addStruct("slice");
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
  const NamedOperatorMap &namedOperators) const
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
  for(auto &data : m_singleOperatorData)
  {
    composite->addOperator(convertOperator(data, composite->getEndProperties(), namedOperators));
  }
  return composite;
}

void NamedOperatorData::defineSchema(inlet::Container &container)
{
  container.addString("name").required();
  defineDimensionsField(container, "start_dimensions", "The initial dimensions of the operator");
  defineUnitsSchema(container,
                    "The units (both start and end) of the operator",
                    "The start units of the operator",
                    "The end units of the operator");
  GeometryOperatorData::defineSchema(container, "value",
                                     "The operation to apply");  //.required();
}

NamedOperatorMapData::NamedOperatorMapData(std::vector<NamedOperatorData> &&operatorData)
  : m_operatorData {operatorData}
{ }

void NamedOperatorMapData::defineSchema(inlet::Container &parent, const std::string &name)
{
  auto &container = parent.addStructArray(name);
  NamedOperatorData::defineSchema(container);
}

NamedOperatorMap NamedOperatorMapData::makeNamedOperatorMap(Dimensions fileDimensions) const
{
  NamedOperatorMap namedOperators;

  for(auto &opData : m_operatorData)
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
  axom::klee::internal::SingleOperatorData operator()(const axom::inlet::Container &base)
  {
    return axom::klee::internal::SingleOperatorData {&base};
  }
};

axom::klee::internal::GeometryOperatorData
FromInlet<axom::klee::internal::GeometryOperatorData>::operator()(const axom::inlet::Container &base)
{
  std::vector<axom::klee::internal::SingleOperatorData> v =
    base.get<std::vector<axom::klee::internal::SingleOperatorData>>();
  return axom::klee::internal::GeometryOperatorData {base.name(), std::move(v)};
}

axom::klee::internal::NamedOperatorData FromInlet<axom::klee::internal::NamedOperatorData>::operator()(
  const axom::inlet::Container &base)
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
FromInlet<axom::klee::internal::NamedOperatorMapData>::operator()(const axom::inlet::Container &base)
{
  return axom::klee::internal::NamedOperatorMapData {
    base.get<std::vector<axom::klee::internal::NamedOperatorData>>()};
}
