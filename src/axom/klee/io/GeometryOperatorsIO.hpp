// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEE_GEOMETRYOPERATORSIO_HPP_
#define AXOM_KLEE_GEOMETRYOPERATORSIO_HPP_

#include "axom/inlet.hpp"
#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Units.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace axom
{
namespace klee
{
class GeometryOperator;
struct TransformableGeometryProperties;

namespace internal
{
using NamedOperatorMap = std::unordered_map<std::string, std::shared_ptr<const GeometryOperator>>;

/// The data for a single operator.
struct SingleOperatorData
{
  const inlet::Container *m_container;
};

/// The data for the "operator" component of "geometry" objects.
class GeometryOperatorData
{
public:
  /// Construct a GeometryOperatorData with no operators.
  GeometryOperatorData() = default;

  /**
   * Construct a GeometryOperatorData with no operators.
   * @param path the path where the operators were defined
   */
  explicit GeometryOperatorData(const Path &path);

  /**
   * Construct a GeometryOperatorData for the given list of operators
   * @param path the path where the operators were defined
   * @param singleOperatorData the data for the individual operators
   */
  explicit GeometryOperatorData(const Path &path,
                                std::vector<SingleOperatorData> &&singleOperatorData);

  /**
   * Define the schema for geometry operators
   * @param parent the parent container
   * @param fieldName the name of the field
   * @param description a description of the field
   * @return the Container for the new item
   */
  static inlet::Container &defineSchema(inlet::Container &parent,
                                        const std::string &fieldName,
                                        const std::string &description);

  /**
   * Make a (possibly null) operator describing the transformation to apply to the geometry
   *
   * @param startProperties properties of the geometry before the first operator
   * @param namedOperators a map of any named operators
   * @return the (possibly null) operator
   */
  std::shared_ptr<GeometryOperator> makeOperator(const TransformableGeometryProperties &startProperties,
                                                 const NamedOperatorMap &namedOperators) const;

  /**
   * Get the path of this operator in the source document
   * @return the operator's path
   */
  const Path &getPath() const { return m_path; }

private:
  Path m_path;
  std::vector<SingleOperatorData> m_singleOperatorData;
};

/// Data for a named operator.
struct NamedOperatorData
{
  std::string name;
  LengthUnit startUnits;
  LengthUnit endUnits;
  bool startDimsSet;
  Dimensions startDims;
  GeometryOperatorData value;

  /**
   * Define the schema for a named operator.
   *
   * @param container the container in which to describe a single named operator
   */
  static void defineSchema(inlet::Container &container);
};

/// Data for all a collection of named operators
struct NamedOperatorMapData
{
  /// Create a NamedOperatorMapData with no operators.
  NamedOperatorMapData() = default;

  /**
   * Create a NamedOperatorMapData with the given list of operators.
   *
   * @param operatorData the data for all the named operators in this map
   */
  explicit NamedOperatorMapData(std::vector<NamedOperatorData> &&operatorData);

  /**
   * Convert the data to a NamedOperatorMap.
   *
   * @param fileDimensions the dimensions that shapes should be in in this file.
   * @return the name of converted operators
   */
  NamedOperatorMap makeNamedOperatorMap(Dimensions fileDimensions) const;

  /**
   * Define the schema for a collection of named operators.
   *
   * @param parent the parent object in which to define the operator map
   * @param name the name of the map
   */
  static void defineSchema(inlet::Container &parent, const std::string &name);

private:
  std::vector<NamedOperatorData> m_operatorData;
};

}  // namespace internal
}  // namespace klee
}  // namespace axom

template <>
struct FromInlet<axom::klee::internal::GeometryOperatorData>
{
  axom::klee::internal::GeometryOperatorData operator()(const axom::inlet::Container &base);
};

template <>
struct FromInlet<axom::klee::internal::NamedOperatorData>
{
  axom::klee::internal::NamedOperatorData operator()(const axom::inlet::Container &base);
};

template <>
struct FromInlet<axom::klee::internal::NamedOperatorMapData>
{
  axom::klee::internal::NamedOperatorMapData operator()(const axom::inlet::Container &base);
};

#endif  // AXOM_KLEE_GEOMETRYOPERATORSIO_HPP_
