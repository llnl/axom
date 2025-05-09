// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IO.hpp"

#include <fstream>
#include <functional>
#include <iterator>
#include <string>
#include <tuple>

#include "axom/inlet.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/GeometryOperatorsIO.hpp"
#include "axom/klee/IOUtil.hpp"
#include "axom/klee/KleeError.hpp"

namespace axom
{
namespace klee
{
namespace
{
// Because we can't have context-aware validation when extracting the
// data from Inlet, we need a set of structs that parallels the real
// classes. These are used to do some basic validation, and then we convert
// the to the real classes, doing more thorough validation.

struct GeometryData
{
  std::string format;
  std::string path;
  LengthUnit startUnits;
  LengthUnit endUnits;
  Dimensions startDimensions;
  bool dimensionsSet;
  internal::GeometryOperatorData operatorData;
  Path pathInFile;
};

struct ShapeData
{
  std::string name;
  std::string material;
  std::vector<std::string> materialsReplaced;
  std::vector<std::string> materialsNotReplaced;
  GeometryData geometry;
};

}  // namespace
}  // namespace klee
}  // namespace axom

template <>
struct FromInlet<axom::klee::ShapeData>
{
  axom::klee::ShapeData operator()(const axom::inlet::Container &base)
  {
    return axom::klee::ShapeData {base.get<std::string>("name"),
                                  base.get<std::string>("material"),
                                  base["replaces"].get<std::vector<std::string>>(),
                                  base["does_not_replace"].get<std::vector<std::string>>(),
                                  base.get<axom::klee::GeometryData>("geometry")};
  }
};

template <>
struct FromInlet<axom::klee::GeometryData>
{
  axom::klee::GeometryData operator()(const axom::inlet::Container &base)
  {
    axom::klee::GeometryData data;
    data.format = base.contains("format") ? base.get<std::string>("format") : "";
    data.path = base.contains("path") ? base.get<std::string>("path") : "";
    data.operatorData = base["operators"].get<axom::klee::internal::GeometryOperatorData>();

    if(base.contains("start_dimensions"))
    {
      data.startDimensions = axom::klee::internal::toDimensions(base["start_dimensions"]);
      data.dimensionsSet = true;
    }
    else
    {
      data.dimensionsSet = false;
    }

    std::tie(data.startUnits, data.endUnits) =
      axom::klee::internal::getOptionalStartAndEndUnits(base);

    data.pathInFile = base.name();

    return data;
  }
};

namespace axom
{
namespace klee
{
namespace
{
using inlet::Container;
using inlet::Field;
using inlet::Inlet;

/**
 * Define the schema for the "geometry" member of shapes
 *
 * @param geometry the Container representing a "geometry" object.
 */
void defineGeometry(Container &geometry)
{
  geometry.addString("format", "The format of the input file").required();
  geometry.addString("path",
                     "The path of the input file, relative to the yaml file."
                     "Required unless 'format' is 'none'");
  internal::defineDimensionsField(geometry,
                                  "start_dimensions",
                                  "The initial dimensions of the geometry file");
  internal::defineUnitsSchema(geometry,
                              "The units in which the geometry file is expressed if the units "
                              "are not embedded. These will also be the units of the operators "
                              "until they are explicitly changed.",
                              "The start units of the shape",
                              "The end units of the shape");
  internal::GeometryOperatorData::defineSchema(geometry,
                                               "operators",
                                               "Operators to apply to this object");
}

/**
 * Define the schema for the list of shapes
 *
 * @param document the Inlet document for which to define the schema
 */
void defineShapeList(Inlet &document)
{
  Container &shapeList = document.addStructArray("shapes", "The list of shapes");
  shapeList.addString("name", "The shape's name").required();
  shapeList.addString("material", "The shape's material").required();
  shapeList.addStringArray("replaces", "The list of materials this shape replaces");
  shapeList.addStringArray("does_not_replace", "The list of materials this shape does not replace");
  auto &geometry =
    shapeList.addStruct("geometry", "Contains information about the shape's geometry");
  defineGeometry(geometry);

  // Verify syntax here, semantics later!!!
  shapeList.registerVerifier(
    [](const Container &shape, std::vector<inlet::VerificationError> *errors) -> bool {
      if(shape.contains("replaces") && shape.contains("does_not_replace"))
      {
        INLET_VERIFICATION_WARNING(shape.name(),
                                   "Can't specify both 'replaces' and 'does_not_replace'",
                                   errors);
        return false;
      }

      if(shape.contains("geometry"))
      {
        const auto geom = shape.get<GeometryData>("geometry");
        if(geom.path.empty() && geom.format != "none")
        {
          INLET_VERIFICATION_WARNING(shape.name(),
                                     fmt::format("'geometry/path' field required unless "
                                                 "'geometry/format' is 'none'. "
                                                 "Provided format was '{}'",
                                                 geom.format),
                                     errors);
          return false;
        }
      }

      return true;
    });
}

/**
 * Define the schema for Klee documents.
 *
 * @param document the Inlet document for which to define the schema
 */
void defineKleeSchema(Inlet &document)
{
  internal::defineDimensionsField(document.getGlobalContainer(), "dimensions").required();
  defineShapeList(document);
  internal::NamedOperatorMapData::defineSchema(document.getGlobalContainer(), "named_operators");
}

/**
 * Create a Shape's Geometry from its raw data
 *
 * \param data the data read from inlet
 * \param fileDimensions the number of dimensions the file expects shapes to
 * have
 * \param namedOperators any named operators that were parsed from the file
 * \return the geometry description for the shape
 */
Geometry convert(GeometryData const &data,
                 Dimensions fileDimensions,
                 internal::NamedOperatorMap const &namedOperators)
{
  TransformableGeometryProperties startProperties;
  startProperties.units = data.startUnits;
  if(data.dimensionsSet)
  {
    startProperties.dimensions = data.startDimensions;
  }
  else
  {
    startProperties.dimensions = fileDimensions;
  }

  Geometry geometry {startProperties,
                     data.format,
                     data.path,
                     data.operatorData.makeOperator(startProperties, namedOperators)};

  if(geometry.getEndProperties().dimensions != fileDimensions)
  {
    throw KleeError(
      {data.pathInFile, "Did not end up in the number of dimensions specified by the file"});
  }
  return geometry;
}

/**
 * Create a Shape from its raw data representation
 *
 * \param data the data read from Inlet
 * \param fileDimensions the number of dimensions the file expects shapes to
 * have
 * \param namedOperators any named operators that were parsed from the file
 * \return the shape as a Shape object
 */
Shape convert(ShapeData const &data,
              Dimensions fileDimensions,
              internal::NamedOperatorMap const &namedOperators)
{
  return Shape {data.name,
                data.material,
                data.materialsReplaced,
                data.materialsNotReplaced,
                convert(data.geometry, fileDimensions, namedOperators)};
}

/**
 * Create a list of Shapes from their raw data representation
 *
 * \param shapeData the data read from Inlet
 * \param fileDimensions the number of dimensions the file expects shapes to
 * have
 * \param namedOperators any named operators that were parsed from the file
 * \return the shape as a Shape object
 */
std::vector<Shape> convert(std::vector<ShapeData> const &shapeData,
                           Dimensions const &fileDimensions,
                           internal::NamedOperatorMap const &namedOperators)
{
  std::vector<Shape> converted;
  converted.reserve(shapeData.size());
  for(auto &data : shapeData)
  {
    converted.emplace_back(convert(data, fileDimensions, namedOperators));
  }
  return converted;
}

/**
 * Get all named geometry operators from the file
 * \param doc the inlet document containing the file
 * \param startDimensions the number of dimensions that operators should
 * start at unless otherwise specified
 * \return all named operators read from the document
 */
internal::NamedOperatorMap getNamedOperators(const inlet::Inlet &doc, Dimensions startDimensions)
{
  if(doc.contains("named_operators"))
  {
    auto opData = doc["named_operators"].get<internal::NamedOperatorMapData>();
    return opData.makeNamedOperatorMap(startDimensions);
  }
  return internal::NamedOperatorMap {};
}
}  // namespace

ShapeSet readShapeSet(std::istream &stream)
{
  std::string contents {std::istreambuf_iterator<char>(stream), {}};

  auto reader = std::unique_ptr<inlet::YAMLReader>(new inlet::YAMLReader());
  reader->parseString(contents);

  sidre::DataStore dataStore;
  Inlet doc(std::move(reader), dataStore.getRoot());
  defineKleeSchema(doc);
  std::vector<inlet::VerificationError> errors;
  if(!doc.verify(&errors))
  {
    if(errors.empty())
    {
      throw KleeError(
        {Path {"<unknown path>"}, "Invalid Klee file given. Check the log for details."});
    }
    throw KleeError(errors);
  }

  auto shapeData = doc["shapes"].get<std::vector<ShapeData>>();
  Dimensions dimensions = internal::toDimensions(doc["dimensions"]);
  auto namedOperators = getNamedOperators(doc, dimensions);
  ShapeSet shapeSet;
  shapeSet.setDimensions(dimensions);
  shapeSet.setShapes(convert(shapeData, dimensions, namedOperators));
  return shapeSet;
}

ShapeSet readShapeSet(const std::string &filePath)
{
  std::ifstream fin {filePath};
  auto shapeSet = readShapeSet(fin);
  fin.close();
  shapeSet.setPath(filePath);
  return shapeSet;
}
}  // namespace klee
}  // namespace axom
