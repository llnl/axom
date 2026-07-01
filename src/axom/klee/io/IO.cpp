// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "IO.hpp"
#include "IOUtil.hpp"
#include "GeometryOperatorsIO.hpp"

#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/klee/KleeError.hpp"

#include "axom/config.hpp"
#include "axom/inlet.hpp"
#ifdef AXOM_USE_LUA
  #include "axom/inlet/LuaReader.hpp"
  #include "axom/sol.hpp"
#endif

#include <algorithm>
#include <cctype>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_set>

namespace axom
{
namespace klee
{
namespace
{
#ifdef AXOM_USE_LUA
class KleeLuaReader : public inlet::LuaReader
{
public:
  void setInputVariable(const std::string& name, const InputVariableValue& value)
  {
    auto lua = solState();
    std::visit([&](const auto& typedValue) { (*lua)[name] = typedValue; }, value);
  }
};
#endif

// Because we can't have context-aware validation when extracting the
// data from Inlet, we need a set of structs that parallels the real
// classes. These are used to do some basic validation, and then we convert
// them to the real classes, doing more thorough validation.

struct GeometryData
{
  std::string format;
  std::string path;
  LengthUnit startUnits {LengthUnit::unspecified};
  LengthUnit endUnits {LengthUnit::unspecified};
  Dimensions startDimensions {Dimensions::Unspecified};
  Dimensions explicitDimensions {Dimensions::Unspecified};
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
  axom::klee::ShapeData operator()(const axom::inlet::Container& base)
  {
    axom::klee::ShapeData data {base.get<std::string>("name"),
                                base.get<std::string>("material"),
                                base["replaces"].get<std::vector<std::string>>(),
                                base["does_not_replace"].get<std::vector<std::string>>(),
                                base.get<axom::klee::GeometryData>("geometry")};
    data.geometry.operatorData.setShapeName(data.name);
    return data;
  }
};

template <>
struct FromInlet<axom::klee::GeometryData>
{
  axom::klee::GeometryData operator()(const axom::inlet::Container& base)
  {
    axom::klee::GeometryData data;
    data.format = base.contains("format") ? base.get<std::string>("format") : "";
    data.path = base.contains("path") ? base.get<std::string>("path") : "";
    data.operatorData = base["operators"].get<axom::klee::internal::GeometryOperatorData>();

    data.startDimensions = base.contains("start_dimensions")
      ? axom::klee::internal::toDimensions(base["start_dimensions"])
      : axom::klee::Dimensions::Unspecified;

    data.explicitDimensions = base.contains("dimensions")
      ? axom::klee::internal::toDimensions(base["dimensions"])
      : axom::klee::Dimensions::Unspecified;

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

/**
 * Define the schema for the "geometry" member of shapes
 *
 * @param geometry the Container representing a "geometry" object.
 */
void defineGeometry(inlet::Container& geometry, bool enableLuaCallbacks)
{
  geometry.addString("format", "The format of the input file").required();
  geometry.addString("path",
                     "The path of the input file, relative to the Klee input deck."
                     "Required unless 'format' is 'none'");
  internal::defineDimensionsField(geometry,
                                  "start_dimensions",
                                  "The initial dimensions of the geometry file");
  internal::defineDimensionsField(geometry,
                                  "dimensions",
                                  "An explicit (final) dimension for the shape."
                                  "This overrides the global 'dimensions' field for this shape.");
  internal::defineUnitsSchema(geometry,
                              "The units in which the geometry file is expressed if the units "
                              "are not embedded. These will also be the units of the operators "
                              "until they are explicitly changed.",
                              "The start units of the shape",
                              "The end units of the shape");
  internal::GeometryOperatorData::defineSchema(geometry,
                                               "operators",
                                               "Operators to apply to this object",
                                               enableLuaCallbacks);
}

/**
 * Define the schema for the list of shapes
 *
 * @param document the Inlet document for which to define the schema
 */
void defineShapeList(inlet::Inlet& document, bool enableLuaCallbacks)
{
  inlet::Container& shapeList = document.addStructArray("shapes", "The list of shapes");

  shapeList.addString("name", "The shape's name").required();
  shapeList.addString("material", "The shape's material").required();
  shapeList.addStringArray("replaces", "The list of materials this shape replaces");
  shapeList.addStringArray("does_not_replace", "The list of materials this shape does not replace");
  auto& geometry =
    shapeList.addStruct("geometry", "Contains information about the shape's geometry");

  defineGeometry(geometry, enableLuaCallbacks);

  // Verify syntax here, semantics later!!!
  shapeList.registerVerifier(
    [](const inlet::Container& shape, std::vector<inlet::VerificationError>* errors) -> bool {
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
          INLET_VERIFICATION_WARNING(  //
            shape.name(),
            axom::fmt::format("'geometry/path' field required unless 'geometry/format' is 'none'. "
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
void defineKleeSchema(inlet::Inlet& document, bool enableLuaCallbacks)
{
  internal::defineDimensionsField(document.getGlobalContainer(), "dimensions").required();
  defineShapeList(document, enableLuaCallbacks);
  internal::NamedOperatorMapData::defineSchema(document.getGlobalContainer(),
                                               "named_operators",
                                               enableLuaCallbacks);
}

/**
 * Create a Shape's Geometry from its raw data
 *
 * \param data the data read from inlet
 * \param fileDimensions the number of dimensions the file expects shapes to have
 * \param namedOperators any named operators that were parsed from the file
 * \return the geometry description for the shape
 */
Geometry convert(GeometryData const& data,
                 Dimensions fileDimensions,
                 internal::NamedOperatorMap const& namedOperators)
{
  const bool has_start_dims = data.startDimensions != Dimensions::Unspecified;
  const bool has_explicit_dims = data.explicitDimensions != Dimensions::Unspecified;

  TransformableGeometryProperties startProperties;
  startProperties.units = data.startUnits;
  if(has_start_dims)
  {
    startProperties.dimensions = data.startDimensions;
  }
  else if(has_explicit_dims)
  {
    startProperties.dimensions = data.explicitDimensions;
  }
  else
  {
    startProperties.dimensions = fileDimensions;
  }

  Geometry geometry {startProperties,
                     data.format,
                     data.path,
                     data.operatorData.makeOperator(startProperties, namedOperators)};

  const auto computed_end_dims = geometry.getEndProperties().dimensions;
  const auto expected_end_dims = has_explicit_dims ? data.explicitDimensions : fileDimensions;
  if(computed_end_dims != expected_end_dims)
  {
    throw KleeError({data.pathInFile,
                     axom::fmt::format("Did not end up with the expected number of dimensions. "
                                       "Expected: {}, but got: {}",
                                       expected_end_dims,
                                       computed_end_dims)});
  }

  return geometry;
}

/**
 * Create a Shape from its raw data representation
 *
 * \param data the data read from Inlet
 * \param fileDimensions the number of dimensions the file expects shapes to have
 * \param namedOperators any named operators that were parsed from the file
 * \return the shape as a Shape object
 */
Shape convert(ShapeData const& data,
              Dimensions fileDimensions,
              internal::NamedOperatorMap const& namedOperators)
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
 * \param fileDimensions the number of dimensions the file expects shapes to have
 * \param namedOperators any named operators that were parsed from the file
 * \return the shape as a Shape object
 */
std::vector<Shape> convert(std::vector<ShapeData> const& shapeData,
                           Dimensions const& fileDimensions,
                           internal::NamedOperatorMap const& namedOperators)
{
  std::vector<Shape> converted;
  converted.reserve(shapeData.size());
  for(auto& data : shapeData)
  {
    converted.emplace_back(convert(data, fileDimensions, namedOperators));
  }
  return converted;
}

/**
 * Get all named geometry operators from the file
 * 
 * \param doc the inlet document containing the file
 * \param startDimensions the number of dimensions that operators should
 * start at unless otherwise specified
 * \return all named operators read from the document
 */
internal::NamedOperatorMap getNamedOperators(const inlet::Inlet& doc, Dimensions startDimensions)
{
  if(doc.contains("named_operators"))
  {
    auto opData = doc["named_operators"].get<internal::NamedOperatorMapData>();
    return opData.makeNamedOperatorMap(startDimensions);
  }
  return internal::NamedOperatorMap {};
}

std::string extensionOf(const std::string& filePath)
{
  const auto slash = filePath.find_last_of("/\\");
  const auto dot = filePath.find_last_of('.');
  if(dot == std::string::npos || (slash != std::string::npos && dot < slash))
  {
    return "";
  }

  std::string extension = filePath.substr(dot);
  axom::utilities::string::toLower(extension);
  return extension;
}

InputFormat inferInputFormat(const std::string& filePath)
{
  const auto extension = extensionOf(filePath);
  if(extension == ".yaml" || extension == ".yml")
  {
    return InputFormat::YAML;
  }
  if(extension == ".lua")
  {
    return InputFormat::Lua;
  }

  throw KleeError(
    {Path {filePath},
     axom::fmt::format("Unsupported Klee input file extension '{}'. Supported extensions are "
                       ".yaml, .yml, and .lua.",
                       extension.empty() ? "<none>" : extension)});
}

bool isLuaIdentifier(const std::string& name)
{
  if(name.empty())
  {
    return false;
  }

  auto isNameStart = [](unsigned char ch) { return std::isalpha(ch) || ch == '_'; };
  auto isNameChar = [](unsigned char ch) { return std::isalnum(ch) || ch == '_'; };

  if(!isNameStart(static_cast<unsigned char>(name.front())))
  {
    return false;
  }
  return std::all_of(name.begin() + 1, name.end(), [&](char ch) {
    return isNameChar(static_cast<unsigned char>(ch));
  });
}

std::unordered_set<std::string> inputVariableNames(const InputVariables& variables)
{
  std::unordered_set<std::string> names;
  for(const auto& entry : variables)
  {
    names.insert(entry.first);
  }
  return names;
}

void validateInputVariables(const InputVariables& variables)
{
  for(const auto& entry : variables)
  {
    if(!isLuaIdentifier(entry.first))
    {
      throw KleeError({Path {entry.first.empty() ? "<empty>" : entry.first},
                       axom::fmt::format("Invalid Klee Lua input variable name '{}'. Input "
                                         "variable names must be Lua identifiers.",
                                         entry.first)});
    }
  }
}

std::unique_ptr<inlet::Reader> createReader(InputFormat format, const InputVariables& variables)
{
  if(format != InputFormat::Lua && !variables.empty())
  {
    throw KleeError(
      {Path {"<unknown path>"}, "Klee input variables are only supported for Lua input decks."});
  }
  validateInputVariables(variables);

  switch(format)
  {
  case InputFormat::YAML:
    return std::make_unique<inlet::YAMLReader>();
  case InputFormat::Lua:
#ifdef AXOM_USE_LUA
  {
    auto reader = std::make_unique<KleeLuaReader>();
    for(const auto& entry : variables)
    {
      reader->setInputVariable(entry.first, entry.second);
    }
    return reader;
  }
#else
    throw KleeError(
      {Path {"<unknown path>"},
       "Lua input decks require Axom configured with AXOM_ENABLE_LUA=ON and Sol library support. "
       "Rebuild Axom with Lua enabled or convert deck to YAML."});
#endif
  }

  throw KleeError({Path {"<unknown path>"}, "Unsupported Klee input format."});
}

void appendUnexpectedGlobalErrors(const inlet::Inlet& doc,
                                  std::vector<inlet::VerificationError>& errors,
                                  const std::unordered_set<std::string>& allowedGlobals)
{
  std::unordered_set<std::string> unexpectedGlobals;
  for(const auto& name : doc.unexpectedNames())
  {
    const auto tokens = axom::utilities::string::split(name, '/');
    unexpectedGlobals.insert(tokens.empty() ? name : tokens.front());
  }

  for(const auto& name : unexpectedGlobals)
  {
    if(allowedGlobals.find(name) != allowedGlobals.end())
    {
      continue;
    }
    errors.push_back({Path {name},
                      axom::fmt::format("Unexpected global variable '{}' in Lua input deck. Use "
                                        "'local' for helper values and functions.",
                                        name)});
  }
}

ShapeSet readShapeSetFromReader(std::unique_ptr<inlet::Reader> reader,
                                bool rejectUnexpectedGlobals,
                                bool enableLuaCallbacks,
                                const std::unordered_set<std::string>& allowedGlobals)
{
  sidre::DataStore dataStore;
  inlet::Inlet doc(std::move(reader), dataStore.getRoot());
  defineKleeSchema(doc, enableLuaCallbacks);
  std::vector<inlet::VerificationError> errors;
  bool verified = doc.verify(&errors);
  if(rejectUnexpectedGlobals)
  {
    appendUnexpectedGlobalErrors(doc, errors, allowedGlobals);
    verified = verified && errors.empty();
  }
  if(!verified)
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
}  // namespace

ShapeSet readShapeSet(std::istream& stream) { return readShapeSet(stream, InputFormat::YAML); }

ShapeSet readShapeSet(std::istream& stream, InputFormat format)
{
  return readShapeSet(stream, format, InputVariables {});
}

ShapeSet readShapeSet(std::istream& stream, InputFormat format, const InputVariables& variables)
{
  std::string contents {std::istreambuf_iterator<char>(stream), {}};

  auto reader = createReader(format, variables);
  reader->parseString(contents);
  return readShapeSetFromReader(std::move(reader),
                                format == InputFormat::Lua,
                                format == InputFormat::Lua,
                                inputVariableNames(variables));
}

ShapeSet readShapeSet(const std::string& filePath)
{
  return readShapeSet(filePath, InputVariables {});
}

ShapeSet readShapeSet(const std::string& filePath, const InputVariables& variables)
{
  const auto format = inferInputFormat(filePath);
  auto reader = createReader(format, variables);
  reader->parseFile(filePath);
  auto shapeSet = readShapeSetFromReader(std::move(reader),
                                         format == InputFormat::Lua,
                                         format == InputFormat::Lua,
                                         inputVariableNames(variables));
  shapeSet.setPath(filePath);
  return shapeSet;
}
}  // namespace klee
}  // namespace axom
