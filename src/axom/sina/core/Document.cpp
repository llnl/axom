// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Document.cpp
 *
 * \brief   Implementation file for Sina Document class
 *
 ******************************************************************************
 */
#include "axom/sina/core/Document.hpp"
#include "axom/sina/core/CurveSet.hpp"
#include "axom/sina/core/Curve.hpp"
#include "axom/sina/core/Record.hpp"
#include "axom/config.hpp"
#include "axom/core/Path.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

#include "conduit.hpp"
#ifdef AXOM_USE_HDF5
  #include "conduit_relay.hpp"
  #include "conduit_relay_io.hpp"
  #include "conduit_relay_io_hdf5.hpp"
#endif

#include <functional>
#include <set>
#include <string>
#include <cstdio>
#include <fstream>
#include <ios>
#include <iostream>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace axom
{
namespace sina
{

namespace
{
char const RECORDS_KEY[] = "records";
char const RELATIONSHIPS_KEY[] = "relationships";
char const SAVE_TMP_FILE_EXTENSION[] = ".sina.tmp";
// These two are used in switch statements--could be cleaner
enum AppendFields {
    ID_FIELD,
    TYPE_FIELD,
    APPLICATION_FIELD,
    USER_FIELD,
    VERSION_FIELD,
    DATA_FIELD,
    USER_DEFINED_FIELD,
    FILES_FIELD,
    CURVE_SETS_FIELD,
    LIBRARY_DATA_FIELD,
    UNKNOWN_FIELD
};
// Really we should have one single source of truth on this one--toplevel Sina, perhaps?
static const std::map<std::string, AppendFields> appendFieldStrings {
        { "id", AppendFields::ID_FIELD }, { "type", AppendFields::TYPE_FIELD },
        { "application", AppendFields::APPLICATION_FIELD }, { "user", AppendFields::USER_FIELD },
        { "version", AppendFields::VERSION_FIELD }, { "data", AppendFields::DATA_FIELD },
        { "user_defined", AppendFields::USER_DEFINED_FIELD }, { "files", AppendFields::FILES_FIELD },
        { "curve_sets", AppendFields::CURVE_SETS_FIELD }, { "library_data", AppendFields::LIBRARY_DATA_FIELD },
    };
}  // namespace

void protocolWarn(std::string const protocol, std::string const &name)
{
  std::unordered_map<std::string, std::string> protocolMessages = {
    {".json", ".json extension not found, did you mean to save to this format?"},
    {".hdf5",
     ".hdf5 extension not found, did you use one of its other supported types? "
     "(h5, hdf, ...)"}};

  Path path(name, '.');

  if(protocol != '.' + path.baseName())
  {
    auto messageIt = protocolMessages.find(protocol);
    if(messageIt != protocolMessages.end())
    {
      std::cerr << messageIt->second;
    }
  }
}

std::string get_supported_file_types()
{
  std::string types = "[";
  for(size_t i = 0; i < supported_types.size(); ++i)
  {
    types += supported_types[i];
    if(i < supported_types.size() - 1)
    {
      types += ", ";
    }
  }
  types += "]";
  return types;
}

void Document::add(std::unique_ptr<Record> record) { records.emplace_back(std::move(record)); }

void Document::add(Relationship relationship)
{
  relationships.emplace_back(std::move(relationship));
}

conduit::Node Document::toNode() const
{
  conduit::Node document(conduit::DataType::object());
  document[RECORDS_KEY] = conduit::Node(conduit::DataType::list());
  document[RELATIONSHIPS_KEY] = conduit::Node(conduit::DataType::list());
  for(auto &record : records)
  {
    auto &list_entry = document[RECORDS_KEY].append();
    list_entry.set_node(record->toNode());
  }
  for(auto &relationship : relationships)
  {
    auto &list_entry = document[RELATIONSHIPS_KEY].append();
    list_entry = relationship.toNode();
  }
  return document;
}

void Document::createFromNode(const conduit::Node &asNode, const RecordLoader &recordLoader)
{
  conduit::Node nodeCopy = asNode;

  auto processChildNodes = [&](const char *key, std::function<void(conduit::Node &)> addFunc) {
    if(nodeCopy.has_child(key))
    {
      conduit::Node &childNodes = nodeCopy[key];

      // -- 1. Check if this node is a primitive leaf (throw immediately if so)
      // Customize these checks to match exactly what you consider "primitive."
      if(childNodes.dtype().is_number() || childNodes.dtype().is_char8_str() ||
         childNodes.dtype().is_string())
      {
        std::ostringstream message;
        message << "The '" << key << "' element of a document cannot be a primitive value.";
        throw std::invalid_argument(message.str());
      }

      // -- 2. Not a primitive. Check if it has no children.
      if(childNodes.number_of_children() == 0)
      {
        // Turn it into an empty list
        childNodes.set(conduit::DataType::list());
      }

      // -- 3. If it's still not a list, throw
      if(!childNodes.dtype().is_list())
      {
        std::ostringstream message;
        message << "The '" << key << "' element of a document must be an array/list.";
        throw std::invalid_argument(message.str());
      }

      // -- 4. Now it's guaranteed to be a list, so iterate
      auto childIter = childNodes.children();
      while(childIter.has_next())
      {
        conduit::Node child = childIter.next();
        addFunc(child);
      }
    }
  };
  processChildNodes(RECORDS_KEY, [&](conduit::Node &record) { add(recordLoader.load(record)); });

  processChildNodes(RELATIONSHIPS_KEY,
                    [&](conduit::Node &relationship) { add(Relationship {relationship}); });
}

Document::Document(conduit::Node const &asNode, RecordLoader const &recordLoader)
{
  this->createFromNode(asNode, recordLoader);
}

Document::Document(std::string const &asJson, RecordLoader const &recordLoader)
{
  conduit::Node asNode;
  asNode.parse(asJson, "json");
  this->createFromNode(asNode, recordLoader);
}

#ifdef AXOM_USE_HDF5
void removeSlashes(const conduit::Node &originalNode, conduit::Node &modifiedNode)
{
  for(auto it = originalNode.children(); it.has_next();)
  {
    it.next();
    std::string key = it.name();
    std::string modifiedKey = axom::utilities::string::replaceAllInstances(key, "/", slashSubstitute);

    modifiedNode[modifiedKey] = it.node();

    if(it.node().dtype().is_object())
    {
      conduit::Node nestedNode;
      removeSlashes(it.node(), nestedNode);
      modifiedNode[modifiedKey].set(nestedNode);
    }
  }
}

void restoreSlashes(const conduit::Node &modifiedNode, conduit::Node &restoredNode)
{
  // Check if List or Object, if its a list the else statement would turn it into an object
  // which breaks the Document

  if(modifiedNode.dtype().is_list())
  {
    // If its empty with no children it's the end of a tree

    for(auto it = modifiedNode.children(); it.has_next();)
    {
      it.next();
      conduit::Node &newChild = restoredNode.append();
      auto data_type = it.node().dtype();

      // Leaves empty nodes empty, if null data is set the
      // Document breaks

      if(data_type.is_string() || data_type.is_number())
      {
        newChild.set(it.node());  // Lists need .set
      }

      // Recursive Call
      if(it.node().number_of_children() > 0)
      {
        restoreSlashes(it.node(), newChild);
      }
    }
  }
  else
  {
    for(auto it = modifiedNode.children(); it.has_next();)
    {
      it.next();
      std::string key = it.name();
      std::string restoredKey =
        axom::utilities::string::replaceAllInstances(key, slashSubstitute, "/");

      // Initialize a new node for the restored key
      conduit::Node &newChild = restoredNode.add_child(restoredKey);
      auto data_type = it.node().dtype();

      // Leaves empty keys empty but continues recursive call if its a list
      if(data_type.is_string() || data_type.is_number() || data_type.is_object())
      {
        newChild.set(it.node());
      }
      else if(data_type.is_list())
      {
        restoreSlashes(it.node(), newChild);  // Handle nested lists
      }

      // If the node has children, recursively restore them
      if(it.node().number_of_children() > 0)
      {
        conduit::Node nestedNode;
        restoreSlashes(it.node(), nestedNode);
        newChild.set(nestedNode);
      }
    }
  }
}

conduit::Node Document::toHDF5Node() const
{
  conduit::Node node;
  conduit::Node &recordsNode = node["records"];
  conduit::Node &relationshipsNode = node["relationships"];

  for(const auto &record : getRecords())
  {
    conduit::Node recordNode = record->toNode();

    removeSlashes(recordNode, recordsNode.append());
  }

  // Process relationships
  for(const auto &relationship : getRelationships())
  {
    conduit::Node relationshipNode = relationship.toNode();

    removeSlashes(relationshipNode, relationshipsNode.append());
  }
  return node;
}

void Document::toHDF5(const std::string &filename) const
{

  conduit::relay::io::save(this->toHDF5Node(), filename, "hdf5");
}
#endif

//

std::string Document::toJson(conduit::index_t indent,
                             conduit::index_t depth,
                             const std::string &pad,
                             const std::string &eoe) const
{
  return this->toNode().to_json("json", indent, depth, pad, eoe);
}

void saveDocument(Document const &document, std::string const &fileName, Protocol protocol)
{
  // It is a common use case for users to want to overwrite their files as
  // the simulation progresses. However, this operation should be atomic so
  // that if a write fails, the old file is left intact. For this reason,
  // we write to a temporary file first and then move the file. The temporary
  // file is in the same directory to ensure that it is part of the same
  // file system as the destination file so that the move operation is
  // atomic.

  std::string tmpFileName = fileName + SAVE_TMP_FILE_EXTENSION;

  switch(protocol)
  {
  case Protocol::JSON:
  {
    protocolWarn(".json", fileName);
    auto asJson = document.toJson();
    std::ofstream fout {tmpFileName};
    fout.exceptions(std::ostream::failbit | std::ostream::badbit);
    fout << asJson;
    fout.close();
  }
  break;
#ifdef AXOM_USE_HDF5
  case Protocol::HDF5:
    protocolWarn(".hdf5", fileName);
    document.toHDF5(tmpFileName);
    break;
#endif
  default:
  {
    std::ostringstream message;
    message << "Invalid format choice. Please choose from one of the supported "
               "protocols: "
            << get_supported_file_types();
    throw std::invalid_argument(message.str());
  }
  }

  if(rename(tmpFileName.c_str(), fileName.c_str()) != 0)
  {
    std::string message {"Could not save to '"};
    message += fileName;
    message += "'";
    throw std::ios::failure {message};
  }
}

Document loadDocument(std::string const &path, Protocol protocol)
{
  return loadDocument(path, createRecordLoaderWithAllKnownTypes(), protocol);
}

Document loadDocument(std::string const &path, RecordLoader const &recordLoader, Protocol protocol)
{
  conduit::Node node, modifiedNode;
  std::ostringstream file_contents;
  std::ifstream file_in {path};

  // Load the file depending on the protocol
  switch(protocol)
  {
  case Protocol::JSON:
    file_contents << file_in.rdbuf();
    file_in.close();
    node.parse(file_contents.str(), "json");
    return Document {node, recordLoader};
#ifdef AXOM_USE_HDF5
  case Protocol::HDF5:
    file_in.close();
    conduit::relay::io::load(path, "hdf5", node);
    restoreSlashes(node, modifiedNode);
    return Document {modifiedNode, recordLoader};
#endif
  default:
    std::ostringstream message;
    message << "Invalid format choice. Please choose from one of the supported "
               "protocols: "
            << get_supported_file_types();
    throw std::invalid_argument(message.str());
    break;
  }
}

// Unified helper function that validates a set of curves (either dependent or independent).
// Parameters:
//   new_curves     : a map of new curves (e.g. std::map<std::string, Curve>).
//   existing_keys  : a set of keys that are present in the existing data.
//   getExistingSize: a callable that takes a curve key and returns the size (int)
//                    from the existing data or -1 if not found.
//   curveType      : "dependent" or "independent" (for erro messages).
//   recordId       : identifier for the record (for error messages).
//   curveSetId     : identifier for the curve set (for error messages).
//   baseline       : reference to an int that will hold the computed baseline value
//                    (initialize to -1 to have the function set baseline).
//
// Templated on the container type to avoid conversion issues (e.g., unordered_map vs. map).
template <typename CurveMap>
bool validate_curves_unified(const CurveMap &new_curves,
                             const std::set<std::string> &existing_keys,
                             std::function<int(const std::string &)> getExistingSize,
                             const std::string &curveType,
                             const std::string &recordId,
                             const std::string &curveSetId,
                             int baseline)
{
  std::set<std::string> unionKeys = existing_keys;
  for(const auto &pair : new_curves)
  {
    unionKeys.insert(pair.first);
  }

  for(const auto &key : unionKeys)
  {
    int newSize = 0;
    auto newItr = new_curves.find(key);
    if(newItr != new_curves.end())
    {
      newSize = static_cast<int>(newItr->second.getValues().size());
    }
    int existingSize = getExistingSize(key);

    // Get total size but ignore -1 returns
    int total = (existingSize >= 0 ? newSize + existingSize : newSize);

    if(baseline < 0)
    {
      baseline = total;
    }
    else if(baseline != total)
    {
      std::cerr << "Error validating " << curveType << ": Record " << recordId << ", Curve Set "
                << curveSetId << ", Curve " << key << " size mismatch (expected " << baseline
                << ", got " << total << ")." << std::endl;
      return false;
    }
  }
  return true;
}

bool validate_curve_sets_json(const DataHolder::CurveSetMap new_curve_sets,
                              const conduit::Node &existing_curve_sets,
                              const std::string recordId)
{
  // Iterate over the keys in the existing curve sets node.
  for(const auto &curveSetId : existing_curve_sets.child_names())
  {
    const conduit::Node &ecs = existing_curve_sets[curveSetId];

    if(new_curve_sets.find(curveSetId) != new_curve_sets.end())
    {
      const auto &new_curve_set = new_curve_sets.at(curveSetId);

      // --- Validate Dependent Curves ---
      std::set<std::string> existingDepKeys;
      if(ecs.has_child("dependent"))
      {
        for(const auto &dep_key : ecs["dependent"].child_names())
        {
          existingDepKeys.insert(dep_key);
        }
      }
      // Use a lambda to get the current size of the dependent curve array.
      auto getExistingDepSize = [&ecs](const std::string &key) -> int {
        if(ecs.has_child("dependent") && ecs["dependent"].has_child(key))
        {
          return static_cast<int>(ecs["dependent"][key]["value"].dtype().number_of_elements());
        }
        return -1;
      };

      if(!new_curve_set.getDependentCurves().empty() &&
         !validate_curves_unified(new_curve_set.getDependentCurves(),
                                  existingDepKeys,
                                  getExistingDepSize,
                                  "dependent",
                                  recordId,
                                  curveSetId,
                                  -1))
      {
        return false;
      }

      // --- Validate Independent Curves ---
      std::set<std::string> existingIndepKeys;
      if(ecs.has_child("independent"))
      {
        for(const auto &indep_key : ecs["independent"].child_names())
        {
          existingIndepKeys.insert(indep_key);
        }
      }
      auto getExistingIndepSize = [&ecs](const std::string &key) -> int {
        if(ecs.has_child("independent") && ecs["independent"].has_child(key))
        {
          return static_cast<int>(ecs["independent"][key]["value"].dtype().number_of_elements());
        }
        return -1;
      };

      if(!new_curve_set.getIndependentCurves().empty() &&
         !validate_curves_unified(new_curve_set.getIndependentCurves(),
                                  existingIndepKeys,
                                  getExistingIndepSize,
                                  "independent",
                                  recordId,
                                  curveSetId,
                                  -1))
      {
        return false;
      }
    }
    else
    {
      std::cerr << "Curve set " << curveSetId << " not found in new data for record " << recordId
                << std::endl;
      return false;
    }
  }
  return true;
}

/*bool append_to_json(const std::string &jsonFilePath,
                    Document const &newData,
                    const int data_protocol,
                    const int udc_protocol)
{
  conduit::Node root;
  try
  {
    root.load(jsonFilePath, "json");
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error loading JSON file: " << e.what() << std::endl;
    return false;
  }

  if((unsigned long)root["records"].number_of_children() != newData.getRecords().size())
  {
    std::cerr << "Mismatch in the number of records." << std::endl;
    return false;
  }

  std::vector<std::function<void()>> write_queue;

  for(auto &new_record : newData.getRecords())
  {
    for(conduit::index_t i = 0; i < root["records"].number_of_children(); ++i)
    {
      conduit::Node &existing_record = root["records"].child(i);
      if(!validate_curve_sets_json(new_record->getCurveSets(),
                                   existing_record["curve_sets"],
                                   existing_record["id"].as_string()))
      {
        return false;
      }
    }
  }

  for(auto &new_record : newData.getRecords())
  {
    bool found = false;
    for(int i = 0; i < root["records"].number_of_children(); ++i)
    {
      conduit::Node &existing_record = root["records"].child(i);

      if(new_record->getId().getId() == existing_record["id"].as_string())
      {
        found = true;
        // ----------------------
        // Queue update of CURVE SETS.
        // ----------------------
        // Example: processing the curve_sets section for a given record.
        if((new_record->getCurveSets().size() > 0) && existing_record.has_child("curve_sets"))
        {
          conduit::Node &existing_curve_sets = existing_record["curve_sets"];
          auto &new_curve_sets = new_record->getCurveSets();
          for(auto &new_curve_set : new_curve_sets)
          {
            auto cs_key = new_curve_set.first;
            auto &new_cs_values = new_curve_set.second;

            // Queue merging of dependent curves.
            for(auto &new_dependent : new_cs_values.getDependentCurves())
            {
              // Get the key and make a copy of the new dependent values.
              auto dep_key = new_dependent.first;
              auto new_dep_vals = new_dependent.second.getValues();  // copy values for lambda capture

              write_queue.push_back([cs_key, dep_key, new_dep_vals, &existing_curve_sets]() {
                // Access the existing values node for this dependent curve.
                conduit::Node &existing_vals_node =
                  existing_curve_sets[cs_key]["dependent"][dep_key]["value"];
                std::vector<double> merged_values;

                conduit::index_t num_elems = existing_vals_node.dtype().number_of_elements();
                for(conduit::index_t i = 0; i < num_elems; i++)
                {
                  merged_values.push_back(existing_vals_node.as_double_array()[i]);
                }

                for(const auto &val : new_dep_vals)
                {
                  merged_values.push_back(static_cast<double>(val));
                }
                // Reset and update the node with merged values.
                existing_vals_node.reset();
                existing_vals_node.set(merged_values);
              });
            }

            // Queue merging of independent curves.
            for(auto &new_independent : new_cs_values.getIndependentCurves())
            {
              auto indep_key = new_independent.first;
              auto new_indep_vals = new_independent.second.getValues();  // copy for lambda capture

              write_queue.push_back([cs_key, indep_key, new_indep_vals, &existing_curve_sets]() {
                conduit::Node &existing_vals_node =
                  existing_curve_sets[cs_key]["independent"][indep_key]["value"];
                std::vector<double> merged_values;

                conduit::index_t num_elems = existing_vals_node.dtype().number_of_elements();
                for(conduit::index_t i = 0; i < num_elems; i++)
                {
                  merged_values.push_back(existing_vals_node.as_double_array()[i]);
                }

                for(const auto &val : new_indep_vals)
                {
                  merged_values.push_back(static_cast<double>(val));
                }
                existing_vals_node.reset();
                existing_vals_node.set(merged_values);
              });
            }
          }
        }

        // ----------------------
        // Queue update of DATA VALUES.
        // ----------------------
        if((new_record->getData().size() > 0) && existing_record.has_child("data"))
        {
          conduit::Node &existing_data_sets = existing_record["data"];
          auto &new_data_sets = new_record->getData();
          for(auto &new_data : new_data_sets)
          {
            auto &new_data_key = new_data.first;
            auto &new_data_pair = new_data.second;
            auto data_key = new_data_key;  // local copy for capture
            conduit::Node obj;
            obj.parse(new_data_pair.toNode().to_json(), "json");
            if(existing_data_sets.has_child(data_key))
            {
              // Duplicate Handling
              switch(data_protocol)
              {
              case 1:
                write_queue.push_back([&existing_data_sets, data_key, obj]() {
                  existing_data_sets[data_key].update(obj);
                });
                break;
              case 2:
                break;
              case 3:
                std::cerr
                  << "Found a duplicate data entry, protocol 3 dictates append cancellation."
                  << std::endl;
                return false;
              default:
                std::cerr << "Invalid Data Protocol Entry. Append cancelled." << std::endl;
                return false;
              }
            }
            else
            {
              write_queue.push_back([&existing_data_sets, data_key, obj]() {
                existing_data_sets[data_key].update(obj);
              });
            }
          }
        }

        // ----------------------
        // Queue update of USER DEFINED CONTENT.
        // ----------------------
        if((!new_record->getUserDefinedContent().dtype().is_empty()) &&
           existing_record.has_child("user_defined"))
        {
          conduit::Node &existing_udc = existing_record["user_defined"];
          auto &new_udc = new_record->getUserDefinedContent();
          for(auto &udc : new_udc.children())
          {
            std::string udc_name = udc.name();
            if(existing_record["user_defined"].has_child(udc_name))
            {
              switch(udc_protocol)
              {
              case 1:
                write_queue.push_back([&existing_udc, udc_name, udc]() {
                  existing_udc[udc_name].set(udc.as_string());
                });
                break;
              case 2:
                break;
              case 3:
                std::cerr << "Found duplicate UDC, protocol 3 dictates append cancellation."
                          << std::endl;
                return false;
              default:
                std::cerr << "Invalid UDC Protocol Entry. Append cancelled." << std::endl;
                return false;
              }
            }
            else
            {
              write_queue.push_back(
                [&existing_udc, udc_name, udc]() { existing_udc[udc_name].set(udc.as_string()); });
            }
          }
        }
      }
    }
    if(!found)
    {
      conduit::Node obj;
      obj.parse(new_record->toNode().to_json(), "json");
      write_queue.push_back([&root, obj]() { root["records"].append().update(obj); });
    }
  }

  // Execute all queued operations.
  for(auto &op : write_queue)
  {
    try
    {
      op();
    }
    catch(const std::exception &e)
    {
      std::cerr << "Error executing queued operation: " << e.what() << std::endl;
      return false;
    }
  }

  try
  {
    root.save(jsonFilePath, "json");
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error saving JSON file: " << e.what() << std::endl;
    return false;
  }

  return true;
}*/

// Specifically validate ONE curve set for ONE DataHolder (record, library_data...) for appending,
// appendTo is notionally const, but the has_path() etc. methods aren't const.
conduit::Node validate_curve_sets(conduit::relay::io::IOHandle &appendTo,
                                 const conduit::Node &appendFrom,
                                 const std::string &endpoint) {
    int baseline = -1;  // baseline is shared across dependent and independent
    conduit::Node msgNode = conduit::Node(conduit::DataType::list());
    const std::vector<std::string> curve_categories = {"dependent", "independent"};
    for(const std::string& curve_cat : curve_categories){
        // First, we need to make sure we didn't "forget" an existing curve: make sure either EVERY
        // curve in the HDF5 is being appended to, or NONE of them are (all curves are new). 
        std::string curves_endpoint = endpoint + "/" + curve_cat;
        std::vector<std::string> curve_names;
        appendTo.list_child_names(curves_endpoint, curve_names);
        u_int curvesWritten = 0;
        for(const auto &cname : curve_names){
            if(appendFrom[curves_endpoint].has_child(cname) && appendFrom[curves_endpoint][cname].number_of_children() > 0){
                curvesWritten ++;
            }
        }
        if(curvesWritten != 0 && curvesWritten != curve_names.size()){
              msgNode.append() = "Failed to append curve sets: ";//parent had " + curve_names.size() + " curves, but append only addressed " + curvesWritten);
        }
        // Now loop through what we've actually got. Once we find something, use it to set the baseline.
        auto curvesIter = appendFrom[curve_cat].children();
        while(curvesIter.has_next()) {
          auto &testCurve = curvesIter.next();
          std::string sub_endpoint;
          int post_append_size = testCurve.number_of_children();
          if(appendTo.has_path(endpoint + "/" + curvesIter.name())){
              conduit::Node n;
              sub_endpoint = endpoint + "/" + curvesIter.name() + "/num_elements";
              appendTo.read(sub_endpoint, n);
              post_append_size += n.to_int();
              if(baseline == -1){baseline = post_append_size;}
              if(post_append_size != baseline){
                msgNode.append() = "Failed to append curve sets: ";// + curvesIter.name() + " had " + post_append_size + " values after appending, but expected " + baseline);
              }
          }
        }
    }
    return msgNode;
}

// Top-level append validation function. Works recursively on library data (hence endpoint)
conduit::Node validate_append( conduit::relay::io::IOHandle &appendTo, const conduit::Node &appendFrom,
                      const std::string &endpoint, const int mergeProtocol){

    conduit::Node msgNode = conduit::Node(conduit::DataType::list());
    // Case one: die if the types disagree. A pingpong_game shouldn't become a billiards_game
    // Validation note: we allow for no type here beecause of library_data. If someone
    // forgot the type for a top-level record, it should've died already.
    if(appendFrom.has_child("type")){
      conduit::Node typeNode;
      appendTo.read(endpoint+"/type", typeNode);
      if(typeNode.to_string() != appendFrom["type"].to_string()){
        msgNode.append() = "Failed to append record: type mismatch"; // TODO: better errors
      }
    }

    // Case two: merge protocol is 3, we need to die if certain fields appear in both places
    if(mergeProtocol == 3){
      const std::vector<std::string> prot3Fields = {"data", "user_defined", "files"};
      for(auto &field : prot3Fields){
        auto dataIter = appendFrom[field].children();
        while(dataIter.has_next()){
          dataIter.next();
          if(appendTo.has_path(endpoint + "/" + field + "/" + dataIter.name())){
            msgNode.append() = "Failed to overlay (append protocol 3) record: conflicting data";
          }
        }
      }
    }

    // Case three: curve sets. Go into each curve set and, if it already exists, make sure the to-be-appended curves are valid.
    if(appendFrom.has_child("curve_sets"))
    {
      std::string subEndpoint;
      auto curveSetsIter = appendFrom["curve_sets"].children();
      while(curveSetsIter.has_next())
      {
        curveSetsIter.next();
        subEndpoint = endpoint + "/curve_sets/" + curveSetsIter.name();
        // We only have to validate if the hdf5 already has a curve set with that name.
        if(appendTo.has_path(subEndpoint)){
          // TODO: merge lists cleanly
          msgNode.append() = validate_curve_sets(appendTo, curveSetsIter.next(), subEndpoint);
        }
      }
    }

    // Case four: library data. Recurse on it if it's already in the hdf5.
    auto libraryIter = appendFrom["library_data"].children();
    std::string subEndpoint;
    while(libraryIter.has_next()){
      libraryIter.next();
      subEndpoint = endpoint + "/library_data/" + libraryIter.name();
      // We only have to validate if the hdf5 already has a library with that name.
      if(appendTo.has_path(endpoint)){
        msgNode.append() = validate_append(appendTo, libraryIter.next(), subEndpoint, mergeProtocol);
      }
    }
    return msgNode;
}

// Avoiding a terrible if/else chunk in append_recordlike_fields and friends.
AppendFields field_lookup(const std::string &input){
    auto itr = appendFieldStrings.find(input);
    if( itr != appendFieldStrings.end() ) {
        return itr->second;
    }
    return AppendFields::UNKNOWN_FIELD; 
}

void append_recordlike_fields(conduit::relay::io::IOHandle &appendTo,
                              conduit::Node &appendFrom,
                              const std::string &endpoint,
                              const int mergeProtocol,
                              bool canAppendCurves){
  auto fieldsIter = appendFrom.children();
  while(fieldsIter.has_next()){
    conduit::Node &recField = fieldsIter.next();
    std::string appendAtEndpoint = endpoint+"/"+fieldsIter.name();
    switch (field_lookup(fieldsIter.name())) {
    case AppendFields::ID_FIELD:  // we already have it, else we couldn't find this
    case AppendFields::TYPE_FIELD: // we already have it, else we exploded above
    case AppendFields::USER_FIELD: // special run fields, ignore? explode?
    case AppendFields::VERSION_FIELD:
    case AppendFields::APPLICATION_FIELD:
      break;
    // For simpler/arbitrary structures, we just overwrite.
    case AppendFields::DATA_FIELD:
    case AppendFields::USER_DEFINED_FIELD:
    case AppendFields::FILES_FIELD:
      // We already exploded for protocol 3 if anything overwrote, so we handle 1 and 3
      if(mergeProtocol==1 || mergeProtocol==3){
        appendTo.write(recField, appendAtEndpoint);
      } else if(mergeProtocol==2){
        auto subFieldIter = appendFrom[fieldsIter.name()].children();
        while(subFieldIter.has_next()){
          conduit::Node &subField = subFieldIter.next();
          if(!appendTo.has_path(appendAtEndpoint)){
            appendTo.write(subField, appendAtEndpoint+"/"+fieldsIter.name()+"/"+subFieldIter.name());
          }
        }
      }
      break;
    case AppendFields::LIBRARY_DATA_FIELD: {
      // We recurse
      auto libraryIter = appendFrom[fieldsIter.name()].children();
      while(libraryIter.has_next()){
        libraryIter.next();
        if(appendTo.has_path(endpoint)){
          append_recordlike_fields(appendTo, libraryIter.next(),
                                   appendAtEndpoint+"/"+fieldsIter.name()+"/"+libraryIter.name(),
                                   mergeProtocol, canAppendCurves);
        }
      }
      break;
    }
    case AppendFields::CURVE_SETS_FIELD:
      // Conduit is kind: the "opts" node, required for HDF5 appending,
      // is simply ignored for JSON.
      
      break;
    default:
        std::cout << "oh no" << std::endl;
    }
  }
}

conduit::Node append(conduit::relay::io::IOHandle &appendTo,
            conduit::Node &appendFrom,
            const int mergeProtocol,
            bool canAppendCurves){

  // We do all of our validation up-front, including/especially library data.
  // There's no validation to perform on things like relationships (right..?)
  auto recordsIter = appendFrom["records"].children();
  conduit::Node msgNode = conduit::Node(conduit::DataType::list());
  while(recordsIter.has_next())
  {
    conduit::Node &n = recordsIter.next();
    msgNode.append() = validate_append(appendTo, n, recordsIter.name(), mergeProtocol);

  }
  // Return with our error list
  if(msgNode.number_of_children() > 0){ return msgNode; }

  // Our validation passed, time to throw it all in!
  recordsIter = appendFrom["records"].children();
  while(recordsIter.has_next())
  {
    conduit::Node &rec = recordsIter.next();
    std::string rec_endpoint = "records/" + rec["id"].to_string() + "/";
    // Easiest case, the record doesn't exist yet. Add it.
    if(!appendTo.has_path(rec_endpoint)){
      appendTo.write(rec, rec_endpoint);
    } else {
      append_recordlike_fields(appendTo, rec, "records/"+recordsIter.name(), mergeProtocol, canAppendCurves);
    }
  }
  return msgNode;
}

conduit::Node append_to_json(const std::string &jsonFilePath,
                    const Document &newData,
                    const int mergeProtocol){
  conduit::relay::io::IOHandle appendTo;
  appendTo.open(jsonFilePath);
  conduit::Node appendFrom = newData.toNode();
  conduit::Node success = append(appendTo, appendFrom, mergeProtocol, false);
  // More conduit kindness: it looks like the writes only resolve once we close.
  // In that case, we don't need to worry about re-re-re-dumping this file with writes.
  appendTo.close();
  return success;
}

conduit::Node append_to_hdf5(const std::string &hdf5FilePath,
                    const Document &newData,
                    const int mergeProtocol){ 
#ifdef AXOM_USE_HDF5
  conduit::relay::io::IOHandle appendTo;
  appendTo.open(hdf5FilePath);
  conduit::Node appendFrom = newData.toHDF5Node();
  conduit::Node msgNode = append(appendTo, appendFrom, mergeProtocol, true);
  appendTo.close();
  return msgNode;
#endif
#ifndef AXOM_USE_HDF5
  conduit::Node msgNode = conduit::Node(conduit::DataType::list());
  msgNode.append("Failed to append Sina HDF5: Axom wasn't build with HDF5");
  return msgNode;
#endif
}


    /*bool found = false;
    for(const auto &rec_name : record_list)
    {
      std::string id_path = "records/" + rec_name + "/id/";
      conduit::Node id;
      existing_file.read(id_path, id);
      if(id.to_string() == "\"" + record->getId().getId() + "\"")
      {
        found = true;
        std::string record_path = "records/" + rec_name;

        // ----------------------
        // Queue update of DATA VALUES.
        // ----------------------
        try
        {
          data_path = "records/" + rec_name + "/data/";
          std::vector<std::string> data_keys;
          existing_file.list_child_names(data_path, data_keys);
          auto &new_data_sets = record->getData();
          for(auto &new_data : new_data_sets)
          {
            auto new_data_key = new_data.first;
            auto new_data_pair = new_data.second;
            if(std::find(data_keys.begin(), data_keys.end(), new_data_key) != data_keys.end())
            {
              // Duplicate Handling
              switch(data_protocol)
              {
              case 1:
                write_queue.push_back([&existing_file, data_path, new_data_key, new_data_pair]() {
                  existing_file.remove(data_path + new_data_key);
                  existing_file.write(new_data_pair.toNode(), data_path + new_data_key);
                });
                break;
              case 2:
                break;
              case 3:
                std::cerr
                  << "Found a duplicate data entry, protocol 3 dictates append cancellation."
                  << std::endl;
                return false;
              default:
                std::cerr << "Invalid Data Protocol Entry. Append cancelled." << std::endl;
                return false;
              }
            }
            else
            {
              write_queue.push_back([&existing_file, data_path, new_data_key, new_data_pair]() {
                existing_file.write(new_data_pair.toNode(), data_path + new_data_key);
              });
            }
          }
        }
        catch(const std::exception &e)
        {
          std::cerr << "Error updating data values: " << e.what() << std::endl;
        }

        // ----------------------
        // Queue update of USER DEFINED CONTENT.
        // ----------------------
        try
        {
          udc_path = "records/" + rec_name + "/user_defined/";
          std::vector<std::string> udc_list;
          existing_file.list_child_names(udc_path, udc_list);
          auto &new_udc_sets = record->getUserDefinedContent();
          for(auto &new_udc : new_udc_sets.children())
          {
            std::string udc_name = new_udc.name();
            if(std::find(udc_list.begin(), udc_list.end(), udc_name) != udc_list.end())
            {
              // Duplicate Handling
              switch(udc_protocol)
              {
              case 1:
                write_queue.push_back([&existing_file, udc_path, udc_name, new_udc]() {
                  existing_file.remove(udc_path + udc_name);
                  existing_file.write(new_udc, udc_path + udc_name);
                });
                break;
              case 2:
                break;
              case 3:
                std::cerr << "Found duplicate UDC, protocol 3 dictates append cancellation."
                          << std::endl;
                return false;
              default:
                std::cerr << "Invalid UDC Protocol Entry. Append cancelled." << std::endl;
                return false;
              }
            }
            else
            {
              write_queue.push_back([&existing_file, udc_path, udc_name, new_udc]() {
                existing_file.write(new_udc, udc_path + udc_name);
              });
            }
          }
        }
        catch(const std::exception &e)
        {
          std::cerr << "Error updating user defined content: " << e.what() << std::endl;
        }

        // ----------------------
        // Queue update of CURVE SETS.
        // ----------------------
        try
        {
          cs_path = "records/" + rec_name + "/curve_sets/";
          std::vector<std::string> curve_sets_list;
          existing_file.list_child_names("records/" + rec_name + "/curve_sets", curve_sets_list);
          auto &new_curve_sets = record->getCurveSets();

          for(const auto &curve_set_pair : new_curve_sets)
          {
            const std::string &cs_name = curve_set_pair.first;
            const auto &cs_data = curve_set_pair.second;

            bool cs_exists = (std::find(curve_sets_list.begin(), curve_sets_list.end(), cs_name) !=
                              curve_sets_list.end());

            if(cs_exists)
            {
              const auto &new_dependents = cs_data.getDependentCurves();
              for(const auto &dep_pair : new_dependents)
              {
                const std::string &dep_name = dep_pair.first;
                std::string curve_set_path =
                  record_path + "/curve_sets/" + cs_name + "/dependent/" + dep_name + "/value";

                if(info.has_path(curve_set_path) && info[curve_set_path].has_child("num_elements"))
                {
                  int current_size = info[curve_set_path]["num_elements"].to_int();
                  opts["offset"] = current_size;
                }
                else
                {
                  opts["offset"] = 0;
                }

                const auto &new_dep_values = dep_pair.second.getValues();
                std::vector<double> current_array(new_dep_values.begin(), new_dep_values.end());
                to_set.set(current_array);

                write_queue.push_back([curve_set_path, hdf5FilePath, opts, to_set]() mutable {
                  conduit::relay::io::hdf5_write(to_set, hdf5FilePath, curve_set_path, opts, true);
                });
              }

              const auto &new_independents = cs_data.getIndependentCurves();
              for(const auto &indep_pair : new_independents)
              {
                const std::string &indep_name = indep_pair.first;
                std::string curve_set_path =
                  record_path + "/curve_sets/" + cs_name + "/independent/" + indep_name + "/value";

                if(info.has_path(curve_set_path) && info[curve_set_path].has_child("num_elements"))
                {
                  int current_size = info[curve_set_path]["num_elements"].to_int();
                  opts["offset"] = current_size;
                }
                else
                {
                  opts["offset"] = 0;
                }

                const auto &new_indep_values = indep_pair.second.getValues();
                std::vector<double> current_array(new_indep_values.begin(), new_indep_values.end());
                to_set.set(current_array);

                write_queue.push_back([curve_set_path, hdf5FilePath, opts, to_set]() mutable {
                  conduit::relay::io::hdf5_write(to_set, hdf5FilePath, curve_set_path, opts, true);
                });
              }
            }
            else
            {
              conduit::Node cs_node;
              const auto &new_dependents = cs_data.getDependentCurves();
              for(const auto &dep_pair : new_dependents)
              {
                const std::string &dep_name = dep_pair.first;
                const auto &new_dep_values = dep_pair.second.getValues();
                std::vector<double> dep_array(new_dep_values.begin(), new_dep_values.end());
                cs_node["dependent"][dep_name]["value"].set(dep_array);
              }

              // Process independent curves.
              const auto &new_independents = cs_data.getIndependentCurves();
              for(const auto &indep_pair : new_independents)
              {
                const std::string &indep_name = indep_pair.first;
                const auto &new_indep_values = indep_pair.second.getValues();
                std::vector<double> indep_array(new_indep_values.begin(), new_indep_values.end());
                cs_node["independent"][indep_name]["value"].set(indep_array);
              }

              write_queue.push_back([record_path, cs_name, hdf5FilePath, cs_node]() mutable {
                conduit::relay::io::hdf5_write(cs_node,
                                               hdf5FilePath,
                                               record_path + "/curve_sets/" + cs_name,
                                               conduit::Node(),
                                               true);
              });
            }
          }
        }
        catch(const std::exception &e)
        {
          std::cerr << "Error updating curve sets: " << e.what() << std::endl;
        }
      }
    }
    if(!found)
    {
      std::string record_path = "records/" + std::to_string(extension++) + "/";
      write_queue.push_back([&existing_file, record_path, rec_ptr = record.get()]() {
        existing_file.write(rec_ptr->toNode(), record_path);
      });
    }
  }

  // Execute all queued operations.
  for(auto &write_op : write_queue)
  {
    try
    {
      write_op();
    }
    catch(const std::exception &e)
    {
      std::cerr << "Error executing queued write operation: " << e.what() << std::endl;
      return false;
    }
  }

  return true;
}*/

}  // namespace sina
}  // namespace axom
