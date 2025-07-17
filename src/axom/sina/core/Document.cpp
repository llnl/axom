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
#include "axom/core.hpp"
#include "axom/fmt.hpp"

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
enum AppendFields
{
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
  {"id", AppendFields::ID_FIELD},
  {"type", AppendFields::TYPE_FIELD},
  {"application", AppendFields::APPLICATION_FIELD},
  {"user", AppendFields::USER_FIELD},
  {"version", AppendFields::VERSION_FIELD},
  {"data", AppendFields::DATA_FIELD},
  {"user_defined", AppendFields::USER_DEFINED_FIELD},
  {"files", AppendFields::FILES_FIELD},
  {"curve_sets", AppendFields::CURVE_SETS_FIELD},
  {"library_data", AppendFields::LIBRARY_DATA_FIELD},
};
}  // namespace

std::vector<std::string> const CURVE_CATEGORIES = {"dependent", "independent"};
#define UNUSED(x) \
  (void)(x)  // TODO: axom likely has a better preference on how to suppress the warning

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
  return axom::fmt::format("[{}]", axom::fmt::join(supported_types, ", "));
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
  // file system as the destination file so that the move operation is atomic.

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

  // windows doesn't let you rename to a destination that already exists
  axom::utilities::filesystem::removeFile(fileName);

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

// Helper function--if a node has children, add them to the provided. Assumes both nodes are conduit lists
void concat_list_node(conduit::Node &concatTo, const conduit::Node &concatFrom)
{
  auto itr = concatFrom.children();
  while(itr.has_next())
  {
    concatTo.append() = itr.next();
  }
}

///////////////////// TEMPLATE HELPER BLOCK -- smooth differences between JSON and HDF5 /////////////////////
// This section exits because the hdf5 uses the format records/<some_num>/{actual_record}, whereas the JSON
// uses a list of records, which the pathlike relay interface doesn't seem to be abe to access children within a list.
conduit::Node relayLikeRead(conduit::Node &appendTo, const std::string &endpoint, int record_num)
{
  return appendTo["records"].child(record_num)[endpoint];
}

conduit::Node relayLikeRead(conduit::relay::io::IOHandle &appendTo,
                            const std::string &endpoint,
                            int record_num)
{
  UNUSED(record_num);
  conduit::Node holderNode;
  appendTo.read(endpoint, holderNode);
  return holderNode;
}

bool relayLikeHasPath(conduit::Node &appendTo, const std::string &endpoint, int record_num)
{
  return appendTo["records"].child(record_num).has_path(endpoint);
}

bool relayLikeHasPath(conduit::relay::io::IOHandle &appendTo,
                      const std::string &endpoint,
                      int record_num)
{
  UNUSED(record_num);
  return appendTo.has_path(endpoint);
}

bool nodeWorkaroundHasChildSlashes(conduit::Node &appendTo,
                                   const std::string &endpoint,
                                   const std::string &child_name,
                                   int record_num)
{
  return appendTo["records"].child(record_num)[endpoint].has_child(child_name);
}

// HDF5 already escapes the slashes, so we don't have to worry.
bool nodeWorkaroundHasChildSlashes(conduit::relay::io::IOHandle &appendTo,
                                   const std::string &endpoint,
                                   const std::string &child_name,
                                   int record_num)
{
  UNUSED(record_num);
  UNUSED(child_name);
  return appendTo.has_path(endpoint);
}

std::vector<std::string> relayLikeListChildNames(conduit::Node &appendTo,
                                                 const std::string &endpoint,
                                                 int record_num)
{
  return appendTo["records"].child(record_num)[endpoint].child_names();
}

std::vector<std::string> relayLikeListChildNames(conduit::relay::io::IOHandle &appendTo,
                                                 const std::string &endpoint,
                                                 int record_num)
{
  UNUSED(record_num);
  std::vector<std::string> nameHolder;
  appendTo.list_child_names(endpoint, nameHolder);
  return nameHolder;
}

int relayLikeNumChildren(conduit::Node &appendTo, const std::string &endpoint, int record_num)
{
  return appendTo["records"].child(record_num)[endpoint].number_of_children();
}

int relayLikeNumChildren(conduit::relay::io::IOHandle &appendTo,
                         const std::string &endpoint,
                         int record_num)
{
  UNUSED(record_num);
  std::vector<std::string> child_name_holder;
  appendTo.list_child_names(endpoint, child_name_holder);
  return child_name_holder.size();  // Feels incorrect, didn't find anything in conduit though
}

void relayLikeWrite(conduit::relay::io::IOHandle &appendTo, conduit::Node &appendFrom, const std::string &endpoint, int record_num){
  UNUSED(record_num);
  //conduit::Node sacrificial;
  //appendTo.read(endpoint, sacrificial);
  //std::cerr << "appending from: " << appendFrom.to_string() << std::endl;
  //std::cerr << "appending to: " << sacrificial.to_string() << std::endl;
  if(relayLikeHasPath(appendTo, endpoint, record_num)){
    appendTo.remove(endpoint);
  }
  appendTo.write(appendFrom, endpoint);
}

void relayLikeWrite(conduit::Node &appendTo, conduit::Node &appendFrom, const std::string &endpoint, int record_num){
    appendTo["records"].child(record_num)[endpoint].update(appendFrom);
}

void relayLikeAddNewRecord(conduit::relay::io::IOHandle &appendTo,
                           conduit::Node &new_record,
                           int new_record_num)
{
  relayLikeWrite(appendTo, new_record, "records/" + std::to_string(new_record_num), new_record_num);
}

void relayLikeAddNewRecord(conduit::Node &appendTo, conduit::Node &new_record, int new_record_num)
{
  UNUSED(new_record_num);
  appendTo["records"].append() = new_record;
}

std::unordered_map<std::string, int> relayLikeRecordOrderMap(conduit::relay::io::IOHandle &appendTo)
{
  std::unordered_map<std::string, int> order_map;
  conduit::Node n;
  std::vector<std::string> child_names;
  appendTo.list_child_names("records/", child_names);
  for(const std::string &child_name : child_names)
  {
    appendTo.read("records/" + child_name + "/id", n);
    order_map.insert(std::make_pair(n.to_string(), std::stoi(child_name)));
  }
  return order_map;
}

uint64_t relayLikeArrayNumElements(conduit::Node &appendTo, const std::string &endpoint, int record_num, const std::string &original_file_path)
{
  UNUSED(original_file_path);
  return appendTo["records"].child(record_num)[endpoint].dtype().number_of_elements();
}

uint64_t relayLikeArrayNumElements(conduit::relay::io::IOHandle &appendTo,
                                const std::string &endpoint,
                                int record_num,
                                const std::string &original_file_path)
{
  // This is the only reason why original_file_path has to be passed all the way down from append()
  // Couldn't find a quick way to get path from IOHandle, even though it probably holds it.
  UNUSED(appendTo);
  UNUSED(record_num);
  conduit::Node metadata_only;
  conduit::relay::io::hdf5_read_info(original_file_path, endpoint, metadata_only);
  return metadata_only["num_elements"].value();
}

void relayLikeAppendCurve(conduit::relay::io::IOHandle &appendTo, conduit::Node &appendFrom, const std::string &endpoint, int record_num,
        const std::string &original_file_path){
  UNUSED(record_num);
  conduit::Node OPTS_NODE;  //are there any static defaults we want to set?
  OPTS_NODE["offset"] = relayLikeArrayNumElements(appendTo, endpoint, record_num, original_file_path);
  appendTo.write(appendFrom, endpoint, OPTS_NODE);
}

void relayLikeAppendCurve(conduit::Node &appendTo, conduit::Node &appendFrom, const std::string &endpoint, int record_num,
    const std::string &original_file_path){
  UNUSED(original_file_path);
  conduit::Node &append_at = appendTo["records"].child(record_num)[endpoint];
  std::vector<double> merged_values(append_at.as_double_ptr(), append_at.as_double_ptr() + append_at.dtype().number_of_elements());
  merged_values.insert(merged_values.end(), appendFrom.as_double_ptr(), appendFrom.as_double_ptr() + appendFrom.dtype().number_of_elements());
  append_at.set(merged_values);
}

std::unordered_map<std::string, int> relayLikeRecordOrderMap(conduit::Node &appendTo){
  std::unordered_map<std::string, int> order_map;
  int num_children = appendTo["records"].number_of_children();
  for(int i = 0; i < num_children; i++)
  {
    order_map.insert(std::make_pair(appendTo["records"].child(i)["id"].to_string(), i));
  }
  return order_map;
}

///////////////////// TEMPLATE HELPER BLOCK END /////////////////////

// Specifically validate ONE curve set for ONE DataHolder (record, library_data...) for appending,
// appendTo is notionally const, but the has_path() etc. methods aren't const.
template <typename ConduitRelayLike>
conduit::Node validateCurveSets(ConduitRelayLike &appendTo,
                                const conduit::Node &appendFrom,
                                const std::string &endpoint,
                                int rec_num,
                                const std::string &original_file_path)
{
  int baseline = -1;  // baseline is shared across dependent and independent
  conduit::Node msgNode = conduit::Node(conduit::DataType::list());
  // First, we need to make sure we didn't "forget" an existing curve: make sure either EVERY
  // curve in the HDF5 is being appended to, or NONE of them are (all curves are new).
  u_int curves_written = 0;
  u_int existing_curves = 0;
  int unappended_baseline = -1;  // Find length of anything we don't append to, for later.
  for(const std::string &curve_cat : CURVE_CATEGORIES)
  {
    std::string curves_endpoint = endpoint + "/" + curve_cat;
    std::vector<std::string> curve_names =
      relayLikeListChildNames(appendTo, curves_endpoint, rec_num);
    for(const std::string &cname : curve_names)
    {
      if(cname == "value" || cname == "tags" || cname == "units")
      {
        // The way we list child names is recursive. thus, if someone
        // has a curve named "value" (or tags or units...) we can't properly verify. TODO: revisit with more conduit
        continue;
      }
      if(appendFrom.has_path(curve_cat + "/" + cname) &&
         appendFrom[curve_cat][cname]["value"].dtype().number_of_elements() > 0)
      {
        curves_written++;
      }
      else if(unappended_baseline == -1)
      { 
        unappended_baseline = relayLikeArrayNumElements(appendTo, curves_endpoint + "/" + cname +"/value", rec_num, original_file_path);
      }
      existing_curves++;  // instead of a .size() to account for value/tags/etc. above
    }
  }
  if(curves_written != 0 && curves_written != existing_curves)
  {
    msgNode.append() = "Failed to append record " + std::to_string(rec_num) +
      ": did not append ALL or NO pre-existing curves (causing append element count mismatch)";
  }

  for(const std::string &curve_cat : CURVE_CATEGORIES)
  {
    // Now loop through what we've actually got. Once we find something, use it to set the baseline.
    if(appendFrom.has_child(curve_cat))
    {
      auto curvesIter = appendFrom[curve_cat].children();
      while(curvesIter.has_next())
      {
        const conduit::Node &testCurve = curvesIter.next()["value"];
        int post_append_size = testCurve.dtype().number_of_elements();
        std::string sub_endpoint = endpoint + "/" + curve_cat + "/" + curvesIter.name() + "/value";
        if(relayLikeHasPath(appendTo, sub_endpoint, rec_num))
        {
          int num_elements = relayLikeArrayNumElements(appendTo, sub_endpoint, rec_num, original_file_path);
          post_append_size += num_elements;
        }
        if(baseline == -1)
        {
          baseline = post_append_size;
        }
        if(post_append_size != baseline ||
           (unappended_baseline != -1 && baseline != unappended_baseline))
        {
          msgNode.append() = "Failed to append record " + std::to_string(rec_num) + "'s curve '" +
            curvesIter.name() + "': count of appended elements would differ between series";
        }
      }
    }
  }
  return msgNode;
}

// Top-level append validation function. Works recursively on library data (hence endpoint)
template <typename ConduitRelayLike>
conduit::Node validateAppendDocument(ConduitRelayLike &appendTo,
                                     const conduit::Node &appendFrom,
                                     const std::string &endpoint,
                                     const int mergeProtocol,
                                     const int record_num,
                                     const std::string &original_file_path)
{
  conduit::Node msgNode = conduit::Node(conduit::DataType::list());
  // Case one: die if the types disagree. A pingpong_game shouldn't become a billiards_game
  // Validation note: we allow for no type here beecause of library_data. If someone
  // forgot the type for a top-level record, it should've died already.
  if(appendFrom.has_child("type"))
  {
    conduit::Node typeNode;
    typeNode = relayLikeRead(appendTo, endpoint + "/type", record_num);
    if(typeNode.to_string() != appendFrom["type"].to_string())
    {
      msgNode.append() = "Failed to append record " + std::to_string(record_num) +
        ": type mismatch " + typeNode.to_string() + "vs " + appendFrom["type"].to_string();
    }
  }
    // Case two: merge protocol is 3, we need to die if certain fields appear in both places
    if(mergeProtocol == 3){
      const std::vector<std::string> prot3Fields = {"data", "user_defined", "files"};
      for(auto &field : prot3Fields){
        if(appendFrom.has_child(field) && relayLikeHasPath(appendTo, endpoint + "/" + field + "/", record_num)){
          auto dataIter = appendFrom[field].children();
          while(dataIter.has_next()){
            dataIter.next();
            if(nodeWorkaroundHasChildSlashes(appendTo, endpoint + "/" + field + "/", dataIter.name(), record_num)){
              msgNode.append() = "Failed to append record " + std::to_string(record_num) + endpoint + " (protocol 3): conflicting "+field+": "+dataIter.name();
            }
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
      const conduit::Node &n = curveSetsIter.next();
      subEndpoint = endpoint + "/curve_sets/" + curveSetsIter.name();
      // We only have to validate if the hdf5 already has a curve set with that name.
      if(relayLikeHasPath(appendTo, subEndpoint, record_num))
      {
        concat_list_node(msgNode, validateCurveSets(appendTo, n, subEndpoint, record_num, original_file_path));
      }
    }
  }

    // Case four: library data. Recurse on it if it's already in the hdf5.
    if(appendFrom.has_child("library_data")){
      auto libraryIter = appendFrom["library_data"].children();
      std::string subEndpoint;
      while(libraryIter.has_next()){
        const conduit::Node &n = libraryIter.next();
        subEndpoint = endpoint + "/library_data/" + libraryIter.name();
        // We only have to validate if the target already has a library with that name.
        if(relayLikeHasPath(appendTo, subEndpoint, record_num)){
          concat_list_node(msgNode,
                           validateAppendDocument(appendTo, n, subEndpoint, mergeProtocol, record_num));
        }
      }
    }
    return msgNode;
}

// Avoiding a terrible if/else chunk in append_recordlike_fields and friends.
AppendFields field_lookup(const std::string &input)
{
  auto itr = appendFieldStrings.find(input);
  if(itr != appendFieldStrings.end())
  {
    return itr->second;
  }
  return AppendFields::UNKNOWN_FIELD;
}

template <typename ConduitRelayLike>
void append_curveset(ConduitRelayLike &appendTo,
                     conduit::Node &appendFrom,
                     const std::string &endpoint,
                     int record_num,
                     const std::string &original_file_path){
    for(const std::string& curve_cat : CURVE_CATEGORIES){
      auto curveIter = appendFrom[curve_cat].children();
      while(curveIter.has_next()){{
        conduit::Node &n = curveIter.next();
        std::string curve_endpoint = endpoint + "/" + curve_cat + "/" + curveIter.name() + "/value";
        if(relayLikeHasPath(appendTo, curve_endpoint, record_num)){
          relayLikeAppendCurve(appendTo, n["value"], curve_endpoint, record_num, original_file_path);
        } else {
          relayLikeWrite(appendTo, n["value"], curve_endpoint, record_num);
        }
      }
    }
  }
}

template <typename ConduitRelayLike>
void append_recordlike_fields(ConduitRelayLike &appendTo,
                              conduit::Node &appendFrom,
                              const std::string &endpoint,
                              const int mergeProtocol,
                              int record_num,
                              const std::string &original_file_path,
                              bool isHDF5)
{
  auto fieldsIter = appendFrom.children();
  while(fieldsIter.has_next())
  {
    conduit::Node &recField = fieldsIter.next();
    std::string appendAtEndpoint = endpoint+"/"+fieldsIter.name() + "/";
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
      if(mergeProtocol == 1 || mergeProtocol == 3)
      {
        if(isHDF5){ // We have to go bit by bit removing each individually, lest we delete something unupdated
          auto subFieldIter = appendFrom[fieldsIter.name()].children();
          while(subFieldIter.has_next())
          {
            conduit::Node &subField = subFieldIter.next();
            relayLikeWrite(appendTo, subField, appendAtEndpoint+subFieldIter.name(), record_num);
          }
        } else { // We can just update everything at a go
          relayLikeWrite(appendTo, recField, appendAtEndpoint, record_num);
        }
      }
      else if(mergeProtocol == 2)
      {
        auto subFieldIter = appendFrom[fieldsIter.name()].children();
        while(subFieldIter.has_next())
        {
          conduit::Node &subField = subFieldIter.next();
          if(!relayLikeHasPath(appendTo, appendAtEndpoint, record_num)){
            relayLikeWrite(appendTo, subField, appendAtEndpoint+subFieldIter.name(), record_num);
          }
        }
      }
      break;
    case AppendFields::LIBRARY_DATA_FIELD:
    {
      // We recurse
      auto libraryIter = appendFrom[fieldsIter.name()].children();
      std::string appendAtEndpoint = endpoint + "/library_data/";
      while(libraryIter.has_next())
      {
        conduit::Node &libraryField = libraryIter.next();
        if(relayLikeHasPath(appendTo, appendAtEndpoint+libraryIter.name(), record_num)){
          append_recordlike_fields(appendTo, libraryField,
                                   appendAtEndpoint+libraryIter.name(),
                                   mergeProtocol, record_num, original_file_path, isHDF5);
        } else {
          relayLikeWrite(appendTo, libraryField, appendAtEndpoint+libraryIter.name(), record_num);
        }
      }
      break;
    }
    case AppendFields::CURVE_SETS_FIELD:
    {
      auto curveSetIter = appendFrom[fieldsIter.name()].children();
      while(curveSetIter.has_next())
      {
        conduit::Node &curveSetField = curveSetIter.next();
        append_curveset(appendTo, curveSetField,
                        appendAtEndpoint+curveSetIter.name(),
                        record_num, original_file_path);
      }
    }
    break;
    default:
      std::cout << "oh no" << std::endl;
    }
  }
}

template <typename ConduitRelayLike>
conduit::Node append(ConduitRelayLike &appendTo,
                     conduit::Node &appendFrom,
                     const int mergeProtocol,
                     bool isHDF5,
                     const std::string &original_file_path)
{
  conduit::Node msgNode = conduit::Node(conduit::DataType::list());
  // We need to figure out where each record is in appendTo, since there's no guarantee in the order
  // of appendFrom.
  std::unordered_map<std::string, int> rec_order = relayLikeRecordOrderMap(appendTo);
  // We do all of our validation up-front, including/especially library data.
  // There's no validation to perform on things like relationships (right..?)
  auto recordsIter = appendFrom["records"].children();
  while(recordsIter.has_next())
  {
    conduit::Node &n = recordsIter.next();
    auto rec_num = rec_order.find(n["id"].to_string());
    // We only validate records we're appending (not just adding). This does mean someone could insert a malformed record,
    // but that's always been the case; validation is meant to assist with catching bad curves, mostly
    if(rec_num != rec_order.end()){
      std::string endpoint = isHDF5 ? "records/" + std::to_string(rec_num->second): "";
      concat_list_node(msgNode, validateAppendDocument(appendTo, n, endpoint, mergeProtocol, rec_num->second, original_file_path));
    }
  }
  // Return with our error list if we errored.
  if(msgNode.number_of_children() > 0)
  {
    return msgNode;
  }

  // Our validation passed, time to throw it all in!
  recordsIter = appendFrom["records"].children();
  int offset = rec_order.size();
  while(recordsIter.has_next())
  {
    conduit::Node &rec = recordsIter.next();
    // Easiest case, the record doesn't exist yet. Add it.
    auto rec_num = rec_order.find(rec["id"].to_string());
    std::string endpoint = isHDF5 ? "records/" + std::to_string(rec_num->second): "";
    if(rec_num == rec_order.end()){
      relayLikeAddNewRecord(appendTo, rec, offset);
      offset ++;
    } else {
      append_recordlike_fields(appendTo, rec, endpoint, mergeProtocol, rec_num->second, original_file_path, isHDF5);
    }
  }
  return msgNode;
}

conduit::Node appendDocumentToJson(const std::string &jsonFilePath,
                    const Document &newData,
                    const int mergeProtocol){
  conduit::Node appendTo;
  appendTo.load(jsonFilePath, "json");
  conduit::Node appendFrom = newData.toNode();
  conduit::Node msgNode = append(appendTo, appendFrom, mergeProtocol, false, jsonFilePath);
  conduit::relay::io::save(appendTo, jsonFilePath);
  return msgNode;
}

conduit::Node appendDocumentToHDF5(const std::string &hdf5FilePath,
                                   const Document &newData,
                                   const int mergeProtocol)
{
#ifdef AXOM_USE_HDF5
  conduit::relay::io::IOHandle appendTo;
  appendTo.open(hdf5FilePath);
  conduit::Node appendFrom = newData.toHDF5Node();
  conduit::Node msgNode = append(appendTo, appendFrom, mergeProtocol, true, hdf5FilePath);
  appendTo.close();
  return msgNode;
#endif
#ifndef AXOM_USE_HDF5
  conduit::Node msgNode = conduit::Node(conduit::DataType::list());
  msgNode.append("Failed to append Sina HDF5: Axom wasn't build with HDF5");
  return msgNode;
#endif
}

}  // namespace sina
}  // namespace axom
