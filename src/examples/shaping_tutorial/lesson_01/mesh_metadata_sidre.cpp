// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
/*!
 * \file mesh_metadata_sidre.cpp
 * \brief An example that uses Sidre to describe a 2D Cartesian mesh.
 *
 * This example demonstrates creating a Sidre DataStore and using Groups and Views
 * to represent mesh metadata including bounding box coordinates and resolution.
 * It also optionally generates the mesh following the conduit mesh blueprint convention.
 * 
 * Example run:
 * ./mesh_metadata_sidre --min_x 0.0 --min_y 0.0 --max_x 2.0 --max_y 3.0 --res_x 20 --res_y 30 [-o]
 */
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/sidre.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

// Axom's config.hpp has defines for the dependencies it was configured against
#if defined(AXOM_USE_CONDUIT)
  #include "conduit.hpp"
  #include "conduit_blueprint.hpp"
  #include "conduit_relay.hpp"
  #include "conduit_relay_io_blueprint.hpp"
#else
  #error Bad configuration: Conduit is a hard dependency for Sidre
#endif

#include <iostream>

/*!
 * \struct Input
 * \brief Struct representing user input parameters for mesh bounding box and resolution.
 *
 * This struct holds the minimum and maximum x and y coordinates defining the bounding box,
 * as well as the resolution (number of cells) in both x and y directions.
 */
struct Input
{
  double min_x = 0.0;
  double min_y = 0.0;
  double max_x = 1.0;
  double max_y = 1.0;
  int res_x = 10;
  int res_y = 20;
};

/*!
 * \brief Sets up the Sidre hierarchy for mesh metadata using the provided input parameters.
 *
 * This function creates a hierarchical structure within the Sidre DataStore to represent
 * mesh metadata, including bounding box coordinates (min and max) and resolution in both
 * x and y directions. It stores these parameters in Groups and Views for easy access and
 * manipulation.
 *
 * \param datastore Reference to the Sidre DataStore to populate.
 * \param input Struct containing bounding box coordinates and resolution data.
 */
void setup_mesh_metadata(axom::sidre::DataStore& datastore, const Input& input)
{
  // Create a root group for the mesh metadata
  axom::sidre::Group* meshGroup = datastore.getRoot()->createGroup("mesh");

  // Create bounding box groups and views
  axom::sidre::Group* minGroup = meshGroup->createGroup("bounding_box/min");
  minGroup->createViewScalar("x", input.min_x);
  minGroup->createViewScalar("y", input.min_y);

  axom::sidre::Group* maxGroup = meshGroup->createGroup("bounding_box/max");
  maxGroup->createViewScalar("x", input.max_x);
  maxGroup->createViewScalar("y", input.max_y);

  // Create resolution group and views
  axom::sidre::Group* resGroup = meshGroup->createGroup("resolution");
  resGroup->createViewScalar("x", input.res_x);
  resGroup->createViewScalar("y", input.res_y);
}

/*!

 * \brief Verifies mesh metadata stored in Sidre DataStore
 *
 * This function accesses the mesh group and verifies that all required metadata
 * (bounding box coordinates and resolution) exists.
 *
 * \param meshGroup Pointer to the Sidre Group representing the mesh metadata.
 * \return True if all metadata is present and valid, false otherwise.
 */
bool verify_mesh_metadata(axom::sidre::Group* meshGroup)
{
  if(!meshGroup)
  {
    SLIC_WARNING("Missing mesh group");
    return false;
  }

  bool valid = true;

  // check bounding_box group
  if(meshGroup->hasGroup("bounding_box"))
  {
    axom::sidre::Group* bbGroup = meshGroup->getGroup("bounding_box");
    for(const std::string path : {"min/x", "min/y", "max/x", "max/y"})
    {
      if(!bbGroup->hasView(path))
      {
        SLIC_WARNING(axom::fmt::format("Missing '{}' view in 'bounding_box' group", path));
        valid = false;
      }
    }
  }
  else
  {
    valid = true;
    SLIC_WARNING("Missing 'bounding_box' group in mesh metadata");
  }

  // check resolution group
  if(meshGroup->hasGroup("resolution"))
  {
    axom::sidre::Group* resGroup = meshGroup->getGroup("resolution");
    for(const std::string path : {"x", "y"})
    {
      if(!resGroup->hasView("x"))
      {
        SLIC_WARNING(axom::fmt::format("Missing '{}' view in 'resolution' group", path));
        valid = false;
      }
    }
  }
  else
  {
    SLIC_WARNING("Missing 'resolution' group in mesh metadata");
    return false;
  }

  return valid;
}

/*!
 * \brief Prints mesh metadata stored in Sidre DataStore
 *
 * This function accesses the mesh group and prints bounding box coordinates and resolution.
 *
 * \param meshGroup Pointer to the Sidre Group representing the mesh metadata.
 */
void print_mesh_metadata(axom::sidre::Group* meshGroup)
{
  if(!meshGroup) return;

  axom::sidre::Group* bbGroup = meshGroup->getGroup("bounding_box");
  if(!bbGroup) return;

  axom::sidre::Group* resGroup = meshGroup->getGroup("resolution");
  if(!resGroup) return;

  SLIC_INFO(axom::fmt::format("Bounding Box Min: ({}, {})",
                              bbGroup->getView("min/x")->getData<double>(),
                              bbGroup->getView("min/y")->getData<double>()));

  SLIC_INFO(axom::fmt::format("Bounding Box Max: ({}, {})",
                              bbGroup->getView("max/x")->getData<double>(),
                              bbGroup->getView("max/y")->getData<double>()));

  SLIC_INFO(axom::fmt::format("Resolution: ({}, {})",
                              resGroup->getView("x")->getData<int>(),
                              resGroup->getView("y")->getData<int>()));
}

/*!
 * \brief Converts Sidre mesh metadata to Conduit mesh blueprint
 *
 * This function extracts mesh metadata from Sidre groups and views and creates
 * a Conduit Node that represents a uniform 2D Cartesian mesh following
 * the Conduit mesh blueprint conventions.
 *
 * \param meshGroup Pointer to the Sidre Group containing mesh metadata
 * \return A Conduit Node containing the mesh blueprint representation
 */
conduit::Node create_mesh_blueprint(axom::sidre::Group* meshGroup)
{
  conduit::Node blueprint;

  if(!meshGroup)
  {
    SLIC_ERROR("Invalid mesh group provided");
    return blueprint;
  }

  // Get bounding box information
  axom::sidre::Group* bbGroup = meshGroup->getGroup("bounding_box");
  if(!bbGroup)
  {
    SLIC_ERROR("Missing bounding_box group in mesh metadata");
    return blueprint;
  }

  double x_min = bbGroup->getView("min/x")->getData<double>();
  double y_min = bbGroup->getView("min/y")->getData<double>();
  double x_max = bbGroup->getView("max/x")->getData<double>();
  double y_max = bbGroup->getView("max/y")->getData<double>();

  // Get resolution information
  axom::sidre::Group* resGroup = meshGroup->getGroup("resolution");
  if(!resGroup)
  {
    SLIC_ERROR("Missing resolution group in mesh metadata");
    return blueprint;
  }

  const int res_x = resGroup->getView("x")->getData<int>();
  const int res_y = resGroup->getView("y")->getData<int>();
  SLIC_ERROR_IF(res_x < 1,
                axom::fmt::format("Resolution in x-coordinate ({}) must be positive", res_x));
  SLIC_ERROR_IF(res_y < 1,
                axom::fmt::format("Resolution in y-coordinate ({}) must be positive", res_y));

  // Create coordset
  blueprint["coordsets/coords/type"] = "uniform";
  blueprint["coordsets/coords/dims/i"] = res_x + 1;
  blueprint["coordsets/coords/dims/j"] = res_y + 1;
  blueprint["coordsets/coords/origin/x"] = x_min;
  blueprint["coordsets/coords/origin/y"] = y_min;
  blueprint["coordsets/coords/spacing/dx"] = (x_max - x_min) / res_x;
  blueprint["coordsets/coords/spacing/dy"] = (y_max - y_min) / res_y;

  // Create topology
  blueprint["topologies/mesh/type"] = "uniform";
  blueprint["topologies/mesh/coordset"] = "coords";

  return blueprint;
}

/*!
 * \brief Saves a Conduit Node blueprint to a file
 *
 * This function verifies that the provided Node conforms to the Conduit mesh blueprint
 * specification and saves it to the specified file format if valid.
 *
 * \param blueprint Conduit Node containing the mesh blueprint to save
 * \param file_name Base name for the output file (without extension)
 * \param file_format Format to save the file in (yaml, json, etc.)
 * \return True if save was successful, false otherwise
 */
bool save_blueprint(const conduit::Node& blueprint,
                    const std::string& file_name,
                    const std::string& file_format = "yaml")
{
  if(blueprint.number_of_children() > 0)
  {
    SLIC_INFO(axom::fmt::format("Saving mesh blueprint to '{}.{}'", file_name, file_format));
    conduit::Node info;
    if(conduit::blueprint::mesh::verify(blueprint, info))
    {
      conduit::relay::io::blueprint::save_mesh(blueprint, file_name, file_format);
      return true;
    }
    else
    {
      SLIC_ERROR("Blueprint verification failed \n" << info.to_string());
      return false;
    }
  }
  else
  {
    SLIC_ERROR("Unable to create mesh blueprint, cannot save file");
    return false;
  }
}

int main(int argc, char** argv)
{
  // initialize SLIC logger
  axom::slic::SimpleLogger logger;

  // parse input using CLI11
  Input input;
  bool output_blueprint = false;
  axom::CLI::App app {"Mesh Metadata Setup"};
  app.add_option("--min_x", input.min_x, "Minimum x coordinate of bounding box");
  app.add_option("--min_y", input.min_y, "Minimum y coordinate of bounding box");
  app.add_option("--max_x", input.max_x, "Maximum x coordinate of bounding box");
  app.add_option("--max_y", input.max_y, "Maximum y coordinate of bounding box");
  app.add_option("--res_x", input.res_x, "Resolution in x direction");
  app.add_option("--res_y", input.res_y, "Resolution in y direction");
  app.add_flag("-o,--output_blueprint", output_blueprint, "Output mesh blueprint to file");
  CLI11_PARSE(app, argc, argv);

  // load parsed data into sidre datastore
  axom::sidre::DataStore datastore;
  setup_mesh_metadata(datastore, input);

  // validate and print results
  auto* root = datastore.getRoot();
  SLIC_ASSERT_MSG(root->hasGroup("mesh"), "Missing expected 'mesh' group");
  const bool is_valid = verify_mesh_metadata(root->getGroup("mesh"));
  if(is_valid)
  {
    SLIC_INFO("Sidre hierarchy was properly set up");
    print_mesh_metadata(root->getGroup("mesh"));
  }
  else
  {
    SLIC_WARNING("Sidre hierarchy was not properly set up");
    std::stringstream sstr;
    root->printTree(0, sstr);
    SLIC_INFO("Sidre hierarchy: \n" << sstr.str());
  }

  // Optionally, create mesh blueprint and save as a yaml file
  if(output_blueprint)
  {
    conduit::Node blueprint = create_mesh_blueprint(root->getGroup("mesh"));
    save_blueprint(blueprint, "uniform_bp");
  }

  return is_valid ? 0 : 1;
}
