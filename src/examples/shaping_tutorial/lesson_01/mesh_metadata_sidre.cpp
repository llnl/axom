#include "axom/config.hpp"
#include "axom/sidre.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include <iostream>

/*!
 * \file basic_sidre_example.cpp
 * \brief A basic example of how to use Sidre to describe a mesh object.
 *
 * This example demonstrates creating a Sidre DataStore and using Groups and Views
 * to represent mesh metadata including bounding box coordinates and resolution.
 * 
 * 
 * Example run:
 * ./basic_sidre_example --min_x 0.0 --min_y 0.0 --max_x 2.0 --max_y 3.0 --res_x 20 --res_y 30
 */

/*!
 * \struct Input
 * \brief Struct representing user input parameters for mesh bounding box and resolution.
 *
 * This struct holds the minimum and maximum x and y coordinates defining the bounding box,
 * as well as the resolution (number of divisions) in both x and y directions.
 * These parameters configure the mesh metadata within the application.
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
void setup_mesh(axom::sidre::DataStore& datastore, const Input& input)
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
                              bbGroup->getView("max/x")->getData<double>()));

  SLIC_INFO(axom::fmt::format("Resolution: ({}, {})",
                              resGroup->getView("x")->getData<int>(),
                              resGroup->getView("y")->getData<int>()));
}

int main(int argc, char** argv)
{
  // initialize SLIC logger
  axom::slic::SimpleLogger logger;

  // parse input
  Input input;
  axom::CLI::App app {"Mesh Metadata Setup"};
  app.add_option("--min_x", input.min_x, "Minimum x coordinate of bounding box");
  app.add_option("--min_y", input.min_y, "Minimum y coordinate of bounding box");
  app.add_option("--max_x", input.max_x, "Maximum x coordinate of bounding box");
  app.add_option("--max_y", input.max_y, "Maximum y coordinate of bounding box");
  app.add_option("--res_x", input.res_x, "Resolution in x direction");
  app.add_option("--res_y", input.res_y, "Resolution in y direction");
  CLI11_PARSE(app, argc, argv);

  // load parsed data into sidre datastore
  axom::sidre::DataStore datastore;
  setup_mesh(datastore, input);

  // print results
  auto* root = datastore.getRoot();
  SLIC_ASSERT_MSG(root->hasGroup("mesh"), "Missing expected 'mesh' group");
  print_mesh_metadata(root->getGroup("mesh"));

  return 0;
}
