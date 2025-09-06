#include "axom/config.hpp"
#include "axom/sidre.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include <iostream>

struct Input
{
  double min_x = 0.0;
  double min_y = 0.0;
  double max_x = 1.0;
  double max_y = 1.0;
  int res_x = 10;
  int res_y = 20;
};

void setup_mesh(axom::sidre::DataStore& datastore, const Input& input)
{
  // Create a root group for the mesh metadata
  axom::sidre::Group* meshGroup = datastore.getRoot()->createGroup("mesh_metadata");

  // Create groups for bounding box min and max
  axom::sidre::Group* minGroup = meshGroup->createGroup("bounding_box/min");
  axom::sidre::Group* maxGroup = meshGroup->createGroup("bounding_box/max");

  // Create views for min and max coordinates (x and y as doubles)
  axom::sidre::View* minXView = minGroup->createViewScalar("x", input.min_x);
  axom::sidre::View* minYView = minGroup->createViewScalar("y", input.min_y);
  axom::sidre::View* maxXView = maxGroup->createViewScalar("x", input.max_x);
  axom::sidre::View* maxYView = maxGroup->createViewScalar("y", input.max_y);

  // Create group for resolution
  axom::sidre::Group* resGroup = meshGroup->createGroup("resolution");

  // Create views for resolution (x and y as integers)
  axom::sidre::View* resXView = resGroup->createViewScalar("x", input.res_x);
  axom::sidre::View* resYView = resGroup->createViewScalar("y", input.res_y);

  SLIC_INFO(axom::fmt::format("Bounding Box Min: ({}, {})",
                              minXView->getData<double>(),
                              maxYView->getData<double>()));

  SLIC_INFO(axom::fmt::format("Bounding Box Max: ({}, {})",

                              minXView->getData<double>(),
                              maxYView->getData<double>()));

  SLIC_INFO(
    axom::fmt::format("Resolution: ({}, {})", resXView->getData<int>(), resYView->getData<int>()));
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;  // Install and initialize the SLIC logging system

  axom::CLI::App app {"Mesh Metadata Setup"};

  Input input;

  app.add_option("--min_x", input.min_x, "Minimum x coordinate of bounding box");
  app.add_option("--min_y", input.min_y, "Minimum y coordinate of bounding box");
  app.add_option("--max_x", input.max_x, "Maximum x coordinate of bounding box");
  app.add_option("--max_y", input.max_y, "Maximum y coordinate of bounding box");
  app.add_option("--res_x", input.res_x, "Resolution in x direction");
  app.add_option("--res_y", input.res_y, "Resolution in y direction");

  CLI11_PARSE(app, argc, argv);

  axom::sidre::DataStore datastore;

  setup_mesh(datastore, input);

  return 0;
}
