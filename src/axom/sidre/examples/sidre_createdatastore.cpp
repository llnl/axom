// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file sidre_createdatastore.cpp
 *  \brief This example code is a basic demonstration of how to use the
 *  Sidre classes to construct a hierarchical data store and save it to
 *  a file.  (This example does not use the parallel I/O facility.)
 *  This also shows how to use Sidre with Conduit to generate and save a
 *  data store using the Mesh Blueprint, for easy data exchange and visualization.
 */

/*
 * Summary of the functions called by main this example:
 *
 * create_datastore           -- Builds a DataStore containing various Groups,
 *                               Views, Buffers, and Attributes
 *
 * create_tiny_datastore      -- Builds a DataStore that holds a small 3D mesh
 *                               conforming to the conduit mesh blueprint
 *
 * access_datastore           -- Exercises access methods that get data from
 *                               the DataStore created by create_datastore
 *
 * iterate_datastore          -- Exercises iterator methods for views, groups,
 *                               buffers and attributes in the DataStore
 *                               created by create_datastore
 *
 * save_as_blueprint          -- Generates a blueprint index and writes data
 *                               to file using direct calls to conduit
 *
 * generate_blueprint         -- Generates a blueprint index and writes data
 *                               to file using Sidre DataStore and Group calls
 *
 * generate_blueprint_to_path -- Similar to above, but uses path arguments to
 *                               find data deeper in Sidre hierarchy
 *
 * generate_spio_blueprint    -- Uses SPIO/IOManager calls to write data to
 *                               file and generate blueprint index
 *
 * generate_blueprint_to_path -- Similar to above, but uses path arguments to
 *                               find data deeper in Sidre hierarchy
 */

/* This example code contains snippets used in the Sidre Sphinx documentation.
 * They begin and end with comments
 *
 * first_example_creategroups_start
 * first_example_creategroups_end
 * first_example_state_start
 * first_example_state_end
 * first_example_nodes_start
 * first_example_nodes_end
 * first_example_fields_start
 * first_example_fields_end
 * first_example_access_start
 * first_example_access_end
 * serial_io_save_start
 * serial_io_save_end
 * tiny_create_start
 * tiny_create_end
 * blueprint_restructure_save_start
 * blueprint_restructure_save_end
 * blueprint_save_start
 * blueprint_save_end
 *
 * each prepended with an underscore.
 */

// Axom headers
#include "axom/core.hpp"
#include "axom/sidre.hpp"
#include "axom/fmt.hpp"

// Conduit headers
#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"

#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
#endif

// C++ headers
#include <cstring>

// "using" directives to simplify code
namespace sidre = axom::sidre;

std::unique_ptr<sidre::DataStore> create_datastore(int* region)
{
  // _first_example_creategroups_start
  // Create Sidre datastore object and get root group
  auto ds = std::unique_ptr<sidre::DataStore> {new sidre::DataStore()};
  sidre::Group* root = ds->getRoot();

  // Create two attributes
  ds->createAttributeScalar("vis", 0);
  ds->createAttributeScalar("restart", 1);

  // Create group children of root group
  sidre::Group* state = root->createGroup("state");
  sidre::Group* nodes = root->createGroup("nodes");
  sidre::Group* fields = root->createGroup("fields");
  // _first_example_creategroups_end

  // _first_example_state_start
  // Populate "state" group
  state->createViewScalar("cycle", 25);
  state->createViewScalar("time", 1.2562e-2);
  state->createViewString("name", "sample_20171206_a");
  // _first_example_state_end

  // _first_example_nodes_start
  int N = 16;
  int nodecount = N * N * N;
  int eltcount = (N - 1) * (N - 1) * (N - 1);

  // Populate "nodes" group
  //
  // "x", "y", and "z" are three views into a shared Sidre buffer object that
  // holds 3 * nodecount doubles.  These views might describe the location of
  // each node in a 16 x 16 x 16 hexahedron mesh.  Each view is described by
  // number of elements, offset, and stride into that data.
  sidre::Buffer* buff = ds->createBuffer(sidre::DOUBLE_ID, 3 * nodecount)->allocate();
  nodes->createView("x", buff)->apply(sidre::DOUBLE_ID, nodecount, 0, 3);
  nodes->createView("y", buff)->apply(sidre::DOUBLE_ID, nodecount, 1, 3);
  nodes->createView("z", buff)->apply(sidre::DOUBLE_ID, nodecount, 2, 3);
  // _first_example_nodes_end

  // _first_example_fields_start
  // Populate "fields" group
  //
  // "temp" is a view into a buffer that is not shared with another View.
  // In this case, the data Buffer is allocated directly through the View
  // object.  Likewise with "rho."  Both Views have the default offset (0)
  // and stride (1).  These Views could point to data associated with
  // each of the 15 x 15 x 15 hexahedron elements defined by the nodes above.
  sidre::View* temp = fields->createViewAndAllocate("temp", sidre::DOUBLE_ID, eltcount);
  sidre::View* rho = fields->createViewAndAllocate("rho", sidre::DOUBLE_ID, eltcount);

  // Explicitly set values for the "vis" Attribute on the "temp" and "rho"
  // buffers.
  temp->setAttributeScalar("vis", 1);
  rho->setAttributeScalar("vis", 1);

  // The "fields" Group also contains a child Group "ext" which holds a pointer
  // to an externally owned integer array.  Although Sidre does not own the
  // data, the data can still be described to Sidre.
  sidre::Group* ext = fields->createGroup("ext");
  // int * region has been passed in as a function argument.  As with "temp"
  // and "rho", view "region" has default offset and stride.
  ext->createView("region", region)->apply(sidre::INT_ID, eltcount);
  // _first_example_fields_end

  return ds;
}

void access_datastore(sidre::DataStore* ds)
{
  // _first_example_access_start
  // Retrieve Group pointers
  sidre::Group* root = ds->getRoot();
  sidre::Group* state = root->getGroup("state");
  sidre::Group* nodes = root->getGroup("nodes");
  sidre::Group* fields = root->getGroup("fields");

  // Accessing a Group that is not there gives a null pointer
  // Requesting a nonexistent View also gives a null pointer
  sidre::Group* goofy = root->getGroup("goofy");
  if(goofy == nullptr)
  {
    std::cout << "no such group: goofy" << std::endl;
  }
  else
  {
    std::cout << "Something is very wrong!" << std::endl;
  }

  // Access items in "state" group
  int cycle = state->getView("cycle")->getScalar();
  double time = state->getView("time")->getScalar();
  const char* name = state->getView("name")->getString();

  // Access some items in "nodes" and "fields" groups
  double* y = nodes->getView("y")->getArray();
  int ystride = nodes->getView("y")->getStride();
  double* temp = fields->getView("temp")->getArray();
  int* region = fields->getView("ext/region")->getArray();

  // Nudge the 3rd node, adjust temp and region of the 3rd element
  y[2 * ystride] += 0.0032;
  temp[2] *= 1.0021;
  region[2] = 6;
  // _first_example_access_end

  // Deal with unused variables
  AXOM_UNUSED_VAR(cycle);
  AXOM_UNUSED_VAR(time);
  AXOM_UNUSED_VAR(name);
  AXOM_UNUSED_VAR(y);
  AXOM_UNUSED_VAR(temp);
  AXOM_UNUSED_VAR(region);
}

void iterate_datastore(sidre::DataStore* ds)
{
  const std::string fill_line = axom::fmt::format("{:=^80}", "");

  std::cout << fill_line << std::endl;

  // iterate through the attributes in ds
  std::cout << "The datastore has the following attributes:\n";
  for(auto& attr : ds->attributes())
  {
    std::cout << axom::fmt::format("  * [{}] '{}' of type {} and default value: {}\n",
                                   attr.getIndex(),
                                   attr.getName(),
                                   conduit::DataType::id_to_name(attr.getTypeID()),
                                   attr.getDefaultNodeRef().to_yaml());
  }

  std::cout << fill_line << std::endl;

  // iterate through the buffers in ds
  std::cout << "The datastore has the following buffers:\n";
  for(auto& buff : ds->buffers())
  {
    std::cout << axom::fmt::format("  * [{}] {} buffer with {} elements of type {} with {} views\n",
                                   buff.getIndex(),
                                   buff.isAllocated() ? "Allocated" : "Unallocated",
                                   buff.getNumElements(),
                                   conduit::DataType::id_to_name(buff.getTypeID()),
                                   buff.getNumViews());
  }

  std::cout << fill_line << std::endl;

  // iterate through the groups of the root group
  std::cout << "The root group has the following groups:\n";
  for(auto& grp : ds->getRoot()->groups())
  {
    std::cout << axom::fmt::format("  * [{}] '{}' with {} groups and {} views\n",
                                   grp.getIndex(),
                                   grp.getName(),
                                   grp.getNumGroups(),
                                   grp.getNumViews());
  }

  std::cout << fill_line << std::endl;

  // iterate through the views of the 'state' group
  std::cout << "The 'state' group has the following views:\n";
  for(auto& view : ds->getRoot()->getGroup("state")->views())
  {
    std::cout << axom::fmt::format("  * [{}] '{}' -- {} view of type {} and {} elements\n",
                                   view.getIndex(),
                                   view.getName(),
                                   view.isAllocated() ? "Allocated" : "Unallocated",
                                   conduit::DataType::id_to_name(view.getTypeID()),
                                   view.getNumElements());
  }

  std::cout << fill_line << std::endl;
}

std::unique_ptr<sidre::DataStore> create_tiny_datastore()
{
  // _tiny_create_start
  auto ds = std::unique_ptr<sidre::DataStore> {new sidre::DataStore()};

  int nodecount = 12;
  int elementcount = 2;

  // Create views and buffers to hold node positions and field values
  sidre::Group* nodes = ds->getRoot()->createGroup("nodes");
  sidre::View* xs = nodes->createViewAndAllocate("xs", sidre::DOUBLE_ID, nodecount);
  sidre::View* ys = nodes->createViewAndAllocate("ys", sidre::DOUBLE_ID, nodecount);
  sidre::View* zs = nodes->createViewAndAllocate("zs", sidre::DOUBLE_ID, nodecount);

  sidre::Group* fields = ds->getRoot()->createGroup("fields");
  sidre::View* nodefield = fields->createViewAndAllocate("nodefield", sidre::INT_ID, nodecount);
  sidre::View* eltfield = fields->createViewAndAllocate("eltfield", sidre::DOUBLE_ID, elementcount);

  // Set node position for two adjacent hexahedrons
  double* xptr = xs->getArray();
  double* yptr = ys->getArray();
  double* zptr = zs->getArray();
  for(int pos = 0; pos < nodecount; ++pos)
  {
    xptr[pos] = ((pos + 1) / 2) % 2;
    yptr[pos] = (pos / 2) % 2;
    zptr[pos] = pos / 4;
  }

  // Assign a value to the node field
  int* nf = nodefield->getArray();
  for(int pos = 0; pos < nodecount; ++pos)
  {
    nf[pos] = static_cast<int>(xptr[pos] + yptr[pos] + zptr[pos]);
  }
  // and to the element field.
  double* ef = eltfield->getArray();
  // There are only two elements.
  ef[0] = 2.65;
  ef[1] = 1.96;

  return ds;
  // _tiny_create_end
}

void setup_blueprint_coords(sidre::DataStore* ds, sidre::Group* coords)
{
  // _blueprint_restructure_coords_start
  // Set up the coordinates as Mesh Blueprint requires
  coords->createViewString("type", "explicit");
  // We use prior knowledge of the layout of the original datastore
  sidre::View* origv = ds->getRoot()->getView("nodes/xs");
  sidre::Group* conduitval = coords->createGroup("values");
  conduitval->createView("x", sidre::DOUBLE_ID, origv->getNumElements(), origv->getBuffer());
  origv = ds->getRoot()->getView("nodes/ys");
  conduitval->createView("y", sidre::DOUBLE_ID, origv->getNumElements(), origv->getBuffer());
  origv = ds->getRoot()->getView("nodes/zs");
  conduitval->createView("z", sidre::DOUBLE_ID, origv->getNumElements(), origv->getBuffer());
  // _blueprint_restructure_coords_end
}

void setup_blueprint_external_coords(sidre::DataStore* ds, sidre::Group* coords)
{
  // _blueprint_external_coords_start
  // Set up the coordinates as Mesh Blueprint requires
  coords->createViewString("type", "explicit");
  // We use prior knowledge of the layout of the original datastore
  sidre::View* origv = ds->getRoot()->getView("nodes/xs");
  sidre::Group* conduitval = coords->createGroup("values");
  conduitval->createView("x",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         static_cast<double*>(origv->getArray()));
  origv = ds->getRoot()->getView("nodes/ys");
  conduitval->createView("y",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         static_cast<double*>(origv->getArray()));
  origv = ds->getRoot()->getView("nodes/zs");
  conduitval->createView("z",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         static_cast<double*>(origv->getArray()));
  // _blueprint_external_coords_end
}

void setup_blueprint_topos(sidre::DataStore* ds, sidre::Group* topos)
{
  // _blueprint_restructure_topo_start
  // Sew the nodes together into the two hexahedra, using prior knowledge.
  sidre::Group* connmesh = topos->createGroup("mesh");
  connmesh->createViewString("type", "unstructured");
  connmesh->createViewString("coordset", "coords");
  sidre::Group* elts = connmesh->createGroup("elements");
  elts->createViewString("shape", "hex");

  // We have two eight-node hex elements, so we need 2 * 8 = 16 ints.
  sidre::View* connectivity = elts->createViewAndAllocate("connectivity", sidre::INT_ID, 16);

  // The Mesh Blueprint connectivity array for a hexahedron lists four nodes on
  // one face arranged by right-hand rule to indicate a normal pointing into
  // the element, then the four nodes of the opposite face arranged to point
  // the normal the same way (out of the element).  This is the same as for
  // a VTK_HEXAHEDRON.  See
  // https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf.

  int* c = connectivity->getArray();

  // First hex.  In this example, the Blueprint node ordering matches the
  // dataset layout.  This is fortuitous but not required.
  c[0] = 0;
  c[1] = 1;
  c[2] = 2;
  c[3] = 3;
  c[4] = 4;
  c[5] = 5;
  c[6] = 6;
  c[7] = 7;

  // Second and last hex
  c[8] = 4;
  c[9] = 5;
  c[10] = 6;
  c[11] = 7;
  c[12] = 8;
  c[13] = 9;
  c[14] = 10;
  c[15] = 11;
  // _blueprint_restructure_topo_end

  // Deal with unused variable
  AXOM_UNUSED_VAR(ds);
}

void setup_blueprint_fields(sidre::DataStore* ds, sidre::Group* fields)
{
  // _blueprint_restructure_field_start
  // Set up the node-centered field
  // Get the original data
  sidre::View* origv = ds->getRoot()->getView("fields/nodefield");
  sidre::Group* nodefield = fields->createGroup("nodefield");
  nodefield->createViewString("association", "vertex");
  nodefield->createViewString("type", "scalar");
  nodefield->createViewString("topology", "mesh");
  nodefield->createView("values", sidre::INT_ID, origv->getNumElements(), origv->getBuffer());

  // Set up the element-centered field
  // Get the original data
  origv = ds->getRoot()->getView("fields/eltfield");
  sidre::Group* eltfield = fields->createGroup("eltfield");
  eltfield->createViewString("association", "element");
  eltfield->createViewString("type", "scalar");
  eltfield->createViewString("topology", "mesh");
  eltfield->createView("values", sidre::DOUBLE_ID, origv->getNumElements(), origv->getBuffer());
  // _blueprint_restructure_field_end
}

void setup_blueprint_external_fields(sidre::DataStore* ds, sidre::Group* fields)
{
  // _blueprint_external_field_start
  // Set up the node-centered field
  // Get the original data
  sidre::View* origv = ds->getRoot()->getView("fields/nodefield");
  sidre::Group* nodefield = fields->createGroup("nodefield");
  nodefield->createViewString("association", "vertex");
  nodefield->createViewString("type", "scalar");
  nodefield->createViewString("topology", "mesh");
  nodefield->createView("values",
                        sidre::INT_ID,
                        origv->getNumElements(),
                        static_cast<int*>(origv->getArray()));

  // Set up the element-centered field
  // Get the original data
  origv = ds->getRoot()->getView("fields/eltfield");
  sidre::Group* eltfield = fields->createGroup("eltfield");
  eltfield->createViewString("association", "element");
  eltfield->createViewString("type", "scalar");
  eltfield->createViewString("topology", "mesh");
  eltfield->createView("values",
                       sidre::DOUBLE_ID,
                       origv->getNumElements(),
                       static_cast<double*>(origv->getArray()));
  // _blueprint_external_field_end
}

void save_as_blueprint(sidre::DataStore* ds)
{
  // _blueprint_restructure_toplevel_start
  // Conduit needs a specific hierarchy.
  // We'll make a new DataStore with that hierarchy, pointing at the
  // application's data.
  std::string mesh_name = "tinymesh";

  // The Conduit specifies top-level groups:
  sidre::Group* mroot = ds->getRoot()->createGroup(mesh_name);
  sidre::Group* coords = mroot->createGroup("coordsets/coords");
  sidre::Group* topos = mroot->createGroup("topologies");
  // no material sets in this example
  sidre::Group* fields = mroot->createGroup("fields");
  // no adjacency sets in this (single-domain) example
  // _blueprint_restructure_toplevel_end

  setup_blueprint_coords(ds, coords);

  setup_blueprint_topos(ds, topos);

  setup_blueprint_fields(ds, fields);

  // _blueprint_restructure_save_start
  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::verify(bp_protocol, mesh_node[mesh_name], info))
  {
    // Generate the Conduit index
    conduit::Node& index = root_node["blueprint_index"];
    conduit::blueprint::mesh::generate_index(mesh_node[mesh_name], mesh_name, 1, index[mesh_name]);

    std::string root_output_path = mesh_name + ".root";
    std::string output_path = mesh_name + ".json";

    root_node["protocol/name"] = "json";
    root_node["protocol/version"] = "0.1";
    root_node["number_of_files"] = 1;
    root_node["number_of_trees"] = 1;
    root_node["file_pattern"] = output_path;
    root_node["tree_pattern"] = "/";

    // Now save both the index and the data set
    conduit::relay::io::save(root_node, root_output_path, "json");
    conduit::relay::io::save(mesh_node, output_path, "json");
  }
  else
  {
    std::cout << "does not conform to Mesh Blueprint: ";
    info.print();
    std::cout << std::endl;
  }
  // _blueprint_restructure_save_end
}

void generate_blueprint(sidre::DataStore* ds)
{
  // _blueprint_generate_toplevel_start
  // Conduit needs a specific hierarchy.
  // We'll make a new Group with that hierarchy, pointing at the
  // application's data.
  std::string domain_name = "domain";
  std::string domain_location = "domain_data/" + domain_name;
  std::string mesh_name = "mesh";
  std::string domain_mesh = domain_location + "/" + mesh_name;

  sidre::Group* mroot = ds->getRoot()->createGroup(domain_location);
  sidre::Group* coords = mroot->createGroup(mesh_name + "/coordsets/coords");
  sidre::Group* topos = mroot->createGroup(mesh_name + "/topologies");
  // no material sets in this example
  sidre::Group* fields = mroot->createGroup(mesh_name + "/fields");
  // no adjacency sets in this (single-domain) example
  // _blueprint_generate_toplevel_end

  setup_blueprint_coords(ds, coords);

  setup_blueprint_topos(ds, topos);

  setup_blueprint_fields(ds, fields);

  // _blueprint_generate_save_start
  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::verify(bp_protocol, mesh_node[domain_mesh], info))
  {
    std::string bp("rootfile_data/blueprint_index/automesh");

    ds->generateBlueprintIndex(domain_mesh, mesh_name, bp, 1);

    sidre::Group* rootfile_grp = ds->getRoot()->getGroup("rootfile_data");
    rootfile_grp->createViewString("protocol/name", "json");
    rootfile_grp->createViewString("protocol/version", "0.1");
    rootfile_grp->createViewScalar("number_of_files", 1);
    rootfile_grp->createViewScalar("number_of_trees", 1);
    rootfile_grp->createViewScalar("file_pattern", "bpgen.json");
    rootfile_grp->createViewScalar("tree_pattern", "/domain");
    rootfile_grp->save("bpgen.root", "json");

    ds->getRoot()->getGroup("domain_data")->save("bpgen.json", "json");
  }
  else
  {
    std::cout << "does not conform to Mesh Blueprint: ";
    info.print();
    std::cout << std::endl;
  }
  // _blueprint_generate_save_end
}

void generate_blueprint_to_path(sidre::DataStore* ds)
{
  // _blueprint_generate_path_start
  // Conduit needs a specific hierarchy.
  // For this example, we locate the domain data at a path
  // that is not the top level of its Sidre heirarchy.
  std::string domain_name = "domain";
  std::string domain_location = "domain_data/level/domains/" + domain_name;
  std::string mesh_name = "pathmesh";
  std::string domain_mesh = domain_location + "/" + mesh_name;

  sidre::Group* mroot = ds->getRoot()->createGroup(domain_location);
  sidre::Group* coords = mroot->createGroup(mesh_name + "/coordsets/coords");
  sidre::Group* topos = mroot->createGroup(mesh_name + "/topologies");
  // no material sets in this example
  sidre::Group* fields = mroot->createGroup(mesh_name + "/fields");
  // no adjacency sets in this (single-domain) example
  // _blueprint_generate_path_end

  setup_blueprint_coords(ds, coords);

  setup_blueprint_topos(ds, topos);

  setup_blueprint_fields(ds, fields);

  // _blueprint_generate_path_save_start
  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::verify(bp_protocol, mesh_node[domain_mesh], info))
  {
    std::string bp("rootfile_data/blueprint_index/automesh");

    ds->generateBlueprintIndex(domain_mesh, mesh_name, bp, 1);

    sidre::Group* rootfile_grp = ds->getRoot()->getGroup("rootfile_data");
    rootfile_grp->createViewString("protocol/name", "json");
    rootfile_grp->createViewString("protocol/version", "0.1");
    rootfile_grp->createViewScalar("number_of_files", 1);
    rootfile_grp->createViewScalar("number_of_trees", 1);
    rootfile_grp->createViewScalar("file_pattern", "pathbp.json");
    rootfile_grp->createViewScalar("tree_pattern", "level/domains/domain");
    rootfile_grp->save("pathbp.root", "json");

    ds->getRoot()->getGroup("domain_data")->save("pathbp.json", "json");
  }
  else
  {
    std::cout << "does not conform to Mesh Blueprint: ";
    info.print();
    std::cout << std::endl;
  }
  // _blueprint_generate_path_save_end
}

#ifdef AXOM_USE_MPI
void generate_spio_blueprint(sidre::DataStore* ds)
{
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // _blueprint_spio_toplevel_start
  std::string domain_name = "domain";
  std::string domain_location = "domain_data/" + domain_name;
  std::string mesh_name = "mesh";
  std::string domain_mesh = domain_location + "/" + mesh_name;

  sidre::Group* mroot = ds->getRoot()->createGroup(domain_location);
  sidre::Group* coords = mroot->createGroup(mesh_name + "/coordsets/coords");
  sidre::Group* topos = mroot->createGroup(mesh_name + "/topologies");
  // no material sets in this example
  sidre::Group* fields = mroot->createGroup(mesh_name + "/fields");
  // no adjacency sets in this (single-domain) example
  // _blueprint_spio_toplevel_end

  setup_blueprint_coords(ds, coords);

  setup_blueprint_topos(ds, topos);

  setup_blueprint_fields(ds, fields);

  // _blueprint_generate_spio_start
  sidre::IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol, mesh_node[domain_mesh], info, MPI_COMM_WORLD))
  {
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif
    std::string output_name = "bpspio";
    if(comm_size > 1)
    {
      output_name = output_name + "_par";
    }

    std::string bp_rootfile = output_name + ".root";

    writer.write(ds->getRoot()->getGroup(domain_location), 1, output_name, protocol);

    writer.writeBlueprintIndexToRootFile(ds, domain_mesh, bp_rootfile, mesh_name);
  }
  // _blueprint_generate_spio_end
}

void generate_spio_blueprint_to_path(sidre::DataStore* ds)
{
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // _blueprint_spio_path_start
  std::string domain_name = "domain";
  std::string domain_location = "domain_data/level/domains/" + domain_name;
  std::string mesh_name = "spiopathmesh";
  std::string domain_mesh = domain_location + "/" + mesh_name;

  sidre::Group* mroot = ds->getRoot()->createGroup(domain_location);
  sidre::Group* coords = mroot->createGroup(mesh_name + "/coordsets/coords");
  sidre::Group* topos = mroot->createGroup(mesh_name + "/topologies");
  // no material sets in this example
  sidre::Group* fields = mroot->createGroup(mesh_name + "/fields");
  // no adjacency sets in this (single-domain) example
  // _blueprint_spio_path_end

  setup_blueprint_external_coords(ds, coords);

  setup_blueprint_topos(ds, topos);

  setup_blueprint_external_fields(ds, fields);

  // _blueprint_generate_spio_path_start
  sidre::IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol, mesh_node[domain_mesh], info, MPI_COMM_WORLD))
  {
    std::string output_name = "pathbpspio";
    if(comm_size > 1)
    {
      output_name = output_name + "_par";
    }

    std::string bp_rootfile = output_name + ".root";
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif

    writer.write(ds->getRoot()->getGroup("domain_data"), 1, output_name, protocol);

    writer.writeBlueprintIndexToRootFile(ds,
                                         domain_mesh,
                                         bp_rootfile,
                                         "level/domains/domain/" + mesh_name);
  }
  // _blueprint_generate_spio_path_end
}
#endif

void serial_save_datastore_and_load_copy_lower(sidre::DataStore* ds)
{
#if defined(AXOM_USE_HDF5)
  std::string protocol = "sidre_hdf5";
  std::string filename = "example.hdf5";
#else
  std::string protocol = "sidre_json";
  std::string filename = "example.json";
#endif

  // _serial_io_save_start
  // Save the data store to a file,
  // saving all Views
  ds->getRoot()->save(filename, protocol);
  // Delete the data hierarchy under the root, then load it from the file
  ds->getRoot()->load(filename, protocol);
  sidre::Group* additional = ds->getRoot()->createGroup("additional");
  additional->createGroup("yetanother");
  // Load another copy of the data store into the "additional" group
  // without first clearing all its contents
  std::string groupname;
  additional->load(filename, protocol, true, groupname);
  // _serial_io_save_end
}

int main(int argc, char** argv)
{
  // Deal with unused variables
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  int region[3375];

#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
#endif

  auto ds = create_datastore(region);
  access_datastore(ds.get());
  iterate_datastore(ds.get());

  auto tds = create_tiny_datastore();
  save_as_blueprint(tds.get());

  auto bds = create_tiny_datastore();
  generate_blueprint(bds.get());

  auto pds = create_tiny_datastore();
  generate_blueprint_to_path(pds.get());

#ifdef AXOM_USE_MPI
  auto sds = create_tiny_datastore();
  auto spds = create_tiny_datastore();
  generate_spio_blueprint(sds.get());
  generate_spio_blueprint_to_path(spds.get());
  MPI_Finalize();
#endif

  return 0;
}
