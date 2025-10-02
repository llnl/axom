// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file sidre_generateindex.cpp
 *  \brief This example code is a basic demonstration of the generation
 *  of a Sidre DataStore containing a mesh conforming to the Mesh Blueprint,
 *  including the automatic generation of a Blueprint index and saving it
 *  to a file.
 */

/*
 * Summary of the functions called by main this example:
 *
 * create_tiny_datastore      -- Builds a DataStore that holds a small 3D mesh
 *                               conforming to the conduit mesh blueprint
 *
 * generate_spio_blueprint    -- Uses SPIO/IOManager calls to write data to
 *                               file and generate blueprint index
 */

// Axom headers
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

// Conduit headers
#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"

#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
#endif

// fmt format header
#include "axom/fmt.hpp"

// C++ headers
#include <cstring>

// "using" directives to simplify code
namespace sidre = axom::sidre;

static int cart_nx = 5;  // nodal domain x size for cartesian mesh example
static int cart_ny = 5;  // nodal domain y size for cartesian mesh example

std::unique_ptr<sidre::DataStore> create_tiny_datastore()
{
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
}

void setup_blueprint_coords(sidre::DataStore* ds, sidre::Group* coords)
{
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
}

void setup_cartesian_coords(sidre::Group* coords, int domain_id)
{
  // Set up uniform cartesian coordinates
  coords->createViewString("type", "uniform");
  coords->createViewScalar("dims/i", cart_nx);
  coords->createViewScalar("dims/j", cart_ny);

  double x_lo = static_cast<double>(domain_id);
  double x_hi = x_lo + 1.0;
  double dx = (x_hi - x_lo) / static_cast<double>(cart_nx - 1);
  double y_lo = 0.0;
  double y_hi = 1.5;
  double dy = (y_hi - y_lo) / static_cast<double>(cart_ny - 1);
  coords->createViewScalar("origin/x", x_lo);
  coords->createViewScalar("origin/y", y_lo);
  coords->createViewScalar("spacing/dx", dx);
  coords->createViewScalar("spacing/dy", dy);
}

void setup_blueprint_topos(sidre::DataStore* ds, sidre::Group* topos)
{
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

  // Deal with unused variable
  AXOM_UNUSED_VAR(ds);
}

void setup_cartesian_topos(sidre::Group* topos)
{
  sidre::Group* structmesh = topos->createGroup("cartesian");
  structmesh->createViewString("type", "structured");
  structmesh->createViewString("coordset", "coords");
  structmesh->createViewScalar("elements/dims/i", cart_nx - 1);
  structmesh->createViewScalar("elements/dims/j", cart_ny - 1);
}

void setup_blueprint_fields(sidre::DataStore* ds, sidre::Group* fields)
{
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
}

void setup_cartesian_fields(sidre::Group* fields)
{
  // Match known values from setup_cartesian_coords
  int ex = cart_nx - 1;
  int ey = cart_ny - 1;
  int eltcount = ex * ey;
  int nodecount = cart_nx * cart_ny;

  // Create a node-centered field and an element-centered field
  sidre::Group* nodefield = fields->createGroup("nodefield");
  nodefield->createViewString("association", "vertex");
  nodefield->createViewString("type", "scalar");
  nodefield->createViewString("topology", "cartesian");
  sidre::View* nodes = nodefield->createViewAndAllocate("values", sidre::DOUBLE_ID, nodecount);

  double* nptr = nodes->getArray();

  for(int j = 0; j < cart_ny; ++j)
  {
    for(int i = 0; i < cart_nx; ++i)
    {
      nptr[j * cart_nx + i] = static_cast<double>(i * j);
    }
  }

  sidre::Group* eltfield = fields->createGroup("eltfield");
  eltfield->createViewString("association", "element");
  eltfield->createViewString("type", "scalar");
  eltfield->createViewString("topology", "cartesian");
  sidre::View* elts = eltfield->createViewAndAllocate("values", sidre::DOUBLE_ID, eltcount);

  double* eptr = elts->getArray();

  for(int j = 0; j < ey; ++j)
  {
    for(int i = 0; i < ex; ++i)
    {
      eptr[j * ex + i] = static_cast<double>(i * j);
    }
  }
}

#ifdef AXOM_USE_MPI
void generate_spio_blueprint(sidre::DataStore* ds, bool dense)
{
  int my_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  std::string domain_name = "domain";
  std::string holder_name = "domain_data";
  std::string domain_location = holder_name + "/" + domain_name;
  std::string mesh_name = "mesh";

  sidre::Group* holder = ds->getRoot()->createGroup(holder_name);

  // dense places domains on all ranks when true. When not true, it
  // leaves off rank zero to show that the index can still be generated
  // with some ranks empty.
  if(dense || my_rank > 0 || comm_size == 1)
  {
    sidre::Group* mroot = holder->createGroup(domain_name);
    sidre::Group* coords = mroot->createGroup("coordsets/coords");
    sidre::Group* topos = mroot->createGroup("topologies");
    // no material sets in this example
    sidre::Group* fields = mroot->createGroup("fields");
    // no adjacency sets in this (single-domain) example

    setup_blueprint_coords(ds, coords);

    setup_blueprint_topos(ds, topos);

    setup_blueprint_fields(ds, fields);
  }

  sidre::IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol, mesh_node[domain_location], info, MPI_COMM_WORLD))
  {
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif

    std::string output_name;
    if(dense)
    {
      output_name = "bpdense";
    }
    else
    {
      output_name = "bpsparse";
    }

    if(comm_size > 1)
    {
      output_name = output_name + "_par";
    }

    std::string bp_rootfile = output_name + ".root";

    writer.write(ds->getRoot()->getGroup("domain_data"), 1, output_name, protocol);

    writer.writeBlueprintIndexToRootFile(ds, domain_location, bp_rootfile, mesh_name);
  }
}

void generate_multidomain_blueprint(sidre::DataStore* ds,
                                    const std::string& filename,
                                    int num_files,
                                    bool use_list)
{
  int my_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // 3 domains on even ranks, 2 domains on odd ranks
  int64_t domain_begin = (5 * (my_rank / 2)) + (3 * (my_rank % 2));
  int64_t domain_end = domain_begin + (3 - (my_rank % 2));

  std::string holder_name = "domain_data";
  std::string mesh_name = "mesh";
  sidre::Group* holder = ds->getRoot()->createGroup(holder_name, use_list);

  std::string domain_pattern = "domain_{domain:06d}";
  ds->getRoot()->createViewString("domain_pattern", domain_pattern);

  for(int64_t i = domain_begin; i < domain_end; ++i)
  {
    sidre::Group* mroot = nullptr;

    // The list format does not use string names to identify the domains
    // held by group 'holder`, while the map format in the else block
    // must give each domain a unique name.
    if(use_list)
    {
      mroot = holder->createUnnamedGroup();
    }
    else
    {
      std::string domain_name = axom::fmt::format("domain_{:06d}", i);
      mroot = holder->createGroup(domain_name);
    }
    sidre::Group* coords = mroot->createGroup("coordsets/coords");
    sidre::Group* topos = mroot->createGroup("topologies");
    // no material sets in this example
    sidre::Group* fields = mroot->createGroup("fields");
    // no adjacency sets in this example
    mroot->createViewScalar("state/domain_id", i);

    setup_blueprint_coords(ds, coords);

    setup_blueprint_topos(ds, topos);

    setup_blueprint_fields(ds, fields);
  }

  sidre::IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol, mesh_node[holder_name], info, MPI_COMM_WORLD))
  {
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif

    std::string output_name = filename;

    if(comm_size > 1)
    {
      output_name = output_name + "_par";
    }

    std::string bp_rootfile = output_name + ".root";

    writer.write(ds->getRoot()->getGroup("domain_data"),
                 std::min(num_files, comm_size),
                 output_name,
                 protocol);

    writer.writeBlueprintIndexToRootFile(ds, holder_name, bp_rootfile, mesh_name);
  }
}

void generate_cartesian_blueprint(sidre::DataStore* ds, const std::string& filename, int num_files)
{
  int my_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // 1 domain on even ranks, 2 domains on odd ranks
  int64_t domain_begin = (3 * (my_rank / 2)) + (1 * (my_rank % 2));
  int64_t domain_end = domain_begin;
  if(my_rank % 2)
  {
    domain_end += 2;
  }
  else
  {
    domain_end++;
  }

  std::string holder_name = "domain_data";
  std::string mesh_name = "";
  sidre::Group* holder = ds->getRoot()->createGroup(holder_name);

  std::string domain_pattern = "domain_{domain:06d}";
  ds->getRoot()->createViewString("domain_pattern", domain_pattern);

  for(int64_t i = domain_begin; i < domain_end; ++i)
  {
    std::string domain_name = axom::fmt::format("domain_{:06d}", i);

    sidre::Group* mroot = holder->createGroup(domain_name);
    sidre::Group* coords = mroot->createGroup("coordsets/coords");
    sidre::Group* topos = mroot->createGroup("topologies");
    // no material sets in this example
    sidre::Group* fields = mroot->createGroup("fields");
    // no adjacency sets in this example
    mroot->createViewScalar("state/domain_id", i);

    setup_cartesian_coords(coords, i);

    setup_cartesian_topos(topos);

    setup_cartesian_fields(fields);
  }

  sidre::IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol, mesh_node[holder_name], info, MPI_COMM_WORLD))
  {
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif

    std::string output_name = filename;

    if(comm_size > 1)
    {
      output_name = output_name + "_par";
    }

    std::string bp_rootfile = output_name + ".root";

    writer.write(ds->getRoot()->getGroup("domain_data"),
                 std::min(num_files, comm_size),
                 output_name,
                 protocol);

    writer.writeBlueprintIndexToRootFile(ds, holder_name, bp_rootfile, mesh_name);
  }
}
#endif

int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);
#ifdef AXOM_USE_MPI
  int num_ranks = 1;

  axom::slic::SimpleLogger logger;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  auto sds = create_tiny_datastore();
  generate_spio_blueprint(sds.get(), true);  //dense

  auto spds = create_tiny_datastore();
  generate_spio_blueprint(spds.get(), false);  //sparse
  spds->getRoot()->destroyGroups();
  spds->getRoot()->destroyViews();

  spds = create_tiny_datastore();
  generate_multidomain_blueprint(spds.get(), "multi1", 1, false);
  spds->getRoot()->destroyGroups();
  spds->getRoot()->destroyViews();

  spds = create_tiny_datastore();
  generate_multidomain_blueprint(spds.get(), "multi1_list", 1, true);
  if(num_ranks > 1)
  {
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();

    spds = create_tiny_datastore();
    generate_multidomain_blueprint(spds.get(), "multi2", 2, false);
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();

    spds = create_tiny_datastore();
    generate_multidomain_blueprint(spds.get(), "multi2_list", 2, true);
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();

    spds = create_tiny_datastore();
    generate_multidomain_blueprint(spds.get(), "multi_all", num_ranks, false);
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();

    spds = create_tiny_datastore();
    generate_multidomain_blueprint(spds.get(), "multi_all_list", num_ranks, true);
  }
  spds->getRoot()->destroyGroups();
  spds->getRoot()->destroyViews();

  generate_cartesian_blueprint(spds.get(), "cart1", 1);
  if(num_ranks > 1)
  {
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();

    spds = create_tiny_datastore();
    generate_cartesian_blueprint(spds.get(), "cart2", 2);
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();

    spds = create_tiny_datastore();
    generate_cartesian_blueprint(spds.get(), "cart_all", num_ranks);
  }

  MPI_Finalize();
#endif

  return 0;
}
