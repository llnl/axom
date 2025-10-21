// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file quest_proe_bbox.cpp
 *  \brief This example code demonstrates how to read a subset of a ProE (Creo)
 *  text-format tetrahedron mesh using a bounding box.  It reads a file
 *  specified with the -f argument and writes the mesh contained therein
 *  into a file specified with the -o argument.
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

// _read_proe_include2_start
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/utils/vtk_utils.hpp"  // for write_vtk
// _read_proe_include2_end

// _read_proe_include1_start
#include "axom/quest/io/ProEReader.hpp"
// _read_proe_include1_end

// _read_proe_typealiases_start
namespace mint = axom::mint;
namespace primal = axom::primal;
namespace slic = axom::slic;

using IndexType = axom::IndexType;
using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
// _read_proe_typealiases_end

//------------------------------------------------------------------------------
void initialize_logger()
{
  // initialize logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  // setup the logstreams
  std::string fmt = "";
  slic::LogStream* logStream = nullptr;

  fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);

  // register stream objects with the logger
  slic::addStreamToAllMsgLevels(logStream);
}

//------------------------------------------------------------------------------
void finalize_logger()
{
  slic::flushStreams();
  slic::finalize();
}

struct Arguments
{
  std::string file_name;
  std::string outfile_name;
  std::vector<double> bbox_min;
  std::vector<double> bbox_max;

  Arguments()
  {
    bbox_min.resize(3);
    bbox_max.resize(3);
  }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-i,--input", this->file_name, "specifies the input mesh file")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_option("-o,--outfile", this->outfile_name, "specifies the output mesh file")->required();

    app.add_option("--min", this->bbox_min, "specifies the minimum of the bounding box")
      ->expected(3)
      ->required();

    app.add_option("--max", this->bbox_max, "specifies the maximum of the bounding box")
      ->expected(3)
      ->required();

    app.get_formatter()->column_width(40);

    // could throw an exception
    app.parse(argc, argv);

    slic::flushStreams();
  }
};

int main(int argc, char** argv)
{
  initialize_logger();
  Arguments args;
  axom::CLI::App app {"Example showing how to read a Pro/E mesh"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    retval = app.exit(e);
    finalize_logger();
    return retval;
  }

  // The constructor of args resizes bbox_min and bbox_max to three elements.
  double* bbox_min = args.bbox_min.data();
  double* bbox_max = args.bbox_max.data();

  SLIC_INFO("Reading file: '" << args.file_name << "'...\n");
  // _read_proe_file_start
  // Read file
  axom::quest::ProEReader reader;
  reader.setFileName(args.file_name);

  // Set up a bounding box to keep only certain tets.
  // bbox_min and bbox_max are pointers to double.
  axom::quest::ProEReader::BBox3D bbox;
  bbox.addPoint(axom::quest::ProEReader::Point3D {bbox_min, 3});
  bbox.addPoint(axom::quest::ProEReader::Point3D {bbox_max, 3});
  // Keep only tets with all four nodes inside the bounding box.
  reader.setTetPredFromBoundingBox(bbox, false);
  // Pass true as the second argument of setTetPredFromBoundingBox() to
  // keep tets with at least one node inside the bounding box.
  // To keep all tets, do not set a TetPred.

  // Read in the file.
  reader.read();

  // Get surface mesh
  UMesh mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);
  // _read_proe_file_end

  SLIC_INFO("Mesh has " << mesh.getNumberOfNodes() << " vertices and " << mesh.getNumberOfCells()
                        << " triangles.");

  axom::mint::write_vtk(&mesh, args.outfile_name);

  finalize_logger();
}
