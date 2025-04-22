// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Mint includes
#include "axom/mint/mesh/Mesh.hpp"       /* for mesh base class  */
#include "axom/mint/utils/su2_utils.hpp" /* for su2 i/o routines */
#include "axom/mint/utils/vtk_utils.hpp" /* for vtk i/o routines */

// Slic includes
#include "axom/slic/interface/slic.hpp"              /* for slic macros     */
#include "axom/slic/streams/GenericOutputStream.hpp" /* logging to terminal */

// C/C++ includes
#include <string>  /* for C++ string */
#include <cstring> /* for strcmp() */

/*!
 * \file
 *
 * \brief A simple example that illustrates how to read in an su2 mesh.
 */

//------------------------------------------------------------------------------
// GLOBALS
//------------------------------------------------------------------------------
static struct
{
  bool dumpVtk;
  std::string file;
} Parameters;

namespace slic = axom::slic;
namespace mint = axom::mint;

//------------------------------------------------------------------------------
// FUNCTION PROTOTYPES
//------------------------------------------------------------------------------
void parse_args(int argc, char** argv);

int main(int argc, char** argv)
{
  // STEP 0: initialize the slic logging environment
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::addStreamToAllMsgLevels(new slic::GenericOutputStream(&std::cout, "[<LEVEL>]: <MESSAGE>\n"));

  // STEP 1: parse command line arguemnts
  parse_args(argc, argv);

  // STEP 2: read SU2 mesh
  SLIC_INFO("reading file [" << Parameters.file << "]\n");

  mint::Mesh* mesh = nullptr;
  int rc = mint::read_su2(Parameters.file, mesh);
  SLIC_ERROR_IF((rc != 0), "Failed to read SU2 file!");
  SLIC_ASSERT(mesh != nullptr);

  // STEP 3: print some mesh information
  SLIC_INFO("MESH DIMENSION:  " << mesh->getDimension());
  SLIC_INFO("NUMBER OF NODES: " << mesh->getNumberOfNodes());
  SLIC_INFO("NUMBER OF CELLS: " << mesh->getNumberOfCells());
  if(mesh->hasMixedCellTypes())
  {
    SLIC_INFO("MIXED_CELL_TOPOLOGY");
  }
  else
  {
    SLIC_INFO("SINGLE_CELL_TOPOLOGY");
  }

  // STEP 4: dump a corresponding VTK file
  if(Parameters.dumpVtk)
  {
    std::string vtkFile = Parameters.file + ".vtk";
    SLIC_INFO("generating vtk file [" << vtkFile << "]");

    mint::write_vtk(mesh, vtkFile);
  }

  // STEP 4: finalize
  delete mesh;
  slic::flushStreams();
  slic::finalize();
  return 0;
}

//------------------------------------------------------------------------------
// FUNCTION PROTOTYPE IMPLEMENTATION
//------------------------------------------------------------------------------
void parse_args(int argc, char** argv)
{
  Parameters.file = "";
  Parameters.dumpVtk = false;

  for(int i = 1; i < argc; ++i)
  {
    if(strcmp(argv[i], "--file") == 0)
    {
      Parameters.file = std::string(argv[++i]);
    }
    else if(strcmp(argv[i], "--vtk") == 0)
    {
      Parameters.dumpVtk = true;
    }
  }

  SLIC_ERROR_IF(Parameters.file.length() == 0, "must specify a valid SU2 file.");
}
