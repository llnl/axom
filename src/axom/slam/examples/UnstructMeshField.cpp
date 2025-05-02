// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file UnstructMeshField.cpp
 *
 * \brief Example that uses Slam for generating and processing a simple 3D mesh
 *
 * \details Loads a hex mesh from a VTK file, generates the Node to Zone
 *  relation and does simple mesh processing.
 *
 * \author T. Brunner (original)
 * \author K. Weiss (modified to use axom's Slam component)
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include "axom/fmt.hpp"

#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

namespace slam = axom::slam;

namespace slamUnstructuredHex
{
using DataType = double;
using Point3 = slam::util::Point3<DataType>;

/**
 * \brief Simple hex mesh for this example.
 *  Contains the necessary types, for a set of nodes and zones,
 *  the incidence relations between them and some fields defined on them.
 */
struct HexMesh
{
public:
  enum
  {
    COORDS_PER_NODE = 3,
    NODES_PER_ZONE = 8
  };

  /// types for sets
  using PositionType = slam::DefaultPositionType;
  using ElementType = slam::DefaultElementType;

  using NodeSet = slam::PositionSet<PositionType, ElementType>;
  using ZoneSet = slam::PositionSet<PositionType, ElementType>;

  /// types for relations
  using STLIndirection = slam::policies::STLVectorIndirection<PositionType, PositionType>;
  using VariableCardinality = slam::policies::VariableCardinality<PositionType, STLIndirection>;
  using NodeToZoneRelation =
    slam::StaticRelation<PositionType, ElementType, VariableCardinality, STLIndirection, NodeSet, ZoneSet>;
  using NodeZoneIterator = NodeToZoneRelation::RelationConstIterator;

  using ZNStride = slam::policies::CompileTimeStride<PositionType, NODES_PER_ZONE>;
  using ConstantCardinality = slam::policies::ConstantCardinality<PositionType, ZNStride>;
  using ZoneToNodeRelation =
    slam::StaticRelation<PositionType, ElementType, ConstantCardinality, STLIndirection, ZoneSet, NodeSet>;
  using ZoneNodeIterator = ZoneToNodeRelation::RelationConstIterator;

  /// types for maps
  using BaseSet = axom::slam::Set<PositionType, ElementType>;
  using NodalPositions = slam::Map<Point3>;
  using ZonalPositions = slam::Map<Point3>;
  using NodeField = slam::Map<DataType>;
  using ZoneField = slam::Map<DataType>;

public:
  /** \brief Simple accessor for the number of nodes in the mesh  */
  PositionType numNodes() const { return nodes.size(); }

  /** \brief Simple accessor for the number of zones in the mesh */
  PositionType numZones() const { return zones.size(); }

public:
  /// Sets in the mesh
  NodeSet nodes;
  ZoneSet zones;

  /// Relations in the mesh
  ZoneToNodeRelation zoneToNodeRelation;  // storage for relation_(3,0): zones -> nodes
  NodeToZoneRelation nodeToZoneRelation;  // storage for relation_(0,3): nodes -> zones

  /// Maps (fields) defined on the mesh -- nodal and zonal positions and scalar fields
  NodalPositions nodePosition;
  ZonalPositions zonePosition;
  ZoneField zoneField;
  NodeField nodeFieldExact;
  NodeField nodeFieldAvg;
};

/// The repository is a proxy for a data allocator/manager
struct Repository
{
  // Define the explicit instances of our local (key/value) datastore for int and double
  using SetType = axom::slam::Set<>;
  using IntsRegistry = slam::FieldRegistry<SetType, SetType::ElementType>;
  using RealsRegistry = slam::FieldRegistry<SetType, double>;
  using IntField = slam::Map<int>;
  using RealField = slam::Map<double>;

  static IntsRegistry intsRegistry;
  static RealsRegistry realsRegistry;
};

Repository::IntsRegistry Repository::intsRegistry;
Repository::RealsRegistry Repository::realsRegistry;

/** A simple class to read a VTK hex mesh */
class SimpleVTKHexMeshReader
{
public:
  // uses RAII to open/close the file
  SimpleVTKHexMeshReader(const std::string& fileName) : vtkMesh(fileName.c_str())
  {
    // Handle errors opening file, if any
    if(!vtkMesh)
    {
      bool fExists = axom::utilities::filesystem::pathExists(fileName);
      SLIC_ERROR(axom::fmt::format("Attempted to open file '{}'. It {}.",
                                   fileName,
                                   fExists ? "exists" : "does not exist"));
    }
  }
  ~SimpleVTKHexMeshReader()
  {
    vtkMesh.close();  // Close the file.
  }

  void parseMeshFile()
  {
    using RealBuf = Repository::RealsRegistry::BufferType;
    using IndexBuf = Repository::IntsRegistry::BufferType;
    using PositionType = HexMesh::PositionType;

    // Read some initial header stuff.  Note: this is not a robust vtkreader
    std::string junk;

    PositionType numNodes;

    /// Read in POINT data, (which we'll call the nodes)
    while(junk != "POINTS")
    {
      vtkMesh >> junk;
    }
    vtkMesh >> numNodes >> junk;

    Repository::intsRegistry.addScalar("num_nodes", numNodes);
    SLIC_INFO("-- Number of nodes: " << numNodes);

    const PositionType numCoords = HexMesh::COORDS_PER_NODE * numNodes;
    RealBuf& pointData = Repository::realsRegistry.addBuffer("node_positions", numCoords);
    for(PositionType idx = 0; idx < numCoords; ++idx)
    {
      vtkMesh >> pointData[idx];
    }

    /// Read in the CELL data (which we'll call zones).
    /// We're going to assume hexahedra (VTK type 12)
    PositionType numZones, listSize, nodeCount;
    vtkMesh >> junk >> numZones >> listSize;
    Repository::intsRegistry.addScalar("num_zones", numZones);
    const PositionType numNodeZoneIndices = HexMesh::NODES_PER_ZONE * numZones;

    SLIC_INFO("-- Number of zones: " << numZones);

    // Note: The VTK format has an extra value per zone for the number of indices
    // This is constant since we're assuming a Hex mesh.  General meshes can be different.
    SLIC_ASSERT_MSG((listSize - numZones) == numNodeZoneIndices,
                    axom::fmt::format("Error while reading mesh!\n "
                                      "numZones = {0}; numZones*{1} = {2}; indices in file = {3}",
                                      numZones,
                                      static_cast<int>(HexMesh::NODES_PER_ZONE),
                                      numNodeZoneIndices,
                                      listSize - numZones));

    IndexBuf& zn_indices =
      Repository::intsRegistry.addBuffer("zone_node_indices", numNodeZoneIndices);
    PositionType idx = 0;
    for(PositionType zoneIdx = 0; zoneIdx < numZones; ++zoneIdx)
    {
      vtkMesh >> nodeCount;
      SLIC_ASSERT(nodeCount == HexMesh::NODES_PER_ZONE);

      for(PositionType n = 0; n < (HexMesh::NODES_PER_ZONE); ++n)
      {
        vtkMesh >> zn_indices[idx++];
      }
    }
  }

private:
  std::ifstream vtkMesh;
};

void readHexMesh(std::string fileName, HexMesh* mesh)
{
  {
    SimpleVTKHexMeshReader vtkMeshReader(fileName);
    vtkMeshReader.parseMeshFile();
  }
  using RealBuf = Repository::RealsRegistry::BufferType;
  using IndexBuf = Repository::IntsRegistry::BufferType;
  using PositionType = HexMesh::PositionType;

  /// Check that the mesh has been loaded properly
  if(!(Repository::intsRegistry.hasScalar("num_nodes") &&
       Repository::intsRegistry.hasScalar("num_zones") &&
       Repository::realsRegistry.hasBuffer("node_positions") &&
       Repository::intsRegistry.hasBuffer("zone_node_indices")))
  {
    SLIC_ERROR("Hex mesh not loaded properly from file " << fileName);
  }

  /// Create the sets of nodes and zones in the mesh
  PositionType numNodes = Repository::intsRegistry.getScalar("num_nodes");
  mesh->nodes = HexMesh::NodeSet(numNodes);

  PositionType numZones = Repository::intsRegistry.getScalar("num_zones");
  mesh->zones = HexMesh::ZoneSet(numZones);

  /// Create the nodal position field
  mesh->nodePosition = HexMesh::NodalPositions(&mesh->nodes);
  RealBuf::iterator ptIt = Repository::realsRegistry.getBuffer("node_positions").begin();
  for(PositionType idx = 0; idx < mesh->numNodes(); ++idx)
  {
    mesh->nodePosition[idx] = Point3(*ptIt++, *ptIt++, *ptIt++);
  }

  /// Create the topological incidence relation from zones to nodes
  IndexBuf& zn_indices = Repository::intsRegistry.getBuffer("zone_node_indices");
  mesh->zoneToNodeRelation = HexMesh::ZoneToNodeRelation(&mesh->zones, &mesh->nodes);
  mesh->zoneToNodeRelation.bindIndices(zn_indices.size(), &zn_indices);

  // Check that the relation is valid
  SLIC_ASSERT_MSG(mesh->zoneToNodeRelation.isValid(),
                  "Error creating (static) relation from zones to nodes!");
  SLIC_INFO("-- numNodesOfZones: " << zn_indices.size());
}

void generateNodeZoneRelation(HexMesh* mesh)
{
  // Create NodeToZone relation by inverting the ZoneToZone relation
  // TODO: This function to invert a relation should be moved into Slam

  using IndexBuf = Repository::IntsRegistry::BufferType;
  using RelationSubset = HexMesh::ZoneToNodeRelation::RelationSubset;
  using PositionType = HexMesh::PositionType;

  /// Step 1: Compute the cardinalities of each node by looping through zone to node relation
  IndexBuf& nzBegins = Repository::intsRegistry.addBuffer("node_zone_begins", mesh->nodes.size() + 1);
  for(PositionType zIdx = 0; zIdx < mesh->numZones(); ++zIdx)
  {
    RelationSubset nSet = mesh->zoneToNodeRelation[zIdx];
    for(PositionType idx = 0; idx < nSet.size(); ++idx)
    {
      ++nzBegins[nSet[idx]];
    }
  }

  /// Step 2: Compute begin offsets for each node based on cardinalities
  // Strategy: perform (inplace) exclusive prefix sum of cardinalities in nzBegins
  PositionType prevVal = nzBegins[0];
  nzBegins[0] = 0;
  for(int i = 1; i <= mesh->numNodes(); ++i)
  {
    PositionType nextVal = nzBegins[i];
    nzBegins[i] = nzBegins[i - 1] + prevVal;
    prevVal = nextVal;
  }

  /// Step 3: Invert the zone_node relation, use nzBegins[node_index] as offset for next zone
  IndexBuf& zIndices =
    Repository::intsRegistry.addBuffer("node_zone_indices", nzBegins[mesh->numNodes()]);
  for(PositionType zIdx = 0; zIdx < mesh->numZones(); ++zIdx)
  {
    RelationSubset nSet = mesh->zoneToNodeRelation[zIdx];
    for(PositionType idx = 0; idx < nSet.size(); ++idx)
    {
      const PositionType nIdx = nSet[idx];
      const PositionType offset = nzBegins[nIdx]++;
      zIndices[offset] = zIdx;
    }
  }

  /// Step 4: Fix begin offsets by shifting back by one index
  for(int i = mesh->numNodes(); i > 0; --i)
  {
    nzBegins[i] = nzBegins[i - 1];
  }
  nzBegins[0] = 0;

  /// We can finally create the node to zone relation
  mesh->nodeToZoneRelation = HexMesh::NodeToZoneRelation(&mesh->nodes, &mesh->zones);
  mesh->nodeToZoneRelation.bindBeginOffsets(mesh->nodes.size(), &nzBegins);
  mesh->nodeToZoneRelation.bindIndices(zIndices.size(), &zIndices);

  SLIC_ASSERT_MSG(mesh->nodeToZoneRelation.isValid(true),
                  "Error creating (static) relation from nodes to zones!\n");

  SLIC_INFO("-- numZonesOfNode: " << zIndices.size());
}

void computeZoneBarycenters(HexMesh* mesh)
{
  using NodeSet = HexMesh::ZoneToNodeRelation::RelationSubset;
  using PositionType = HexMesh::PositionType;

  // Compute the zone positions as the the averages
  // of the positions of the nodes around each zone
  mesh->zonePosition = HexMesh::ZonalPositions(&mesh->zones);

  // Outer loop over each zone in the mesh
  for(PositionType zIdx = 0; zIdx < mesh->numZones(); ++zIdx)
  {
    Point3 zonePos;

    // Inner loop over each node of the zone
    const NodeSet& nodeSet = mesh->zoneToNodeRelation[zIdx];
    for(PositionType idx = 0; idx < nodeSet.size(); ++idx)
    {
      zonePos += mesh->nodePosition[nodeSet[idx]];
    }
    zonePos /= nodeSet.size();

    mesh->zonePosition[zIdx] = zonePos;
  }
}

void createZoneRadiusField(HexMesh* mesh)
{
  using PositionType = HexMesh::PositionType;

  // Compute a zone field based on the L2 norm of their position vectors
  mesh->zoneField = HexMesh::ZoneField(&mesh->zones);

  for(PositionType zIdx = 0; zIdx < mesh->numZones(); ++zIdx)
  {
    mesh->zoneField[zIdx] = mesh->zonePosition[zIdx].radius();
  }
}

DataType computeNodalErrors(HexMesh* mesh)
{
  // Compute the node average version, and the approximation error
  using ZoneSet = HexMesh::NodeToZoneRelation::RelationSubset;
  using PositionType = HexMesh::PositionType;

  mesh->nodeFieldAvg = HexMesh::NodeField(&mesh->nodes);
  mesh->nodeFieldExact = HexMesh::NodeField(&mesh->nodes);
  double errSqSum = 0.0;

  // Outer loop over each node
  for(PositionType nIdx = 0; nIdx < mesh->numNodes(); ++nIdx)
  {
    // What's the radius?
    const DataType nodalValExact = mesh->nodePosition[nIdx].radius();

    // Inner loop over each zone of the node to find average value
    const ZoneSet& zoneSet = mesh->nodeToZoneRelation[nIdx];
    DataType nodalAvg = DataType();
    for(PositionType idx = 0; idx < zoneSet.size(); ++idx)
    {
      nodalAvg += mesh->zoneField[zoneSet[idx]];
    }
    nodalAvg /= zoneSet.size();

    // Update fields and compute error
    mesh->nodeFieldAvg[nIdx] = nodalAvg;
    mesh->nodeFieldExact[nIdx] = nodalValExact;

    const double err = nodalAvg - nodalValExact;
    errSqSum += err * err;
  }

  DataType err = std::sqrt(errSqSum / mesh->numNodes());
  SLIC_INFO("-> The L2-ish error in the node average radius was " << err);

  return err;
}

}  // end namespace slamUnstructuredHex

int main(int argc, char** argv)
{
  using namespace slamUnstructuredHex;

  axom::slic::SimpleLogger logger;

#ifndef USE_ONE
  int const NUM_RESOLUTIONS = 4;
#else
  int const NUM_RESOLUTIONS = 1;
#endif

  int fileResolutions[] = {1, 2, 4, 8};
  DataType expectedResults[] = {0.10736689892, 0.037977237476, 0.013251067479, 0.0046357167735};

  std::string dataDir;
  if(argc > 1)
  {
    dataDir = std::string(argv[1]);
  }
  else
  {
    SLIC_WARNING("Usage: 'slam_unstructMesh_ex <data_directory>'");

#ifdef AXOM_DATA_DIR
    namespace fs = axom::utilities::filesystem;
    dataDir = fs::joinPath(AXOM_DATA_DIR, "slam");

    SLIC_INFO("Using default data directory '" << dataDir << "'");
#else
    axom::utilities::processAbort();
#endif
  }

  int numFailedTests = 0;
  for(int res = 0; res < NUM_RESOLUTIONS; ++res)
  {
    std::string meshName = axom::fmt::format("{}/ball_{}.vtk", dataDir, fileResolutions[res]);

    SLIC_INFO("Loading mesh file '" << meshName << "' and generating zone-> node relation");

    HexMesh hexMesh;
    readHexMesh(meshName, &hexMesh);
    SLIC_ASSERT(hexMesh.zoneToNodeRelation.isValid());
    //--------------------------------------------------------------

    // Now build the node to zone relation
    SLIC_INFO("Generating node->zone relation");
    generateNodeZoneRelation(&hexMesh);
    SLIC_ASSERT(hexMesh.nodeToZoneRelation.isValid());

    //--------------------------------------------------------------
    // Now that we have the mesh in memory, we can start to do things with it.
    SLIC_INFO("Computing zone barycenters using zone->node relation");
    computeZoneBarycenters(&hexMesh);

    SLIC_INFO("Generating a zone-centered radius field");
    createZoneRadiusField(&hexMesh);

    DataType errVal = computeNodalErrors(&hexMesh);

    // Some error checking based on precomputed values
    if(!axom::utilities::isNearlyEqual(errVal, expectedResults[res]))
    {
      SLIC_WARNING("Error differed from expected value -- "
                   << axom::fmt::format("Expected {}, but got {} (difference: {}",
                                        expectedResults[res],
                                        errVal,
                                        errVal - expectedResults[res]));

      ++numFailedTests;
    }

    SLIC_INFO("-- done.\n");
  }

  //--------------------------------------------------------------
  SLIC_INFO(axom::fmt::format("-- {} tests out of {} passed",
                              NUM_RESOLUTIONS - numFailedTests,
                              NUM_RESOLUTIONS));

  return (numFailedTests == 0) ? 0 : 1;
}
