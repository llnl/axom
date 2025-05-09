// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// arbitrary mesh of polygons in XY geom
// Thu Mar 26 09:38:50 PDT 2015

#ifndef __PolygonMeshXY_hh__
#define __PolygonMeshXY_hh__

#include "VectorXY.hpp"
#include "TinyHydroTypes.hpp"

#include <string>

namespace tinyHydro {

  class PolygonMeshXY
  {

  public:
    // c'tor for doing quads
    PolygonMeshXY(int kmax, int lmax,
        double xmin = 0.0, double xmax = 1.0,
        double ymin = 0.0, double ymax = 1.0);

    // c'tor for reading in a mesh file
    enum MeshFileType {MT0, SILO, PDB, TRIANGLE};
    PolygonMeshXY(std::string /*filename*/, MeshFileType /*meshfiletype = MT0*/)
    {
      SLIC_ERROR("not implemented yet\n");
    }


    // functions that access the mesh topology
    int       numNodes()                           const { return nodes.size(); }
    int       numZones()                           const { return zones.size(); }
    int       zoneNumNodes(int i)                  const { return zoneToNodes.size(i); }
    int       zoneNode(int iz, int in)             const { return zoneToNodes[iz][in]; }

    // functions that access and modify the mesh geometry
    void      moveNodesToPosition(const NodalVectorField& newPos);
    void      computeNewGeometry(void);

    // accessors, mostly used from Python
    double    zoneVol(int i)                    const { return zoneVolume[i]; }
    VectorXY  getPos(int i)                   const { return nodePos[i]; }
    VectorXY  getZonePos(int i)               const { return zonePos[i]; }
    VectorXY  zoneNodePos(int iz, int in)     const { return nodePos[ zoneNode(iz,in) ]; }

    void      setPos(int i, const VectorXY & v)         { nodePos[i] = v; }

    VectorXY  meshAverageKLZMemOrderA();

    void      dumpMesh();

  public:
    // Mesh entities -- sets
    ZoneSet zones;                  // Set of zones in the mesh -- PositionSet (wrapper around an int)
    NodeSet nodes;                  // Set of nodes in the mesh -- PositionSet (wrapper around an int)
    FaceSet faces;                  // Set of faces in the mesh -- added in SLAM version (wrapper around an int)
    CornerSet corners;              // Set of corners in the mesh -- added in SLAM version (wrapper around an int)

    // Mesh topology -- relations
    ZoneToFaceRelation zoneToFaces;
    ZoneToNodeRelation zoneToNodes;

    // Mesh geometry -- fields
    NodalVectorField nodePos;         // node positions
    ZonalVectorField zonePos;         // zone positions (barycenters)
    FaceVectorField faceArea;         // face areas -- scaled outward-facing normals
    ZonalScalarField zoneVolume;      // zone volumes (scalar)

    double timeElapsed;
  };

} // end namespace tinyHydro

#endif      // __PolygonMeshXY_hh__
