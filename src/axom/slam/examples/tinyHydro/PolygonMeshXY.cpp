// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "PolygonMeshXY.hpp"

#include "axom/core/utilities/Timer.hpp"
#include "axom/slic.hpp"

#include <sstream>

namespace tinyHydro {


//----------------------------------------------
  PolygonMeshXY::PolygonMeshXY(int kmax, int lmax,
      double xmin, double xmax,
      double ymin, double ymax)
  {
    // This c'tor constructs a mesh that looks like a logically regular
    // mesh of quads, but in fact it is arbitrary polygons.  It's
    // intended for making quick test meshes and to do direct speed
    // comparisons between the polygon and regular quad mesh.

    SLIC_ASSERT(kmax > 0 && lmax > 0);

    nodes = NodeSet(kmax * lmax);
    zones = ZoneSet((kmax - 1) * (lmax - 1));
    faces = FaceSet( zones.size() * 4);
    corners = CornerSet( zones.size() * 4);

  #ifndef TINY_HYDRO_REGULAR_MESH
    /// Setup mesh topology
    zoneToNodes = ZoneToNodeRelation(&zones, &nodes);
    zoneToFaces = ZoneToFaceRelation(&zones, &faces);

    IndexBuffer& nodeIndices = DataRegistry::setRegistry.addBuffer("zone_node_relation_data", corners.size());
    IndexBuffer& faceIndices = DataRegistry::setRegistry.addBuffer("zone_face_relation_data", faces.size());


    // set up connectivity data, knowing it's logical orthogonal quads
    const int nZones = numZones();
    for (int i = 0; i < nZones; i++)
    {
      int k = i % (kmax - 1);
      int l = i / (kmax - 1);

      nodeIndices[4 * i + 0] = k     + l * kmax;
      nodeIndices[4 * i + 1] = k + 1 + l * kmax;
      nodeIndices[4 * i + 2] = k + 1 + (l + 1) * kmax;
      nodeIndices[4 * i + 3] = k     + (l + 1) * kmax;

      faceIndices[4 * i]   = 4 * i;
      faceIndices[4 * i + 1] = 4 * i + 1;
      faceIndices[4 * i + 2] = 4 * i + 2;
      faceIndices[4 * i + 3] = 4 * i + 3;
    }
    zoneToNodes.bindIndices( nodeIndices.size(), &nodeIndices );
    zoneToFaces.bindIndices( faceIndices.size(), &faceIndices );
  #endif

    /// Geometric fields on the mesh
    nodePos = NodalVectorField(&nodes);
    zoneVolume = ZonalScalarField(&zones);
    zonePos = ZonalVectorField(&zones);
    faceArea = FaceVectorField(&faces);

    // set up initial mesh as regular orthonormal grid
    double dx = (xmax - xmin) / (kmax - 1);
    double dy = (ymax - ymin) / (lmax - 1);

    for (int l = 0; l < lmax; l++)
    {
      for (int k = 0; k < kmax; k++)
      {
        nodePos[k + l * kmax].x = xmin + k * dx;
        nodePos[k + l * kmax].y = ymin + l * dy;
      }
    }

    computeNewGeometry();

    //dumpMesh();
  }


//----------------------------------------------
  void PolygonMeshXY::computeNewGeometry(void)
  {
    SLIC_ASSERT(numZones() > 0 && numNodes() > 2);

    // zone volumes and positions
    const int nZones = numZones();
    for (int iz = 0; iz < nZones; iz++)
    {
      ZNodeSet zNodes = zoneToNodes[iz];
      ZFaceSet zFaces = zoneToFaces[iz];

      double zVol = 0.;
      VectorXY zPos(0.,0.);

      const int numZNodes = zNodes.size();
      for (int in = 0; in < numZNodes; in++)
      {
        // assumption is that zNodes lists nodes in counter-clockwise order around the zone
        ZNodeSet::ModularIntType modIdx(in, numZNodes);
        const VectorXY& n0Pos = nodePos[ zNodes[modIdx  ] ];
        const VectorXY& n1Pos = nodePos[ zNodes[modIdx + 1] ];

        zVol += n0Pos.cross(n1Pos);     // sum volume
        zPos += n1Pos;                  // sum position

        // Face areas are outward normal vectors
        faceArea[zFaces[modIdx] ] = (n0Pos - n1Pos).perp();
      }
      zoneVolume[iz] = 0.5 * zVol; // correct 2x volume from cross product

      // normalize position; if we did it by volume it would be more
      // correct (centroid), but arithmetic average usually works fine
      zonePos[iz] =  (1.0 / numZNodes) * zPos;
    }

    // check our final volumes are positive
    for (int iz = 0; iz < nZones; iz++)
      SLIC_ASSERT(zoneVolume[iz] > 0.0);
  }
//----------------------------------------------
  void PolygonMeshXY::moveNodesToPosition(const NodalVectorField& newPos)
  {
    nodePos.copy(newPos);
  }


//----------------------------------------------

  VectorXY PolygonMeshXY::meshAverageKLZMemOrderA()
  {
    VectorXY ret;

    axom::utilities::Timer timer;

    timer.start();

    const int nZones = numZones();
    for (int i = 0; i < nZones; i++)
    {
      ret.accum(zonePos[i]);
    }
    ret *= (1. / nZones);

    timer.stop();

    SLIC_INFO("Elapsed time for averaging was " << timer.elapsedTimeInSec() << " seconds.");

    return ret;
  }
//----------------------------------------------



  void PolygonMeshXY::dumpMesh()
  {
    SLIC_INFO( axom::fmt::format("Mesh has {} nodes and {} zones", nodes.size(), zones.size()) );

    std::stringstream nodesStr;
    nodesStr << "Nodes";
    for(int i = 0; i< nodes.size(); ++i)
    {
      VectorXY p = getPos(i);
      nodesStr << axom::fmt::format("\n\t Node {} -- pos ({},{})", i, p.x, p.y);
    }
    SLIC_INFO(nodesStr.str() << "\n----\n");

    std::stringstream zonesStr;
    zonesStr << "Zones";
    for(int i = 0; i< zones.size(); ++i)
    {
      VectorXY p = getZonePos(i);
      ZoneToNodeRelation::RelationSubset zNodes = zoneToNodes[i];
      ZoneToNodeRelation::RelationSubset zFaces = zoneToFaces[i];
      zonesStr  << "\n\t Zone " << i
                << axom::fmt::format("-- pos ({},{})", p.x, p.y)
                << " -- vol " << zoneVol(i)
                << "-- zNumNodes " << zNodes.size()
                << axom::fmt::format("-- zoneNodes {} {} {} {}", zNodes[0], zNodes[1], zNodes[2], zNodes[3])
                << axom::fmt::format("-- zoneFaces {} {} {} {}", zFaces[0], zFaces[1], zFaces[2], zFaces[3]);
    }
    SLIC_INFO(zonesStr.str() << "\n----\n");

  }

} // end namespace tinyHydro
