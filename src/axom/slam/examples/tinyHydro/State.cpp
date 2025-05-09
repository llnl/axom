// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#include <sstream>

#include "State.hpp"

#include "PolygonMeshXY.hpp"
#include "Part.hpp"

#include "axom/slic.hpp"

namespace tinyHydro {


//----------------------------------------------
  State::State(PolygonMeshXY & theMesh)
      : mesh(&theMesh),
        nParts(0),
        maxNParts(100)
  {
    SLIC_DEBUG("in State c'tor");

    parts = new Part[maxNParts];

    velocity = NodalVectorField( &mesh->nodes);
    position = NodalVectorField( &mesh->nodes);
  }
//----------------------------------------------
  State::~State(void)
  {
    delete [] parts;
  }

//----------------------------------------------
  Part * State::getPart(int i)
  {
    return &(parts[i]);
  }
//----------------------------------------------
// addPart
  void State::addPart(Part * partPtr)
  {
    SLIC_DEBUG("in State::addPart");

    SLIC_ASSERT_MSG( nParts < maxNParts,
        "tried to add more than " << maxNParts << " parts, limit is hard-wired allowed, craaashing!");

    parts[nParts] = *partPtr; // copy in part data to the list
    nParts++; // now is next unused Part ptr
  }

//----------------------------------------------
// Copy operator
  State::State(const State & arg)
      : mesh(arg.mesh),
        nParts(arg.nParts),
        maxNParts(arg.maxNParts)
  {
    // printf("in State::copy \n");

    // copy over the parts
    parts = new Part[maxNParts];
    for (int i = 0; i < nParts; i++)
    {
      parts[i] = arg.parts[i];
    }

    // copy other data
    velocity = NodalVectorField(arg.velocity);
    position = NodalVectorField(arg.position);
  }

//----------------------------------------------
// Assignment operator
  State & State::operator=(const State & rhs)
  {
    // printf("State::assignment \n");

    // check for self-assignment
    if (this == &rhs)
      return *this;

    // assignment can only happen when states
    // are defined on the same mesh
    SLIC_ASSERT(mesh == rhs.mesh);
    SLIC_ASSERT(maxNParts > rhs.nParts);

    // free up old parts data and copy rhs's parts
    delete [] parts;
    parts = new Part[maxNParts];
    nParts = rhs.nParts;
    for (int i = 0; i < nParts; i++)
    {
      parts[i] = rhs.parts[i];
    }

    // copy other data from rhs
    velocity = NodalVectorField( rhs.velocity);
    position = NodalVectorField( rhs.position);

    return *this;
  }

//----------------------------------------------
// addition in place
  State & State::operator+=(const State & rhs)
  {
    // printf("in State::operator+= \n");
    SLIC_ASSERT(mesh == rhs.mesh);
    SLIC_ASSERT(nParts == rhs.nParts);

    // part data
    for (int i = 0; i < nParts; i++)
    {
      parts[i] += rhs.parts[i];
    }

    // mesh-based data
    const int numNodes = mesh->numNodes();
    for (int i = 0; i < numNodes; i++)
    {
      velocity[i] += rhs.velocity[i];
      position[i] += rhs.position[i];
    }

    return *this;
  }
//----------------------------------------------
// scalar multiply in place
  State & State::operator*=(const double s)
  {
    // printf("in State::operator*= \n");

    // part data
    for (int i = 0; i < nParts; i++)
    {
      parts[i] *= s;
    }

    // mesh-based data
    const int numNodes = mesh->numNodes();
    for (int i = 0; i < numNodes; i++)
    {
      velocity[i] *= s;
      position[i] *= s;
    }

    return *this;
  }
//----------------------------------------------
//----------------------------------------------

  void State::dumpState()
  {
    SLIC_INFO(  "\nState has " << nParts << " parts");

    SLIC_INFO(  "\nParts");
    for(int i = 0; i< nParts; ++i)
    {
      parts[i].dumpPart();
    }

    std::stringstream velStr;
    velStr << "\nVelocities";
    for(int i = 0; i< mesh->numNodes(); ++i)
    {
      VectorXY v = velocity[i];
      VectorXY p = position[i];
      velStr  << "\n\t "
              << "Node " << i
              << "-- pos (" << p.x << ", " << p.y << ") "
              << "-- vel (" << v.x << ", " << v.y << ") ";
    }

    SLIC_INFO( velStr.str() << "\n----\n");

  }


} // end namespace tinyHydro
