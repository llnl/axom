// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// State class that holds all the material data.  In the
// multi-material world, this will mostly be a vector of Part structs
// that hold the data, plus a few functions to do mesh sums and
// averages of quantities.

#include "axom/slam.hpp"

#include "VectorXY.hpp"
#include "TinyHydroTypes.hpp"
#include "Part.hpp"


namespace tinyHydro {


  class PolygonMeshXY;


  class State
  {
  public:
    State(PolygonMeshXY & theMesh);
    ~State(void);

    // Copy
    State( const State & arg);
    // Assignment
    State & operator=(const State & rhs);

    // addition in place
    State & operator+=(const State & rhs);
    // scalar multiply in place
    State & operator*=(const double s);

    // python accessors
    /* double averageRho(int i) const; */
    /* double totalE(int i) const; */
    VectorXY  u(int i)            const {return velocity[i]; }

    // Related to the part of a state
    Part *    partBegin(void)                {return parts; }
    Part *    getPart(int i);
    void      addPart(Part * newPart);

    // Access positions and velocities
    void      setU(int i, const VectorXY& val) { velocity[i] = val; }
    void      setX(int i, const VectorXY& val) { position[i] = val; }

    void      dumpState();

  public:
    PolygonMeshXY * const mesh;

    Part * parts;
    int nParts;
    const int maxNParts;

    NodalVectorField velocity;
    NodalVectorField position;
  };


} // end namespace tinyHydro
