// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


// Hydro driver
// Fri Mar 27 14:14:28 PDT 2015

// HYDRO C: Physics-wise, this hydro is still pretty dumb: multiple
// materials, gamma-law gas.  The mesh is an arbitrary collection of
// polygons.  Architecturally, we use as few temporaries as possible.
// No sweet returning of objects when you can instead pass in an array
// and modify it in place.  No STL containers.  Explicit memory
// allocation of all the memory we ever need at problem start, and
// explicit deletes when the hydro object is destroyed.

#include "State.hpp"

#include "TinyHydroTypes.hpp"

namespace tinyHydro {

  class Hydro
  {
  public:

  public:
    Hydro(State * s); // using a pointer here with PyBindGen prevents
                      // State being garbage collected by Python
    ~Hydro();

    // problem control and advance methods
    void      step(double dt);
    double    newDT();
    void      steps(int numSteps);
    void      advance(double stopTime);
    void      initialize(void); // call once before taking steps
    void      setBC(const char * boundary, double xVel, double yVel);

    // compute time derivatives from first arg, store rates in 2nd arg
    void      calcDerivs(const State & state, State & dState, double dt);

    // physics methods
    void      calcQ(const State & state);
    void      calcForce(const State & state);
    void      calcDivU(const NodalVectorField& u);
    void      applyVelocityBC(NodalVectorField& u);
    double    totalEnergy(const State & s) const;
    double    totalEnergy(void) const;
    void      calcMaxCs(void);

    // update all the algebraically determinable stuff, like density
    void      updateConstitutives(State & state);

    // setter methods (used from Python)
    void      setState(const State s) {state = s; }

    // getter methods
    State *   getState(void) {return &state; }
    double    getQ(int i) {return Q[i]; }
    VectorXY  getForce(int i) {return force[i]; }
    int       numBCnodes(int bc);
    int       bcNode(int bc, int node);

    bool      fuzzyEqual(double a, double b, double tol = 1.0e-15)
    {
      return (fabs(a - b) < tol);
    }

  public:
    PolygonMeshXY&       mesh;

    State&               state;
    State halfStep;                         // re-used temporary for RK2
    State dState;                           // re-used temporary for RK2

    ZonalScalarField Q;                     // re-used for Q
    ZonalScalarField divu;                  // for dt, Q, and some PdV work schemes

    NodalVectorField force;                 // re-used for force calc
    NodalVectorField halfStepVelocity;      // re-used for work calc

    BoundaryEdgeSet boundaryEdgeSet;                                  // SLAM -- new set for domain boundaries (of size 4)
    BoundaryEdgeVectorField bcVelocity;     // velocity BC values     // SLAM -- map to velocity boundary conditions for each edge
    NodeSubset bcNodes[NUM_DOMAIN_BOUNDARIES];                        // SLAM -- indirection sets of nodes on each boundary edge


    ZonalScalarField*    initialMass;   // per material initial mass of zones

    ZonalScalarField totalMass;         // overall mass in zone
    NodalScalarField nodeMass;          // also keep node masses

    ZonalScalarField*    pressure;      // per material zone pressure

    ZonalScalarField totalPressure;     // sum of pressures in a zone
    ZonalScalarField maxCs;             //


    int cycle;
    int dtZone;
    double time;
    double Cq, Cl;                      // Q coefficients
    double cfl;
    int printCycle;

  private:
    // Unused variables...
    // int numBCs;                                                      // <UNUSED>?

#ifdef TINY_HYDRO_USE_CORNER_DATA
    VectorXY * cornerforce; // re-used for force calc                // <UNUSED>?
#endif
  };

} // end namespace tinyHydro
