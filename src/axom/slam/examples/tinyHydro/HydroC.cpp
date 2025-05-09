// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <algorithm>
#include <string.h>

#include "HydroC.hpp"
#include "PolygonMeshXY.hpp"
#include "Part.hpp"
// #include "/usr/include/gperftools/profiler.h"

#include "axom/slic.hpp"

namespace tinyHydro {

//----------------------------------------------
// c'tor
  Hydro::Hydro(State * s)
      : mesh(*(s->mesh)),
        state(*s), halfStep(*s), dState(*s),
        cycle(0), dtZone(-1), time(0.0),
        Cq(1.0), Cl(1.0), cfl(0.7)
  {
    SLIC_INFO("using HydroC");

    SLIC_ASSERT(state.nParts > 0);
    SLIC_ASSERT(mesh.numZones() > 0);
    SLIC_ASSERT(mesh.numNodes() > 0);

    // make memory for Q calculation
    Q     = ZonalScalarField(&mesh.zones);
    force = NodalVectorField(&mesh.nodes);

#ifdef TINY_HYDRO_USE_CORNER_DATA
    cornerforce = new VectorXY * [nParts];
    for (int p = 0; p < nParts; p++)
      cornerforce[p] = new VectorXY[4 * state.parts[p].nzones];
#endif

    const int nParts = state.nParts;

    // make memory for holding zone and node masses
    initialMass = new ZonalScalarField[nParts];
    for (int p = 0; p < nParts; p++)
      initialMass[p] = ZonalScalarField(&state.parts[p].zones);

    totalMass        = ZonalScalarField( &mesh.zones);
    totalPressure    = ZonalScalarField( &mesh.zones);
    maxCs            = ZonalScalarField( &mesh.zones);

    nodeMass         = NodalScalarField( &mesh.nodes);

    // space for half-step velocity
    halfStepVelocity = NodalVectorField(&mesh.nodes);

    // for dt, Q, and some PdV work schemes
    divu = ZonalScalarField( &mesh.zones);


    // zone pressure
    pressure = new ZonalScalarField[nParts];
    for (int p = 0; p < nParts; p++)
      pressure[p] = ZonalScalarField(&state.parts[p].zones);


    // velocity BC space.  BC's are implicitly ordered as 'bottom',
    // 'right', 'top', 'left', numbered 0 through 3.
    // initial BC's are 0xdeadbeef, which we take as a special value
    // meaning free boundary condition.
    bcVelocity = BoundaryEdgeVectorField(&boundaryEdgeSet, VectorXY(0xdeadbeef,0xdeadbeef) );

  #ifndef TINY_HYDRO_REGULAR_MESH
    // compute bounding box of domain
    VectorXY bbLower( mesh.nodePos[0]);
    VectorXY bbUpper( mesh.nodePos[0]);
    for (int n = 1; n < mesh.numNodes(); n++)
    {
      const VectorXY & pos = mesh.nodePos[n];
      if (pos.x < bbLower.x)
        bbLower.x = pos.x;
      if (pos.x > bbUpper.x)
        bbUpper.x = pos.x;
      if (pos.y < bbLower.y)
        bbLower.y = pos.y;
      if (pos.y > bbUpper.y)
        bbUpper.y = pos.y;
    }

    // Find nodes on each boundary
    IndexBuffer* bcNodeList[NUM_DOMAIN_BOUNDARIES];
    for(int i = 0; i< boundaryEdgeSet.size(); ++i)
    {
      bcNodeList[i] = &DataRegistry::setRegistry.addBuffer(axom::fmt::format("bc_nodes_{}",i));
    }

    for (int n = 0; n < mesh.numNodes(); n++)
    {
      const VectorXY & pos = mesh.nodePos[n];
      if (fuzzyEqual(pos.y, bbLower.y))
        bcNodeList[0]->push_back(n);
      if (fuzzyEqual(pos.x, bbUpper.x))
        bcNodeList[1]->push_back(n);
      if (fuzzyEqual(pos.y, bbUpper.y))
        bcNodeList[2]->push_back(n);
      if (fuzzyEqual(pos.x, bbLower.x))
        bcNodeList[3]->push_back(n);
    }
    // Setup the node lists as indirection sets
    for(int i = 0; i< boundaryEdgeSet.size(); ++i)
    {
      bcNodes[i] =  NodeSubset::SetBuilder()
          .size(bcNodeList[i]->size())
          .data(bcNodeList[i]);
    }
  #endif
  }

//----------------------------------------------
// d'tor
  Hydro::~Hydro()
  {
    // delete [] cornerforce;
    delete [] initialMass;
    delete [] pressure;
  }
//----------------------------------------------
// set hydro ready to begin taking steps
  void Hydro::initialize(void)
  {
    // calc new geometry data in mesh, including new zone volumes
    mesh.computeNewGeometry();

    // get our state ready with initial mesh positions
    state.position.copy( mesh.nodePos );

    // tuck away initial zone and node masses and never forget them
    totalMass.clear();
    nodeMass.clear();
    for (int p = 0; p < state.nParts; p++) // loop over parts
    {
      Part& part = state.parts[p];
      ZoneSubset& zones = part.zones;
      const ZonalScalarField& rho = part.density;
      for (int iz = 0; iz < part.numZones(); ++iz) // loop over zones
      {
        const IndexType meshZoneId = zones[iz];
        double massInit = rho[iz] * mesh.zoneVolume[ meshZoneId ];
        initialMass[p][iz] = massInit;
        totalMass[meshZoneId] += massInit;

        ZNodeSet zNodes = mesh.zoneToNodes[ meshZoneId ];
        const int znSize = zNodes.size();
        for (int in = 0; in < znSize; ++in)  // loop over nodes
        {
          const int node = zNodes[in];
          nodeMass[node] += 0.25 * massInit;
        }
      } // end loop over zones
    } // end loop over parts
  }

//----------------------------------------------
// computes the time derivatives, stores them in 2nd arg
  void Hydro::calcDerivs(const State & s, State & dState, double dt)
  {
    // printf("Hydro:: blago calcDerivs 0\n");

    const int nnodes = mesh.numNodes();

    // calc forces
    calcForce(s);

    // Accelerations are the delta velocities for dState, so let's
    // label that explicitly.  Re-naming makes it less confusing, at
    // least to me.
    NodalVectorField& accel = dState.velocity;

    // calc accelerations:
    for (int n = 0; n < nnodes; n++)
    {
      accel[n] = (1. / nodeMass[n]) * force[n];
    }

    // halfStepVelocity = state.u + 0.5*dt*accel

    // Velocities must be anchored at the half-way point in time to get
    // the conservation properties we like. Note that state.u is the
    // beginning of timestep velocity, which we use even during the
    // second partial RK2 step.
    for (int i = 0; i < nnodes; i++)
    {
      halfStepVelocity[i] = accel[i] * 0.5 * dt + state.velocity[i];
    }

    // apply velocity BCs to half-step velocity to prevent
    // spurious position changes and PdV work
    applyVelocityBC(halfStepVelocity);

    // velocities are the delta of position for dState.
    dState.position.copy(halfStepVelocity);

    // calc rate of (P+Q)dV work: to be energy conservative we want
    // the velocity divergence computed on the half-step velocity
    calcDivU(halfStepVelocity);

#ifdef TINY_HYDRO_USE_CORNER_DATA
    ////// interestingly, in the finite volume approach, the
    ////// force.dot(velocity) energy update is algebraically identical
    ////// to the simpler P*div(u) calculation.  So we do the simpler one.
    const int kmax = mesh.kmax;
    for (int z = 0; z < mesh.nzones; z++)
      * /
      {
        dState.energyPerMass[z] =
            -( cornerforce[4 * z].dot(halfStepVelocity[LL(z)])
            + cornerforce[4 * z + 1].dot(halfStepVelocity[LR(z)])
            + cornerforce[4 * z + 2].dot(halfStepVelocity[UR(z)])
            + cornerforce[4 * z + 3].dot(halfStepVelocity[UL(z)])
            ) /
            initialMass[z];
      }
#endif

    // the energy production rate is simply -(P+Q) div(u)
    for (int p = 0; p < s.nParts; p++)
    {
      const ZonalScalarField& rho = s.parts[p].density;
      const ZonalScalarField& P = pressure[p];
      ZonalScalarField& de = dState.parts[p].energyPerMass;
      const int nzones = state.parts[p].numZones();
      const ZoneSubset& zones = state.parts[p].zones;
      for (int z = 0; z < nzones; z++)
      {
        const IndexType meshZoneId = zones[z];
        de[z] = -(P[z] + Q[meshZoneId]) * divu[meshZoneId] / rho[z];
      }
    }
    // printf("Hydro:: blago calcDerivs 9\n");

  }
//----------------------------------------------
// compute the total force that pushes the nodes
  void Hydro::calcForce(const State & s)
  {
    SLIC_ASSERT(s.mesh == &mesh);

    // compute grad of zone-centered (P+Q), return pointer to array of
    // node-centered VectorXY forces.
    force.clear();

    // zero out total before summing
    totalPressure.clear();

    // get pressure from EOS: gamma law gas, P = rho*e*(gamma-1)
    for (int p = 0; p < s.nParts; p++)
    {
      Part& part = s.parts[p];
      const ZoneSubset& zones = part.zones;
      const int nzones = zones.size();
      const double gammaMinusOne = part.gamma - 1.0;
      const ZonalScalarField& rho = part.density;
      const ZonalScalarField& e = part.energyPerMass;
      ZonalScalarField& P = pressure[p];
      for (int z = 0; z < nzones; z++)
      {
        const IndexType meshZoneId = zones[z];
        double zonePr = gammaMinusOne * rho[z] * e[z];
        P[z] = zonePr;
        totalPressure[meshZoneId] += zonePr;
      }
    }

    // get Q
    calcQ(state);

    // compute grad via finite volume method
    // grad(P) = (1/Vn) sum(P*A) for P,A around node
    // for each zone, accumulate P*A to its nodes
    const int nZones = mesh.numZones();
    for (int iz = 0; iz < nZones; iz++)
    {
      const double totP = 0.5 * (totalPressure[iz] + Q[iz]);

      ZNodeSet zNodes = mesh.zoneToNodes[ iz];
      ZFaceSet zFaces = mesh.zoneToFaces[ iz];

      const int relnSize = zNodes.size();
      ZNodeSet::ModularIntType modIdx(-1, relnSize);    // start at one-before first face
      for (int in = 0; in < relnSize; ++in)
      {
        const int face0 = zFaces[ modIdx++];    // get idx of previous face and increment
        const int face1 = zFaces[ modIdx  ];
        const int node  = zNodes[ modIdx  ];

        // multiply areas by pressure + Q (and div 2 for area factor)
        force[node] += totP * (mesh.faceArea[face0] + mesh.faceArea[face1]);
      }
    }
  }
//----------------------------------------------
// compute the artificial viscosity
  void Hydro::calcQ(const State & s)
  {
    // Q = rho((L*div(u))**2 - cs*L*div(u))

    // calculate divergence of velocity
    calcDivU(s.velocity);

    // set max sound speed
    calcMaxCs();

    // compute Q for each zone
    const int nZones = mesh.numZones();
    for (int z = 0; z < nZones; ++z)
    {
      if (divu[z] < 0.0)
      {
        const double vol = mesh.zoneVolume[z];
        const double L = sqrt(vol);
        const double deltaU = L * divu[z];
        const double rho = totalMass[z] / vol;
        Q[z] = rho * (Cq * deltaU * deltaU - Cl * maxCs[z] * deltaU);
      }
      else
        Q[z] = 0.0;
    }
  }

//----------------------------------------------
// update all the algebraically determinable stuff
  void Hydro::updateConstitutives(State & s)
  {
    // move mesh to new node positions
    mesh.moveNodesToPosition(s.position);

    // calc new geometry data in mesh, including new zone volumes
    mesh.computeNewGeometry();

    for (int p = 0; p < s.nParts; p++)
    {
      Part& part = s.parts[p];
      const ZoneSubset& zones = part.zones;
      const int nzones = part.numZones();
      ZonalScalarField& rho = part.density;
      const ZonalScalarField& massInit = initialMass[p];
      const ZonalScalarField& zoneVol = mesh.zoneVolume;
      for (int z = 0; z < nzones; z++)
      {
        rho[z] = massInit[z] / zoneVol[ zones[z] ];
      }
    }

    // apply BC's to newly computed velocity
    applyVelocityBC(s.velocity);
  }

//----------------------------------------------
  void Hydro::step(double dt)
  {
    SLIC_INFO("cycle " << cycle << " time = " << time << " dt = " << dt);
    // printf("Hydro::step blago 1\n");

    // ----------------------------------
    // update physics variables using RK2
    // ----------------------------------
    // printf("Hydro::step blago 2\n");
    time += 0.5 * dt;
    calcDerivs(state, dState, 0.5 * dt); // dState = time derivatives
    // halfstep = state + dt/2 * dState
    halfStep = state;

    dState *= 0.5 * dt;
    halfStep += dState;

    // update all the algebraically determinable stuff
    updateConstitutives(halfStep);

    // second step of our RK2 integration
    calcDerivs(halfStep, dState, dt);
    // final state = state + dt * dState(halfStep)
    dState *= dt;
    state += dState;

    // update all the algebraically determinable stuff
    updateConstitutives(state);

    // update control info
    time += 0.5 * dt; // full time update
    cycle++;

    // printf("Hydro::step blago 9\n");

  }
//----------------------------------------------
// calc max sound speed in each zone
  void Hydro::calcMaxCs(void)
  {
    // zero to start with
    maxCs.clear();

    for (int p = 0; p < state.nParts; p++)
    {
      const Part& part = state.parts[p];
      const ZoneSubset& zones = part.zones;
      const int nzones = part.numZones();
      const double gammaGammaMinusOne = part.gamma * (part.gamma - 1.0);
      const ZonalScalarField& e = part.energyPerMass;
      for (int z = 0; z < nzones; z++)
      {
        const IndexType meshZoneID = zones[z];
        double cs2 = gammaGammaMinusOne * e[z];
        maxCs[meshZoneID] = std::max(maxCs[meshZoneID], cs2);
      }
    }
    // we saved all the sqrts for the very end
    const int nZones = mesh.numZones();
    for (int z = 0; z < nZones; ++z)
    {
      maxCs[z] = std::sqrt(maxCs[z]);
    }
  }

//----------------------------------------------
  double Hydro::newDT(void)
  {
    double minDT = 1.0e99; // impossibly high dt

    dtZone = -1; // impossible zone index

    calcDivU(state.velocity); // we need velocity divergence
    calcMaxCs(); // need max sound speed

    const int nZones = mesh.numZones();
    for (int z = 0; z < nZones; z++)
    {
      double L2 = mesh.zoneVolume[z];
      double mc = maxCs[z];
      double du = divu[z];
      double dt2 = L2 / (mc * mc + L2 * du * du);
      if (minDT >= dt2)
      {
        minDT = dt2;
        dtZone = z;
      }
    }

    //printf("\n\t\tMin DT is %g -- computed in zone %i \n\n", minDT, dtZone);

    return cfl * sqrt(minDT);
  }
//----------------------------------------------
// compute the divergence of velocity
  void Hydro::calcDivU(const NodalVectorField& velocity)
  {
    // div u = 1/V * sum(u dot dA)
    const int nZones = mesh.numZones();

    for (int iz = 0; iz < nZones; iz++)
    {
      ZNodeSet zNodes = mesh.zoneToNodes[ iz];
      ZFaceSet zFaces = mesh.zoneToFaces[ iz];
      const int relnSize = zNodes.size();
      double du = 0.;
      ZNodeSet::ModularIntType modIdx(-1, relnSize);    // start at one-before first face
      for (int in = 0; in < relnSize; in++)
      {
        const int face0 = zFaces[modIdx++];
        const int node  = zNodes[modIdx  ];
        const int face1 = zFaces[modIdx  ];

        VectorXY A = mesh.faceArea[face0] + mesh.faceArea[face1];
        du += velocity[node].dot(A);
      }
      divu[iz] = du / (2.0 * mesh.zoneVolume[iz]); // the 2 normalizes the areas
    }
  }

//----------------------------------------------
  void Hydro::setBC(const char * boundary, double xVel, double yVel)
  {
    switch (boundary[0] )
    {
    case 'b':
      bcVelocity[0] = VectorXY(xVel, yVel);
      break;
    case 'r':
      bcVelocity[1] = VectorXY(xVel, yVel);
      break;
    case 't':
      bcVelocity[2] = VectorXY(xVel, yVel);
      break;
    case 'l':
      bcVelocity[3] = VectorXY(xVel, yVel);
      break;
    default:
      SLIC_DEBUG("ERROR, bad value specified for BC\n");
      break;
    }
  }
//----------------------------------------------
  void Hydro::applyVelocityBC(NodalVectorField& u)
  {
    const int bSize = boundaryEdgeSet.size();

    for(int i = 0; i < bSize; ++i)
    {
      NodeSubset& nodes = bcNodes[i];
      if (bcVelocity[i].x != 0xdeadbeef)
      {
        for(int j = 0; j < nodes.size(); ++j)
          u[ nodes[j] ].x = bcVelocity[i].x;
      }
      if (bcVelocity[i].y != 0xdeadbeef)
      {
        for(int j = 0; j < nodes.size(); ++j)
          u[ nodes[j] ].y = bcVelocity[i].y;
      }
    }
  }
//----------------------------------------------
  void Hydro::steps(int numSteps)
  {
    for (int i = 0; i < numSteps; i++)
    {
      double dt = newDT();
      step(dt);
    }
  }
//----------------------------------------------
  void Hydro::advance(double stopTime)
  {
    /* ProfilerStart("profile.out"); */
    /* struct ProfilerState ps; */

    while (time < stopTime)
    {
      double dt = newDT();
      step(dt);
    }

    /* ProfilerGetCurrentState(&ps); */
    /* printf("profiler: enabled = %d, samples gathered = %d\n", ps.enabled, ps.samples_gathered); */
    /* ProfilerStop(); */

  }
//----------------------------------------------
  double Hydro::totalEnergy(const State & s) const
  {
    double totE = 0.0;

    for (int p = 0; p < s.nParts; p++)
    {
      const int nzones = s.parts[p].numZones();
      const ZonalScalarField& e = s.parts[p].energyPerMass;
      const ZonalScalarField& m = initialMass[p];
      for (int z = 0; z < nzones; z++)
      {
        totE += e[z] * m[z];
      }
    }

    const int nNodes = mesh.numNodes();
    for (int n = 0; n < nNodes; n++)
    {
      const VectorXY& vel = s.velocity[n];
      totE += 0.5 * nodeMass[n] * vel.dot(vel);
    }
    return totE;
  }
//----------------------------------------------
  int Hydro::numBCnodes(int bc)
  {
    return bcNodes[bc].size();
  }
//----------------------------------------------
  int Hydro::bcNode(int bc, int node)
  {
    return bcNodes[bc][node];
  }
//----------------------------------------------
  double Hydro::totalEnergy(void) const
  {
    return totalEnergy(state);
  }
//----------------------------------------------


} // end namespace tinyHydro
