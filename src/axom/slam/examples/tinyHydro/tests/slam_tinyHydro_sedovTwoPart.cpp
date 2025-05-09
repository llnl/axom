// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#include "axom/config.hpp"
#include "axom/core/utilities/Timer.hpp"
#include "axom/slic.hpp"

#include "../PolygonMeshXY.hpp"
#include "../HydroC.hpp"

#include <cmath>
#include <vector>

using namespace tinyHydro;

std::ostream& operator<<(std::ostream& os, const VectorXY & vec)
{
  os << "(" << vec.x << "," << vec.y << ")";
  return os;
}


void tinyHydroSedov_2part()
{
  SLIC_INFO("\n\t***********************************\n"
      << "\tSedov cylindrical square grid test.\n"
      << "\t***********************************");

  int steps = 0;
  int n = 80;
  //double dx = 1.0 / n;
  int kz = n;
  int lz = n;
  PolygonMeshXY mesh(kz + 1,lz + 1);

  // Create a part
  std::vector<int> reg1, reg2;

  for(int i = 0; i< mesh.numZones(); ++i )
  {
    VectorXY pos = mesh.getZonePos(i);
    if(pos.y > pos.x)
      reg1.push_back(i);
    else
      reg2.push_back(i);
  }

  State s(mesh);

  double gamma = 7.0 / 5.0;
  Part p1(reg1, gamma);
  Part p2(reg2, gamma);

  SLIC_INFO("**setting initial density, energy, velocity for sedov problem");
  double rho0 = 1.0;
  for (unsigned int i = 0; i < reg1.size(); i++)
  {
    p1.setRho(i,rho0);
    p1.setE(i, 0.);
  }
  for (unsigned int i = 0; i < reg2.size(); i++)
  {
    p2.setRho(i,rho0);
    p2.setE(i, 0.);
  }
  for (int i = 0; i < mesh.numNodes(); i++)
    s.setU(i, VectorXY(0.,0.));


  // set energy spike at the origin
  double zonemass = rho0 * mesh.zoneVol(0);
  double E = 4.48333289847;
  p2.setE(0, E / (4.0 * zonemass)); // puts shock at r=0.8 at t=0.3, in XY coords

  s.addPart(&p1);
  s.addPart(&p2);


  SLIC_INFO("**making hydro");
  Hydro h(&s);

  SLIC_INFO("**Set boundary conditions!");
  h.setBC("bottom", 0xdeadbeef, 0);     // x can slide, y fixed
  h.setBC("left",   0,          0xdeadbeef); // x fixed, y can slide
  h.setBC("top",    0,          0);     // fixed x,y
  h.setBC("right",  0,          0);     // fixed x,y


  SLIC_INFO( "**Initialize" );
  h.initialize();

  h.cfl = 0.1;

  double E0 = h.totalEnergy();

  axom::utilities::Timer timer;
  timer.start();
  if(steps > 0)
  {
    h.steps(steps);
  }
  else
  {
    h.advance(0.3);
  }
  timer.stop();

  SLIC_INFO("--");
  SLIC_INFO("Elapsed time after advancing was " << timer.elapsedTimeInSec() << " seconds.");


  double E1 = h.totalEnergy();

  SLIC_INFO("\t Total final energy: " << E1 );
  SLIC_INFO("\t (E1-E0)/E0 = : " << (E1 - E0) / E0 );

  SLIC_INFO(" *** PASS ***");
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
int main()
{
  int result = 0;
  axom::slic::SimpleLogger logger;

  tinyHydroSedov_2part();

  return result;
}
