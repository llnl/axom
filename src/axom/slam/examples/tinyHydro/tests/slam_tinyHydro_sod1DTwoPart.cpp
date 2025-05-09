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


void tinyHydroSod1D_2part()
{
  SLIC_INFO("\n\t***********************************\n"
      << "\t      Sod 1D test, two parts.\n"
      << "\t***********************************");


  double goalTime = 0.171;   // puts shock right on 0.80
  double gamma = 7.0 / 5.0;
  int steps = 0;

  // zoning
  int n = 200;
  double dx = 1.0 / static_cast<double>(n);
  int kz = n;
  int lz = 1;

  // initial parameters
  double Pp = 1.0;   // left pressure
  double dp = 1.0;   // left density
  double ep = Pp / (dp * (gamma - 1)); //left specific energy

  //double c0 = sqrt(gamma*Pp/dp); // left sound speed

  double Pm = 1.0 / 10.0; // right pressure
  double dm = 1.0 / 8.0; // right density
  double em = Pm / (dm * (gamma - 1)); // right specific energy

  PolygonMeshXY mesh(kz + 1,lz + 1);

  for(int i = mesh.numNodes() / 2; i< mesh.numNodes(); ++i )
  {
    VectorXY pos = mesh.getPos(i);
    mesh.setPos(i, VectorXY(pos.x, dx));
  }

  // Create a part
  std::vector<int> reg1, reg2;
  int halfZones = mesh.numZones() / 2;
  for(int i = 0; i< halfZones; ++i )
  {
    reg1.push_back(i);
    reg2.push_back(i + halfZones);
  }

  SLIC_ASSERT(  reg1.size() == reg2.size());
  SLIC_ASSERT(  reg1.size() + reg2.size() == static_cast<size_t>( mesh.numZones() ));

  Part p1(reg1, gamma);
  Part p2(reg2, gamma);
  State s(mesh);

  SLIC_INFO("**setting initial density, energy, velocity for sod problem");
  for (unsigned int i = 0; i < reg1.size(); i++)
  {
    p1.setRho(i,dp);
    p1.setE(i, ep);
  }
  for (unsigned int i = 0; i < reg2.size(); i++)
  {
    p2.setRho(i,dm);
    p2.setE(i, em);
  }
  for (int i = 0; i < mesh.numNodes(); i++)
    s.setU(i, VectorXY(0.,0.));

  s.addPart(&p1);
  s.addPart(&p2);

  SLIC_INFO("**making hydro");
  Hydro h(&s);


  SLIC_INFO("**Set boundary conditions!");
  h.setBC("bottom", 0xdeadbeef, 0);
  h.setBC("top",    0xdeadbeef, 0);
  h.setBC("left",   0,          0);
  h.setBC("right",  0,          0);


  SLIC_INFO("**Initialize");
  h.initialize();

  double dt = h.newDT();
  SLIC_INFO("  newDT is " << dt);
  h.cfl = 0.2;


  //bool testQ = false;
  //bool testP = false;
  bool special = false;
  double E0 = 0.;   // starting energy

  axom::utilities::Timer timer;
  timer.start();

  if (special)
  {
    dt = 0.2;
    h.Cq = 0.0;
    h.Cl = 0.0;
    E0 = h.totalEnergy();
    SLIC_INFO("  total starting energy = " << E0 );
    SLIC_INFO("**stepping..." );
    h.step(dt);
  }
  else if(steps > 0)
  {
    E0 = h.totalEnergy();
    SLIC_INFO("  total starting energy = " << E0 );
    SLIC_INFO("**stepping..." );
    h.steps(steps);
  }
  else
  {
    E0 = h.totalEnergy();
    SLIC_INFO("  total starting energy = " << E0 );
    SLIC_INFO("**advancing...");
    h.advance(goalTime);     // puts shock at 0.8 (p.57, Triangular Hydro Schemes notebook)
  }
  double E1 = h.totalEnergy();

  timer.stop();

  SLIC_INFO("--");
  SLIC_INFO("Elapsed time after advancing was " << timer.elapsedTimeInSec() << " seconds.");

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
  
  tinyHydroSod1D_2part();

  return result;
}
