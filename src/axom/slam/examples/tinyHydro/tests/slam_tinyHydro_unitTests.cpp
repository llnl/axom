// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#include "gtest/gtest.h"

#include "axom/core.hpp"
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



TEST(slam_tinyHydro,test_02_density_with_prescribed_velocity)
{
  SLIC_INFO("\t******************************************************\n"
      << "\tmove mesh with prescribed velocity, test density change.\n"
      << "\t*******************************************************");

  int n = 100;
  int kz = n;
  int lz = n;
  PolygonMeshXY mesh(kz + 1,lz + 1);

  // Create a part
  std::vector<int> zones( mesh.numZones() );

  for(int i = 0; i< mesh.numZones(); ++i )
    zones[i] = i;

  Part p(zones);
  State s(mesh);

  SLIC_INFO("**setting initial densities");
  double rho0 = 1.0;
  for (int i = 0; i < mesh.numZones(); i++)
  {
    p.setRho(i,rho0);
    p.setE(i,0.0);
  }

  // set energy spike at the origin
  //double zonemass = rho0 * mesh.zoneVolume[0];
  //double E = 4.48333289847;
  //p.setE(0, E/(4.0*zonemass)); // puts shock at r=0.8 at t=0.3, in XY coords


  SLIC_INFO("**setting initial velocities");
  for (int i = 0; i < mesh.numNodes(); i++)
  {
    VectorXY pos = mesh.getPos(i);
    s.setU(i, VectorXY(0.5 - pos.x, 0.5 - pos.y));
  }

  s.addPart(&p);
  Part* pp = s.getPart(0);
  AXOM_UNUSED_VAR(pp);

  SLIC_INFO("**making hydro");
  Hydro h(&s);

  SLIC_INFO("**Before we begin:");
  h.Cl = 0.;
  h.Cq = 0.;

  double dt = 0.1;
  h.initialize();


  h.step(dt);

  SLIC_INFO("**Taking steps:");
  while( h.time < 0.499999 )
    h.step(dt);

  SLIC_INFO("\t hydro cycle = " << h.cycle
                                << " hydro time  = " << h.time);

  double rhoTheory = (1.0 / (1.0 - h.time));
  rhoTheory = rhoTheory * rhoTheory * rho0;

  double tol = 1.0e-10;
  AXOM_UNUSED_VAR(tol);
  SLIC_ASSERT_MSG(
    std::fabs(pp->rho(0) - rhoTheory) < tol
    && std::fabs(pp->rho(n - 1) - rhoTheory) < tol
    , "FAIL -- densities are not correct\n");

  SLIC_INFO(" *** PASS *** " );
}


TEST(slam_tinyHydro,test_03_gradAndForce)
{
  SLIC_INFO("\t******************************************************\n"
      << "\t     test gradient operator and force calculation.\n"
      << "\t*******************************************************");


  using IntVec = std::vector<int>;
  IntVec vals;
  vals. push_back(0);
  vals. push_back(1);
  vals. push_back(10);
  vals. push_back(11);
  vals. push_back(12);
  vals. push_back(13);

  int n = 10;
  int kz = n;
  int lz = n;
  PolygonMeshXY mesh(kz + 1,lz + 1);

  // Create a part
  std::vector<int> zones( mesh.numZones() );
  for(int i = 0; i< mesh.numZones(); ++i )
    zones[i] = i;

  Part p(zones);
  State s(mesh);

  SLIC_INFO("**setting initial density, energy, to make a linear pressure profile");
  double rho0 = 1.0;
  for (int i = 0; i < mesh.numZones(); i++)
  {
    p.setRho(i,rho0);
    VectorXY zPos = mesh.getZonePos(i);
    p.setE(i, zPos.x);
  }

  s.addPart(&p);
  double dt = 0.1;

  SLIC_INFO("**making hydro");
  Hydro h(&s);
  h.setBC("bottom", 0xdeadbeef, 0);   // fix bottom nodes in y


  SLIC_INFO("**Before we begin:");
  h.initialize();

  for(IntVec::iterator it = vals.begin(); it != vals.end(); ++it)
  {
    SLIC_INFO("\tdensity[ " << *it << "] = " <<  p.rho( *it)
                            << " energy[ " << *it << "] = " << p.e(*it));
  }

  SLIC_INFO("\te0, e1, e2, e3, e4 = "
      << p.e(0) << " "
      << p.e(1) << " "
      << p.e(2) << " "
      << p.e(3) << " "
      << p.e(4) << " "
  );

  SLIC_INFO("**Computing forces:");

  h.calcForce(s);

  for(IntVec::iterator it = vals.begin(); it != vals.end(); ++it)
  {
    VectorXY f = h.getForce(*it);
    SLIC_INFO("\tForce on node " << *it << " = " << f );
  }


  double tol = 1e-12;
  AXOM_UNUSED_VAR(tol);
  VectorXY f12 = h.getForce(12);
  VectorXY f13 = h.getForce(13);
  SLIC_ASSERT_MSG( std::fabs(f12.x - f13.x) < tol
      && std::fabs(f12.y) < tol
      && std::fabs(f13.y) < tol,
      "force calculation FAILS --"
      << "  f12: " << f12
      << ", f13: " << f13
      << ", diff in x: " << std::fabs(f12.x - f13.x)
      << ", tolerance: " << tol
  );

  SLIC_INFO("**Testing acceleration:");

  h.Cq = 0.0;   // turn off Q
  h.Cl = 0.0;   // turn off Q
  SLIC_INFO("\ttaking a step..");
  h.step(dt);


  for(IntVec::iterator it = vals.begin(); it != vals.end(); ++it)
  {
    VectorXY f = s.u(*it);
    SLIC_INFO("\tNew velocity on node " << *it << ": " << f << ")");
  }

  VectorXY u0 = s.u(0);
  VectorXY u11 = s.u(11);
  VectorXY u12 = s.u(12);
  SLIC_ASSERT_MSG( std::fabs(u0.x - u12.x) < tol
      && std::fabs(u11.y) < tol
      && std::fabs(u12.y) < tol,
      "acceleration/velocity calculation FAILS");


  SLIC_INFO(" *** PASS ***");

  // Deal with unused variables
  AXOM_UNUSED_VAR(f12);
  AXOM_UNUSED_VAR(f13);
  AXOM_UNUSED_VAR(u0);
  AXOM_UNUSED_VAR(u11);
  AXOM_UNUSED_VAR(u12);
}



TEST(slam_tinyHydro,test_04_BC)
{
  SLIC_INFO("\t******************************************************\n"
      << "\t           testing boundary conditions.\n"
      << "\t*******************************************************");

  int n = 10;
  int kz = n;
  int lz = n;
  PolygonMeshXY mesh(kz + 1,lz + 1);

  // Create a part
  std::vector<int> zones( mesh.numZones() );

  for(int i = 0; i< mesh.numZones(); ++i )
    zones[i] = i;

  Part p(zones);
  State s(mesh);

  SLIC_INFO("**setting initial density, energy, to make a linear pressure profile");
  double rho0 = 1.0;
  for (int i = 0; i < mesh.numZones(); i++)
  {
    p.setRho(i,rho0);
    VectorXY zPos = mesh.getZonePos(i);
    p.setE(i, zPos.x);
  }

  s.addPart(&p);
  double dt = 0.1;

  SLIC_INFO("**making hydro");
  Hydro h(&s);


  SLIC_INFO("**setting boundary conditions");
  h.setBC("bottom", 0xdeadbeef, 0);   // fix bottom nodes in y
  SLIC_INFO("\tbottom free in x, y is fixed");

  h.setBC("top", 0xdeadbeef, 0);
  SLIC_INFO("\ttop free in x, y is fixed");

  //h.setBC("left", 0,0); # leave unset to test free condition default
  SLIC_INFO("\tleft free in x, free in y");

  h.setBC("right", -1,0);
  SLIC_INFO("\tright -1 in x, y is fixed");

  h.initialize();


  SLIC_INFO("**stepping");
  h.step(dt);

  VectorXY u0 = s.u(0);
  VectorXY u11 = s.u(11);
  VectorXY u21 = s.u(21);
  VectorXY u120 = s.u(120);

  VectorXY n0 = mesh.getPos(0);
  VectorXY n11 = mesh.getPos(11);
  VectorXY n21 = mesh.getPos(21);
  VectorXY n120 = mesh.getPos(120);


  SLIC_INFO(" Bottom left -- u: " << u0   << " -- n: " << n0);
  SLIC_INFO(" Left side -- u: " << u11  << " -- n: " << n11);
  SLIC_INFO(" Right side -- u: " << u21  << " -- n: " << n21);
  SLIC_INFO(" Top right -- u: " << u120 << " -- n: " << n120);


  double tol = 1e-21;
  AXOM_UNUSED_VAR(tol);

  SLIC_ASSERT_MSG( std::fabs(u0.y) < tol
      && std::fabs(u11.x) > tol
      && std::fabs(u21.x + 1.0) < tol
      && std::fabs(u120.x + 1.0) < tol
      && std::fabs(u120.y ) < tol
      , "BC test FAILS ");

  SLIC_INFO(" *** PASS *** " );
}



TEST(slam_tinyHydro,test_05_newDT_Noh)
{
  SLIC_INFO("\t******************************************************\n"
      << "\t            testing newDT method -- Noh.\n"
      << "\t*******************************************************");

  int n = 10;
  int kz = n;
  int lz = n;
  PolygonMeshXY mesh(kz + 1,lz + 1);

  // Create a part
  std::vector<int> zones( mesh.numZones() );

  for(int i = 0; i< mesh.numZones(); ++i )
    zones[i] = i;

  Part p(zones);
  State s(mesh);

  SLIC_INFO("**setting initial density, energy, velocity for Noh problem");
  double rho0 = 1.0;
  for (int i = 0; i < mesh.numZones(); i++)
  {
    p.setRho(i,rho0);
    //VectorXY zPos = mesh.getZonePos(i);
    p.setE(i, 0.);
  }
  for (int i = 0; i < mesh.numNodes(); i++)
  {
    s.setU(i,VectorXY(-1.,0.));
  }

  s.addPart(&p);


  SLIC_INFO("**making hydro");
  Hydro h(&s);

  for (int i = 0; i < h.numBCnodes(3); i++)     //set left nodes (left=3) to zero velocity
  {
    h.getState()->setU(h.bcNode(3,i),VectorXY(0,0));
  }

  SLIC_INFO("**Set boundary conditions!");
  h.setBC("bottom", 0xdeadbeef, 0);
  h.setBC("top",    0xdeadbeef, 0);
  h.setBC("left",   0,          0xdeadbeef);
  h.setBC("right",  -1,         0xdeadbeef);


  SLIC_INFO("**Printing positions of boundary condition nodes:");
  for(int j = 0; j< 4; ++j)
  {
    switch (j)
    {
    case 0:
      SLIC_INFO("bottom boundary nodes");
      break;
    case 1:
      SLIC_INFO("right boundary nodes");
      break;
    case 2:
      SLIC_INFO("top boundary nodes");
      break;
    case 3:
      SLIC_INFO("left boundary nodes");
      break;
    }
    for (int i = 0; i < h.numBCnodes(j); ++i)       //set left nodes (left=3) to zero velocity
    {
      int nodeNum = h.bcNode(j,i);
      SLIC_INFO("\t" << nodeNum << " -- " <<   s.u(nodeNum) );

    }
  }

  SLIC_INFO("");



  h.initialize();

  SLIC_INFO("**Get DT!");
  double dt = h.newDT();
  SLIC_INFO("\tnewDT = " << dt );

  double tol = 1.0e-16;
  AXOM_UNUSED_VAR( tol);

  double expDT = 0.1 * h.cfl;
  AXOM_UNUSED_VAR( expDT);

  SLIC_ASSERT_MSG( std::fabs(dt - expDT) < tol,
      " newDT calculation FAILS -- expected dt = " << expDT << " but got " << dt << " instead, leaving " << dt - expDT);

  SLIC_INFO(" *** PASS *** " );
}



TEST(slam_tinyHydro,test_05_newDT_Sedov)
{
  SLIC_INFO("\t******************************************************\n"
      << "\t              testing newDT method -- Sedov.\n"
      << "\t*******************************************************");

  int n = 10;
  int kz = n;
  int lz = n;
  PolygonMeshXY mesh(kz + 1,lz + 1);

  // Create a part
  std::vector<int> zones( mesh.numZones() );

  for(int i = 0; i< mesh.numZones(); ++i )
    zones[i] = i;

  Part p(zones);
  State s(mesh);

  SLIC_INFO("**setting initial density, energy, velocity for Sedov problem");
  double rho0 = 1.0;
  for (int i = 0; i < mesh.numZones(); i++)
  {
    p.setRho(i,rho0);
    //kz = i%n;
    //lz = i / n;
    p.setE(i, 0.);
  }

  // set energy spike at the origin
  double zonemass = rho0 * mesh.zoneVol(0);
  double E = 4.48333289847;
  p.setE(0, E / (4.0 * zonemass)); // puts shock at r=0.8 at t=0.3, in XY coords


  s.addPart(&p);


  SLIC_INFO("**making hydro");
  Hydro h(&s);

  for (int i = 0; i < mesh.numNodes(); i++)
    s.setU(i,VectorXY(0,0));

  SLIC_INFO("**set boundary conditions");
  h.setBC("bottom", 0xdeadbeef, 0);
  h.setBC("top",    0xdeadbeef, 0xdeadbeef);
  h.setBC("left",   0,          0xdeadbeef);
  h.setBC("right",  0xdeadbeef, 0xdeadbeef);

  h.initialize();

  SLIC_INFO("**Get DT!");
  double dt = h.newDT();
  SLIC_INFO("\tnewDT = " << dt );

  double L = 0.1;
  double cfl = 0.7;
  double cs = sqrt(10 * E / (4 * zonemass * 9));
  double theoryDT = cfl * L / cs;
  AXOM_UNUSED_VAR( theoryDT);

  double tol = 1.0e-16;
  AXOM_UNUSED_VAR( tol);

  SLIC_ASSERT_MSG( std::fabs(dt - theoryDT) < tol,
      " newDT calculation FAILS -- expected dt = " << theoryDT << " but code got " << dt << ". Diff:" <<  dt - theoryDT);

  SLIC_INFO(" *** PASS ***");
}


TEST(slam_tinyHydro,test_06_PdV_work)
{
  SLIC_INFO("\t******************************************************\n"
      << "\t              testing PdV work.\n"
      << "\t*******************************************************");

  int steps = 100;

  int n = 1;
  int kz = n;
  int lz = n;
  PolygonMeshXY mesh(kz + 1,lz + 1);

  // Create a part
  std::vector<int> zones( mesh.numZones() );

  for(int i = 0; i< mesh.numZones(); ++i )
    zones[i] = i;

  Part p(zones);
  State s(mesh);

  SLIC_INFO("**setting initial density, energy, velocity for 1x1 compression problem");
  double rho0 = 1.0;
  double e0 = 1.0;
  double u0 = -0.1;

  for (int i = 0; i < mesh.numZones(); i++)
  {
    p.setRho(i,rho0);
    p.setE(i, e0);
  }

  s.addPart(&p);

  // right nodes compress at speed u0
  s.setU( 1,  VectorXY(u0,0));
  s.setU( 3,  VectorXY(u0,0));
  // set left nodes to zero velocity
  s.setU( 0,  VectorXY(0,0));
  s.setU( 2,  VectorXY(0,0));



  SLIC_INFO("**making hydro");
  Hydro h(&s);

  SLIC_INFO("**Set boundary conditions!");
  h.setBC("bottom", 0xdeadbeef, 0);
  h.setBC("top",    0xdeadbeef, 0);
  h.setBC("left",   0,          0xdeadbeef);
  h.setBC("right",  u0,         0xdeadbeef);


  SLIC_INFO("**Turn off Q!");
  h.Cq = 0.;
  h.Cl = 0.;

  SLIC_INFO("**Initialize");
  h.initialize();

  double dt = 0.01;

  SLIC_INFO("**take " << steps << " steps with dt = " << dt );
  for (int i = 0; i < mesh.numZones(); i++)
    h.step(dt);

  SLIC_INFO("**done stepping. ");

  double tol = 1.0e-6;
  AXOM_UNUSED_VAR(tol);

  double theoryRho = 1.0 / (1.0 + h.time * u0);
  double rhoCode = h.getState()->getPart(0)->rho(0);
  AXOM_UNUSED_VAR(rhoCode);

  SLIC_ASSERT_MSG( std::fabs(theoryRho - rhoCode) < tol,
      "density calculation FAILS -- rhoCode = " << rhoCode << " but should be " << theoryRho );


  double theoryE = std::pow(theoryRho / rho0, 2.0 / 3.0); // e = e0*(rho/rho0)**(gamma-1)
  double eCode = h.getState()->getPart(0)->e(0);

  SLIC_INFO("\t h.getState().getPart(0).e(0) " << eCode );
  SLIC_INFO("\t theory e " << theoryE );

  SLIC_ASSERT_MSG( std::fabs(theoryE - eCode) < tol,
      "PdV calculation FAILS");



  SLIC_INFO(" *** PASS *** ");
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
