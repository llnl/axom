// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


// Part class holds the material data for a single material.

#ifndef __Part__
#define __Part__

#include <vector>

#include "TinyHydroTypes.hpp"

namespace tinyHydro {

  class Part
  {
  public:
    Part()
        : zones(0)
          , gamma(0.)
    {}

    Part(const int * zoneList, int nzones, double userGamma = 5.0 / 3.0);
    Part(const std::vector<int>& zoneList, double userGamma = 5.0 / 3.0);

    // Copy
    Part( const Part & arg);
    // Assignment
    Part & operator =(const Part & rhs);
    // addition in place
    Part & operator +=(const Part & rhs);
    // scalar multiply in place
    Part & operator *=(const double s);

    // python accessors -- access Part's fields
    double          rho(int i)                const { return density[i]; }
    double          e(int i)                  const { return energyPerMass[i]; }
    void            setRho(int i, double val)       { density[i] = val; }
    void            setE(int i, double val)         { energyPerMass[i] = val; }

    // access Part's topology
    int             meshZone(int i)           const { return zones[i]; }
    int             numZones()                const { return zones.size(); }

    void            dumpPart();

  public:
    ZoneSubset zones;                       // Indirection set of zones in a part/region
    ZonalScalarField density;
    ZonalScalarField energyPerMass;
    ZonalScalarField volumeFraction;

    double gamma;

  };

} // end namespace tinyHydro

#endif
