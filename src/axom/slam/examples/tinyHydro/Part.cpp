// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Part.hpp"

#include <sstream>
#include <string>
#include <cstring>
#include <algorithm>  // for std::copy

#include "axom/slic.hpp"

namespace tinyHydro {


//----------------------------------------------
  Part::Part(const int * zoneList, int nzones, double gamma)
      : gamma(gamma)
  {
    SLIC_DEBUG("in Part c'tor (from ptr)");

    IndexBuffer& zoneBuf = DataRegistry::setRegistry.addNamelessBuffer(nzones);
    std::copy(zoneList, zoneList + nzones, zoneBuf.begin());
    zones            = ZoneSubset::SetBuilder().size(zoneBuf.size()).data(&zoneBuf);

    density          = ZonalScalarField(&zones);
    energyPerMass    = ZonalScalarField(&zones);
    volumeFraction   = ZonalScalarField(&zones);
  }
//----------------------------------------------
  Part::Part(const std::vector<int>& zoneList, double gamma)
      : gamma(gamma)
  {
    SLIC_DEBUG("in Part c'tor (from vector)");

    IndexBuffer& zoneBuf = DataRegistry::setRegistry.addNamelessBuffer(zoneList.size());
    std::copy(zoneList.begin(), zoneList.end(), zoneBuf.begin());
    zones            = ZoneSubset::SetBuilder().size(zoneBuf.size()).data(&zoneBuf);

    density          = ZonalScalarField(&zones);
    energyPerMass    = ZonalScalarField(&zones);
    volumeFraction   = ZonalScalarField(&zones);
  }

//----------------------------------------------
// Copy
  Part::Part(const Part & arg)
      :    gamma(arg.gamma)
  {
    SLIC_DEBUG("in Part copy .ctor");

    IndexBuffer& zoneBuf = DataRegistry::setRegistry.addNamelessBuffer(arg.zones.size());
    std::copy(&arg.zones[0], &arg.zones[0] + arg.zones.size(), zoneBuf.begin());
    zones            = ZoneSubset::SetBuilder().size(zoneBuf.size()).data(&zoneBuf);

    // copy field data
    density         = ZonalScalarField(arg.density);
    energyPerMass   = ZonalScalarField(arg.energyPerMass);
    volumeFraction  = ZonalScalarField(arg.volumeFraction);
  }

//----------------------------------------------
// Assignment
  Part & Part::operator=(const Part & rhs)
  {
    //SLIC_DEBUG("in Part assignment");

    // check for self-assignment
    if (this == &rhs)
      return *this;

    // copy over the gamma value
    gamma = rhs.gamma;


    IndexBuffer& zoneBuf = DataRegistry::setRegistry.addNamelessBuffer(rhs.zones.size());
    std::copy(&rhs.zones[0], &rhs.zones[0] + rhs.zones.size(), zoneBuf.begin());
    zones            = ZoneSubset::SetBuilder().size(zoneBuf.size()).data(&zoneBuf);

    // copy field data
    density          = ZonalScalarField(rhs.density);
    energyPerMass    = ZonalScalarField(rhs.energyPerMass);
    volumeFraction   = ZonalScalarField(rhs.volumeFraction);

    return *this;
  }

//----------------------------------------------
// addition in place
  Part & Part::operator+=(const Part & rhs)
  {
    // printf("in Part::operator+= \n");
    SLIC_ASSERT(numZones() == rhs.numZones() );

    const int nZones = numZones();
    for (int i = 0; i < nZones; i++)
    {
      density[i]        += rhs.density[i];
      energyPerMass[i]  += rhs.energyPerMass[i];
      volumeFraction[i] += rhs.volumeFraction[i];
    }

    return *this;
  }
//----------------------------------------------
// scalar multiply in place
  Part & Part::operator*=(const double s)
  {
    // printf("in Part::operator*= \n");

    const int nZones = numZones();

    for (int i = 0; i < nZones; i++)
    {
      density[i]        *= s;
      energyPerMass[i]  *= s;
      volumeFraction[i] *= s;
    }

    return *this;
  }
//----------------------------------------------
//----------------------------------------------


  void Part::dumpPart()
  {
    SLIC_INFO("Part has " << numZones() << " zones" );

    std::stringstream zoneStr;
    zoneStr << "\nZones:";
    for(int i = 0; i< numZones(); ++i)
    {
      zoneStr << "\n\t Zone " << i << " -- idx " << zones[i];
    }
    SLIC_INFO( zoneStr.str() );
  }

} // end namespace tinyHydro
