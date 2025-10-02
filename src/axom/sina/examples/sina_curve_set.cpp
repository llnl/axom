// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

#include <iostream>
#include <cmath>
#include <vector>

struct BounceData
{
  std::vector<double> time;
  std::vector<double> xPosition;
  std::vector<double> yPosition;
  std::vector<double> xVelocity;
  std::vector<double> yVelocity;
};

BounceData generateBounceData(double initialY,
                              double initialXVelocity,
                              double coefficientOfRestitution,
                              double timeStep,
                              double simulationTime,
                              double airResistanceCoefficient)
{
  BounceData data;
  double time = 0.0;
  double yPosition = initialY;
  double xPosition = 0.0;
  double xVelocity = initialXVelocity;
  double yVelocity = 0.0;

  while(time < simulationTime)
  {
    // Update position and velocity based on time step
    xPosition += xVelocity * timeStep;
    yVelocity += -9.81 * timeStep;  // Acceleration due to gravity

    // Apply air resistance on x-velocity
    xVelocity -= airResistanceCoefficient * xVelocity * timeStep;

    // Update height (y position)
    yPosition += yVelocity * timeStep;

    // Check for bounce (when y position becomes negative)
    if(yPosition < 0)
    {
      yPosition = 0;                           // Set y position to 0 at the ground
      yVelocity *= -coefficientOfRestitution;  // Reverse y velocity with energy loss
    }

    // Add data points for this time step
    data.time.push_back(time);
    data.xPosition.push_back(xPosition);
    data.yPosition.push_back(yPosition);
    data.xVelocity.push_back(xVelocity);
    data.yVelocity.push_back(yVelocity);

    time += timeStep;
  }

  return data;
}

void addCurveSet(axom::sina::Record &record, BounceData bounceData, std::string curveName)
{
  // Create the curve set object
  axom::sina::CurveSet bounceCurveSet {curveName};

  // Add the independent variable
  axom::sina::Curve timeCurve {"time", bounceData.time};
  timeCurve.setUnits("seconds");
  bounceCurveSet.addIndependentCurve(timeCurve);

  // Add the dependent variables
  bounceCurveSet.addDependentCurve(axom::sina::Curve {"x_position", bounceData.xPosition});
  bounceCurveSet.addDependentCurve(axom::sina::Curve {"y_position", bounceData.yPosition});
  axom::sina::Curve xVelCurve {"x_velocity", bounceData.xVelocity};
  xVelCurve.setUnits("m/s");
  bounceCurveSet.addDependentCurve(xVelCurve);
  axom::sina::Curve yVelCurve {"y_velocity", bounceData.yVelocity};
  yVelCurve.setUnits("m/s");
  bounceCurveSet.addDependentCurve(yVelCurve);

  // Add the curve set to the record
  record.add(bounceCurveSet);
}

int main()
{
  double initialY = 10.0;
  double initialXVelocity = 2.0;
  double coefficientOfRestitution = 0.7;
  double timeStep = 0.2;
  double simulationTime = 5.0;
  double airResistanceCoefficient = 0.01;

  BounceData bounceData = generateBounceData(initialY,
                                             initialXVelocity,
                                             coefficientOfRestitution,
                                             timeStep,
                                             simulationTime,
                                             airResistanceCoefficient);

  axom::sina::Document doc;

  axom::sina::ID id {"ball_bounce_run", axom::sina::IDType::Global};
  std::unique_ptr<axom::sina::Record> study {new axom::sina::Record {id, "ball bounce study"}};

  addCurveSet(*study, bounceData, "ball_bounce");
  doc.add(std::move(study));
  axom::sina::saveDocument(doc, "ball_bounce.json");

  return 0;
}