// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet.hpp"

#include <iostream>
#include <unordered_map>
#include "axom/slic/core/SimpleLogger.hpp"

using axom::inlet::FunctionType;
using axom::inlet::Inlet;
using axom::inlet::LuaReader;

namespace inlet = axom::inlet;
// _inlet_userdef_simple_start
struct Car
{
  std::string make;
  std::string color;
  int seats;
  int horsepower;

  // A function can be used to define a schema for a particular struct
  // For convenience, it can be implemented as a static method of the struct
  static void defineSchema(inlet::Container& car_schema)
  {
    car_schema.addString("make", "Make of car");
    car_schema.addString("color", "Color of car").defaultValue("red");
    car_schema.addInt("seats", "Number of seats").range(2, 7);
    car_schema.addInt("horsepower", "Amount of horsepower");
  }
};
// _inlet_userdef_simple_end

// Additionally, each class should specialize this struct as follows
// in the global namespace so that Inlet can access it
// _inlet_userdef_simple_frominlet_start
template <>
struct FromInlet<Car>
{
  Car operator()(const inlet::Container& input_data)
  {
    Car result;
    result.make = input_data["make"];
    result.color = input_data["color"];
    result.seats = input_data["seats"];
    result.horsepower = input_data["horsepower"];
    return result;
  }
};
// _inlet_userdef_simple_frominlet_end

const std::string input = R"(
  -- A single car
  car = {
      make = "BestCompany",
      seats = 2,
      horsepower = 200,
      color = "blue"
  }

  -- A collection of cars
  fleet = {
    {
      make = "Globex Corp",
      seats = 3,
      horsepower = 155,
      color = "green"
    },
    {
      make = "Initech",
      seats = 4,
      horsepower = 370
    },
    {
      make = "Cyberdyne",
      seats = 2,
      horsepower = 101,
      color = "silver"
    }
  }
)";

int main()
{
  // Inlet requires a SLIC logger to be initialized to output runtime information
  axom::slic::SimpleLogger logger;

  auto lr = std::make_unique<LuaReader>();
  lr->parseString(input);
  Inlet inlet(std::move(lr));

  // _inlet_userdef_simple_usage_start
  // Create a container off the global container for the car object
  // then define its schema
  auto& car_schema = inlet.addStruct("car", "Vehicle description");
  Car::defineSchema(car_schema);
  // _inlet_userdef_simple_usage_end

  // _inlet_userdef_collection_usage_start
  // Create a fleet of cars with the same Car::defineSchema
  auto& fleet_schema = inlet.addStructArray("fleet", "A collection of cars");
  Car::defineSchema(fleet_schema);
  // _inlet_userdef_collection_usage_end

  if(!inlet.verify())
  {
    SLIC_ERROR("Inlet failed to verify against provided schema");
  }

  // Extract the car object
  Car car = inlet["car"].get<Car>();
  std::cout << "Info on single car:" << std::endl;
  std::cout << "make = " << car.make << std::endl;
  std::cout << "color = " << car.color << std::endl;
  std::cout << "seats = " << car.seats << std::endl;
  std::cout << "horsepower = " << car.horsepower << std::endl;

  // Extract the cars in the fleet
  auto fleet = inlet["fleet"].get<std::vector<Car>>();
  std::cout << "\nInfo on cars in fleet:" << std::endl;
  for(const Car& car : fleet)
  {
    std::cout << "make = " << car.make << std::endl;
    std::cout << "color = " << car.color << std::endl;
    std::cout << "seats = " << car.seats << std::endl;
    std::cout << "horsepower = " << car.horsepower << std::endl;
  }
}
