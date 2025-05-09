// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"
#include "axom/slic.hpp"

using ValueTypeUnderlying = typename std::underlying_type<axom::sina::ValueType>::type;

void printType(axom::sina::Datum datum, std::string datumName, const std::string& errMsg)
{
  auto datumType = static_cast<ValueTypeUnderlying>(datum.getType());
  SLIC_ASSERT_MSG(static_cast<bool>(std::is_same<decltype(datumType), ValueTypeUnderlying>::value),
                  errMsg);
  AXOM_UNUSED_VAR(errMsg);

  std::cout << datumName << " type: " << datumType << std::endl;
}

int main(void)
{
  // Initialize slic
  axom::slic::initialize();

  // Define 3 different datums
  axom::sina::Datum myDatum {12.34};
  std::string value = "foobar";
  axom::sina::Datum myOtherDatum {value};
  std::vector<double> scalars = {1, 2, 20.0};
  axom::sina::Datum myArrayDatum {scalars};

  // Prints 1, corresponding to Scalar
  printType(myDatum,
            "myDatum",
            "myDatumType did not match the expected type 'Scalar' (numerically "
            "represented as 1).");

  // Prints 0, corresponding to String
  printType(myOtherDatum,
            "myOtherDatum",
            "myDatumType did not match the expected type 'String' (numerically "
            "represented as 0).");

  // Prints 3, corresponding to ScalarArray
  printType(myArrayDatum,
            "myArrayDatum",
            "myArrayDatum did not match the expected type 'ScalarArray' "
            "(numerically represented as 3).");

  // Finalize slic
  axom::slic::finalize();
}