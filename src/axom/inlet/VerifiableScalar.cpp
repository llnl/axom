// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/VerifiableScalar.hpp"

namespace axom
{
namespace inlet
{
VerifiableScalar& VerifiableScalar::registerVerifier(std::function<bool(const Field&)> verifier)
{
  return registerVerifier(
    [verifier](const Field& field, std::vector<VerificationError>*) { return verifier(field); });
}
}  // namespace inlet
}  // namespace axom
