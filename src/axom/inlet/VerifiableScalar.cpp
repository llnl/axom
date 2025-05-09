// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
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
