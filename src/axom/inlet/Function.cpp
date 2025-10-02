// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Function.hpp"

namespace axom
{
namespace inlet
{
std::string Function::name() const
{
  return utilities::string::removePrefix(m_sidreRootGroup->getPathName(),
                                         m_sidreGroup->getPathName());
}

Function& Function::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Function specific Sidre Datastore Group not set");
  setFlag(*m_sidreGroup, *m_sidreRootGroup, detail::REQUIRED_FLAG, isRequired);
  return *this;
}

bool Function::isRequired() const
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Function specific Sidre Datastore Group not set");
  return checkFlag(*m_sidreGroup, *m_sidreRootGroup, detail::REQUIRED_FLAG);
}

Function& Function::registerVerifier(Verifier lambda)
{
  SLIC_WARNING_IF(m_verifier,
                  fmt::format("[Inlet] Verifier for Function "
                              "already set: {0}",
                              name()));
  m_verifier = lambda;
  return *this;
}

bool Function::verify(std::vector<VerificationError>* errors) const
{
  const bool this_function_exists = static_cast<bool>(m_func);
  // If this function was required, make sure something was defined in it
  bool verified = verifyRequired(*m_sidreGroup, this_function_exists, "Function", errors);
  // Verify this Function if a lambda was configured
  if(this_function_exists && m_verifier && !m_verifier(*this, errors))
  {
    verified = false;
    const std::string msg =
      fmt::format("[Inlet] Function failed verification: {0}", m_sidreGroup->getPathName());
    INLET_VERIFICATION_WARNING(m_sidreGroup->getPathName(), msg, errors);
  }

  return verified;
}

}  // end namespace inlet
}  // end namespace axom
