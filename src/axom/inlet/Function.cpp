// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Function.hpp"

namespace axom
{
namespace inlet
{
namespace detail
{
void addSignatureToGroup(const FunctionTag ret_type,
                         const std::vector<FunctionTag>& arg_types,
                         axom::sidre::Group* group)
{
  group->createViewScalar("return_type", static_cast<int>(ret_type));
  auto args_view =
    group->createViewAndAllocate("function_arguments", sidre::INT_ID, arg_types.size());
  int* args_array = args_view->getArray();
  for(const auto arg_type : arg_types)
  {
    *args_array++ = static_cast<int>(arg_type);
  }
}

bool hasSignature(const FunctionTag ret_type,
                  const std::vector<FunctionTag>& arg_types,
                  const axom::sidre::Group* group)
{
  if(!group->hasView("return_type") || !group->hasView("function_arguments"))
  {
    return false;
  }
  const int stored_ret = group->getView("return_type")->getData();
  const auto args_view = group->getView("function_arguments");
  if(stored_ret != static_cast<int>(ret_type) ||
     args_view->getNumElements() != static_cast<axom::IndexType>(arg_types.size()))
  {
    return false;
  }
  const int* stored_args = args_view->getData();
  for(std::size_t i = 0; i < arg_types.size(); ++i)
  {
    if(stored_args[i] != static_cast<int>(arg_types[i]))
    {
      return false;
    }
  }
  return true;
}

}  // namespace detail

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
