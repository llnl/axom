// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_ERROR_HPP_
#define AXOM_KLEE_ERROR_HPP_

#include "axom/inlet/inlet_utils.hpp"

#include <exception>
#include <vector>

namespace axom
{
namespace klee
{
/**
 * Describes an error that occurred while parsing a Klee file.
 */
class KleeError : public std::exception
{
public:
  /**
   * Create a KleeError from a single verification error.
   * @param error the VerificationError describing the failure
   */
  explicit KleeError(const inlet::VerificationError &error);

  /**
   * Create a KleeError from a vector of verification errors. There must
   * be at least one error provided.
   * @param errors the list VerificationError describing the failures. Must
   * have at least one.
   */
  explicit KleeError(const std::vector<inlet::VerificationError> &errors);

  /**
   * A description of the first error.
   * @return the message of the first error
   */
  const char *what() const noexcept override;

  /**
   * Get the list of all the errors.
   * @return all the errors which caused this exception
   */
  const std::vector<inlet::VerificationError> &getErrors() const { return m_errors; }

private:
  std::vector<inlet::VerificationError> m_errors;
};

}  // namespace klee
}  // namespace axom

#endif  // AXOM_KLEE_ERROR_HPP_
