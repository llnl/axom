// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file JSONReader.hpp
 *
 * \brief This file contains the class definition of the JSONReader.
 *******************************************************************************
 */

#ifndef INLET_JSONREADER_HPP
#define INLET_JSONREADER_HPP

#include "axom/inlet/ConduitReader.hpp"

#include "conduit.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class JSONReader
 *
 * \brief A Reader that is able to read variables from a JSON file.
 *
 * \see Reader
 *******************************************************************************
 */
class JSONReader : public ConduitReader
{
public:
  JSONReader() : ConduitReader("json") { }
};

}  // end namespace inlet
}  // end namespace axom

#endif
