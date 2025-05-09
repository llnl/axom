// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file YAMLReader.hpp
 *
 * \brief This file contains the class definition of the YAMLReader.
 *******************************************************************************
 */

#ifndef INLET_YAMLREADER_HPP
#define INLET_YAMLREADER_HPP

#include "axom/inlet/ConduitReader.hpp"

#include "conduit.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class YAMLReader
 *
 * \brief A Reader that is able to read variables from a YAML file.
 *
 * \see Reader
 *******************************************************************************
 */
class YAMLReader : public ConduitReader
{
public:
  YAMLReader() : ConduitReader("yaml") { }
};

}  // end namespace inlet
}  // end namespace axom

#endif
