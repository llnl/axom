// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Attribute.cpp
 *
 * \brief   Implementation file for Attribute class.
 *
 ******************************************************************************
 */

// Associated header file
#include "Attribute.hpp"

// Other axom headers
//#include "axom/core/Types.hpp"
//#include "axom/slic/interface/slic.hpp"

// Sidre component headers
//#include "Buffer.hpp"
//#include "Group.hpp"
//#include "DataStore.hpp"

namespace axom
{
namespace sidre
{
/*
 *************************************************************************
 *
 * PRIVATE ctor for Attribute not associated with any data.
 *
 *************************************************************************
 */
Attribute::Attribute(const std::string& name) : m_name(name), m_index(InvalidIndex) { }

/*
 *************************************************************************
 *
 * PRIVATE dtor.
 *
 *************************************************************************
 */
Attribute::~Attribute() { }

} /* end namespace sidre */
} /* end namespace axom */
