// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_ADIAK_HPP
#define SINA_ADIAK_HPP

/*!
 ******************************************************************************
 *
 * \file AdiakWriter.hpp
 *
 * \brief   Header file for the Adiak Sina callback function.
 *
 ******************************************************************************
 */

#include "axom/config.hpp"
#ifdef AXOM_USE_ADIAK

  #include <string>
  #include <type_traits>

  #include "axom/sina/core/ConduitUtil.hpp"
  #include "axom/sina/core/Record.hpp"
  #include "axom/sina/core/Run.hpp"

extern "C" {
  #include "adiak_tool.h"
}

namespace axom
{
namespace sina
{

/**
 * \brief The callback function to pass to Adiak in order to write collected data to a Sina Record.
 *
 * To register this with Adiak, you should call adiak_register_cb, passing in
 * a pointer to the record as the last value. Your code should look something
 * like this:
 *
 * \code
 *     #include "axom/sina.hpp"
 *
 *     axom::sina::Record record{axom::sina::ID{"my_id", axom::sina::IDType::Local}, "my_record_type"};
 *     adiak_register_cb(1, adiak_category_all, axom::sina::adiakSinaCallback, 0, &record);
 * \endcode
 *
 * \attention Not everything that Sina can capture an be captured through the
 *            Adiak API. For example, there is currently no support in Adiak to capture
 *            anything like a CurveSet. As a result, to do that, you must hold on to
 *            the Record object passed here as the opaque value and manipulate it directly.
 **/
void adiakSinaCallback(const char *name,
                       adiak_category_t category,
                       const char *subcategory,
                       adiak_value_t *value,
                       adiak_datatype_t *t,
                       void *opaque_value);

}  // namespace sina
}  // namespace axom

#endif  // AXOM_USE_ADIAK

#endif  // SINA_ADIAK_HPP
