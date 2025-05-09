// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 *  \file SidreTypes.h
 *
 *  \brief File containing C types and constants used in
 *         SiDRe toolkit component.
 *
 */

/*
 * Note: Use only C code in this file.
 *       It is part of the C wrapper.
 */

#ifndef SIDRETYPES_H
#define SIDRETYPES_H

// Axom includes
#include "axom/config.hpp"
#include "axom/sidre/core/SidreDataTypeIds.h"

// C includes
#include <stdint.h> /* for int64_t */

#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
typedef int64_t SIDRE_IndexType;
#else
typedef int32_t SIDRE_IndexType;
#endif

typedef short SIDRE_TypeID;
typedef int SIDRE_TypeIDint;

#define SIDRE_InvalidName NULL

#endif  // SIDRETYPES_H
