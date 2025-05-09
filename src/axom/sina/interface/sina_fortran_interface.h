// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

extern "C" char *Get_File_Extension(char *);
extern "C" void create_document_and_run_(char *);
extern "C" axom::sina::Record *Sina_Get_Run();
extern "C" void sina_add_file_to_record_(char *);
extern "C" void sina_add_file_with_mimetype_to_record_(char *, char *);
extern "C" void write_sina_document_protocol_(char *, int *);
extern "C" void write_sina_document_noprotocol_(char *);
extern "C" void sina_add_long_(char *, long long int *, char *, char *);
extern "C" void sina_add_int_(char *, int *, char *, char *);
extern "C" void sina_add_float_(char *, float *, char *, char *);
extern "C" void sina_add_double_(char *, double *, char *, char *);
extern "C" void sina_add_logical_(char *, bool *, char *, char *);
extern "C" void sina_add_string_(char *, char *, char *, char *);
extern "C" void sina_add_curveset_(char *);
extern "C" void sina_add_curve_double_(char *, char *, double *, int *, bool *);
extern "C" void sina_add_curve_float_(char *, char *, float *, int *, bool *);
extern "C" void sina_add_curve_int_(char *, char *, int *, int *, bool *);
extern "C" void sina_add_curve_long_(char *, char *, long long int *, int *, bool *);
