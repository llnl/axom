// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

extern "C" void sina_set_default_record_type_(char *);
extern "C" void sina_create_record_(char *, char *, int, int);
extern "C" void sina_add_file_(char *, char *, char *, int, int, int);
extern "C" axom::sina::Record *Sina_Get_Record(char *);
extern "C" char *Get_File_Extension(char *);
extern "C" axom::sina::Record *Sina_Get_Run();
extern "C" void sina_add_file_to_record_(char *, char *, int, int);
extern "C" void sina_add_file_with_mimetype_to_record_(char *, char *, char *, int, int, int);
extern "C" void sina_write_document_noprotocol_nopreserve_nomerge(char *);
extern "C" void sina_write_document_protocol_nopreserve_nomerge(char *, int *);
extern "C" void sina_write_document_protocol_preserve_nomerge(char *, int *, int *);
extern "C" void sina_write_document_all_args(char *, int *, int *, int *);
// extern "C" void write_sina_document_noprotocol_(char *);
extern "C" void sina_add_long_(char *, long long int *, char *, char *, char *, int, int, int, int);
extern "C" void sina_add_int_(char *, int *, char *, char *, char *, int, int, int, int);
extern "C" void sina_add_float_(char *, float *, char *, char *, char *, int, int, int, int);
extern "C" void sina_add_double_(char *, double *, char *, char *, char *, int, int, int, int);
extern "C" void sina_add_logical_(char *, bool *, char *, char *, char *, int, int, int, int);
extern "C" void sina_add_string_(char *, char *, char *, char *, char *, int, int, int, int);
extern "C" void sina_add_curveset_(char *, char *, int, int);
extern "C" void sina_add_curve_double_(char *, char *, double *, int *, int *, char *);
extern "C" void sina_add_curve_float_(char *, char *, float *, int *, int *, char *);
extern "C" void sina_add_curve_int_(char *, char *, int *, int *, int *, char *);
extern "C" void sina_add_curve_long_(char *, char *, long long int *, int *, int *, char *);
// Curve Ordering Functions
extern "C" void sina_set_curves_order_(int *);
extern "C" void sina_set_record_curves_order_(char *, int *);
