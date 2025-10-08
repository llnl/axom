// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

extern "C" char *Get_File_Extension(char *);
extern "C" axom::sina::Record *Sina_Get_Run();
extern "C" void sina_add_file_to_record_(char *);
extern "C" void sina_add_file_with_mimetype_to_record_(char *, char *);
extern "C" void write_sina_document_protocol_(char *, int *);
extern "C" void write_sina_document_noprotocol_(char *);
extern "C" void sina_add_long_(char *, long long int *, char *, char *, char *);
extern "C" void sina_add_int_(char *, int *, char *, char *, char *);
extern "C" void sina_add_float_(char *, float *, char *, char *, char *);
extern "C" void sina_add_double_(char *, double *, char *, char *, char *);
extern "C" void sina_add_logical_(char *, bool *, char *, char *, char *);
extern "C" void sina_add_string_(char *, char *, char *, char *, char *);
extern "C" void sina_add_curveset_(char *, char *);
extern "C" void sina_add_curve_double_(char *, char *, double *, int *, bool *, char *);
extern "C" void sina_add_curve_float_(char *, char *, float *, int *, bool *, char *);
extern "C" void sina_add_curve_int_(char *, char *, int *, int *, bool *, char *);
extern "C" void sina_add_curve_long_(char *, char *, long long int *, int *, bool *, char *);
// Save/Output Functions
extern "C" int sina_save_document_fortran(void* , const char* , int);
extern "C" int sina_output_document_to_json_fortran(void* , const char*);
extern "C" int sina_output_document_to_hdf5_fortran(void* , const char*);
// Append Functions
extern "C" int sina_append_document_fortran(void* , const char* , int , int);
extern "C" int sina_append_document_to_json_fortran(void* , const char* , int);
extern "C" int sina_append_document_to_hdf5_fortran(void* , const char* , int);
// Curve Ordering Functions
extern "C" int sina_record_set_curve_order_fortran(void* , int);
extern "C" int sina_record_get_curve_order_fortran(void* , int*);
