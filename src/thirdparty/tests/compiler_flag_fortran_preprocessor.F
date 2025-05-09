! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

!-------------------------------------------------------------------------------
! Simple test program to determine if the C preprocessor is
! automatically invoked for files with extension .F (note: capital F)
!
! Note: Our project policy is for Fortran source files that require the
!       C preprocessor to have the file extension *.F  (note: uppercase F)
!       and for files that do not require the C preprocessor
!       to have the extension *.f  (note: lowercase f).
!
! Note: If you are unable (or unwilling) to change the file's extension,
!       the following flags might be helpful:
!       * On gfortran (and iFort on linux), "-cpp" invokes the C preprocessor
!       * IBM's xlf does not appear to have a flag to directly invoke the C
!         preprocessor, but does have a way to remap a single file type
!         to invoke the C preprocessor on a compilation line.
!         If you have a *.f file, '-qsuffix=cpp=f' should work on xlf.
!         If you have a *.f90 file, '-qsuffix=cpp=f90' should work on xlf.
!-------------------------------------------------------------------------------
program preprocessor_tester

#if 0
  Random stuff that should be discarded if file is preprocessed
#endif

  print *, "Fortran's C preprocessor is working!"

end program
