! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

module slic_interface_test
  use iso_c_binding
  use fruit
  use axom_slic
  implicit none

contains
!------------------------------------------------------------------------------

  subroutine initialize_finalize
    call assert_false(slic_is_initialized(), "slic initialized() not called")

    call slic_initialize()
    call assert_true(slic_is_initialized(), "slic initialized() is called")

    call slic_finalize()
    call assert_false(slic_is_initialized(), "slic finalize() is called")

    ! Check that we can initialize and finalize slic a second time
    call slic_initialize()
    call assert_true(slic_is_initialized(), "slic initialized() again")

    call slic_finalize()
    call assert_false(slic_is_initialized(), "slic finalize() again")

  end subroutine initialize_finalize

!------------------------------------------------------------------------------

  subroutine logging_level
    integer n

    call slic_initialize()
    
    call slic_set_logging_msg_level( message_error )
    call assert_true(slic_get_logging_msg_level() == message_error, "Error level")

    call slic_set_logging_msg_level( message_warning )
    call assert_true(slic_get_logging_msg_level() == message_warning, "Warning level")

    call slic_set_logging_msg_level( message_info )
    call assert_true(slic_get_logging_msg_level() == message_info, "Info level")

    call slic_set_logging_msg_level( message_debug )
    call assert_true(slic_get_logging_msg_level() == message_debug, "Debug level")

    call slic_finalize()
  end subroutine logging_level

!------------------------------------------------------------------------------
end module slic_interface_test
!------------------------------------------------------------------------------

program fortran_test
  use iso_c_binding
  use fruit
  use slic_interface_test
  implicit none
  logical ok

  call init_fruit

  call initialize_finalize
  call logging_level

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test
