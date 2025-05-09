! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

program slic_f_logger
  ! Use slic module
  use iso_c_binding
  use axom_slic
  implicit none

  type(SlicGenericOutputStream) fileStream
  type(SlicGenericOutputStream) consoleStream
  integer n
  real(C_DOUBLE) :: u
  integer(C_INT) :: random_msg_level

  ! STEP 0: Initialize logger
  if (.not. slic_is_initialized()) then
    call slic_initialize()
  endif
  call slic_set_logging_msg_level( message_debug )
  call slic_disable_abort_on_error()

  ! STEP 1: Create log streams
  ! setup log stream for ALL messages, including ERROR and WARNING
  fileStream = SlicGenericOutputStream("flogger.dat", "<MESSAGE>\n\t<TIMESTAMP>\n\tLEVEL=<LEVEL>\n\n")

  consoleStream = SlicGenericOutputStream("cerr", "[<LEVEL>]: <MESSAGE>\n")

  ! STEP 2: add streams to logger
  call slic_add_stream_to_msg_level(fileStream, message_error)
  call slic_add_stream_to_msg_level(fileStream, message_warning)
  call slic_add_stream_to_all_msg_levels(consoleStream)

  ! STEP 3: Loop 10 times and generate random logging events
  do n = 1, 10
    call random_number(u)
    random_msg_level = FLOOR(4*u)
    call slic_log_message(random_msg_level,'Here is a random message!')
  enddo

  ! finalize logging, flush streams for file
  call slic_flush_streams()
  call slic_finalize()

end program slic_f_logger
