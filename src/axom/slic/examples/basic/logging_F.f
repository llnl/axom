! Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
! other Axom Project Developers. See the top-level LICENSE file for details.
!
! SPDX-License-Identifier: (BSD-3-Clause)

program slic_logging
  ! Use slic module
  use iso_c_binding
  use axom_slic
  implicit none

  type(SlicGenericOutputStream) stream

  ! initialize slic environment
  call slic_initialize()

  ! set the message level
  call slic_set_logging_msg_level( message_debug )

  call slic_disable_abort_on_error()

  ! register log stream
  stream = SlicGenericOutputStream("cout", "<MESSAGE>\n\t<TIMESTAMP>\n\tLEVEL=<LEVEL>\n\n")
  call slic_add_stream_to_all_msg_levels(stream)

  ! log messages
  call slic_log_message(message_debug,'Here is a debug message!')
  call slic_log_message(message_info, 'Here is an info mesage!')
  call slic_log_message(message_warning, 'Here is a warning!')
  call slic_log_message(message_error, 'Here is an error message!')

  ! finalize slic environment
  call slic_finalize()

end program slic_logging
