program example
  use sina_functions
  use hdf5_config
  implicit none

  ! data types
  integer (KIND=4) :: int_val
  integer (KIND=8) :: long_val
  real :: real_val
  double precision :: double_val
  character :: char_val
  logical :: is_val
  integer :: i
  logical :: independent
  
  ! 1D real Array
  real, dimension(20) :: real_arr
  double precision, dimension(20) :: double_arr
  
  
  ! Strings
  character(:), allocatable :: fle_nme
  character(:), allocatable :: ofle_nme
  character(17) :: wrk_dir
  character(29) :: full_path
  character(36) :: ofull_path
  character(:), allocatable  :: rec_id
  character(:), allocatable :: mime_type
  character(:), allocatable :: tag
  character(:), allocatable :: units 
  character(20) :: json_fn
  character(20) :: hdf5_fn
  character(15) :: name
  character(25) :: curve
  
  ! 1D integer Array
  integer, dimension(20) :: int_arr
  integer (kind=8), dimension(20) :: long_arr
  
  int_val = 10
  long_val = 1000000000
  real_val = 1.234567
  double_val = 1./1.2345678901234567890123456789
  char_val = 'A'
  is_val = .false.
  
  do i = 1, 20
    real_arr(i) = i
    double_arr(i) = i*2.
    int_arr(i) = i*3
    long_arr(i) = i*4
  end do
  
  rec_id = make_cstring('my_rec_id')
  fle_nme = 'my_file.txt'
  ofle_nme = 'my_other_file.txt'
  wrk_dir = '/path/to/my/file/'
  full_path = make_cstring(wrk_dir//''//fle_nme)
  ofull_path = make_cstring(wrk_dir//''//ofle_nme)
  json_fn = make_cstring('sina_dump.json')
  if (use_hdf5) then
    hdf5_fn = make_cstring('sina_dump.hdf5')
  end if
  
  
  mime_type = make_cstring('')
  units = make_cstring('')
  tag = make_cstring('')
  
  print *,rec_id

  ! ========== USAGE ==========
  
  ! create sina record and document
  print *,'Creating the document'
  call create_document_and_record(rec_id)
  
  ! add file to sina record
  print *,'Adding a file to the Sina record'
  call sina_add_file(full_path, mime_type)
  mime_type = make_cstring('png')
  print *,'Adding another file (PNG) to the Sina record'
  call sina_add_file(ofull_path, mime_type)
  print *, "Adding int", int_val
  name = make_cstring('int')
  call sina_add(name, int_val, units, tag)
  print *, "Adding logical"
  name = make_cstring('logical')
  call sina_add(name, is_val, units, tag)
  print *, "Adding long"
  name = make_cstring('long')
  call sina_add(name, long_val, units, tag)
  print *, "Adding real"
  name = make_cstring('real')
  call sina_add(name, real_val, units, tag)
  print *, "Adding double"
  name = make_cstring('double')
  call sina_add(name, double_val, units, tag)
  print *, "Adding char"
  name = make_cstring('char')
  call sina_add(name, trim(char_val)//char(0), units, tag)
  units = make_cstring("kg")
  print *, "Adding int", int_val
  name = make_cstring('u_int')
  call sina_add(name, int_val, units, tag)
  print *, "Adding logical"
  name = make_cstring('u_logical')
  is_val = .true.
  call sina_add(name, is_val, units, tag)
  print *, "Adding long"
  name = make_cstring('u_long')
  call sina_add(name, long_val, units, tag)
  print *, "Adding real"
  name = make_cstring('u_real')
  call sina_add(name, real_val, units, tag)
  print *, "Adding double"
  name = make_cstring('u_double')
  call sina_add(name, double_val, units, tag)
  
  print *, "Adding double with tag"
  name = make_cstring('u_double_w_tag')
  tag = make_cstring('new_fancy_tag')
  call sina_add(name, double_val, units, tag)

  print *, "Adding char"
  name = make_cstring('u_char')
  call sina_add(name, trim(char_val)//char(0), units, tag)

  deallocate(tag)
 
  name = make_cstring('my_curveset')
  call sina_add_curveset(name)

  curve = make_cstring('my_indep_curve_double')
  independent = .TRUE.
  call sina_add_curve(name, curve, double_arr, size(double_arr), independent)
  curve = make_cstring('my_indep_curve_real')
  call sina_add_curve(name, curve, real_arr, size(real_arr), independent)
  curve = make_cstring('my_indep_curve_int')
  call sina_add_curve(name, curve, int_arr, size(int_arr), independent)
  curve = make_cstring('my_indep_curve_long')
  call sina_add_curve(name, curve, long_arr, size(long_arr), independent)
  curve = make_cstring('my_dep_curve_double')
  independent = .false.
  call sina_add_curve(name, curve, double_arr, size(double_arr), independent)
  curve = make_cstring('my_dep_curve_double_2')
  call sina_add_curve(name, curve, double_arr, size(double_arr), independent)
  curve = make_cstring('my_dep_curve_real')
  call sina_add_curve(name, curve, real_arr, size(real_arr), independent)
  curve = make_cstring('my_dep_curve_int')
  call sina_add_curve(name, curve, int_arr, size(int_arr), independent)
  curve = make_cstring('my_dep_curve_long')
  call sina_add_curve(name, curve, long_arr, size(long_arr), independent)
  ! write out the Sina Document
  print *,'Writing out the Sina Document'
  call write_sina_document(json_fn)
  if (use_hdf5) then
    call write_sina_document(hdf5_fn, 1)
  end if

  
contains
  function make_cstring(string) result(cstring)
    character(len=*), intent(in) :: string
    character(len=:), allocatable :: cstring
    cstring = trim(string) // char(0)
  end function make_cstring
  

end program example
