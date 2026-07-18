program example
  use sina_functions
  use sina_hdf5_config
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
  character(:), allocatable  :: rec_id, rec2_id
  character(:), allocatable :: mime_type
  character(:), allocatable :: tag
  character(:), allocatable :: units 
  character(20) :: json_fn
  character(20) :: hdf5_fn
  character(15) :: name
  character(20) :: name2
  character(25) :: curve
  character(11) :: custom_type

  ! 1D integer Array
  integer, dimension(20) :: int_arr
  integer (kind=8), dimension(20) :: long_arr
  
  integer :: num_args
  character(len=100) :: arg
  integer :: status 

  
  int_val = 10
  long_val = 1000000000
  real_val = 1.234567
  double_val = 1./1.2345678901234567890123456789
  char_val = 'A'
  is_val = .false.


  ! Set default value
  arg = "sina_dump"
  
  ! Override if argument was passed
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
  end if

  do i = 1, 20
    real_arr(i) = i
    double_arr(i) = i*2.
    int_arr(i) = i*3
    long_arr(i) = i*4
  end do
  
  rec_id = make_cstring('my_rec_id')
  rec2_id = make_cstring('my_rec_2_id')
  fle_nme = 'my_file.txt'
  ofle_nme = 'my_other_file.txt'
  wrk_dir = '/path/to/my/file/'
  full_path = make_cstring(wrk_dir//''//fle_nme)
  ofull_path = make_cstring(wrk_dir//''//ofle_nme)
  json_fn = make_cstring(trim(arg) // '.json')
  if (use_hdf5) then
    hdf5_fn = make_cstring(trim(arg) // '.hdf5')
  end if
  
  
  mime_type = make_cstring('')
  units = make_cstring('')
  tag = make_cstring('')
  
  print *,rec_id, rec2_id

  ! ========== USAGE ==========
  
  ! create sina record and document
  print *,'Creating the document'
  ! changing default order for my_curveset to reverse alphabetical
  call sina_set_curves_order(3)
  call sina_create_record(rec_id)
  ! set this record curves order to alphabetical
  call sina_set_record_curves_order(rec_id, 2)
  print *,'Creating the document and second record'
  custom_type = make_cstring('custom_type')
  call sina_create_record(rec2_id, custom_type)
  
  ! add file to sina record
  print *,'Adding a file to the Sina record'
  call sina_add_file(full_path, mime_type)
  call sina_add_file(full_path, mime_type, rec2_id)
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
  call sina_add(name, double_val, units, tag, rec2_id)
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

  print *, "Adding char type"
  name = make_cstring('u_char')
  call sina_add(name, trim(char_val)//char(0), units, tag)

  deallocate(tag)
 
  name = make_cstring('my_curveset')
  name2 = make_cstring('my_other_curveset')
  call sina_add_curveset(name)
  call sina_add_curveset(name2, rec2_id)

  curve = make_cstring('my_indep_curve_double')
  independent = .TRUE.
  call sina_add_curve(name, curve, double_arr, size(double_arr), independent)
  call sina_add_curve(name2, curve, double_arr, size(double_arr), independent, rec2_id)
  curve = make_cstring('my_indep_curve_real')
  call sina_add_curve(name, curve, real_arr, size(real_arr), independent)
  call sina_add_curve(name2, curve, real_arr, size(real_arr), independent, rec2_id)
  curve = make_cstring('my_indep_curve_int')
  call sina_add_curve(name, curve, int_arr, size(int_arr), independent)
  call sina_add_curve(name2, curve, int_arr, size(int_arr), independent, rec2_id)
  curve = make_cstring('my_indep_curve_long')
  call sina_add_curve(name, curve, long_arr, size(long_arr), independent)
  call sina_add_curve(name2, curve, long_arr, size(long_arr), independent, rec2_id)
  curve = make_cstring('my_dep_curve_double')
  independent = .false.
  call sina_add_curve(name, curve, double_arr, size(double_arr), independent)
  call sina_add_curve(name2, curve, double_arr, size(double_arr), independent, rec2_id)
  curve = make_cstring('my_dep_curve_double_2')
  call sina_add_curve(name, curve, double_arr, size(double_arr), independent)
  curve = make_cstring('my_dep_curve_real')
  call sina_add_curve(name, curve, real_arr, size(real_arr), independent)
  curve = make_cstring('my_dep_curve_int')
  call sina_add_curve(name, curve, int_arr, size(int_arr), independent)
  curve = make_cstring('my_dep_curve_long')
  call sina_add_curve(name, curve, long_arr, size(long_arr), independent)
  ! write out the Sina Document
  print *,'Writing out the Sina Document as json, preserve records'
  if (use_hdf5) then
    call sina_write_document(json_fn, 0, 1)
    print *,'Writing out the Sina Document as hdf5, yank all records'
    call sina_write_document(hdf5_fn, 1)
  else
    call sina_write_document(json_fn)
  end if

  ! set default record type
  rec_id = make_cstring('fortran_test')
  call sina_set_default_record_type(rec_id)
  
  ! Let's add another record
  rec2_id = make_cstring('my_rec_3_id')
  call sina_create_record(rec2_id)
  curve = make_cstring('my_indep_curve_double')
  independent = .true.
  print*, 'adding curve to rec3 double', independent, curve
  call sina_add_curve(name2, curve, double_arr, size(double_arr), independent, rec2_id)
  
  ! And save the hdf5 only with autodetect
  print*, 'saving to', hdf5_fn
  call sina_write_document(hdf5_fn)
  call sina_add_curve(name2, curve, double_arr, size(double_arr), independent, rec2_id)
  call sina_write_document(hdf5_fn)

  ! ========== CLEANUP - Deallocate only allocatable strings ==========
  if (allocated(rec_id)) deallocate(rec_id)
  if (allocated(rec2_id)) deallocate(rec2_id)
  if (allocated(fle_nme)) deallocate(fle_nme)
  if (allocated(ofle_nme)) deallocate(ofle_nme)
  if (allocated(mime_type)) deallocate(mime_type)
  if (allocated(units)) deallocate(units)
  if (allocated(tag)) deallocate(tag)
  ! Note: full_path, ofull_path, custom_type are fixed-length, not allocatable

  
contains
function make_cstring2(fstr) result(cstr)
  use iso_c_binding
  character(*), intent(in) :: fstr
  character(len=len_trim(fstr)+1, kind=c_char) :: cstr

  if (len_trim(fstr) > 0) then
    cstr(1:len_trim(fstr)) = fstr(1:len_trim(fstr))
  endif
  cstr(len_trim(fstr)+1:len_trim(fstr)+1) = c_null_char
end function make_cstring2 

end program example
