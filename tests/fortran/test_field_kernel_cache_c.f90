program test_field_kernel_cache_c
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_loc, c_ptr
  use bem_field_kernel_c, only: beach_kernel_build, beach_kernel_create, beach_kernel_destroy, &
                                beach_kernel_get_periodic_cache_info, beach_kernel_invalid_argument, &
                                beach_kernel_ok, beach_kernel_set_periodic_cache
  use bem_kinds, only: i32
  use test_support, only: assert_equal_i32, assert_true, delete_file_if_exists, remove_empty_directory, &
                          test_begin, test_summary
  implicit none

  type(c_ptr) :: cold_handle, warm_handle
  integer(c_int) :: status
  integer(c_int), target :: cold_hit, cold_count, cold_fingerprint_len, cold_path_len
  integer(c_int), target :: warm_hit, warm_count, warm_fingerprint_len, warm_path_len
  integer(c_int), target :: periodic_axes(2)
  real(c_double), target :: periodic_len(2), box_min(3), box_max(3)
  real(c_double), target :: v0(3, 4), v1(3, 4), v2(3, 4)
  character(kind=c_char), allocatable, target :: cache_dir_bytes(:)
  character(kind=c_char), target :: cold_fingerprint(17), warm_fingerprint(17)
  character(kind=c_char), target :: cold_path(513), warm_path(513)
  character(kind=c_char), target :: short_fingerprint(1), short_path(1)
  character(len=256) :: cache_dir
  character(len=512) :: cold_path_text, warm_path_text
  integer :: clock_count

  call system_clock(count=clock_count)
  write (cache_dir, '(a,i0)') 'build/test_field_kernel_cache_c_', clock_count
  call to_c_bytes(trim(cache_dir), cache_dir_bytes)
  periodic_axes = [1_c_int, 2_c_int]
  periodic_len = [2.0d0, 2.0d0]
  box_min = [0.0d0, 0.0d0, -1.0d0]
  box_max = [2.0d0, 2.0d0, 1.0d0]
  v0(:, 1) = [0.30d0, 0.40d0, -0.20d0]
  v0(:, 2) = [1.50d0, 0.50d0, 0.16d0]
  v0(:, 3) = [0.60d0, 1.60d0, 0.24d0]
  v0(:, 4) = [1.64d0, 1.44d0, -0.12d0]
  v1 = v0
  v2 = v0
  v1(1, :) = v1(1, :) + 0.05d0
  v2(2, :) = v2(2, :) + 0.05d0

  call test_begin('field_kernel_cache_c_rejects_invalid_configuration')
  status = beach_kernel_create(cold_handle)
  call assert_equal_i32(status, beach_kernel_ok, 'invalid fixture create')
  status = beach_kernel_set_periodic_cache(cold_handle, c_loc(cache_dir_bytes), 0_c_int, 2.5d-9)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'empty cache path')
  status = beach_kernel_set_periodic_cache( &
           cold_handle, c_loc(cache_dir_bytes), int(size(cache_dir_bytes), c_int), 0.0d0 &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'zero generation tolerance')
  status = beach_kernel_destroy(cold_handle)
  call assert_equal_i32(status, beach_kernel_ok, 'invalid fixture destroy')

  call test_begin('field_kernel_cache_c_cold_then_warm')
  call build_cached_handle(cold_handle)
  short_fingerprint = achar(88, kind=c_char)
  short_path = achar(89, kind=c_char)
  cold_fingerprint_len = -1_c_int
  cold_path_len = -1_c_int
  status = beach_kernel_get_periodic_cache_info( &
           cold_handle, c_loc(cold_hit), c_loc(cold_count), c_loc(short_fingerprint), 1_c_int, &
           c_loc(cold_fingerprint_len), c_loc(short_path), 1_c_int, c_loc(cold_path_len) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cold short diagnostics status')
  call assert_true(cold_fingerprint_len > 0_c_int, 'short getter fingerprint required length')
  call assert_true(cold_path_len > 0_c_int, 'short getter path required length')
  call assert_true(short_fingerprint(1) == achar(88, kind=c_char), 'short getter must not truncate fingerprint')
  call assert_true(short_path(1) == achar(89, kind=c_char), 'short getter must not truncate path')
  status = beach_kernel_get_periodic_cache_info( &
           cold_handle, c_loc(cold_hit), c_loc(cold_count), c_loc(cold_fingerprint), 17_c_int, &
           c_loc(cold_fingerprint_len), c_loc(cold_path), 513_c_int, c_loc(cold_path_len) &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'cold diagnostics status')
  call assert_equal_i32(cold_hit, 0_i32, 'cold cache hit')
  call assert_equal_i32(cold_count, 1_i32, 'cold operator build count')
  call assert_true(cold_fingerprint_len > 0_c_int, 'cold fingerprint must be present')
  call assert_true(cold_path_len > 0_c_int, 'cold path must be present')
  call assert_true(cold_fingerprint(int(cold_fingerprint_len) + 1) == achar(0, kind=c_char), 'cold fingerprint NUL')
  call assert_true(cold_path(int(cold_path_len) + 1) == achar(0, kind=c_char), 'cold path NUL')
  cold_path_text = c_text(cold_path, cold_path_len)
  status = beach_kernel_destroy(cold_handle)
  call assert_equal_i32(status, beach_kernel_ok, 'cold destroy status')

  call build_cached_handle(warm_handle)
  status = beach_kernel_get_periodic_cache_info( &
           warm_handle, c_loc(warm_hit), c_loc(warm_count), c_loc(warm_fingerprint), 17_c_int, &
           c_loc(warm_fingerprint_len), c_loc(warm_path), 513_c_int, c_loc(warm_path_len) &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'warm diagnostics status')
  call assert_equal_i32(warm_hit, 1_i32, 'warm cache hit')
  call assert_equal_i32(warm_count, 0_i32, 'warm operator build count')
  call assert_equal_i32(warm_fingerprint_len, cold_fingerprint_len, 'fingerprint length agreement')
  call assert_equal_i32(warm_path_len, cold_path_len, 'path length agreement')
  call assert_true(warm_fingerprint(int(warm_fingerprint_len) + 1) == achar(0, kind=c_char), 'warm fingerprint NUL')
  call assert_true(warm_path(int(warm_path_len) + 1) == achar(0, kind=c_char), 'warm path NUL')
  call assert_true( &
    all(warm_fingerprint(1:int(warm_fingerprint_len)) == cold_fingerprint(1:int(cold_fingerprint_len))), &
    'cold/warm fingerprint agreement' &
    )
  warm_path_text = c_text(warm_path, warm_path_len)
  call assert_true(trim(warm_path_text) == trim(cold_path_text), 'cold/warm path agreement')
  call assert_true(index(trim(warm_path_text), trim(cache_dir)//'/') == 1, 'cache path must use configured directory')
  status = beach_kernel_destroy(warm_handle)
  call assert_equal_i32(status, beach_kernel_ok, 'warm destroy status')

  call delete_file_if_exists(trim(cold_path_text))
  call delete_file_if_exists(trim(cold_path_text)//'.lock')
  call remove_empty_directory(trim(cache_dir))
  call test_summary()

contains

  subroutine build_cached_handle(handle)
    type(c_ptr), intent(out) :: handle

    status = beach_kernel_create(handle)
    call assert_equal_i32(status, beach_kernel_ok, 'cached create status')
    status = beach_kernel_set_periodic_cache( &
             handle, c_loc(cache_dir_bytes), int(size(cache_dir_bytes), c_int), 2.5d-9 &
             )
    call assert_equal_i32(status, beach_kernel_ok, 'cached setter status')
    status = beach_kernel_build( &
             handle, 4_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.45d0, 2_c_int, 2_c_int, 1_c_int, &
             c_loc(periodic_axes), c_loc(periodic_len), 1_c_int, 3_c_int, 0.0d0, 4_c_int, &
             c_loc(box_min), c_loc(box_max) &
             )
    call assert_equal_i32(status, beach_kernel_ok, 'cached build status')
  end subroutine build_cached_handle

  subroutine to_c_bytes(text, bytes)
    character(len=*), intent(in) :: text
    character(kind=c_char), allocatable, target, intent(out) :: bytes(:)
    integer :: i

    allocate (bytes(len(text)))
    do i = 1, len(text)
      bytes(i) = char(iachar(text(i:i)), kind=c_char)
    end do
  end subroutine to_c_bytes

  function c_text(bytes, byte_count) result(text)
    character(kind=c_char), intent(in) :: bytes(:)
    integer(c_int), intent(in) :: byte_count
    character(len=512) :: text
    integer :: i

    text = ''
    do i = 1, int(byte_count)
      text(i:i) = achar(iachar(bytes(i)))
    end do
  end function c_text

end program test_field_kernel_cache_c
