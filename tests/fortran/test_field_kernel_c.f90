program test_field_kernel_c
  use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, ieee_value
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_loc, c_null_char, c_null_ptr, c_ptr
  use bem_field_kernel_c, only: beach_kernel_build, beach_kernel_create, &
                                beach_kernel_abi_major, beach_kernel_abi_minor, &
                                beach_kernel_destroy, beach_kernel_eval_e, beach_kernel_eval_e_direct, &
                                beach_kernel_eval_phi, beach_kernel_eval_phi_direct, &
                                beach_kernel_force_on_charges, beach_kernel_get_abi_version, &
                                beach_kernel_get_build_info, &
                                beach_kernel_get_periodic_cache_info, &
                                beach_kernel_invalid_argument, beach_kernel_invalid_handle, beach_kernel_not_ready, &
                                beach_kernel_ok, beach_kernel_set_periodic_cache, beach_kernel_update_charges
  use bem_kinds, only: dp, i32
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use bem_version, only: beach_build_info
  use test_support, only: assert_allclose_1d, assert_equal_i32, assert_true, test_begin, test_summary
  implicit none

  type(c_ptr) :: handle
  integer(c_int) :: status
  real(c_double), target :: v0(3, 2), v1(3, 2), v2(3, 2), src_q(2), target_pos(3, 2), e(3, 2), e_direct(3, 2)
  real(c_double), target :: target_q(1), origin(3), force(3), torque(3)
  real(c_double), target :: panel_q(2), phi(2), phi_direct(2)
  real(c_double) :: expected_e(3), expected_force(3), expected_torque(3)
  real(dp) :: expected_phi
  type(panel_geometry_type) :: geometry
  integer(i32) :: geometry_status
  integer(c_int), target :: abi_major, abi_minor, build_info_length, cache_hit, cache_build_count
  integer(c_int), target :: fingerprint_length, path_length
  character(kind=c_char), target :: cache_path(8), embedded_nul_path(3), invalid_utf8_path(2), long_path(257)
  character(kind=c_char), target :: valid_utf8_path(5), overlong_path(2), surrogate_path(3), out_of_range_path(4)
  character(kind=c_char), target :: truncated_utf8_path(2), blank_path(1), trailing_blank_path(2)
  character(kind=c_char), target :: build_info_buffer(512), fingerprint_buffer(17), path_buffer(513)
  integer(c_int), target :: cache_periodic_axes(2)
  real(c_double), target :: cache_periodic_len(2), cache_box_min(3), cache_box_max(3)
  integer(i32) :: panel_target_idx
  integer :: i

  call test_begin('field_kernel_c_abi_version')
  abi_major = -1_c_int
  abi_minor = -1_c_int
  status = beach_kernel_get_abi_version(c_null_ptr, c_loc(abi_minor))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'ABI version NULL major output')
  status = beach_kernel_get_abi_version(c_loc(abi_major), c_null_ptr)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'ABI version NULL minor output')
  status = beach_kernel_get_abi_version(c_loc(abi_major), c_loc(abi_minor))
  call assert_equal_i32(status, beach_kernel_ok, 'ABI version getter status')
  call assert_equal_i32(abi_major, beach_kernel_abi_major, 'ABI major version')
  call assert_equal_i32(abi_minor, beach_kernel_abi_minor, 'ABI minor version')
  call assert_equal_i32(abi_major, 2_i32, 'panel-only ABI major version')

  call test_begin('field_kernel_c_build_info')
  build_info_buffer = achar(88, kind=c_char)
  build_info_length = -1_c_int
  status = beach_kernel_get_build_info(c_null_ptr, 512_c_int, c_loc(build_info_length))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'build-info NULL buffer')
  status = beach_kernel_get_build_info(c_loc(build_info_buffer), 1_c_int, c_loc(build_info_length))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'build-info undersized buffer')
  call assert_equal_i32(build_info_length, int(len(beach_build_info), i32), 'build-info required length')
  status = beach_kernel_get_build_info(c_loc(build_info_buffer), 512_c_int, c_loc(build_info_length))
  call assert_equal_i32(status, beach_kernel_ok, 'build-info getter status')
  do i = 1, len(beach_build_info)
    call assert_true(transfer(build_info_buffer(i), ' ') == beach_build_info(i:i), 'build-info payload mismatch')
  end do
  call assert_true(build_info_buffer(len(beach_build_info) + 1) == c_null_char, 'build-info must be NUL terminated')

  call test_begin('field_kernel_c_cache_abi_validation')
  call set_c_bytes(cache_path, [99, 97, 99, 104, 101, 45, 195, 169])
  call set_c_bytes(embedded_nul_path, [97, 0, 98])
  call set_c_bytes(invalid_utf8_path, [195, 40])
  call set_c_bytes(valid_utf8_path, [97, 240, 159, 152, 128])
  call set_c_bytes(overlong_path, [192, 128])
  call set_c_bytes(surrogate_path, [237, 160, 128])
  call set_c_bytes(out_of_range_path, [244, 144, 128, 128])
  call set_c_bytes(truncated_utf8_path, [226, 130])
  call set_c_bytes(blank_path, [32])
  call set_c_bytes(trailing_blank_path, [97, 32])
  long_path = achar(97, kind=c_char)
  fingerprint_buffer = achar(88, kind=c_char)
  path_buffer = achar(89, kind=c_char)

  status = beach_kernel_set_periodic_cache(c_null_ptr, c_loc(cache_path), 8_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_handle, 'cache setter NULL handle')
  status = beach_kernel_get_periodic_cache_info( &
           c_null_ptr, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), 17_c_int, &
           c_loc(fingerprint_length), c_loc(path_buffer), 513_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_handle, 'cache getter NULL handle')

  status = beach_kernel_create(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'cache ABI create status')
  status = beach_kernel_get_periodic_cache_info( &
           handle, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), 17_c_int, &
           c_loc(fingerprint_length), c_loc(path_buffer), 513_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_not_ready, 'cache getter before build')

  status = beach_kernel_set_periodic_cache(handle, c_null_ptr, 8_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter NULL path')
  status = beach_kernel_set_periodic_cache(handle, c_loc(cache_path), 0_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter empty path')
  status = beach_kernel_set_periodic_cache(handle, c_loc(long_path), 257_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter oversized path')
  status = beach_kernel_set_periodic_cache(handle, c_loc(embedded_nul_path), 3_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter embedded NUL')
  status = beach_kernel_set_periodic_cache(handle, c_loc(invalid_utf8_path), 2_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter invalid UTF-8')
  status = beach_kernel_set_periodic_cache(handle, c_loc(overlong_path), 2_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter overlong UTF-8')
  status = beach_kernel_set_periodic_cache(handle, c_loc(surrogate_path), 3_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter surrogate UTF-8')
  status = beach_kernel_set_periodic_cache(handle, c_loc(out_of_range_path), 4_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter out-of-range UTF-8')
  status = beach_kernel_set_periodic_cache(handle, c_loc(truncated_utf8_path), 2_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter truncated UTF-8')
  status = beach_kernel_set_periodic_cache(handle, c_loc(blank_path), 1_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter blank path')
  status = beach_kernel_set_periodic_cache(handle, c_loc(trailing_blank_path), 2_c_int, 1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter trailing blank path')
  status = beach_kernel_set_periodic_cache(handle, c_loc(cache_path), 8_c_int, 0.0d0)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter zero tolerance')
  status = beach_kernel_set_periodic_cache(handle, c_loc(cache_path), 8_c_int, -1.0d-8)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter negative tolerance')
  status = beach_kernel_set_periodic_cache( &
           handle, c_loc(cache_path), 8_c_int, ieee_value(0.0d0, ieee_quiet_nan) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter NaN tolerance')
  status = beach_kernel_set_periodic_cache( &
           handle, c_loc(cache_path), 8_c_int, ieee_value(0.0d0, ieee_positive_inf) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache setter infinite tolerance')
  status = beach_kernel_set_periodic_cache(handle, c_loc(valid_utf8_path), 5_c_int, 2.5d-9)
  call assert_equal_i32(status, beach_kernel_ok, 'cache setter valid four-byte UTF-8 path')
  status = beach_kernel_set_periodic_cache(handle, c_loc(cache_path), 8_c_int, 2.5d-9)
  call assert_equal_i32(status, beach_kernel_ok, 'cache setter valid UTF-8 path')

  v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  v1(:, 1) = [0.2d0, 0.0d0, 0.0d0]
  v2(:, 1) = [0.0d0, 0.2d0, 0.0d0]
  v0(:, 2) = [1.0d0, 0.0d0, 0.0d0]
  v1(:, 2) = [1.2d0, 0.0d0, 0.0d0]
  v2(:, 2) = [1.0d0, 0.2d0, 0.0d0]
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_null_ptr, c_loc(v2), 0.5d0, 8_c_int, 4_c_int, 0_c_int, c_null_ptr, &
           c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'panel build rejects NULL vertex array')
  v2(:, 1) = v0(:, 1)
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 4_c_int, 0_c_int, c_null_ptr, &
           c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'panel build rejects degenerate triangles')
  v2(:, 1) = [0.0d0, 0.2d0, 0.0d0]
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 0_c_int, 0_c_int, c_null_ptr, &
           c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'panel build requires positive FMM order')
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 4_c_int, 0_c_int, c_null_ptr, &
           c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'cache ABI non-cached build')

  fingerprint_length = -1_c_int
  path_length = -1_c_int
  status = beach_kernel_get_periodic_cache_info( &
           handle, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), 0_c_int, &
           c_loc(fingerprint_length), c_loc(path_buffer), 513_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache getter requires fingerprint NUL capacity')
  call assert_equal_i32(fingerprint_length, 0_i32, 'cache getter short buffer fingerprint required length')
  call assert_equal_i32(path_length, 0_i32, 'cache getter short buffer path required length')
  call assert_true(fingerprint_buffer(1) == achar(88, kind=c_char), 'cache getter must not truncate fingerprint')
  status = beach_kernel_get_periodic_cache_info( &
           handle, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), 17_c_int, &
           c_loc(fingerprint_length), c_loc(path_buffer), 0_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache getter requires path NUL capacity')
  call assert_equal_i32(fingerprint_length, 0_i32, 'cache getter zero path capacity fingerprint length')
  call assert_equal_i32(path_length, 0_i32, 'cache getter zero path capacity required length')

  status = beach_kernel_get_periodic_cache_info( &
           handle, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), 17_c_int, &
           c_loc(fingerprint_length), c_loc(path_buffer), 513_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'cache getter non-cached status')
  call assert_equal_i32(cache_hit, 0_i32, 'non-cached cache hit')
  call assert_equal_i32(cache_build_count, 0_i32, 'non-cached operator build count')
  call assert_equal_i32(fingerprint_length, 0_i32, 'non-cached fingerprint length')
  call assert_equal_i32(path_length, 0_i32, 'non-cached path length')
  call assert_true(fingerprint_buffer(1) == c_null_char, 'non-cached fingerprint NUL')
  call assert_true(path_buffer(1) == c_null_char, 'non-cached path NUL')
  status = beach_kernel_get_periodic_cache_info( &
           handle, c_null_ptr, c_loc(cache_build_count), c_loc(fingerprint_buffer), 17_c_int, &
           c_loc(fingerprint_length), c_loc(path_buffer), 513_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache getter NULL hit output')
  status = beach_kernel_get_periodic_cache_info( &
           handle, c_loc(cache_hit), c_null_ptr, c_loc(fingerprint_buffer), 17_c_int, &
           c_loc(fingerprint_length), c_loc(path_buffer), 513_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache getter NULL build count output')
  status = beach_kernel_get_periodic_cache_info( &
           handle, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), 17_c_int, &
           c_null_ptr, c_loc(path_buffer), 513_c_int, c_loc(path_length) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cache getter NULL fingerprint length output')

  cache_periodic_axes = [1_c_int, 2_c_int]
  cache_periodic_len = [2.0d0, 2.0d0]
  cache_box_min = [0.0d0, 0.0d0, -1.0d0]
  cache_box_max = [2.0d0, 2.0d0, 1.0d0]
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 2_c_int, 1_c_int, &
           c_loc(cache_periodic_axes), c_loc(cache_periodic_len), 1_c_int, 2_c_int, 0.0d0, 4_c_int, &
           c_loc(cache_box_min), c_loc(cache_box_max) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'removed far-correction ABI code 2')
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 2_c_int, 1_c_int, &
           c_loc(cache_periodic_axes), c_loc(cache_periodic_len), 0_c_int, 3_c_int, 0.0d0, 4_c_int, &
           c_loc(cache_box_min), c_loc(cache_box_max) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cached build requires an image layer')
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 2_c_int, 1_c_int, &
           c_loc(cache_periodic_axes), c_loc(cache_periodic_len), 1_c_int, 3_c_int, 0.0d0, 0_c_int, &
           c_loc(cache_box_min), c_loc(cache_box_max) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cached build requires an Ewald layer')
  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 0_c_int, 1_c_int, &
           c_loc(cache_periodic_axes), c_loc(cache_periodic_len), 1_c_int, 3_c_int, 0.0d0, 4_c_int, &
           c_loc(cache_box_min), c_loc(cache_box_max) &
           )
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'cached panel build requires positive order')
  status = beach_kernel_destroy(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'cache ABI destroy status')

  call test_begin('field_kernel_c_panel_free_eval_and_force')

  v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  v1(:, 1) = [0.2d0, 0.0d0, 0.0d0]
  v2(:, 1) = [0.0d0, 0.2d0, 0.0d0]
  v0(:, 2) = [1.0d0, 0.0d0, 0.0d0]
  v1(:, 2) = [1.2d0, 0.0d0, 0.0d0]
  v2(:, 2) = [1.0d0, 0.2d0, 0.0d0]
  src_q = [1.0d-9, -2.0d-9]
  target_pos(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  target_q = [3.0d-9]
  origin = [0.0d0, 0.0d0, 0.0d0]

  status = beach_kernel_create(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'create status')

  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 4_c_int, 0_c_int, c_null_ptr, &
           c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'build status')

  status = beach_kernel_update_charges(handle, 2_c_int, c_loc(src_q))
  call assert_equal_i32(status, beach_kernel_ok, 'update status')

  status = beach_kernel_eval_e_direct(handle, -1_c_int, c_loc(target_pos), c_loc(e_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'eval_e_direct negative target count')
  status = beach_kernel_eval_phi_direct(handle, huge(0_c_int), c_loc(target_pos), c_loc(phi_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'eval_phi_direct oversized target count')
  status = beach_kernel_eval_e_direct(handle, 1_c_int, c_null_ptr, c_loc(e_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'eval_e_direct NULL target')
  status = beach_kernel_eval_phi_direct(handle, 1_c_int, c_loc(target_pos), c_null_ptr)
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'eval_phi_direct NULL output')
  target_pos(:, 1) = [ieee_value(0.0d0, ieee_quiet_nan), 1.0d0, 0.0d0]
  status = beach_kernel_eval_e_direct(handle, 1_c_int, c_loc(target_pos), c_loc(e_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'eval_e_direct nonfinite target')
  status = beach_kernel_eval_phi_direct(handle, 1_c_int, c_loc(target_pos), c_loc(phi_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'eval_phi_direct nonfinite target')
  target_pos(:, 1) = [0.0d0, 1.0d0, 0.0d0]

  e = 0.0d0
  status = beach_kernel_eval_e(handle, 1_c_int, c_loc(target_pos), c_loc(e))
  call assert_equal_i32(status, beach_kernel_ok, 'eval_e status')

  call direct_panel_field(v0, v1, v2, src_q, target_pos(:, 1), expected_phi, expected_e)
  call assert_allclose_1d( &
    e(:, 1), expected_e, 1.0d-12*max(1.0d0, maxval(abs(expected_e))), 'eval_e value' &
    )
  status = beach_kernel_eval_e_direct(handle, 1_c_int, c_loc(target_pos), c_loc(e_direct))
  call assert_equal_i32(status, beach_kernel_ok, 'eval_e_direct status')
  status = beach_kernel_eval_phi_direct(handle, 1_c_int, c_loc(target_pos), c_loc(phi_direct))
  call assert_equal_i32(status, beach_kernel_ok, 'eval_phi_direct status')
  call assert_allclose_1d( &
    e_direct(:, 1), expected_e, 1.0d-12*max(1.0d0, maxval(abs(expected_e))), 'eval_e_direct value' &
    )
  call assert_allclose_1d( &
    phi_direct(:1), [expected_phi], 1.0d-12*max(1.0d0, abs(expected_phi)), 'eval_phi_direct value' &
    )

  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 4_c_int, 1_c_int, &
           c_loc(cache_periodic_axes), c_loc(cache_periodic_len), 1_c_int, 1_c_int, 0.0d0, 4_c_int, &
           c_loc(cache_box_min), c_loc(cache_box_max) &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'periodic direct rejection build status')
  status = beach_kernel_eval_e_direct(handle, 1_c_int, c_loc(target_pos), c_loc(e_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'uncharged periodic eval_e_direct rejection')
  status = beach_kernel_eval_phi_direct(handle, 1_c_int, c_loc(target_pos), c_loc(phi_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'uncharged periodic eval_phi_direct rejection')
  status = beach_kernel_update_charges(handle, 2_c_int, c_loc(src_q))
  call assert_equal_i32(status, beach_kernel_ok, 'periodic direct rejection update status')
  status = beach_kernel_eval_e_direct(handle, 1_c_int, c_loc(target_pos), c_loc(e_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'periodic eval_e_direct rejection')
  status = beach_kernel_eval_phi_direct(handle, 1_c_int, c_loc(target_pos), c_loc(phi_direct))
  call assert_equal_i32(status, beach_kernel_invalid_argument, 'periodic eval_phi_direct rejection')

  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 4_c_int, 0_c_int, c_null_ptr, &
           c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'restore free build status')
  status = beach_kernel_update_charges(handle, 2_c_int, c_loc(src_q))
  call assert_equal_i32(status, beach_kernel_ok, 'restore free update status')

  force = 0.0d0
  torque = 0.0d0
  status = beach_kernel_force_on_charges( &
           handle, 1_c_int, c_loc(target_pos), c_loc(target_q), c_loc(origin), c_loc(force), c_loc(torque) &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'force status')
  expected_force = target_q(1)*expected_e
  expected_torque = cross(target_pos(:, 1) - origin, expected_force)
  call assert_allclose_1d(force, expected_force, 1.0d-20, 'force value')
  call assert_allclose_1d(torque, expected_torque, 1.0d-20, 'torque value')

  status = beach_kernel_destroy(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'destroy status')

  call test_begin('field_kernel_c_panel_eval')
  v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  v1(:, 1) = [2.0d0, 0.0d0, 0.0d0]
  v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  panel_q(1) = 2.5d-9
  target_pos(:, 1) = [0.25d0, 0.2d0, 0.4d0]
  target_pos(:, 2) = [0.25d0, 0.2d0, 0.0d0]
  status = beach_kernel_create(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'panel create status')
  status = beach_kernel_build( &
           handle, 1_c_int, c_loc(v0), c_loc(v1), c_loc(v2), 0.5d0, 8_c_int, 4_c_int, &
           0_c_int, c_null_ptr, c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'panel build status')
  status = beach_kernel_update_charges(handle, 1_c_int, c_loc(panel_q))
  call assert_equal_i32(status, beach_kernel_ok, 'panel update status')
  status = beach_kernel_eval_e(handle, 1_c_int, c_loc(target_pos), c_loc(e))
  call assert_equal_i32(status, beach_kernel_ok, 'panel eval_e status')
  status = beach_kernel_eval_phi(handle, 1_c_int, c_loc(target_pos), c_loc(phi))
  call assert_equal_i32(status, beach_kernel_ok, 'panel eval_phi status')
  status = beach_kernel_eval_e_direct(handle, 2_c_int, c_loc(target_pos), c_loc(e_direct))
  call assert_equal_i32(status, beach_kernel_ok, 'panel eval_e_direct status')
  status = beach_kernel_eval_phi_direct(handle, 2_c_int, c_loc(target_pos), c_loc(phi_direct))
  call assert_equal_i32(status, beach_kernel_ok, 'panel eval_phi_direct status')
  do panel_target_idx = 1_i32, 2_i32
    call init_panel_geometry(v0(:, 1), v1(:, 1), v2(:, 1), geometry, geometry_status)
    call assert_equal_i32(geometry_status, panel_geometry_ok, 'panel fixture geometry status')
    call panel_potential_field( &
      geometry, panel_q(1), target_pos(:, panel_target_idx), panel_side_principal_value, expected_phi, expected_e &
      )
    call assert_allclose_1d( &
      e_direct(:, panel_target_idx), expected_e, 1.0d-12*max(1.0d0, maxval(abs(expected_e))), &
      'panel eval_e_direct value' &
      )
    call assert_allclose_1d( &
      phi_direct(panel_target_idx:panel_target_idx), [expected_phi], 1.0d-12*max(1.0d0, abs(expected_phi)), &
      'panel eval_phi_direct value' &
      )
    if (panel_target_idx == 1_i32) then
      call assert_allclose_1d(e(:, 1), expected_e, 1.0d-12*max(1.0d0, maxval(abs(expected_e))), 'panel eval_e value')
      call assert_allclose_1d(phi(:1), [expected_phi], 1.0d-12*max(1.0d0, abs(expected_phi)), 'panel eval_phi value')
    end if
  end do
  status = beach_kernel_destroy(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'panel destroy status')

  call test_summary()

contains

  subroutine set_c_bytes(output, bytes)
    character(kind=c_char), intent(out) :: output(:)
    integer, intent(in) :: bytes(:)
    integer :: i

    do i = 1, size(output)
      output(i) = achar(bytes(i), kind=c_char)
    end do
  end subroutine set_c_bytes

  subroutine direct_panel_field(vertices0, vertices1, vertices2, charges, target, phi_out, field_out)
    real(c_double), intent(in) :: vertices0(:, :), vertices1(:, :), vertices2(:, :)
    real(c_double), intent(in) :: charges(:), target(3)
    real(c_double), intent(out) :: phi_out, field_out(3)
    type(panel_geometry_type) :: panel_geometry
    integer(i32) :: panel_status
    real(dp) :: panel_phi, panel_field(3)
    integer :: i

    phi_out = 0.0d0
    field_out = 0.0d0
    do i = 1, size(charges)
      call init_panel_geometry( &
        vertices0(:, i), vertices1(:, i), vertices2(:, i), panel_geometry, panel_status &
        )
      call assert_equal_i32(panel_status, panel_geometry_ok, 'panel fixture geometry status')
      call panel_potential_field( &
        panel_geometry, charges(i), target, panel_side_principal_value, panel_phi, panel_field &
        )
      phi_out = phi_out + panel_phi
      field_out = field_out + panel_field
    end do
  end subroutine direct_panel_field

  pure function cross(a, b) result(c)
    real(c_double), intent(in) :: a(3), b(3)
    real(c_double) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

end program test_field_kernel_c
