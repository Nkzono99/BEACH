program test_periodic2_operator_cache
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, build_plan, destroy_plan
  use test_support, only: test_begin, test_end, test_summary, assert_true, assert_equal_i32, assert_allclose_1d
  implicit none

  type(fmm_options_type) :: options
  type(fmm_plan_type) :: probe_plan, cold_plan, warm_plan, rebuilt_plan
  real(dp) :: src_pos(3, 4)
  real(dp) :: cold_start, cold_end, warm_start, warm_end
  character(len=512) :: cache_path
  integer :: unit, ios

  src_pos(:, 1) = [0.15_dp, 0.20_dp, -0.10_dp]
  src_pos(:, 2) = [0.75_dp, 0.25_dp, 0.08_dp]
  src_pos(:, 3) = [0.30_dp, 0.80_dp, 0.12_dp]
  src_pos(:, 4) = [0.82_dp, 0.72_dp, -0.06_dp]
  options%theta = 0.45_dp
  options%leaf_max = 2_i32
  options%order = 2_i32
  options%use_periodic2 = .true.
  options%periodic_axes = [1_i32, 2_i32]
  options%periodic_len = [1.0_dp, 1.0_dp]
  options%periodic_image_layers = 1_i32
  options%periodic_far_correction = 'cached_kneq0'
  options%periodic_ewald_layers = 2_i32
  options%periodic_cache_dir = '.'
  options%periodic_generation_tolerance = 1.0e-8_dp
  options%target_box_min = [0.0_dp, 0.0_dp, -0.5_dp]
  options%target_box_max = [1.0_dp, 1.0_dp, 0.5_dp]

  call build_plan(probe_plan, src_pos, options)
  cache_path = probe_plan%periodic_cache_path
  call destroy_plan(probe_plan)
  call delete_cache(cache_path)

  call test_begin('cold_build_and_warm_hit')
  call cpu_time(cold_start)
  call build_plan(cold_plan, src_pos, options)
  call cpu_time(cold_end)
  call assert_true(cold_plan%periodic_root_operator_ready, 'cold cached operator must be ready')
  call assert_true(.not. cold_plan%periodic_cache_hit, 'cold operator must report a cache miss')
  call assert_equal_i32(cold_plan%periodic_operator_build_count, 1_i32, 'cold operator build count')
  call cpu_time(warm_start)
  call build_plan(warm_plan, src_pos, options)
  call cpu_time(warm_end)
  call assert_true(warm_plan%periodic_cache_hit, 'second identical plan must hit the cache')
  call assert_equal_i32(warm_plan%periodic_operator_build_count, 0_i32, 'warm operator build count')
  call assert_true( &
    trim(warm_plan%periodic_cache_fingerprint) == trim(cold_plan%periodic_cache_fingerprint), &
    'warm fingerprint mismatch' &
    )
  call assert_allclose_1d( &
    reshape(warm_plan%periodic_root_operator, [size(warm_plan%periodic_root_operator)]), &
    reshape(cold_plan%periodic_root_operator, [size(cold_plan%periodic_root_operator)]), 0.0_dp, &
    'warm operator differs from cold operator' &
    )
  write (*, '(a,2(es12.5,1x))') 'periodic cache build times(cold,warm)=', cold_end - cold_start, warm_end - warm_start
  call assert_true(warm_end - warm_start < 0.5_dp*(cold_end - cold_start), &
                   'warm cache build must be clearly faster than cold generation')
  call test_end()

  call test_begin('corrupt_cache_rebuilds')
  open (newunit=unit, file=trim(cache_path), access='stream', form='unformatted', status='old', action='write', iostat=ios)
  call assert_equal_i32(int(ios, i32), 0_i32, 'open operator cache for corruption')
  write (unit, pos=1, iostat=ios) 'BROKEN-CACHE!!!'
  close (unit)
  call assert_equal_i32(int(ios, i32), 0_i32, 'corrupt operator cache write')
  call build_plan(rebuilt_plan, src_pos, options)
  call assert_true(.not. rebuilt_plan%periodic_cache_hit, 'corrupt cache must not be reported as a hit')
  call assert_equal_i32(rebuilt_plan%periodic_operator_build_count, 1_i32, 'corrupt cache rebuild count')
  call assert_allclose_1d( &
    reshape(rebuilt_plan%periodic_root_operator, [size(rebuilt_plan%periodic_root_operator)]), &
    reshape(cold_plan%periodic_root_operator, [size(cold_plan%periodic_root_operator)]), 0.0_dp, &
    'rebuilt operator differs from deterministic cold operator' &
    )
  call test_end()

  call destroy_plan(cold_plan)
  call destroy_plan(warm_plan)
  call destroy_plan(rebuilt_plan)
  call delete_cache(cache_path)
  call test_summary()

contains

  subroutine delete_cache(path)
    character(len=*), intent(in) :: path
    call delete_file(path)
    call delete_file(trim(path)//'.lock')
  end subroutine delete_cache

  subroutine delete_file(path)
    character(len=*), intent(in) :: path
    logical :: exists
    integer :: delete_unit
    inquire (file=trim(path), exist=exists)
    if (.not. exists) return
    open (newunit=delete_unit, file=trim(path), status='old')
    close (delete_unit, status='delete')
  end subroutine delete_file
end program test_periodic2_operator_cache
