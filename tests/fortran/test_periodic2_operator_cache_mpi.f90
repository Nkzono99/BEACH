program test_periodic2_operator_cache_mpi
  use bem_kinds, only: dp, i32
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_is_root, &
                     mpi_bcast_i32_array, mpi_allreduce_sum_i32_scalar, mpi_allreduce_min_real_dp_array, &
                     mpi_allreduce_max_real_dp_array, mpi_world_barrier
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, build_plan, destroy_plan
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          assert_allclose_1d, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mpi_context) :: mpi
  type(fmm_options_type) :: options
  type(fmm_plan_type) :: cold_plan, warm_plan, rebuilt_plan
  real(dp) :: src_pos(3, 4), local_target_min(1), local_target_max(1)
  integer(i32) :: global_local_target_count, cache_nonce(1)
  integer :: clock_count, unit, ios
  character(len=512) :: cache_path, cache_dir

  call mpi_initialize(mpi)
  call configure_fixture(src_pos, options)
  cache_nonce = 0_i32
  if (mpi_is_root(mpi)) then
    call system_clock(count=clock_count)
    cache_nonce(1) = int(clock_count, i32)
  end if
  call mpi_bcast_i32_array(mpi, cache_nonce, 0_i32)
  write (cache_dir, '(a,i0)') 'test_periodic2_operator_cache_mpi_tmp_', cache_nonce(1)
  options%periodic_cache_dir = trim(cache_dir)

  call test_init(3)
  call test_begin('mpi_distributed_cold_build')
  call build_plan(cold_plan, src_pos, options)
  cache_path = cold_plan%periodic_cache_path
  call assert_true(.not. cold_plan%periodic_cache_hit, 'cold cache miss must be reported on every rank')
  call assert_equal_i32( &
    cold_plan%periodic_operator_build_count, merge(1_i32, 0_i32, mpi_is_root(mpi)), &
    'only the root rank may publish the generated operator' &
    )
  call assert_true(cold_plan%periodic_root_operator_ready, 'collective cached operator must be ready')
  call assert_true(len_trim(cold_plan%periodic_cache_fingerprint) > 0, 'fingerprint must be present on every rank')
  global_local_target_count = cold_plan%periodic_operator_local_target_count
  call mpi_allreduce_sum_i32_scalar(mpi, global_local_target_count)
  call assert_equal_i32( &
    global_local_target_count, cold_plan%periodic_root_target_count, &
    'cold MPI target work must cover each target exactly once' &
    )
  local_target_min(1) = real(cold_plan%periodic_operator_local_target_count, dp)
  local_target_max = local_target_min
  call mpi_allreduce_min_real_dp_array(mpi, local_target_min)
  call mpi_allreduce_max_real_dp_array(mpi, local_target_max)
  call assert_true(local_target_max(1) - local_target_min(1) <= 1.0_dp, 'cold MPI target work is imbalanced')
  call test_end()

  call test_begin('mpi_warm_cache_broadcast')
  call build_plan(warm_plan, src_pos, options)
  call assert_true(warm_plan%periodic_cache_hit, 'warm cache hit must be reported on every rank')
  call assert_equal_i32(warm_plan%periodic_operator_build_count, 0_i32, 'warm cache must not regenerate the operator')
  call assert_equal_i32(warm_plan%periodic_operator_local_target_count, 0_i32, 'warm cache local target count')
  call assert_true(warm_plan%periodic_root_operator_ready, 'broadcast cached operator must be ready')
  call assert_equal_i32( &
    warm_plan%periodic_root_target_count, cold_plan%periodic_root_target_count, &
    'broadcast target count differs from the generated operator' &
    )
  call assert_true( &
    all(warm_plan%periodic_root_target_nodes == cold_plan%periodic_root_target_nodes), &
    'broadcast target nodes differ from the generated operator' &
    )
  call assert_allclose_1d( &
    reshape(warm_plan%periodic_root_operator, [size(warm_plan%periodic_root_operator)]), &
    reshape(cold_plan%periodic_root_operator, [size(cold_plan%periodic_root_operator)]), 0.0_dp, &
    'broadcast operator differs from the generated operator' &
    )
  call assert_true( &
    trim(warm_plan%periodic_cache_fingerprint) == trim(cold_plan%periodic_cache_fingerprint), &
    'cache fingerprint differs across cold and warm builds' &
    )
  call assert_true(trim(warm_plan%periodic_cache_path) == trim(cold_plan%periodic_cache_path), &
                   'cache path differs across cold and warm builds')
  call test_end()

  call test_begin('mpi_corrupt_cache_collective_rebuild')
  if (mpi_is_root(mpi)) then
    open (newunit=unit, file=trim(cache_path), access='stream', form='unformatted', status='old', &
          action='write', iostat=ios)
    call assert_equal_i32(int(ios, i32), 0_i32, 'open operator cache for corruption')
    write (unit, pos=1, iostat=ios) 'BROKEN-CACHE!!!'
    close (unit)
    call assert_equal_i32(int(ios, i32), 0_i32, 'corrupt operator cache write')
  end if
  call mpi_world_barrier(mpi)
  call build_plan(rebuilt_plan, src_pos, options)
  call assert_true(.not. rebuilt_plan%periodic_cache_hit, 'corrupt cache must be rejected on every rank')
  call assert_equal_i32( &
    rebuilt_plan%periodic_operator_build_count, merge(1_i32, 0_i32, mpi_is_root(mpi)), &
    'only the root rank may publish the rebuilt operator' &
    )
  global_local_target_count = rebuilt_plan%periodic_operator_local_target_count
  call mpi_allreduce_sum_i32_scalar(mpi, global_local_target_count)
  call assert_equal_i32( &
    global_local_target_count, rebuilt_plan%periodic_root_target_count, &
    'corrupt-cache rebuild must cover each target exactly once' &
    )
  call assert_true(rebuilt_plan%periodic_root_operator_ready, 'collectively rebuilt operator must be ready')
  call assert_equal_i32( &
    rebuilt_plan%periodic_root_target_count, cold_plan%periodic_root_target_count, &
    'rebuilt target count differs from the generated operator' &
    )
  call assert_true( &
    all(rebuilt_plan%periodic_root_target_nodes == cold_plan%periodic_root_target_nodes), &
    'rebuilt target nodes differ from the generated operator' &
    )
  call assert_allclose_1d( &
    reshape(rebuilt_plan%periodic_root_operator, [size(rebuilt_plan%periodic_root_operator)]), &
    reshape(cold_plan%periodic_root_operator, [size(cold_plan%periodic_root_operator)]), 0.0_dp, &
    'collectively rebuilt operator differs from the generated operator' &
    )
  call test_end()

  call destroy_plan(cold_plan)
  call destroy_plan(warm_plan)
  call destroy_plan(rebuilt_plan)
  if (mpi_is_root(mpi)) then
    call delete_file_if_exists(cache_path)
    call delete_file_if_exists(trim(cache_path)//'.lock')
    call remove_empty_directory(trim(cache_dir))
  end if
  call mpi_world_barrier(mpi)
  call test_summary()
  call mpi_shutdown(mpi)

contains

  subroutine configure_fixture(positions, config)
    real(dp), intent(out) :: positions(3, 4)
    type(fmm_options_type), intent(out) :: config

    positions(:, 1) = [0.15_dp, 0.20_dp, -0.10_dp]
    positions(:, 2) = [0.75_dp, 0.25_dp, 0.08_dp]
    positions(:, 3) = [0.30_dp, 0.80_dp, 0.12_dp]
    positions(:, 4) = [0.82_dp, 0.72_dp, -0.06_dp]
    config%theta = 0.45_dp
    config%leaf_max = 2_i32
    config%order = 2_i32
    config%use_periodic2 = .true.
    config%periodic_axes = [1_i32, 2_i32]
    config%periodic_len = [1.0_dp, 1.0_dp]
    config%periodic_image_layers = 1_i32
    config%periodic_far_correction = 'cached_kneq0'
    config%periodic_ewald_layers = 2_i32
    config%periodic_generation_tolerance = 1.0e-8_dp
    config%target_box_min = [0.0_dp, 0.0_dp, -0.5_dp]
    config%target_box_max = [1.0_dp, 1.0_dp, 0.5_dp]
  end subroutine configure_fixture

end program test_periodic2_operator_cache_mpi
