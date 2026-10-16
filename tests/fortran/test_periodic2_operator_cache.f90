program test_periodic2_operator_cache
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, build_panel_plan, &
                                  update_state, destroy_state, destroy_plan
  use bem_coulomb_fmm_periodic_root_ops, only: precompute_periodic_root_operator
  use test_support, only: test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          remove_empty_directory
  implicit none

  type(fmm_options_type) :: options
  type(fmm_plan_type) :: cold_plan, warm_plan, rebuilt_plan, recovered_plan
  real(dp) :: src_pos(3, 4)
  character(len=512) :: cache_path, cache_dir
  integer :: unit, ios, clock_count

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
  call system_clock(count=clock_count)
  write (cache_dir, '(a,i0)') 'test_periodic2_operator_cache_tmp_', clock_count
  options%periodic_cache_dir = trim(cache_dir)
  options%periodic_generation_tolerance = 1.0e-8_dp
  options%target_box_min = [0.0_dp, 0.0_dp, -0.5_dp]
  options%target_box_max = [1.0_dp, 1.0_dp, 0.5_dp]

  call test_begin('cold_build_and_warm_hit')
  call build_plan(cold_plan, src_pos, options)
  cache_path = cold_plan%periodic_cache_path
  call assert_true(.not. cold_plan%periodic_cache_hit, 'cold operator must report a cache miss')
  call assert_equal_i32(cold_plan%periodic_operator_build_count, 1_i32, 'cold operator build count')
  call build_plan(warm_plan, src_pos, options)
  call assert_true(warm_plan%periodic_cache_hit, 'second identical plan must hit the cache')
  call assert_equal_i32(warm_plan%periodic_operator_build_count, 0_i32, 'warm operator build count')
  call assert_true( &
    trim(warm_plan%periodic_cache_fingerprint) == trim(cold_plan%periodic_cache_fingerprint), &
    'warm fingerprint mismatch' &
    )
  call test_end()

  call test_begin('corrupt_cache_rebuilds_and_is_reused')
  open (newunit=unit, file=trim(cache_path), access='stream', form='unformatted', status='old', action='write', iostat=ios)
  call assert_equal_i32(int(ios, i32), 0_i32, 'open operator cache for corruption')
  write (unit, pos=1, iostat=ios) 'BROKEN-CACHE!!!'
  close (unit)
  call assert_equal_i32(int(ios, i32), 0_i32, 'corrupt operator cache write')
  call build_plan(rebuilt_plan, src_pos, options)
  call assert_true(.not. rebuilt_plan%periodic_cache_hit, 'corrupt cache must not be reported as a hit')
  call assert_equal_i32(rebuilt_plan%periodic_operator_build_count, 1_i32, 'corrupt cache rebuild count')
  call build_plan(recovered_plan, src_pos, options)
  call assert_true(recovered_plan%periodic_cache_hit, 'rebuilt cache must be reusable')
  call assert_equal_i32(recovered_plan%periodic_operator_build_count, 0_i32, 'recovered warm build count')
  call test_end()

  call destroy_plan(cold_plan)
  call destroy_plan(warm_plan)
  call destroy_plan(rebuilt_plan)
  call destroy_plan(recovered_plan)
  call delete_cache(cache_path)
  call check_box_mesh_reuse_and_translations()
  call check_zero_root_moments()
  call remove_empty_directory(cache_dir)
  call test_summary()

contains

  subroutine check_box_mesh_reuse_and_translations()
    type(fmm_plan_type) :: first, second, third
    type(fmm_options_type) :: config
    real(dp) :: a(3, 3), b(3, 3), c(3, 3), center(3)
    real(dp), allocatable :: shared(:, :, :)
    integer(i32) :: idx, shared_count, nanchor, slot, slot2, cell
    character(len=512) :: path

    config = options
    config%periodic_operator_source_depth = 1_i32
    config%periodic_operator_target_depth = 1_i32
    do idx = 1_i32, 3_i32
      center = src_pos(:, idx)
      a(:, idx) = center + [-0.02_dp, -0.01_dp, 0.0_dp]
      b(:, idx) = center + [0.02_dp, -0.01_dp, 0.0_dp]
      c(:, idx) = center + [0.0_dp, 0.02_dp, 0.0_dp]
    end do
    call test_begin('box_operator_reused_by_different_panel_mesh')
    call build_panel_plan(first, a(:, 1:2), b(:, 1:2), c(:, 1:2), config)
    path = first%periodic_cache_path
    call assert_equal_i32(first%periodic_operator_new_cell_count, first%periodic_operator_built_count, &
                          'cold plan builds exactly its occupied cells')
    config%leaf_max = 1_i32
    call build_panel_plan(second, a, b, c, config)
    call assert_true(first%periodic_cache_fingerprint == second%periodic_cache_fingerprint, &
                     'mesh and adaptive topology must not affect the fingerprint')
    call assert_equal_i32(second%periodic_root_source_count, 8_i32, 'every box cell must stay addressable')
    call assert_true(.not. second%periodic_cache_hit, 'a newly occupied cell is a partial cache miss')
    call assert_equal_i32(second%periodic_operator_new_cell_count, &
                          second%periodic_operator_built_count - first%periodic_operator_built_count, &
                          'only cells missing from the cache must be built')
    call assert_true(second%periodic_operator_new_cell_count > 0_i32, 'the third panel must occupy a new cell')
    nanchor = first%periodic_operator_anchor_count
    do slot = 1_i32, first%periodic_operator_built_count
      cell = first%periodic_operator_cells(slot)
      slot2 = second%periodic_operator_cell_slot(cell)
      call assert_true(slot2 > 0_i32, 'cells of the first mesh must stay available')
      call assert_true(maxval(abs(first%periodic_root_operator(:, :, (slot - 1)*nanchor + 1:slot*nanchor) - &
                                  second%periodic_root_operator(:, :, (slot2 - 1)*nanchor + 1:slot2*nanchor))) == 0.0_dp, &
                       'cached cell blocks must be reused verbatim')
    end do
    do idx = 1_i32, second%periodic_operator_anchor_count
      call assert_equal_i32(second%target_node_depth(second%periodic_operator_anchor_nodes(idx)), 1_i32, &
                            'box target anchors must map to the current adaptive tree')
    end do
    call build_panel_plan(third, a, b, c, config)
    call assert_true(third%periodic_cache_hit, 'merged cache must cover the second mesh')
    call assert_equal_i32(third%periodic_operator_build_count, 0_i32, 'merged cache must not rebuild')
    call destroy_plan(third)
    config%leaf_max = options%leaf_max
    call build_panel_plan(third, a(:, 1:2), b(:, 1:2), c(:, 1:2), config)
    call assert_true(third%periodic_cache_hit, 'merging must keep the cells of the first mesh')
    call destroy_plan(third)
    call test_end()

    call test_begin('shared_translations_and_reflections_match_every_pair')
    shared = first%periodic_root_operator
    shared_count = first%periodic_operator_class_count
    call precompute_periodic_root_operator(first, reuse_translations=.false.)
    call assert_true(shared_count < first%periodic_operator_class_count, 'shared path must reduce Ewald classes')
    call assert_true(maxval(abs(shared - first%periodic_root_operator)) < &
                     5.0e-12_dp*max(1.0_dp, maxval(abs(shared))), 'class sharing must agree within rounding error')
    write (*, '(a,es12.5)') 'shared/direct operator max difference=', &
      maxval(abs(shared - first%periodic_root_operator))
    call test_end()
    call destroy_plan(first)
    call destroy_plan(second)
    call delete_cache(path)
  end subroutine check_box_mesh_reuse_and_translations

  subroutine check_zero_root_moments()
    type(fmm_options_type) :: config
    type(fmm_plan_type) :: neutral
    type(fmm_state_type) :: state
    real(dp) :: positions(3, 4), charge(4)
    character(len=512) :: path

    call test_begin('zero_global_moments_do_not_disable_box_cell_correction')
    config = options
    config%order = 1_i32
    config%leaf_max = 8_i32
    config%periodic_operator_source_depth = 1_i32
    config%periodic_operator_target_depth = 1_i32
    positions(:, 1) = [0.25_dp, 0.25_dp, 0.0_dp]
    positions(:, 2) = [0.75_dp, 0.25_dp, 0.0_dp]
    positions(:, 3) = [0.25_dp, 0.75_dp, 0.0_dp]
    positions(:, 4) = [0.75_dp, 0.75_dp, 0.0_dp]
    charge = [1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp]
    call build_plan(neutral, positions, config)
    path = neutral%periodic_cache_path
    call update_state(neutral, state, charge)
    ! Analytic checkerboard: monopole and dipole vanish exactly; its quadrupole field remains.
    call assert_true(all(state%multipole(:, 1) == 0.0_dp), 'checkerboard global monopole and dipole must vanish')
    call assert_true(any(abs(state%local) > 1.0e-12_dp), 'nonzero cell moments must still produce a far local field')
    call destroy_state(state)
    call destroy_plan(neutral)
    call delete_cache(path)
    call test_end()
  end subroutine check_zero_root_moments

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
