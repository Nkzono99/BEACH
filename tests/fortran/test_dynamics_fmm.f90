!> FMM 精度テスト: free-space / periodic2 / fallback / image layers / M2L キャッシュ再利用 / モード選択。
program test_dynamics_fmm
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_templates, only: make_sphere
  use bem_field, only: electric_field_at
  use bem_field_solver, only: field_solver_type
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_correction_all_sources
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  call test_init(10)

  call test_begin('fmm_field_accuracy')
  call test_fmm_field_accuracy()
  call test_end()

  call test_begin('fmm_free_box_dual_target_accuracy')
  call test_fmm_free_box_dual_target_accuracy()
  call test_end()

  call test_begin('fmm_core_free_box_dual_target_accuracy')
  call test_fmm_core_free_box_dual_target_accuracy()
  call test_end()

  call test_begin('fmm_periodic2_field_accuracy')
  call test_fmm_periodic2_field_accuracy()
  call test_end()

  call test_begin('fmm_core_periodic2_field_accuracy')
  call test_fmm_core_periodic2_field_accuracy()
  call test_end()

  call test_begin('fmm_periodic2_fallback_accuracy')
  call test_fmm_periodic2_fallback_accuracy()
  call test_end()

  call test_begin('fmm_periodic2_image_layers_accuracy')
  call test_fmm_periodic2_image_layers_accuracy()
  call test_end()

  call test_begin('fmm_periodic2_m2l_cache_reuse')
  call test_fmm_periodic2_m2l_cache_reuse()
  call test_end()

  call test_begin('fmm_periodic2_default_m2l_root_oracle_mode')
  call test_fmm_periodic2_default_m2l_root_oracle_mode()
  call test_end()

  call test_begin('fmm_periodic2_none_disables_far_correction')
  call test_fmm_periodic2_none_disables_far_correction()
  call test_end()

  call test_summary()

contains

  subroutine test_fmm_field_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, norm_r

    call make_sphere(mesh_fmm, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%tree_theta = 0.5d0
    sim%tree_leaf_max = 16_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 271828 + 29*i
    end do
    call random_seed(put=seed)
    deallocate (seed)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 200_i32
      call random_number(r)
      r = 100.0d0*(r - 0.5d0)
      norm_r = sqrt(sum(r*r))
      if (norm_r < 20.0d0) then
        if (norm_r > 1.0d-12) then
          r = r*(20.0d0/norm_r)
        else
          r = [20.0d0, 0.0d0, 0.0d0]
        end if
      end if

      call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') 'test_fmm_field_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count >= 100_i32, 'fmm accuracy test has too few valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'fmm E relative error exceeds 5e-3')
  end subroutine test_fmm_field_accuracy

  subroutine test_fmm_free_box_dual_target_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, center_dist

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=18_i32, n_lat=10_i32, center=[0.5d0, 0.5d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call assert_true(solver%target_tree_ready, 'fmm free/use_box should enable dual target tree')
    call assert_allclose_1d( &
      solver%target_node_center(:, 1_i32), 0.5d0*(sim%box_min + sim%box_max), 1.0d-12, &
      'fmm free/use_box target root center mismatch' &
      )
    call assert_allclose_1d( &
      solver%target_node_half_size(:, 1_i32), 0.5d0*(sim%box_max - sim%box_min), 1.0d-12, &
      'fmm free/use_box target root half-size mismatch' &
      )

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 314159 + 31*i
    end do
    call random_seed(put=seed)
    deallocate (seed)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 200_i32
      call random_number(r)
      r = sim%box_min + r*(sim%box_max - sim%box_min)
      center_dist = sqrt((r(1) - 0.5d0)**2 + (r(2) - 0.5d0)**2 + r(3)**2)
      if (center_dist < 0.35d0) cycle

      call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_free_box_dual_target_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count >= 100_i32, 'fmm free/use_box test has too few valid samples')
    call assert_true(max_rel_err <= 2.5d-2, 'fmm free/use_box E relative error exceeds 2.5e-2')
  end subroutine test_fmm_free_box_dual_target_accuracy

  subroutine test_fmm_core_free_box_dual_target_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, center_dist

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=16_i32, n_lat=8_i32, center=[0.5d0, 0.5d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 0.0d0
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call assert_true(solver%fmm_use_core, 'softening=0 free FMM should use the core path')
    call assert_true(solver%target_tree_ready, 'core free/use_box should expose target tree metadata')
    call assert_allclose_1d( &
      solver%target_node_center(:, 1_i32), 0.5d0*(sim%box_min + sim%box_max), 1.0d-12, &
      'core free/use_box target root center mismatch' &
      )
    call assert_allclose_1d( &
      solver%target_node_half_size(:, 1_i32), 0.5d0*(sim%box_max - sim%box_min), 1.0d-12, &
      'core free/use_box target root half-size mismatch' &
      )

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 271829 + 17*i
    end do
    call random_seed(put=seed)
    deallocate (seed)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 120_i32
      call random_number(r)
      r = sim%box_min + r*(sim%box_max - sim%box_min)
      center_dist = sqrt((r(1) - 0.5d0)**2 + (r(2) - 0.5d0)**2 + r(3)**2)
      if (center_dist < 0.32d0) cycle

      call electric_field_at(mesh_fmm, r, 0.0d0, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count >= 80_i32, 'core free/use_box test has too few valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'core free/use_box E relative error exceeds 5e-3')
  end subroutine test_fmm_core_free_box_dual_target_accuracy

  subroutine test_fmm_periodic2_field_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_reference(mesh_fmm, solver, r, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_periodic2_field_accuracy(l2p): valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count == 6_i32, 'fmm periodic2 accuracy test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-1, 'fmm periodic2 E relative error exceeds 5e-1')
  end subroutine test_fmm_periodic2_field_accuracy

  subroutine test_fmm_core_periodic2_field_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=16_i32, n_lat=8_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 0.0d0
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call assert_true(solver%fmm_use_core, 'softening=0 periodic2 FMM should use the core path')
    call assert_true(solver%target_tree_ready, 'core periodic2 path should expose target tree metadata')
    call assert_true(solver%nleaf > 0_i32, 'core periodic2 path should expose leaf metadata')

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_reference(mesh_fmm, solver, r, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'core periodic2 accuracy test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-1, 'core periodic2 E relative error exceeds 5e-1')
  end subroutine test_fmm_core_periodic2_field_accuracy

  subroutine test_fmm_periodic2_fallback_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    queries(:, 1) = [0.15d0, 0.15d0, 1.10d0]
    queries(:, 2) = [0.85d0, 0.20d0, 1.25d0]
    queries(:, 3) = [0.20d0, 0.80d0, 1.40d0]
    queries(:, 4) = [0.75d0, 0.75d0, -1.15d0]
    queries(:, 5) = [0.55d0, 0.35d0, -1.30d0]
    queries(:, 6) = [0.35d0, 0.60d0, 1.55d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_reference(mesh_fmm, solver, r, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_periodic2_fallback_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count == 6_i32, 'fmm periodic2 fallback test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-1, 'fmm periodic2 fallback E relative error exceeds 5e-1')
  end subroutine test_fmm_periodic2_fallback_accuracy

  subroutine test_fmm_periodic2_image_layers_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 2_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_reference(mesh_fmm, solver, r, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_periodic2_image_layers_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count == 6_i32, 'fmm periodic2 image-layer test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-1, 'fmm periodic2 image-layer E relative error exceeds 5e-1')
  end subroutine test_fmm_periodic2_image_layers_accuracy

  subroutine test_fmm_periodic2_m2l_cache_reuse()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)

    mesh_fmm%q_elem = 2.0d0*mesh_fmm%q_elem
    call solver%refresh(mesh_fmm)
  end subroutine test_fmm_periodic2_m2l_cache_reuse

  subroutine test_fmm_periodic2_default_m2l_root_oracle_mode()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver_default = field_solver_type()
    type(field_solver_type) :: solver_root = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i
    real(dp) :: r(3), e_default(3), e_root(3)
    real(dp) :: max_delta_default_root

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver_default%init(mesh_fmm, sim)
    call solver_default%refresh(mesh_fmm)

    sim%field_periodic_far_correction = 'm2l_root_oracle'
    sim%field_periodic_ewald_layers = 4_i32
    call solver_root%init(mesh_fmm, sim)
    call solver_root%refresh(mesh_fmm)

    call assert_true( &
      trim(solver_default%periodic_far_correction) == 'm2l_root_oracle', &
      'periodic2 default should normalize to m2l_root_oracle' &
      )
    call assert_true( &
      trim(solver_root%periodic_far_correction) == 'm2l_root_oracle', &
      'explicit periodic2 m2l_root_oracle should be preserved' &
      )
    call assert_true( &
      trim(solver_default%fmm_core_options%periodic_far_correction) == 'm2l_root_oracle', &
      'normalized periodic2 far correction should reach FMM core options' &
      )

    max_delta_default_root = 0.0d0
    do i = 1_i32, 4_i32
      r = [0.15d0 + 0.2d0*real(i - 1_i32, dp), 0.20d0 + 0.15d0*real(i - 1_i32, dp), -0.6d0 + 0.35d0*real(i - 1_i32, dp)]
      call solver_default%eval_e(mesh_fmm, r, e_default)
      call solver_root%eval_e(mesh_fmm, r, e_root)
      max_delta_default_root = max(max_delta_default_root, sqrt(sum((e_default - e_root)*(e_default - e_root))))
    end do

    call assert_true( &
      max_delta_default_root <= 1.0d-18, &
      'default periodic2 and explicit m2l_root_oracle should agree' &
      )
  end subroutine test_fmm_periodic2_default_m2l_root_oracle_mode

  subroutine test_fmm_periodic2_none_disables_far_correction()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver_none = field_solver_type()
    type(sim_config) :: sim

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_far_correction = 'none'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver_none%init(mesh_fmm, sim)
    call solver_none%refresh(mesh_fmm)

    call assert_true(trim(solver_none%periodic_far_correction) == 'none', 'periodic2 none should be preserved')
    call assert_true( &
      trim(solver_none%fmm_core_options%periodic_far_correction) == 'none', &
      'periodic2 none should reach FMM core options' &
      )
    call assert_true( &
      .not. solver_none%fmm_core_plan%periodic_root_operator_ready, &
      'periodic2 none should not build a root-operator far correction' &
      )
    call assert_true( &
      .not. solver_none%fmm_core_plan%periodic_ewald%ready, &
      'periodic2 none should not precompute oracle Ewald data' &
      )
  end subroutine test_fmm_periodic2_none_disables_far_correction

  subroutine assign_periodic_test_dipole_charges(mesh, scale)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: scale
    integer(i32) :: elem_idx, npos, nneg
    real(dp) :: neg_scale

    npos = 0_i32
    nneg = 0_i32
    do elem_idx = 1_i32, mesh%nelem
      if (mesh%center_z(elem_idx) >= 0.0d0) then
        npos = npos + 1_i32
      else
        nneg = nneg + 1_i32
      end if
    end do

    if (npos <= 0_i32 .or. nneg <= 0_i32) then
      do elem_idx = 1_i32, mesh%nelem
        if (mod(elem_idx, 2_i32) == 0_i32) then
          mesh%q_elem(elem_idx) = -scale
        else
          mesh%q_elem(elem_idx) = scale
        end if
      end do
      return
    end if

    neg_scale = scale*real(npos, dp)/real(nneg, dp)
    do elem_idx = 1_i32, mesh%nelem
      if (mesh%center_z(elem_idx) >= 0.0d0) then
        mesh%q_elem(elem_idx) = scale
      else
        mesh%q_elem(elem_idx) = -neg_scale
      end if
    end do
  end subroutine assign_periodic_test_dipole_charges

  subroutine electric_field_at_periodic2_images(mesh, r, softening, box_min, box_max, periodic_axes, nimg, e)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: softening
    real(dp), intent(in) :: box_min(3), box_max(3)
    integer(i32), intent(in) :: periodic_axes(2)
    integer(i32), intent(in) :: nimg
    real(dp), intent(out) :: e(3)

    integer(i32) :: i, ix, iy
    integer(i32) :: axis1, axis2
    real(dp) :: l1, l2, soft2, shift1, shift2
    real(dp) :: src(3), d(3), r2, inv_r3
    real(dp) :: ex, ey, ez

    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    l1 = box_max(axis1) - box_min(axis1)
    l2 = box_max(axis2) - box_min(axis2)
    soft2 = softening*softening
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0

    do i = 1_i32, mesh%nelem
      do ix = -nimg, nimg
        do iy = -nimg, nimg
          src = [mesh%center_x(i), mesh%center_y(i), mesh%center_z(i)]
          shift1 = real(ix, dp)*l1
          shift2 = real(iy, dp)*l2
          src(axis1) = src(axis1) + shift1
          src(axis2) = src(axis2) + shift2
          d = r - src
          r2 = sum(d*d) + soft2
          inv_r3 = 1.0d0/(sqrt(r2)*r2)
          ex = ex + mesh%q_elem(i)*inv_r3*d(1)
          ey = ey + mesh%q_elem(i)*inv_r3*d(2)
          ez = ez + mesh%q_elem(i)*inv_r3*d(3)
        end do
      end do
    end do

    e(1) = k_coulomb*ex
    e(2) = k_coulomb*ey
    e(3) = k_coulomb*ez
  end subroutine electric_field_at_periodic2_images

  subroutine electric_field_at_periodic2_reference(mesh, solver, r, e)
    type(mesh_type), intent(in) :: mesh
    type(field_solver_type), intent(in) :: solver
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: e(3)
    real(dp) :: wrapped_r(3), e_core(3)
    integer(i32) :: axis

    wrapped_r = r
    do axis = 1_i32, 2_i32
      wrapped_r(solver%periodic_axes(axis)) = solver%target_box_min(solver%periodic_axes(axis)) + &
                                              modulo( &
                                              wrapped_r(solver%periodic_axes(axis)) - &
                                              solver%target_box_min(solver%periodic_axes(axis)), &
                                              solver%periodic_len(axis) &
                                              )
    end do

    call electric_field_at_periodic2_images( &
      mesh, wrapped_r, solver%softening, solver%target_box_min, solver%target_box_max, solver%periodic_axes, &
      solver%periodic_image_layers, e &
      )
    if (trim(solver%periodic_far_correction) /= 'm2l_root_oracle') return

    e_core = e/k_coulomb
    call add_periodic2_exact_ewald_correction_all_sources(solver%fmm_core_plan, solver%fmm_core_state, wrapped_r, e_core)
    e = k_coulomb*e_core
  end subroutine electric_field_at_periodic2_reference

end program test_dynamics_fmm
