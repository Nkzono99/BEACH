!> treecode/FMM ソルバモード自動選択とパラメータチューニングの検証テスト。
program test_dynamics_field_solver
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_templates, only: make_plane, make_sphere
  use bem_field, only: electric_field_at, electric_potential_at
  use bem_field_solver, only: field_solver_type
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, panel_surface_side_ok
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  call test_init(12)

  call test_begin('field_solver_auto_mode')
  call test_field_solver_auto_mode()
  call test_end()

  call test_begin('field_solver_explicit_mode_autotune')
  call test_field_solver_explicit_mode_autotune()
  call test_end()

  call test_begin('treecode_field_accuracy')
  call test_treecode_field_accuracy()
  call test_end()

  call test_begin('treecode_same_sign_root_monopole')
  call test_treecode_same_sign_root_monopole()
  call test_end()

  call test_begin('treecode_mixed_sign_cancellation_sweep')
  call test_treecode_mixed_sign_cancellation_sweep()
  call test_end()

  call test_begin('zero_softening_self_singularity_skip')
  call test_zero_softening_self_singularity_skip()
  call test_end()

  call test_begin('direct_field_length_normalization_accuracy')
  call test_direct_field_length_normalization_accuracy()
  call test_end()

  call test_begin('direct_potential_length_normalization_accuracy')
  call test_direct_potential_length_normalization_accuracy()
  call test_end()

  call test_begin('treecode_field_length_normalization_accuracy')
  call test_treecode_field_length_normalization_accuracy()
  call test_end()

  call test_begin('fmm_field_length_normalization_accuracy')
  call test_fmm_field_length_normalization_accuracy()
  call test_end()

  call test_begin('fmm_field_box_normalization_origin')
  call test_fmm_field_box_normalization_origin()
  call test_end()

  call test_begin('direct_triangle_panel_contract')
  call test_direct_triangle_panel_contract()
  call test_end()

  call test_summary()

contains

  subroutine test_field_solver_auto_mode()
    type(mesh_type) :: mesh_small, mesh_large
    type(mesh_type) :: mesh_mid, mesh_dense, mesh_huge
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim

    call make_plane(mesh_small, size_x=1.0d0, size_y=1.0d0, nx=1_i32, ny=1_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_small, sim)
    call assert_true(trim(solver%mode) == 'direct', 'auto solver should select direct for small meshes')

    call make_sphere(mesh_large, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 64_i32
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 1_i32
    call solver%init(mesh_large, sim)
    call assert_true(trim(solver%mode) == 'treecode', 'auto solver should select treecode for dense meshes')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'auto theta mismatch for nelem<1500')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'auto leaf_max mismatch for nelem<1500')

    call make_plane(mesh_mid, size_x=1.0d0, size_y=1.0d0, nx=30_i32, ny=30_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 1_i32
    call solver%init(mesh_mid, sim)
    call assert_close_dp(solver%theta, 0.50d0, 1.0d-12, 'auto theta mismatch for 1500<=nelem<10000')
    call assert_equal_i32(solver%leaf_max, 16_i32, 'auto leaf_max mismatch for 1500<=nelem<10000')

    call make_plane(mesh_dense, size_x=1.0d0, size_y=1.0d0, nx=80_i32, ny=80_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 1_i32
    call solver%init(mesh_dense, sim)
    call assert_close_dp(solver%theta, 0.58d0, 1.0d-12, 'auto theta mismatch for 10000<=nelem<50000')
    call assert_equal_i32(solver%leaf_max, 20_i32, 'auto leaf_max mismatch for 10000<=nelem<50000')

    call make_plane(mesh_huge, size_x=1.0d0, size_y=1.0d0, nx=160_i32, ny=160_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 1_i32
    call solver%init(mesh_huge, sim)
    call assert_close_dp(solver%theta, 0.65d0, 1.0d-12, 'auto theta mismatch for nelem>=50000')
    call assert_equal_i32(solver%leaf_max, 24_i32, 'auto leaf_max mismatch for nelem>=50000')
  end subroutine test_field_solver_auto_mode

  subroutine test_field_solver_explicit_mode_autotune()
    type(mesh_type) :: mesh_large
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim

    call make_sphere(mesh_large, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 3_i32
    call solver%init(mesh_large, sim)
    call assert_true(trim(solver%mode) == 'treecode', 'treecode mode mismatch')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'treecode autotune theta mismatch')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'treecode autotune leaf_max mismatch')

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 3_i32
    sim%has_tree_theta = .true.
    sim%has_tree_leaf_max = .true.
    call solver%init(mesh_large, sim)
    call assert_close_dp(solver%theta, 0.95d0, 1.0d-12, 'treecode explicit theta override mismatch')
    call assert_equal_i32(solver%leaf_max, 3_i32, 'treecode explicit leaf_max override mismatch')

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%tree_theta = 0.92d0
    sim%tree_leaf_max = 5_i32
    call solver%init(mesh_large, sim)
    call assert_true(trim(solver%mode) == 'fmm', 'fmm mode mismatch')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'fmm autotune theta mismatch')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'fmm autotune leaf_max mismatch')

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%tree_theta = 0.92d0
    sim%tree_leaf_max = 5_i32
    sim%has_tree_theta = .true.
    sim%has_tree_leaf_max = .true.
    call solver%init(mesh_large, sim)
    call assert_close_dp(solver%theta, 0.92d0, 1.0d-12, 'fmm explicit theta override mismatch')
    call assert_equal_i32(solver%leaf_max, 5_i32, 'fmm explicit leaf_max override mismatch')
  end subroutine test_field_solver_explicit_mode_autotune

  subroutine test_treecode_field_accuracy()
    type(mesh_type) :: mesh_tree
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_tree(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, norm_r

    call make_sphere(mesh_tree, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])
    mesh_tree%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.10d0
    sim%tree_leaf_max = 6_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_tree, sim)
    call solver%refresh(mesh_tree)

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 314159 + 17*i
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

      call electric_field_at(mesh_tree, r, sim%softening, e_direct)
      call solver%eval_e(mesh_tree, r, e_tree)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_tree - e_direct)*(e_tree - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count >= 100_i32, 'treecode accuracy test has too few valid samples')
    call assert_true(max_rel_err <= 1.0d-3, 'treecode E relative error exceeds 1e-3')
  end subroutine test_treecode_field_accuracy

  subroutine test_treecode_same_sign_root_monopole()
    type(mesh_type) :: mesh_tree
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp), parameter :: source_charge = 1.0d-12
    real(dp) :: r(3), e_direct(3), e_tree(3), norm_direct, relative_error

    call init_treecode_monopole_fixture(mesh_tree, solver, sim)
    mesh_tree%q_elem = source_charge
    call solver%refresh(mesh_tree)

    r = [0.0d0, 100.0d0, 0.0d0]
    call electric_field_at(mesh_tree, r, sim%softening, e_direct)
    call solver%eval_e(mesh_tree, r, e_tree)

    norm_direct = sqrt(sum(e_direct*e_direct))
    call assert_true(norm_direct > 0.0d0, 'same-sign direct field must be nonzero')
    relative_error = sqrt(sum((e_tree - e_direct)*(e_tree - e_direct)))/norm_direct
    call assert_true(relative_error > 1.0d-4, 'same-sign root should retain the monopole path')
    call assert_true(relative_error < 2.0d-4, 'same-sign root monopole characterization changed')
    call assert_true(relative_error <= 1.0d-3, 'same-sign root monopole exceeds the tree accuracy contract')
  end subroutine test_treecode_same_sign_root_monopole

  subroutine test_treecode_mixed_sign_cancellation_sweep()
    type(mesh_type) :: mesh_tree
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp), parameter :: source_charge = 1.0d-12
    real(dp), parameter :: old_cancellation_tol = 1.0d-10
    real(dp), parameter :: target_ratio(2) = [0.5d0*old_cancellation_tol, 2.0d0*old_cancellation_tol]
    character(len=16), parameter :: sample_name(2) = [character(len=16) :: 'below', 'above']
    character(len=256) :: assertion_message
    integer(i32) :: sample
    real(dp) :: delta, actual_ratio, r(3), e_direct(3), e_tree(3), norm_direct, relative_error

    call init_treecode_monopole_fixture(mesh_tree, solver, sim)
    r = [0.0d0, 100.0d0, 0.0d0]

    do sample = 1_i32, 2_i32
      delta = 2.0d0*target_ratio(sample)/(1.0d0 + target_ratio(sample))
      mesh_tree%q_elem = [source_charge, -source_charge*(1.0d0 - delta)]
      actual_ratio = abs(sum(mesh_tree%q_elem))/sum(abs(mesh_tree%q_elem))
      if (sample == 1_i32) then
        call assert_true(actual_ratio < old_cancellation_tol, 'below-threshold cancellation sample crossed the old limit')
      else
        call assert_true(actual_ratio > old_cancellation_tol, 'above-threshold cancellation sample crossed the old limit')
      end if

      call solver%refresh(mesh_tree)
      call electric_field_at(mesh_tree, r, sim%softening, e_direct)
      call solver%eval_e(mesh_tree, r, e_tree)

      norm_direct = sqrt(sum(e_direct*e_direct))
      call assert_true(norm_direct > 0.0d0, 'mixed-sign direct field must be nonzero')
      relative_error = sqrt(sum((e_tree - e_direct)*(e_tree - e_direct)))/norm_direct
      write (assertion_message, '(a,a,a,es12.4,a,es12.4)') &
        'mixed-sign ', trim(sample_name(sample)), ' old-threshold sample: ratio=', actual_ratio, &
        ', relative error=', relative_error
      call assert_true(relative_error <= 1.0d-3, trim(assertion_message))
    end do
  end subroutine test_treecode_mixed_sign_cancellation_sweep

  subroutine init_treecode_monopole_fixture(mesh_tree, solver, sim)
    type(mesh_type), intent(out) :: mesh_tree
    type(field_solver_type), intent(out) :: solver
    type(sim_config), intent(out) :: sim
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)

    v0(:, 1) = [-1.0d0, -0.01d0, -0.01d0]
    v1(:, 1) = [-1.0d0, 0.02d0, -0.01d0]
    v2(:, 1) = [-1.0d0, -0.01d0, 0.02d0]
    v0(:, 2) = [1.0d0, -0.01d0, -0.01d0]
    v1(:, 2) = [1.0d0, 0.02d0, -0.01d0]
    v2(:, 2) = [1.0d0, -0.01d0, 0.02d0]
    call init_mesh(mesh_tree, v0, v1, v2)

    sim = sim_config()
    sim%softening = 0.0d0
    sim%field_solver = 'treecode'
    sim%field_normalization = 'si'
    sim%tree_theta = 0.5d0
    sim%tree_leaf_max = 1_i32
    sim%has_tree_theta = .true.
    sim%has_tree_leaf_max = .true.
    call solver%init(mesh_tree, sim)

    call assert_true(solver%child_count(1) > 0_i32, 'treecode cancellation fixture root must be internal')
  end subroutine init_treecode_monopole_fixture

  subroutine test_zero_softening_self_singularity_skip()
    type(mesh_type) :: mesh_direct
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_solver(3)

    call make_plane(mesh_direct, size_x=1.0d0, size_y=1.0d0, nx=1_i32, ny=1_i32, center=[0.0d0, 0.0d0, 0.0d0])
    mesh_direct%q_elem = 0.0d0
    mesh_direct%q_elem(1) = 1.0d-12
    r = [mesh_direct%center_x(1), mesh_direct%center_y(1), mesh_direct%center_z(1)]

    call electric_field_at(mesh_direct, r, 0.0d0, e_direct)
    call assert_true(all(e_direct == 0.0d0), 'direct self field should be skipped')

    sim = sim_config()
    sim%softening = 0.0d0
    sim%field_solver = 'direct'
    call solver%init(mesh_direct, sim)
    call solver%refresh(mesh_direct)
    call solver%eval_e(mesh_direct, r, e_solver)
    call assert_true(all(e_solver == 0.0d0), 'solver direct self field should be skipped')
  end subroutine test_zero_softening_self_singularity_skip

  subroutine test_direct_field_length_normalization_accuracy()
    type(mesh_type) :: mesh_direct
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_solver(3), norm_direct, rel_err

    call make_sphere(mesh_direct, radius=2.0d-5, n_lon=8_i32, n_lat=4_i32, center=[5.0d-5, 0.0d0, 0.0d0])
    mesh_direct%q_elem = 1.0d-13

    sim = sim_config()
    sim%softening = 1.0d-7
    sim%field_solver = 'direct'
    sim%field_normalization = 'length'
    sim%field_length_scale = 1.0d-4
    call solver%init(mesh_direct, sim)
    call solver%refresh(mesh_direct)

    r = [8.0d-5, -6.0d-5, 7.0d-5]
    call electric_field_at(mesh_direct, r, sim%softening, e_direct)
    call solver%eval_e(mesh_direct, r, e_solver)

    norm_direct = sqrt(sum(e_direct*e_direct))
    rel_err = sqrt(sum((e_solver - e_direct)*(e_solver - e_direct)))/norm_direct
    call assert_true(rel_err <= 1.0d-12, 'normalized direct E relative error exceeds 1e-12')
    call assert_close_dp( &
      solver%field_output_scale, k_coulomb/(sim%field_length_scale*sim%field_length_scale), 1.0d3, &
      'normalized direct E output scale mismatch' &
      )
  end subroutine test_direct_field_length_normalization_accuracy

  subroutine test_direct_potential_length_normalization_accuracy()
    type(mesh_type) :: mesh_direct
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), phi_direct, phi_solver

    call make_sphere(mesh_direct, radius=2.0d-5, n_lon=8_i32, n_lat=4_i32, center=[5.0d-5, 0.0d0, 0.0d0])
    mesh_direct%q_elem = 1.0d-13

    sim = sim_config()
    sim%softening = 1.0d-7
    sim%field_solver = 'direct'
    sim%field_normalization = 'length'
    sim%field_length_scale = 1.0d-4
    call solver%init(mesh_direct, sim)
    call solver%refresh(mesh_direct)

    r = [8.0d-5, -6.0d-5, 7.0d-5]
    call electric_potential_at(mesh_direct, r, sim%softening, phi_direct)
    call solver%eval_potential(mesh_direct, sim, r, phi_solver)

    call assert_close_dp(phi_solver, phi_direct, max(1.0d-12, abs(phi_direct)*1.0d-12), 'direct potential mismatch')
  end subroutine test_direct_potential_length_normalization_accuracy

  subroutine test_treecode_field_length_normalization_accuracy()
    type(mesh_type) :: mesh_tree
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_tree(3), norm_direct, rel_err

    call make_sphere(mesh_tree, radius=2.0d-5, n_lon=16_i32, n_lat=8_i32, center=[5.0d-5, 0.0d0, 0.0d0])
    mesh_tree%q_elem = 1.0d-13

    sim = sim_config()
    sim%softening = 1.0d-7
    sim%field_solver = 'treecode'
    sim%field_normalization = 'length'
    sim%field_length_scale = 1.0d-4
    sim%tree_theta = 0.10d0
    sim%tree_leaf_max = 10000_i32
    sim%has_tree_leaf_max = .true.
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_tree, sim)
    call solver%refresh(mesh_tree)

    r = [8.0d-5, -6.0d-5, 7.0d-5]
    call electric_field_at(mesh_tree, r, sim%softening, e_direct)
    call solver%eval_e(mesh_tree, r, e_tree)

    norm_direct = sqrt(sum(e_direct*e_direct))
    rel_err = sqrt(sum((e_tree - e_direct)*(e_tree - e_direct)))/norm_direct
    call assert_true(rel_err <= 1.0d-12, 'normalized treecode E relative error exceeds 1e-12')
  end subroutine test_treecode_field_length_normalization_accuracy

  subroutine test_fmm_field_length_normalization_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_fmm(3), norm_direct, rel_err

    call make_sphere(mesh_fmm, radius=2.0d-5, n_lon=16_i32, n_lat=8_i32, center=[5.0d-5, 0.0d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-13

    sim = sim_config()
    sim%softening = 1.0d-7
    sim%field_solver = 'fmm'
    sim%field_normalization = 'length'
    sim%field_length_scale = 1.0d-4
    sim%tree_theta = 0.4d0
    sim%tree_leaf_max = 12_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    r = [8.0d-5, -6.0d-5, 7.0d-5]
    call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
    call solver%eval_e(mesh_fmm, r, e_fmm)

    norm_direct = sqrt(sum(e_direct*e_direct))
    rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_direct
    call assert_true(rel_err <= 5.0d-3, 'normalized FMM E relative error exceeds 5e-3')
    call assert_close_dp( &
      solver%fmm_core_options%softening, sim%softening/sim%field_length_scale, 1.0d-15, &
      'normalized FMM softening mismatch' &
      )
  end subroutine test_fmm_field_length_normalization_accuracy

  subroutine test_fmm_field_box_normalization_origin()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_fmm(3), norm_direct, rel_err

    call make_sphere(mesh_fmm, radius=2.0d-5, n_lon=16_i32, n_lat=8_i32, center=[1.05d-3, 2.05d-3, -0.95d-3])
    mesh_fmm%q_elem = 1.0d-13

    sim = sim_config()
    sim%softening = 1.0d-7
    sim%field_solver = 'fmm'
    sim%field_normalization = 'box'
    sim%use_box = .true.
    sim%box_min = [1.0d-3, 2.0d-3, -1.0d-3]
    sim%box_max = [1.1d-3, 2.1d-3, -0.8d-3]
    sim%tree_theta = 0.4d0
    sim%tree_leaf_max = 12_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    r = [1.08d-3, 1.94d-3, -0.93d-3]
    call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
    call solver%eval_e(mesh_fmm, r, e_fmm)

    norm_direct = sqrt(sum(e_direct*e_direct))
    rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_direct
    call assert_true(rel_err <= 5.0d-3, 'box-normalized FMM E relative error exceeds 5e-3')
    call assert_close_dp(solver%fmm_core_options%target_box_min(1), 0.0d0, 1.0d-15, 'box-normalized xmin mismatch')
    call assert_close_dp(solver%fmm_core_options%target_box_min(2), 0.0d0, 1.0d-15, 'box-normalized ymin mismatch')
    call assert_close_dp(solver%fmm_core_options%target_box_min(3), 0.0d0, 1.0d-15, 'box-normalized zmin mismatch')
  end subroutine test_fmm_field_box_normalization_origin

  subroutine test_direct_triangle_panel_contract()
    type(mesh_type) :: mesh_panel
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    type(field_physics_config) :: field_config
    type(periodic2_physics_config) :: periodic_config
    type(panel_kernel_config) :: panel_config
    type(panel_geometry_type) :: geometry
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1), target(3)
    real(dp) :: field(3), expected_field(3), potential, expected_potential, mesh_potential(1)
    integer(i32) :: status
    character(len=128) :: message

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(mesh_panel, v0, v1, v2)
    mesh_panel%q_elem(1) = 2.0e-12_dp
    call resolve_panel_surface_sides(mesh_panel, 'normal_plus', status, message)
    call assert_equal_i32(status, panel_surface_side_ok, 'panel side setup failed')

    sim = sim_config()
    sim%field_solver = 'direct'
    sim%field_bc_mode = 'free'
    sim%softening = 0.0_dp
    field_config = field_physics_config(backend='direct', normalization='si')
    periodic_config = periodic2_physics_config()
    panel_config = panel_kernel_config( &
                   source_model='triangle_p0', kernel_id='triangle_p0_exact_direct', &
                   surface_side_policy='per_element' &
                   )
    call solver%init(mesh_panel, sim, field_config, periodic_config, panel_config)

    target = [0.2_dp, 0.3_dp, 0.4_dp]
    call solver%eval_e(mesh_panel, target, field)
    call solver%eval_potential(mesh_panel, sim, target, potential)
    call init_panel_geometry(v0(:, 1), v1(:, 1), v2(:, 1), geometry, status)
    call panel_potential_field( &
      geometry, mesh_panel%q_elem(1), target, panel_side_principal_value, expected_potential, expected_field &
      )
    call assert_allclose_1d(field, expected_field, 1.0e-12_dp*max(1.0_dp, maxval(abs(expected_field))), &
                            'panel solver field mismatch')
    call assert_close_dp(potential, expected_potential, 1.0e-12_dp*max(1.0_dp, abs(expected_potential)), &
                         'panel solver potential mismatch')
    call solver%compute_mesh_potential(mesh_panel, sim, mesh_potential)
    call panel_potential_field( &
      geometry, mesh_panel%q_elem(1), geometry%centroid, panel_side_principal_value, expected_potential, expected_field &
      )
    call assert_close_dp(mesh_potential(1), expected_potential, 1.0e-12_dp*max(1.0_dp, abs(expected_potential)), &
                         'panel mesh self potential mismatch')
  end subroutine test_direct_triangle_panel_contract

end program test_dynamics_field_solver
