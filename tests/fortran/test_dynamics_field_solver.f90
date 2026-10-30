!> treecode/FMM ソルバモード自動選択とパラメータチューニングの検証テスト。
program test_dynamics_field_solver
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config
  use bem_mesh, only: init_mesh
  use bem_templates, only: make_plane, make_sphere
  use bem_field_solver, only: field_solver_type, resolve_field_solver_tree_params
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, panel_surface_side_ok
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use test_panel_oracle, only: panel_sum_field, panel_sum_potential
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  call test_init(10)

  call test_begin('field_solver_auto_mode')
  call test_field_solver_auto_mode()
  call test_end()

  call test_begin('field_solver_explicit_mode_autotune')
  call test_field_solver_explicit_mode_autotune()
  call test_end()

  call test_begin('treecode_field_accuracy')
  call test_treecode_field_accuracy()
  call test_end()

  call test_begin('treecode_mixed_sign_cancellation_sweep')
  call test_treecode_mixed_sign_cancellation_sweep()
  call test_end()

  call test_begin('direct_length_normalization_field_and_potential')
  call test_direct_length_normalization_field_and_potential()
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

  call test_begin('treecode_triangle_panel_field_and_potential')
  call test_treecode_triangle_panel_field_and_potential()
  call test_end()

  call test_summary()

contains

  subroutine test_field_solver_auto_mode()
    type(mesh_type) :: mesh_small, mesh_threshold
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32), parameter :: nelem_cases(8) = [1_i32, 1499_i32, 1500_i32, 9999_i32, 10000_i32, 49999_i32, &
                                                 50000_i32, 250000_i32]
    real(dp), parameter :: theta_cases(8) = [0.40d0, 0.40d0, 0.50d0, 0.50d0, 0.58d0, 0.58d0, 0.65d0, 0.65d0]
    integer(i32), parameter :: leaf_cases(8) = [12_i32, 12_i32, 16_i32, 16_i32, 20_i32, 20_i32, 24_i32, 24_i32]
    real(dp) :: theta
    integer(i32) :: leaf_max, k

    call make_plane(mesh_small, size_x=1.0d0, size_y=1.0d0, nx=1_i32, ny=1_i32, center=[0.0d0, 0.0d0, 0.0d0])
    call mark_normal_plus(mesh_small)
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_small, sim)
    call assert_true(trim(solver%mode) == 'direct', 'auto solver should select direct for small meshes')

    call make_plane(mesh_threshold, size_x=1.0d0, size_y=1.0d0, nx=4_i32, ny=4_i32, center=[0.0d0, 0.0d0, 0.0d0])
    call mark_normal_plus(mesh_threshold)
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = mesh_threshold%nelem
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 1_i32
    call solver%init(mesh_threshold, sim)
    call assert_true(trim(solver%mode) == 'fmm', 'auto solver should select FMM at tree_min_nelem')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'auto FMM must use the tuned theta, not the unset default')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'auto FMM must use the tuned leaf_max, not the unset default')

    ! The element-count table itself is the pure resolver used by init and by the summary writer.
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 1_i32
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 1_i32
    do k = 1, size(nelem_cases)
      call resolve_field_solver_tree_params(nelem_cases(k), sim, theta, leaf_max)
      call assert_close_dp(theta, theta_cases(k), 1.0d-12, 'auto theta mismatch in the element-count table')
      call assert_equal_i32(leaf_max, leaf_cases(k), 'auto leaf_max mismatch in the element-count table')
    end do
    sim%tree_min_nelem = 64_i32
    call resolve_field_solver_tree_params(63_i32, sim, theta, leaf_max)
    call assert_close_dp(theta, 0.95d0, 0.0d0, 'direct auto selection keeps the configured theta')
    call assert_equal_i32(leaf_max, 1_i32, 'direct auto selection keeps the configured leaf_max')
  end subroutine test_field_solver_auto_mode

  subroutine test_field_solver_explicit_mode_autotune()
    type(mesh_type) :: mesh_small
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    ! Explicit tree modes autotune from the element count even below tree_min_nelem, and an explicitly
    ! set theta or leaf_max overrides only its own tuned value. The table is checked through the pure
    ! resolver that init uses; two inits check the wiring for each tree backend.
    character(len=8), parameter :: mode_cases(8) = [character(len=8) :: &
                                                    'treecode', 'treecode', 'treecode', 'treecode', &
                                                    'fmm', 'fmm', 'fmm', 'fmm']
    logical, parameter :: has_theta_cases(8) = [.false., .true., .true., .false., .false., .true., .true., .false.]
    logical, parameter :: has_leaf_cases(8) = [.false., .true., .false., .true., .false., .true., .false., .true.]
    real(dp), parameter :: theta_expected(8) = [0.40d0, 0.95d0, 0.95d0, 0.40d0, 0.50d0, 0.95d0, 0.95d0, 0.50d0]
    integer(i32), parameter :: leaf_expected(8) = [12_i32, 3_i32, 12_i32, 3_i32, 16_i32, 3_i32, 16_i32, 3_i32]
    integer(i32), parameter :: nelem_cases(8) = [528_i32, 528_i32, 528_i32, 528_i32, &
                                                 2000_i32, 2000_i32, 2000_i32, 2000_i32]
    real(dp) :: theta
    integer(i32) :: leaf_max, k

    do k = 1, size(mode_cases)
      sim = sim_config()
      sim%field_solver = mode_cases(k)
      sim%tree_min_nelem = 100000_i32
      sim%tree_theta = 0.95d0
      sim%tree_leaf_max = 3_i32
      sim%has_tree_theta = has_theta_cases(k)
      sim%has_tree_leaf_max = has_leaf_cases(k)
      call resolve_field_solver_tree_params(nelem_cases(k), sim, theta, leaf_max)
      call assert_close_dp(theta, theta_expected(k), 1.0d-12, 'explicit tree mode theta resolution mismatch')
      call assert_equal_i32(leaf_max, leaf_expected(k), 'explicit tree mode leaf_max resolution mismatch')
    end do

    call make_plane(mesh_small, size_x=1.0d0, size_y=1.0d0, nx=2_i32, ny=2_i32, center=[0.0d0, 0.0d0, 0.0d0])
    call mark_normal_plus(mesh_small)

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.93d0
    sim%tree_leaf_max = 3_i32
    sim%has_tree_theta = .true.
    call solver%init(mesh_small, sim)
    call assert_true(trim(solver%mode) == 'treecode', 'treecode mode mismatch')
    call assert_close_dp(solver%theta, 0.93d0, 1.0d-12, 'treecode init must keep the explicit theta')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'treecode init must autotune the unset leaf_max')

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%tree_theta = 0.92d0
    sim%tree_leaf_max = 5_i32
    sim%has_tree_leaf_max = .true.
    call solver%init(mesh_small, sim)
    call assert_true(trim(solver%mode) == 'fmm', 'fmm mode mismatch')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'fmm init must autotune the unset theta')
    call assert_equal_i32(solver%leaf_max, 5_i32, 'fmm init must keep the explicit leaf_max')
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
    real(dp) :: e_init(3), phi_init, phi_tree

    call make_sphere(mesh_tree, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])
    call mark_normal_plus(mesh_tree)
    mesh_tree%q_elem = 1.0d-12

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.10d0
    sim%tree_leaf_max = 6_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_tree, sim)

    ! The far target accepts tree nodes, so the comparison sees the refreshed node moments.
    r = [30.0d0, 5.0d0, -4.0d0]
    call solver%eval_e(mesh_tree, r, e_init)
    call solver%eval_potential(mesh_tree, sim, r, phi_init)
    call solver%refresh(mesh_tree)
    call solver%eval_e(mesh_tree, r, e_tree)
    call solver%eval_potential(mesh_tree, sim, r, phi_tree)
    call assert_unchanged_refresh(e_init, e_tree, phi_init, phi_tree, 'treecode')

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

      e_direct = panel_sum_field(mesh_tree, r)
      call solver%eval_e(mesh_tree, r, e_tree)
      ! MAX may drop a NaN argument, so the treecode field is checked before its error is folded in.
      call assert_true(all(ieee_is_finite(e_tree)), 'treecode returned a non-finite field')

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
      e_direct = panel_sum_field(mesh_tree, r)
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

  !> init already prepares every backend; refreshing with the same charges must leave E and phi unchanged.
  subroutine assert_unchanged_refresh(e_before, e_after, phi_before, phi_after, backend)
    real(dp), intent(in) :: e_before(3), e_after(3), phi_before, phi_after
    character(len=*), intent(in) :: backend
    real(dp) :: tolerance

    tolerance = 8.0_dp*epsilon(1.0_dp)*maxval(abs(e_before))
    call assert_allclose_1d(e_after, e_before, tolerance, trim(backend)//' E changed after an unchanged refresh')
    tolerance = 8.0_dp*epsilon(1.0_dp)*abs(phi_before)
    call assert_close_dp(phi_after, phi_before, tolerance, trim(backend)//' phi changed after an unchanged refresh')
  end subroutine assert_unchanged_refresh

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
    call mark_normal_plus(mesh_tree)

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%field_normalization = 'si'
    sim%tree_theta = 0.5d0
    sim%tree_leaf_max = 1_i32
    sim%has_tree_theta = .true.
    sim%has_tree_leaf_max = .true.
    call solver%init(mesh_tree, sim)

    call assert_true(solver%child_count(1) > 0_i32, 'treecode cancellation fixture root must be internal')
  end subroutine init_treecode_monopole_fixture

  subroutine test_direct_length_normalization_field_and_potential()
    type(mesh_type) :: mesh_direct
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_solver(3), e_init(3), norm_direct, rel_err
    real(dp) :: phi_direct, phi_solver, phi_init

    call make_sphere(mesh_direct, radius=2.0d-5, n_lon=8_i32, n_lat=4_i32, center=[5.0d-5, 0.0d0, 0.0d0])
    call mark_normal_plus(mesh_direct)
    mesh_direct%q_elem = 1.0d-13

    sim = sim_config()
    sim%field_solver = 'direct'
    sim%field_normalization = 'length'
    sim%field_length_scale = 1.0d-4
    call solver%init(mesh_direct, sim)

    r = [8.0d-5, -6.0d-5, 7.0d-5]
    call solver%eval_e(mesh_direct, r, e_init)
    call solver%eval_potential(mesh_direct, sim, r, phi_init)
    call solver%refresh(mesh_direct)
    e_direct = panel_sum_field(mesh_direct, r)
    call solver%eval_e(mesh_direct, r, e_solver)
    phi_direct = panel_sum_potential(mesh_direct, r)
    call solver%eval_potential(mesh_direct, sim, r, phi_solver)

    call assert_unchanged_refresh(e_init, e_solver, phi_init, phi_solver, 'direct')
    norm_direct = sqrt(sum(e_direct*e_direct))
    rel_err = sqrt(sum((e_solver - e_direct)*(e_solver - e_direct)))/norm_direct
    call assert_true(rel_err <= 1.0d-12, 'normalized direct E relative error exceeds 1e-12')
    call assert_close_dp(phi_solver, phi_direct, max(1.0d-12, abs(phi_direct)*1.0d-12), 'direct potential mismatch')
  end subroutine test_direct_length_normalization_field_and_potential

  subroutine test_treecode_field_length_normalization_accuracy()
    type(mesh_type) :: mesh_tree
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_tree(3), norm_direct, rel_err

    call make_sphere(mesh_tree, radius=2.0d-5, n_lon=16_i32, n_lat=8_i32, center=[5.0d-5, 0.0d0, 0.0d0])
    call mark_normal_plus(mesh_tree)
    mesh_tree%q_elem = 1.0d-13

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%field_normalization = 'length'
    sim%field_length_scale = 1.0d-4
    sim%tree_theta = 0.10d0
    sim%tree_leaf_max = 10000_i32
    sim%has_tree_leaf_max = .true.
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_tree, sim)

    r = [8.0d-5, -6.0d-5, 7.0d-5]
    e_direct = panel_sum_field(mesh_tree, r)
    call solver%eval_e(mesh_tree, r, e_tree)

    norm_direct = sqrt(sum(e_direct*e_direct))
    rel_err = sqrt(sum((e_tree - e_direct)*(e_tree - e_direct)))/norm_direct
    call assert_true(rel_err <= 1.0d-12, 'normalized treecode E relative error exceeds 1e-12')
  end subroutine test_treecode_field_length_normalization_accuracy

  subroutine test_fmm_field_length_normalization_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_fmm(3), e_init(3), norm_direct, rel_err
    real(dp) :: phi_direct, phi_fmm, phi_init

    call make_sphere(mesh_fmm, radius=2.0d-5, n_lon=16_i32, n_lat=8_i32, center=[5.0d-5, 0.0d0, 0.0d0])
    call mark_normal_plus(mesh_fmm)
    mesh_fmm%q_elem = 1.0d-13

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_normalization = 'length'
    sim%field_length_scale = 1.0d-4
    sim%tree_theta = 0.4d0
    sim%tree_leaf_max = 12_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_fmm, sim)

    r = [8.0d-5, -6.0d-5, 7.0d-5]
    call solver%eval_e(mesh_fmm, r, e_init)
    call solver%eval_potential(mesh_fmm, sim, r, phi_init)
    call solver%refresh(mesh_fmm)
    e_direct = panel_sum_field(mesh_fmm, r)
    call solver%eval_e(mesh_fmm, r, e_fmm)
    phi_direct = panel_sum_potential(mesh_fmm, r)
    call solver%eval_potential(mesh_fmm, sim, r, phi_fmm)
    call assert_unchanged_refresh(e_init, e_fmm, phi_init, phi_fmm, 'fmm')

    norm_direct = sqrt(sum(e_direct*e_direct))
    rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_direct
    call assert_true(rel_err <= 5.0d-3, 'normalized FMM E relative error exceeds 5e-3')
    call assert_close_dp( &
      phi_fmm, phi_direct, max(1.0d-12, abs(phi_direct)*5.0d-3), &
      'normalized FMM potential relative error exceeds 5e-3' &
      )
    call assert_close_dp( &
      solver%fmm_core_options%softening, 0.0_dp, 1.0d-15, &
      'triangle P0 FMM must keep internal kernel regularization disabled' &
      )
  end subroutine test_fmm_field_length_normalization_accuracy

  subroutine test_fmm_field_box_normalization_origin()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: targets(3, 2), r(3), e_direct(3), e_fmm(3), norm_direct, rel_err
    real(dp) :: phi_direct, phi_fmm
    integer(i32) :: k

    call make_sphere(mesh_fmm, radius=2.0d-5, n_lon=16_i32, n_lat=8_i32, center=[1.05d-3, 2.05d-3, -0.95d-3])
    call mark_normal_plus(mesh_fmm)
    mesh_fmm%q_elem = 1.0d-13

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_normalization = 'box'
    sim%use_box = .true.
    sim%box_min = [1.0d-3, 2.0d-3, -1.0d-3]
    sim%box_max = [1.1d-3, 2.1d-3, -0.8d-3]
    sim%tree_theta = 0.4d0
    sim%tree_leaf_max = 12_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_fmm, sim)

    ! The first target lies below box_min(2) and takes the direct fallback. The second lies inside the
    ! box, about five sphere radii from the charge, so it is evaluated through the target tree and the
    ! translated local expansions of the nonzero-origin box.
    targets(:, 1) = [1.08d-3, 1.94d-3, -0.93d-3]
    targets(:, 2) = [1.09d-3, 2.08d-3, -0.86d-3]
    do k = 1_i32, 2_i32
      r = targets(:, k)
      e_direct = panel_sum_field(mesh_fmm, r)
      call solver%eval_e(mesh_fmm, r, e_fmm)
      phi_direct = panel_sum_potential(mesh_fmm, r)
      call solver%eval_potential(mesh_fmm, sim, r, phi_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_direct
      call assert_true(rel_err <= 5.0d-3, 'box-normalized FMM E relative error exceeds 5e-3')
      call assert_close_dp(phi_fmm, phi_direct, abs(phi_direct)*5.0d-3, &
                           'box-normalized FMM potential relative error exceeds 5e-3')
    end do
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
    field_config = field_physics_config(backend='direct', normalization='si')
    periodic_config = periodic2_physics_config()
    panel_config = panel_kernel_config( &
                   kernel_id='triangle_p0_exact_direct', surface_side_policy='per_element' &
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

  subroutine test_treecode_triangle_panel_field_and_potential()
    type(mesh_type) :: mesh_panel
    type(field_solver_type) :: direct_solver = field_solver_type()
    type(field_solver_type) :: tree_solver = field_solver_type()
    type(sim_config) :: direct_sim, tree_sim
    type(field_physics_config) :: direct_field_config, tree_field_config
    type(periodic2_physics_config) :: periodic_config
    type(panel_kernel_config) :: direct_panel_config, tree_panel_config
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
    real(dp) :: near_target(3), far_target(3), charge_center(3), displacement(3)
    real(dp) :: direct_field(3), tree_field(3), expected_field(3)
    real(dp) :: direct_potential, tree_potential, expected_potential, distance, total_charge
    real(dp) :: direct_mesh_potential(2), tree_mesh_potential(2), centroid_radius
    integer(i32) :: status
    character(len=128) :: message

    v0(:, 1) = [-1.4_dp, -0.4_dp, 0.0_dp]
    v1(:, 1) = [-0.6_dp, -0.4_dp, 0.0_dp]
    v2(:, 1) = [-1.0_dp, 0.4_dp, 0.0_dp]
    v0(:, 2) = [0.6_dp, -0.4_dp, 0.0_dp]
    v1(:, 2) = [1.4_dp, -0.4_dp, 0.0_dp]
    v2(:, 2) = [1.0_dp, 0.4_dp, 0.0_dp]
    call init_mesh(mesh_panel, v0, v1, v2)
    mesh_panel%q_elem = [2.0e-12_dp, 1.0e-12_dp]
    call resolve_panel_surface_sides(mesh_panel, 'normal_plus', status, message)
    call assert_equal_i32(status, panel_surface_side_ok, 'tree panel side setup failed')

    direct_sim = sim_config()
    direct_sim%field_solver = 'direct'
    direct_sim%field_bc_mode = 'free'
    direct_sim%field_normalization = 'length'
    direct_sim%field_length_scale = 3.5_dp
    direct_field_config = field_physics_config(backend='direct', normalization='length')
    direct_panel_config = panel_kernel_config( &
                          kernel_id='triangle_p0_exact_direct', surface_side_policy='per_element' &
                          )
    periodic_config = periodic2_physics_config()
    call direct_solver%init(mesh_panel, direct_sim, direct_field_config, periodic_config, direct_panel_config)

    tree_sim = direct_sim
    tree_sim%field_solver = 'treecode'
    tree_sim%tree_theta = 1.0_dp
    tree_sim%tree_leaf_max = 1_i32
    tree_sim%has_tree_theta = .true.
    tree_sim%has_tree_leaf_max = .true.
    tree_field_config = field_physics_config(backend='treecode', normalization='length')
    tree_panel_config = panel_kernel_config( &
                        kernel_id='triangle_p0_exact_tree_near', surface_side_policy='per_element' &
                        )
    call tree_solver%init(mesh_panel, tree_sim, tree_field_config, periodic_config, tree_panel_config)

    call assert_true(trim(tree_solver%mode) == 'treecode', 'triangle panel treecode mode mismatch')
    call assert_true(tree_solver%child_count(1) > 0_i32, 'triangle panel tree root must be internal')
    centroid_radius = sqrt(sum(tree_solver%node_half_size(:, 1)**2))
    call assert_true( &
      tree_solver%node_radius(1) > centroid_radius, &
      'triangle panel tree MAC radius must include full panel vertices' &
      )

    near_target = mesh_panel%centers(:, 1)
    call direct_solver%eval_e(mesh_panel, near_target, direct_field)
    call tree_solver%eval_e(mesh_panel, near_target, tree_field)
    call assert_allclose_1d( &
      tree_field, direct_field, 1.0e-12_dp*max(1.0_dp, maxval(abs(direct_field))), &
      'triangle panel tree near field must use the exact panel kernel' &
      )
    call direct_solver%eval_potential(mesh_panel, direct_sim, near_target, direct_potential)
    call tree_solver%eval_potential(mesh_panel, tree_sim, near_target, tree_potential)
    call assert_close_dp( &
      tree_potential, direct_potential, 1.0e-12_dp*max(1.0_dp, abs(direct_potential)), &
      'triangle panel tree near potential must use the exact panel kernel' &
      )

    far_target = [0.0_dp, 10.0_dp, 5.0_dp]
    total_charge = sum(mesh_panel%q_elem)
    charge_center = matmul(mesh_panel%centers, mesh_panel%q_elem)/total_charge
    displacement = far_target - charge_center
    distance = sqrt(sum(displacement*displacement))
    expected_field = k_coulomb*total_charge*displacement/(distance**3)
    expected_potential = k_coulomb*total_charge/distance
    call tree_solver%eval_e(mesh_panel, far_target, tree_field)
    call tree_solver%eval_potential(mesh_panel, tree_sim, far_target, tree_potential)
    call assert_allclose_1d( &
      tree_field, expected_field, 1.0e-12_dp*max(1.0_dp, maxval(abs(expected_field))), &
      'triangle panel tree far field must use the node monopole' &
      )
    call assert_close_dp( &
      tree_potential, expected_potential, 1.0e-12_dp*max(1.0_dp, abs(expected_potential)), &
      'triangle panel tree far potential must use the node monopole' &
      )

    call direct_solver%compute_mesh_potential(mesh_panel, direct_sim, direct_mesh_potential)
    call tree_solver%compute_mesh_potential(mesh_panel, tree_sim, tree_mesh_potential)
    call assert_allclose_1d( &
      tree_mesh_potential, direct_mesh_potential, &
      1.0e-12_dp*max(1.0_dp, maxval(abs(direct_mesh_potential))), &
      'triangle panel tree mesh potential must preserve exact near and self terms' &
      )
  end subroutine test_treecode_triangle_panel_field_and_potential

  subroutine mark_normal_plus(mesh)
    type(mesh_type), intent(inout) :: mesh

    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals
  end subroutine mark_normal_plus

end program test_dynamics_field_solver
