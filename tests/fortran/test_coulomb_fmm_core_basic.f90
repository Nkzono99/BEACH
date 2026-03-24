!> 疎結合 Coulomb FMM コア基本テスト: P2M/M2M, free-space 精度, softening, dual tree, state reuse, adapter。
program test_coulomb_fmm_core_basic
  use bem_constants, only: k_coulomb
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_templates, only: make_sphere
  use bem_field, only: electric_field_at
  use bem_field_solver, only: field_solver_type
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, update_state, eval_point, &
                                  destroy_plan, destroy_state
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_close_dp, assert_allclose_1d, assert_equal_i32
  implicit none

  call test_init(7)

  call test_begin('p2m_m2m_root_moments')
  call test_p2m_m2m_root_moments()
  call test_end()

  call test_begin('free_field_accuracy')
  call test_free_field_accuracy()
  call test_end()

  call test_begin('softened_free_field_accuracy')
  call test_softened_free_field_accuracy()
  call test_end()

  call test_begin('target_box_dual_tree')
  call test_target_box_dual_tree()
  call test_end()

  call test_begin('state_update_reuse')
  call test_state_update_reuse()
  call test_end()

  call test_begin('field_solver_core_adapter')
  call test_field_solver_core_adapter()
  call test_end()

  call test_begin('field_solver_core_softened_adapter')
  call test_field_solver_core_softened_adapter()
  call test_end()

  call test_summary()

contains

  subroutine test_p2m_m2m_root_moments()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:), expected(:)
    integer(i32) :: alpha_idx

    allocate (plan, state)
    call make_free_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 2_i32
    options%order = 4_i32
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    allocate (expected(plan%ncoef))
    call accumulate_direct_moments(plan, src_pos, q, plan%node_center(:, 1_i32), expected)
    do alpha_idx = 1_i32, plan%ncoef
      call assert_close_dp( &
        state%multipole(alpha_idx, 1_i32), expected(alpha_idx), &
        max(1.0d-12, 1.0d-11*abs(expected(alpha_idx))), 'root multipole moment mismatch' &
        )
    end do

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_p2m_m2m_root_moments

  subroutine test_free_field_accuracy()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 8), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    allocate (plan, state)
    call make_free_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 4_i32
    options%order = 4_i32
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [1.8d0, 1.2d0, 1.5d0]
    queries(:, 2) = [-1.7d0, 1.3d0, 1.1d0]
    queries(:, 3) = [1.6d0, -1.4d0, 1.4d0]
    queries(:, 4) = [-1.5d0, -1.6d0, 1.2d0]
    queries(:, 5) = [1.9d0, 0.8d0, -1.7d0]
    queries(:, 6) = [-1.8d0, 1.0d0, -1.5d0]
    queries(:, 7) = [1.4d0, -1.9d0, -1.3d0]
    queries(:, 8) = [-1.6d0, -1.2d0, -1.8d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_free(src_pos, q, queries(:, i), e_ref)
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 8_i32, 'free core test lost valid samples')
    call assert_true(max_rel_err <= 3.0d-3, 'free core FMM relative error exceeds 3e-3')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_free_field_accuracy

  subroutine test_softened_free_field_accuracy()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 4), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, valid_count

    allocate (plan, state)
    call make_free_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 4_i32
    options%order = 4_i32
    options%softening = 5.0d-2
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    queries(:, 1) = [1.8d0, 1.2d0, 1.5d0]
    queries(:, 2) = [-1.7d0, 1.3d0, 1.1d0]
    queries(:, 3) = [1.9d0, 0.8d0, -1.7d0]
    queries(:, 4) = [-1.6d0, -1.2d0, -1.8d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_free(src_pos, q, queries(:, i), e_ref, options%softening)
      norm_ref = sqrt(sum(e_ref*e_ref))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'softened free core test lost valid samples')
    call assert_true(max_rel_err <= 3.0d-3, 'softened free core FMM relative error exceeds 3e-3')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_softened_free_field_accuracy

  subroutine test_target_box_dual_tree()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)

    allocate (plan)
    call make_periodic_sources(src_pos, q)
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)

    call assert_true( &
      plan%target_tree_ready, &
      'core free/use_box should enable dual target tree' &
      )
    call assert_allclose_1d( &
      plan%target_node_center(:, 1_i32), [0.5d0, 0.5d0, 0.0d0], 1.0d-12, &
      'target root center mismatch' &
      )
    call assert_allclose_1d( &
      plan%target_node_half_size(:, 1_i32), [0.5d0, 0.5d0, 1.0d0], 1.0d-12, &
      'target root half-size mismatch' &
      )

    call destroy_plan(plan)
  end subroutine test_target_box_dual_tree

  subroutine test_state_update_reuse()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q1(:), q2(:)
    real(dp) :: e1(3), e2(3)
    integer(i32) :: pair_count, build_count

    allocate (plan, state)
    call make_periodic_sources(src_pos, q1)
    q2 = 2.0d0*q1
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    pair_count = plan%m2l_pair_count
    build_count = plan%m2l_build_count

    call update_state(plan, state, q1)
    call eval_point(plan, state, [0.35d0, 0.45d0, 0.55d0], e1)
    call update_state(plan, state, q2)
    call eval_point(plan, state, [0.35d0, 0.45d0, 0.55d0], e2)

    call assert_equal_i32(plan%m2l_pair_count, pair_count, 'state update should preserve M2L pair count')
    call assert_equal_i32(plan%m2l_build_count, build_count, 'state update should not rebuild the M2L cache')
    call assert_equal_i32(state%update_count, 2_i32, 'state update count mismatch')
    call assert_true(sqrt(sum((e2 - 2.0d0*e1)*(e2 - 2.0d0*e1))) <= 1.0d-10, 'field should scale linearly with charge')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_state_update_reuse

  subroutine test_field_solver_core_adapter()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_fmm(3), rel_err, norm_ref, max_rel_err
    integer(i32) :: i, valid_count

    call make_sphere(mesh_fmm, radius=0.35d0, n_lon=12_i32, n_lat=6_i32, center=[0.0d0, 0.0d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%softening = 0.0d0
    call solver%init(mesh_fmm, sim)
    call assert_true(solver%fmm_use_core, 'softening=0 free FMM should use the core path')
    call solver%refresh(mesh_fmm)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 4_i32
      r = [1.2d0 + 0.1d0*real(i - 1_i32, dp), -0.8d0 + 0.2d0*real(i - 1_i32, dp), 0.9d0 - 0.1d0*real(i - 1_i32, dp)]
      call electric_field_at(mesh_fmm, r, 0.0d0, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)
      norm_ref = sqrt(sum(e_direct*e_direct))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'core adapter test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'core adapter relative error exceeds 5e-3')
  end subroutine test_field_solver_core_adapter

  subroutine test_field_solver_core_softened_adapter()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    real(dp) :: r(3), e_direct(3), e_fmm(3), rel_err, norm_ref, max_rel_err
    integer(i32) :: i, valid_count

    call make_sphere(mesh_fmm, radius=0.35d0, n_lon=12_i32, n_lat=6_i32, center=[0.0d0, 0.0d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%softening = 1.0d-4
    call solver%init(mesh_fmm, sim)
    call assert_true(solver%fmm_use_core, 'softening>0 free FMM should use the core path')
    call solver%refresh(mesh_fmm)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 4_i32
      r = [1.2d0 + 0.1d0*real(i - 1_i32, dp), -0.8d0 + 0.2d0*real(i - 1_i32, dp), 0.9d0 - 0.1d0*real(i - 1_i32, dp)]
      call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)
      norm_ref = sqrt(sum(e_direct*e_direct))
      if (norm_ref <= 1.0d-16) cycle
      rel_err = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))/norm_ref
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 4_i32, 'softened core adapter test lost valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'softened core adapter relative error exceeds 5e-3')
  end subroutine test_field_solver_core_softened_adapter

  subroutine make_free_sources(src_pos, q)
    real(dp), allocatable, intent(out) :: src_pos(:, :)
    real(dp), allocatable, intent(out) :: q(:)
    real(dp), parameter :: axis_vals(3) = [-0.45d0, 0.0d0, 0.45d0]
    integer(i32) :: ix, iy, iz, idx

    allocate (src_pos(3, 27), q(27))
    idx = 0_i32
    do ix = 1_i32, 3_i32
      do iy = 1_i32, 3_i32
        do iz = 1_i32, 3_i32
          idx = idx + 1_i32
          src_pos(:, idx) = [axis_vals(ix), axis_vals(iy), axis_vals(iz)]
          q(idx) = real((-1_i32)**idx, dp)*1.0d-12*(1.0d0 + 0.05d0*real(mod(idx, 5_i32), dp))
        end do
      end do
    end do
  end subroutine make_free_sources

  subroutine make_periodic_sources(src_pos, q)
    real(dp), allocatable, intent(out) :: src_pos(:, :)
    real(dp), allocatable, intent(out) :: q(:)

    allocate (src_pos(3, 8), q(8))
    src_pos(:, 1) = [0.20d0, 0.20d0, -0.40d0]
    src_pos(:, 2) = [0.80d0, 0.20d0, -0.20d0]
    src_pos(:, 3) = [0.25d0, 0.75d0, 0.10d0]
    src_pos(:, 4) = [0.75d0, 0.70d0, 0.35d0]
    src_pos(:, 5) = [0.40d0, 0.35d0, -0.70d0]
    src_pos(:, 6) = [0.60d0, 0.55d0, 0.55d0]
    src_pos(:, 7) = [0.30d0, 0.60d0, 0.85d0]
    src_pos(:, 8) = [0.65d0, 0.30d0, -0.85d0]
    q = [ &
        -1.0d-12, 1.1d-12, -1.2d-12, 1.3d-12, &
        -1.4d-12, 1.5d-12, -1.6d-12, 1.3d-12 &
        ]
  end subroutine make_periodic_sources

  subroutine accumulate_direct_moments(plan, src_pos, q, center, moments)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_pos(:, :), q(:), center(3)
    real(dp), intent(out) :: moments(:)
    integer(i32) :: idx, alpha_idx
    real(dp) :: d(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    moments = 0.0d0
    do idx = 1_i32, int(size(q), i32)
      d = src_pos(:, idx) - center
      call build_powers(d, plan%options%order, xpow, ypow, zpow)
      do alpha_idx = 1_i32, plan%ncoef
        moments(alpha_idx) = moments(alpha_idx) + q(idx)*xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx)) &
                             *zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
      end do
    end do
  end subroutine accumulate_direct_moments

  subroutine direct_field_free(src_pos, q, r, e, softening)
    real(dp), intent(in) :: src_pos(:, :), q(:), r(3)
    real(dp), intent(out) :: e(3)
    real(dp), intent(in), optional :: softening
    integer(i32) :: idx
    real(dp) :: d(3), r2, inv_r3, soft2

    e = 0.0d0
    soft2 = 0.0d0
    if (present(softening)) soft2 = softening*softening
    do idx = 1_i32, int(size(q), i32)
      d = r - src_pos(:, idx)
      r2 = sum(d*d) + soft2
      if (r2 <= tiny(1.0d0)) cycle
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      e = e + q(idx)*inv_r3*d
    end do
  end subroutine direct_field_free

  subroutine build_powers(d, order, xpow, ypow, zpow)
    real(dp), intent(in) :: d(3)
    integer(i32), intent(in) :: order
    real(dp), intent(out) :: xpow(0:order), ypow(0:order), zpow(0:order)
    integer(i32) :: k

    xpow(0) = 1.0d0
    ypow(0) = 1.0d0
    zpow(0) = 1.0d0
    do k = 1_i32, order
      xpow(k) = xpow(k - 1_i32)*d(1)
      ypow(k) = ypow(k - 1_i32)*d(2)
      zpow(k) = zpow(k - 1_i32)*d(3)
    end do
  end subroutine build_powers

end program test_coulomb_fmm_core_basic
