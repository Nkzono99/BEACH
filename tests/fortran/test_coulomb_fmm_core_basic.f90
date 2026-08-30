!> 疎結合 Coulomb FMM コア基本テスト: P2M/M2M, free-space 精度, softening, dual tree, state reuse。
program test_coulomb_fmm_core_basic
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, update_state, eval_point, &
                                  destroy_plan, destroy_state
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_close_dp, assert_allclose_1d
  implicit none

  call test_init(5)

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

  call test_summary()

contains

  subroutine test_p2m_m2m_root_moments()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: center(3)

    allocate (plan, state)
    call make_free_sources(src_pos, q)
    q = q/1.0d-12
    options%theta = 0.55d0
    options%leaf_max = 2_i32
    options%order = 4_i32
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    center = 0.5d0*(minval(src_pos, dim=2) + maxval(src_pos, dim=2))
    call assert_true(plan%node_max_depth > 0_i32, 'root moment fixture must exercise M2M')
    call assert_root_moment(plan, state, src_pos, q, center, [0_i32, 0_i32, 0_i32], 'monopole')
    call assert_root_moment(plan, state, src_pos, q, center, [1_i32, 0_i32, 0_i32], 'x dipole')
    call assert_root_moment(plan, state, src_pos, q, center, [0_i32, 1_i32, 1_i32], 'yz quadrupole')
    call assert_root_moment(plan, state, src_pos, q, center, [4_i32, 0_i32, 0_i32], 'x quartic moment')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_p2m_m2m_root_moments

  subroutine test_free_field_accuracy()
    call assert_free_field_accuracy(0.0d0, 'free core')
  end subroutine test_free_field_accuracy

  subroutine test_softened_free_field_accuracy()
    ! Keep the softening signal larger than the accepted order-4 FMM error.
    call assert_free_field_accuracy(5.0d-1, 'softened free core')
  end subroutine test_softened_free_field_accuracy

  subroutine assert_free_field_accuracy(softening, label)
    real(dp), intent(in) :: softening
    character(len=*), intent(in) :: label
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: queries(3, 8), e_fmm(3), e_ref(3)
    real(dp) :: norm_ref, rel_err, max_rel_err
    integer(i32) :: i, max_error_index

    allocate (plan, state)
    call make_free_sources(src_pos, q)
    options%theta = 0.55d0
    options%leaf_max = 4_i32
    options%order = 4_i32
    options%softening = softening
    options%target_box_min = [-2.0d0, -2.0d0, -2.0d0]
    options%target_box_max = [2.0d0, 2.0d0, 2.0d0]
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

    call assert_true(plan%target_tree_ready, trim(label)//' must build the target tree')
    call assert_true(plan%m2l_pair_count > 0_i32, trim(label)//' fixture must exercise the M2L path')
    max_rel_err = 0.0d0
    max_error_index = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(plan, state, queries(:, i), e_fmm)
      call direct_field_free(src_pos, q, queries(:, i), options%softening, e_ref)
      call assert_true(all(ieee_is_finite(e_fmm)), trim(label)//' FMM field must be finite')
      call assert_true(all(ieee_is_finite(e_ref)), trim(label)//' reference field must be finite')
      norm_ref = sqrt(sum(e_ref*e_ref))
      call assert_true(norm_ref > 1.0d-16, trim(label)//' reference field must be nonzero')
      rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
      call assert_true(ieee_is_finite(rel_err), trim(label)//' relative error must be finite')
      if (rel_err > max_rel_err) then
        max_rel_err = rel_err
        max_error_index = i
      end if
    end do

    write (*, '(A,A,ES12.5,A,I0)') trim(label), ': max_rel_err=', max_rel_err, ', query=', max_error_index
    call assert_true(max_rel_err <= 5.0d-2, trim(label)//' FMM relative error exceeds 5e-2')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine assert_free_field_accuracy

  subroutine test_target_box_dual_tree()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q(:)
    real(dp) :: target(3), e_fmm(3), e_ref(3), norm_ref, rel_err, tolerance

    allocate (plan, state)
    call make_periodic_sources(src_pos, q)
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    call update_state(plan, state, q)

    call assert_true(plan%target_tree_ready, 'core free/use_box should enable dual target tree')
    call assert_true(plan%m2l_pair_count > 0_i32, 'dual-tree fixture must exercise the M2L path')

    target = [0.0d0, 0.45d0, 0.55d0]
    call eval_point(plan, state, target, e_fmm)
    call direct_field_free(src_pos, q, target, 0.0d0, e_ref)
    norm_ref = sqrt(sum(e_ref*e_ref))
    call assert_true(norm_ref > 1.0d-16, 'target-box boundary reference field must be nonzero')
    rel_err = sqrt(sum((e_fmm - e_ref)*(e_fmm - e_ref)))/norm_ref
    call assert_true(rel_err <= 1.0d-2, 'target-box boundary FMM relative error exceeds 1e-2')

    target = [-0.1d0, 0.45d0, 0.55d0]
    call eval_point(plan, state, target, e_fmm)
    call direct_field_free(src_pos, q, target, 0.0d0, e_ref)
    tolerance = 1.0d-12*max(maxval(abs(e_ref)), tiny(1.0d0))
    call assert_allclose_1d(e_fmm, e_ref, tolerance, 'outside target box should use exact direct fallback')

    call destroy_state(state)
    call destroy_plan(plan)
  end subroutine test_target_box_dual_tree

  subroutine test_state_update_reuse()
    type(fmm_plan_type), allocatable :: plan
    type(fmm_state_type), allocatable :: reused_state, fresh_state
    type(fmm_options_type) :: options
    real(dp), allocatable :: src_pos(:, :), q1(:), q2(:)
    real(dp) :: e_before(3), e_reused(3), e_fresh(3), tolerance

    allocate (plan, reused_state, fresh_state)
    call make_periodic_sources(src_pos, q1)
    q2 = q1
    q2(1) = 0.25d0*q1(1)
    q2(4) = -0.5d0*q1(4)
    q2(7) = 1.75d0*q1(7)
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0d0, 1.0d0]
    options%periodic_image_layers = 1_i32
    options%target_box_min = [0.0d0, 0.0d0, -1.0d0]
    options%target_box_max = [1.0d0, 1.0d0, 1.0d0]
    call build_plan(plan, src_pos, options)
    call assert_true(plan%m2l_pair_count > 0_i32, 'state reuse fixture must exercise the M2L path')

    call update_state(plan, reused_state, q1)
    call eval_point(plan, reused_state, [0.35d0, 0.45d0, 0.55d0], e_before)
    call update_state(plan, reused_state, q2)
    call eval_point(plan, reused_state, [0.35d0, 0.45d0, 0.55d0], e_reused)
    call update_state(plan, fresh_state, q2)
    call eval_point(plan, fresh_state, [0.35d0, 0.45d0, 0.55d0], e_fresh)

    tolerance = 1.0d-12*max(maxval(abs(e_fresh)), tiny(1.0d0))
    call assert_allclose_1d(e_reused, e_fresh, tolerance, 'reused state must match a fresh update for new charges')
    call assert_true( &
      sqrt(sum((e_reused - e_before)*(e_reused - e_before))) > tolerance, &
      'non-proportional charge update must change the evaluated field' &
      )

    call destroy_state(fresh_state)
    call destroy_state(reused_state)
    call destroy_plan(plan)
  end subroutine test_state_update_reuse

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

  subroutine assert_root_moment(plan, state, src_pos, q, center, alpha, label)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: src_pos(:, :), q(:), center(3)
    integer(i32), intent(in) :: alpha(3)
    character(len=*), intent(in) :: label
    integer(i32) :: alpha_idx
    real(dp) :: expected, tolerance

    alpha_idx = plan%alpha_map(alpha(1), alpha(2), alpha(3))
    expected = direct_cartesian_moment(src_pos, q, center, alpha)
    tolerance = 1.0d-11*max(1.0d0, abs(expected))
    call assert_close_dp(state%multipole(alpha_idx, 1_i32), expected, tolerance, trim(label)//' root moment mismatch')
  end subroutine assert_root_moment

  real(dp) function direct_cartesian_moment(src_pos, q, center, alpha) result(moment)
    real(dp), intent(in) :: src_pos(:, :), q(:), center(3)
    integer(i32), intent(in) :: alpha(3)
    integer(i32) :: idx, axis, power
    real(dp) :: d(3), alpha_factorial

    alpha_factorial = 1.0d0
    do axis = 1_i32, 3_i32
      do power = 2_i32, alpha(axis)
        alpha_factorial = alpha_factorial*real(power, dp)
      end do
    end do
    moment = 0.0d0
    do idx = 1_i32, int(size(q), i32)
      d = src_pos(:, idx) - center
      moment = moment + q(idx)*product(d**alpha)/alpha_factorial
    end do
  end function direct_cartesian_moment

  subroutine direct_field_free(src_pos, q, r, softening, e)
    real(dp), intent(in) :: src_pos(:, :), q(:), r(3)
    real(dp), intent(in) :: softening
    real(dp), intent(out) :: e(3)
    integer(i32) :: idx
    real(dp) :: d(3), r2, inv_r3, soft2

    e = 0.0d0
    soft2 = softening*softening
    do idx = 1_i32, int(size(q), i32)
      d = r - src_pos(:, idx)
      r2 = sum(d*d) + soft2
      if (r2 <= tiny(1.0d0)) cycle
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      e = e + q(idx)*inv_r3*d
    end do
  end subroutine direct_field_free

end program test_coulomb_fmm_core_basic
