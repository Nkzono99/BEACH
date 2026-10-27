!> periodic2 FMM コアの有限画像と seam 不変性テスト。
program test_coulomb_fmm_core_periodic
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, update_state, eval_point, &
                                  eval_potential_point, destroy_plan, destroy_state
  use bem_coulomb_fmm_periodic_ewald, only: precompute_periodic2_ewald_data, periodic2_ewald_exp_erfc, &
                                            add_periodic2_exact_ewald_correction_single_source, &
                                            add_periodic2_exact_ewald_potential_correction_single_source, &
                                            add_periodic2_exact_ewald_field_potential_single_source
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true
  implicit none

  call test_init(2)

  call test_begin('periodic2_finite_shell_seam_contract')
  call test_periodic2_finite_shell_seam_contract()
  call test_end()

  call test_begin('periodic2_ewald_teacher_far_free_axis_is_finite')
  call test_periodic2_ewald_teacher_far_free_axis()
  call test_end()

  call test_summary()

contains

  !> free 軸距離が周期長の 20 倍でも Ewald teacher が有限で、格子の一様面電場に一致する（#65）。
  !! Lx=Ly=1、既定 4 層の mode (4,4) では |k dz| = 710.9 が exp の上限を超え、旧式 exp*erfc は inf*0 = NaN だった。
  !! 小さい |dz| では k≠0 項が旧式 exp(kz)*erfc(karg0 + adz) とビット一致することも確認する。
  !! operator cache は作らず、Ewald 事前計算だけを持つ plan で teacher kernel を直接評価する。
  subroutine test_periodic2_ewald_teacher_far_free_axis()
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: alphas(2) = [0.4_dp, 0.6_dp]
    real(dp), parameter :: far_dz(2) = [20.0_dp, -20.0_dp]
    real(dp), parameter :: near_dz(5) = [-2.0_dp, -0.37_dp, 0.0_dp, 0.37_dp, 2.0_dp]
    type(fmm_plan_type), allocatable :: plan
    real(dp) :: src(3, 1), q(1), target(3), e_sep(3), e_comb(3), e_inner(3), lattice(3)
    real(dp) :: phi_sep, phi_comb, phi_inner, sheet_field, alpha, kmag, karg0, kz, adz, term
    integer :: ia, iz, h1, h2
    logical :: near_matches_old

    allocate (plan)
    plan%options%use_periodic2 = .true.
    plan%options%periodic_far_correction = 'cached_kneq0'
    plan%options%periodic_axes = [1_i32, 2_i32]
    plan%options%periodic_len = [1.0_dp, 1.0_dp]
    plan%options%periodic_image_layers = 1_i32
    plan%options%periodic_ewald_layers = 4_i32
    src(:, 1) = [0.3_dp, 0.4_dp, 0.1_dp]
    q = 1.0_dp

    do ia = 1, size(alphas)
      alpha = alphas(ia)
      plan%options%periodic_ewald_alpha = alpha
      call precompute_periodic2_ewald_data(plan)
      call assert_true(plan%periodic_ewald%ready, 'periodic2 Ewald data must be ready')

      ! mode (4, 4): exp(|k dz|) overflows, and the exact term is below the binary64 range.
      kmag = 2.0_dp*pi*sqrt(32.0_dp)
      karg0 = 0.5_dp*kmag/alpha
      kz = kmag*far_dz(1)
      adz = alpha*far_dz(1)
      term = periodic2_ewald_exp_erfc(kz, karg0, adz)
      call assert_true(term >= 0.0_dp .and. term <= exp(-kz), 'overflowing exp*erfc term must stay in [0, exp(-|k dz|)]')

      near_matches_old = .true.
      do h1 = 0, 4
        do h2 = 0, 4
          if (h1 == 0 .and. h2 == 0) cycle
          kmag = 2.0_dp*pi*sqrt(real(h1*h1 + h2*h2, dp))
          karg0 = 0.5_dp*kmag/alpha
          do iz = 1, size(near_dz)
            near_matches_old = near_matches_old .and. &
                               periodic2_ewald_exp_erfc(kmag*near_dz(iz), karg0, alpha*near_dz(iz)) == &
                               exp(kmag*near_dz(iz))*erfc(karg0 + alpha*near_dz(iz)) .and. &
                               periodic2_ewald_exp_erfc(-kmag*near_dz(iz), karg0, -alpha*near_dz(iz)) == &
                               exp(-kmag*near_dz(iz))*erfc(karg0 - alpha*near_dz(iz))
          end do
        end do
      end do
      call assert_true(near_matches_old, 'small |dz| exp*erfc terms must keep the original binary64 values')

      do iz = 1, size(far_dz)
        target = [0.6_dp, 0.2_dp, src(3, 1) + far_dz(iz)]
        e_sep = 0.0_dp
        e_comb = 0.0_dp
        phi_sep = 0.0_dp
        phi_comb = 0.0_dp
        call add_periodic2_exact_ewald_correction_single_source(plan, q(1), src(:, 1), target, e_sep)
        call add_periodic2_exact_ewald_potential_correction_single_source(plan, q(1), src(:, 1), target, phi_sep)
        call add_periodic2_exact_ewald_field_potential_single_source(plan, q(1), src(:, 1), target, e_comb, phi_comb)
        call assert_true( &
          all(ieee_is_finite(e_sep)) .and. ieee_is_finite(phi_sep) .and. &
          all(ieee_is_finite(e_comb)) .and. ieee_is_finite(phi_comb), &
          'periodic2 Ewald teacher must stay finite at |dz| = 20 Lx' &
          )

        ! The correction replaces the directly summed inner images, so adding them back gives the full lattice,
        ! whose field far from a unit cell carrying charge q is the uniform sheet field 2*pi*q/(Lx*Ly).
        call direct_periodic2(src, q, target, [1_i32, 2_i32], [1.0_dp, 1.0_dp], 1_i32, e_inner, phi_inner)
        lattice = e_sep + e_inner
        sheet_field = sign(2.0_dp*pi*q(1), far_dz(iz))
        call assert_true( &
          abs(lattice(1)) <= 1.0e-12_dp .and. abs(lattice(2)) <= 1.0e-12_dp .and. &
          abs(lattice(3) - sheet_field) <= 1.0e-12_dp*abs(sheet_field), &
          'far periodic2 lattice field must equal the uniform sheet field' &
          )
        call assert_true( &
          maxval(abs(e_comb - e_sep)) <= 1.0e-12_dp*abs(sheet_field) .and. &
          abs(phi_comb - phi_sep) <= 1.0e-12_dp*max(1.0_dp, abs(phi_sep)), &
          'combined Ewald teacher must match the separate field and potential kernels' &
          )
      end do
    end do
  end subroutine test_periodic2_ewald_teacher_far_free_axis

  subroutine test_periodic2_finite_shell_seam_contract()
    type(fmm_plan_type), allocatable :: split_plan, unwrapped_plan
    type(fmm_state_type), allocatable :: split_state, unwrapped_state
    type(fmm_options_type) :: options
    real(dp) :: split_pos(3, 6), unwrapped_pos(3, 6), q(6), queries(3, 4)
    real(dp) :: split_field(3), unwrapped_field(3), reference_field(3)
    real(dp) :: split_phi, unwrapped_phi, reference_phi
    real(dp) :: field_scale, potential_scale, field_error, potential_error, field_delta, potential_delta
    real(dp) :: max_field_error, max_potential_error, max_field_delta, max_potential_delta
    integer(i32) :: i

    allocate (split_plan, unwrapped_plan, split_state, unwrapped_state)
    split_pos(:, 1) = [0.02_dp, 0.03_dp, -0.30_dp]
    split_pos(:, 2) = [0.04_dp, 0.97_dp, 0.20_dp]
    split_pos(:, 3) = [0.96_dp, 0.02_dp, -0.10_dp]
    split_pos(:, 4) = [0.98_dp, 0.98_dp, 0.40_dp]
    split_pos(:, 5) = [0.50_dp, 0.50_dp, 0.00_dp]
    split_pos(:, 6) = [0.03_dp, 0.95_dp, -0.50_dp]
    q = [-1.0e-12_dp, 1.1e-12_dp, -0.9e-12_dp, 1.2e-12_dp, -1.3e-12_dp, 1.0e-12_dp]

    ! Same physical lattice as split_pos, explicitly unwrapped at the largest gaps.
    unwrapped_pos = split_pos
    unwrapped_pos(1, 1) = unwrapped_pos(1, 1) + 1.0_dp
    unwrapped_pos(1, 2) = unwrapped_pos(1, 2) + 1.0_dp
    unwrapped_pos(1, 6) = unwrapped_pos(1, 6) + 1.0_dp
    unwrapped_pos(2, 1) = unwrapped_pos(2, 1) + 1.0_dp
    unwrapped_pos(2, 3) = unwrapped_pos(2, 3) + 1.0_dp

    options%theta = 0.55_dp
    options%leaf_max = 2_i32
    options%order = 4_i32
    options%use_periodic2 = .true.
    options%periodic_axes = [1_i32, 2_i32]
    options%periodic_len = [1.0_dp, 1.0_dp]
    options%periodic_image_layers = 1_i32
    options%periodic_far_correction = 'none'
    options%target_box_min = [0.0_dp, 0.0_dp, -1.0_dp]
    options%target_box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    call build_plan(split_plan, split_pos, options)
    call update_state(split_plan, split_state, q)
    call build_plan(unwrapped_plan, unwrapped_pos, options)
    call update_state(unwrapped_plan, unwrapped_state, q)

    queries(:, 1) = [0.15_dp, 0.15_dp, -0.50_dp]
    queries(:, 2) = [0.85_dp, 0.85_dp, 0.30_dp]
    queries(:, 3) = [0.50_dp, 0.50_dp, -0.80_dp]
    queries(:, 4) = [0.30_dp, 0.70_dp, 0.60_dp]

    max_field_error = 0.0_dp
    max_potential_error = 0.0_dp
    max_field_delta = 0.0_dp
    max_potential_delta = 0.0_dp
    do i = 1_i32, int(size(queries, 2), i32)
      call eval_point(split_plan, split_state, queries(:, i), split_field)
      call eval_point(unwrapped_plan, unwrapped_state, queries(:, i), unwrapped_field)
      call eval_potential_point(split_plan, split_state, queries(:, i), split_phi)
      call eval_potential_point(unwrapped_plan, unwrapped_state, queries(:, i), unwrapped_phi)
      call direct_periodic2( &
        unwrapped_pos, q, queries(:, i), options%periodic_axes, options%periodic_len, &
        options%periodic_image_layers, reference_field, reference_phi &
        )
      call assert_true( &
        all(ieee_is_finite(split_field)) .and. all(ieee_is_finite(unwrapped_field)) .and. &
        all(ieee_is_finite(reference_field)) .and. ieee_is_finite(split_phi) .and. &
        ieee_is_finite(unwrapped_phi) .and. ieee_is_finite(reference_phi), &
        'periodic2 finite-shell evaluations must return finite field and potential' &
        )
      field_scale = sqrt(sum(reference_field*reference_field))
      potential_scale = abs(reference_phi)
      call assert_true(field_scale > 1.0e-16_dp, 'periodic2 reference field must be nonzero')
      call assert_true(potential_scale > 1.0e-16_dp, 'periodic2 reference potential must be nonzero')
      field_error = sqrt(sum((split_field - reference_field)**2))/field_scale
      potential_error = abs(split_phi - reference_phi)/potential_scale
      field_delta = sqrt(sum((split_field - unwrapped_field)**2))/field_scale
      potential_delta = abs(split_phi - unwrapped_phi)/potential_scale
      call assert_true( &
        ieee_is_finite(field_error) .and. ieee_is_finite(potential_error) .and. &
        ieee_is_finite(field_delta) .and. ieee_is_finite(potential_delta), &
        'periodic2 finite-shell errors and seam deltas must be finite' &
        )
      max_field_error = max(max_field_error, field_error)
      max_potential_error = max(max_potential_error, potential_error)
      max_field_delta = max(max_field_delta, field_delta)
      max_potential_delta = max(max_potential_delta, potential_delta)
    end do

    write (*, '(A,4(ES12.5,1X))') &
      'periodic2 finite-shell errors and seam deltas(field,potential)=', &
      max_field_error, max_potential_error, max_field_delta, max_potential_delta
    call assert_true(max_field_error <= 8.0e-2_dp, 'periodic2 finite-shell field error exceeds 8e-2')
    call assert_true(max_potential_error <= 8.0e-2_dp, 'periodic2 finite-shell potential error exceeds 8e-2')
    call assert_true(max_field_delta <= 1.0e-12_dp, 'periodic2 seam representation changes the field')
    call assert_true(max_potential_delta <= 1.0e-12_dp, 'periodic2 seam representation changes the potential')

    call destroy_state(unwrapped_state)
    call destroy_state(split_state)
    call destroy_plan(unwrapped_plan)
    call destroy_plan(split_plan)
  end subroutine test_periodic2_finite_shell_seam_contract

  subroutine direct_periodic2(src_pos, q, r, periodic_axes, periodic_len, nimg, e, potential)
    real(dp), intent(in) :: src_pos(:, :), q(:), r(3), periodic_len(2)
    integer(i32), intent(in) :: periodic_axes(2), nimg
    real(dp), intent(out) :: e(3), potential
    integer(i32) :: idx, img_i, img_j, axis1, axis2
    real(dp) :: shifted(3), d(3), r2, inv_r, inv_r3

    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    e = 0.0_dp
    potential = 0.0_dp
    do idx = 1_i32, int(size(q), i32)
      do img_i = -nimg, nimg
        do img_j = -nimg, nimg
          shifted = src_pos(:, idx)
          shifted(axis1) = shifted(axis1) + real(img_i, dp)*periodic_len(1)
          shifted(axis2) = shifted(axis2) + real(img_j, dp)*periodic_len(2)
          d = r - shifted
          r2 = sum(d*d)
          inv_r = 1.0_dp/sqrt(r2)
          inv_r3 = inv_r/r2
          e = e + q(idx)*inv_r3*d
          potential = potential + q(idx)*inv_r
        end do
      end do
    end do
  end subroutine direct_periodic2

end program test_coulomb_fmm_core_periodic
