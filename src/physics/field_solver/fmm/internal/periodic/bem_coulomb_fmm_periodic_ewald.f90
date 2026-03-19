!> periodic2 build-only Ewald oracle と fallback exact correction。
module bem_coulomb_fmm_periodic_ewald
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: inv_sqrt_pi, fmm_plan_type, fmm_state_type, reset_periodic2_ewald_data
  use bem_coulomb_fmm_periodic, only: use_periodic2_m2l_root_oracle
  implicit none
  private

  real(dp), parameter :: pi_dp = acos(-1.0d0)
  real(dp), parameter :: two_pi_dp = 2.0d0*pi_dp

  public :: resolve_periodic2_ewald_alpha
  public :: precompute_periodic2_ewald_data
  public :: add_periodic2_exact_ewald_correction_single_source
  public :: add_periodic2_exact_ewald_correction_all_sources

contains

  !> periodic2 Ewald の減衰係数 alpha を決定する。
  !! @param[in] plan FMM 計画。
  !! @return 有効な alpha。
  real(dp) function resolve_periodic2_ewald_alpha(plan)
    type(fmm_plan_type), intent(in) :: plan
    real(dp) :: min_periodic_len

    resolve_periodic2_ewald_alpha = plan%options%periodic_ewald_alpha
    if (resolve_periodic2_ewald_alpha > 0.0d0) return
    min_periodic_len = min(plan%options%periodic_len(1), plan%options%periodic_len(2))
    if (min_periodic_len <= 0.0d0) then
      resolve_periodic2_ewald_alpha = 0.0d0
      return
    end if
    resolve_periodic2_ewald_alpha = 1.2d0/(real(plan%options%periodic_image_layers + 1_i32, dp)*min_periodic_len)
  end function resolve_periodic2_ewald_alpha

  !> periodic2 Ewald の事前計算データを作成する。
  !! @param[inout] plan FMM 計画。
  subroutine precompute_periodic2_ewald_data(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: img_i, img_j, h1, h2, k_idx
    real(dp) :: alpha, k1, k2, kmag

    call reset_periodic2_ewald_data(plan%periodic_ewald)
    if (.not. use_periodic2_m2l_root_oracle(plan)) return

    alpha = resolve_periodic2_ewald_alpha(plan)
    if (alpha <= 0.0d0) return

    plan%periodic_ewald%axis1 = plan%options%periodic_axes(1)
    plan%periodic_ewald%axis2 = plan%options%periodic_axes(2)
    plan%periodic_ewald%axis_free = 6_i32 - plan%periodic_ewald%axis1 - plan%periodic_ewald%axis2
    if (plan%periodic_ewald%axis_free <= 0_i32 .or. plan%periodic_ewald%axis_free > 3_i32) return

    plan%periodic_ewald%nimg = max(0_i32, plan%options%periodic_image_layers)
    plan%periodic_ewald%kmax = max(1_i32, plan%options%periodic_ewald_layers)
    plan%periodic_ewald%img_outer = plan%periodic_ewald%nimg + plan%periodic_ewald%kmax
    plan%periodic_ewald%alpha = alpha
    plan%periodic_ewald%soft2 = plan%options%softening*plan%options%softening
    plan%periodic_ewald%cell_area = plan%options%periodic_len(1)*plan%options%periodic_len(2)
    if (plan%periodic_ewald%cell_area <= 0.0d0) then
      call reset_periodic2_ewald_data(plan%periodic_ewald)
      return
    end if
    plan%periodic_ewald%k0_pref = 2.0d0*pi_dp/plan%periodic_ewald%cell_area

    plan%periodic_ewald%screen_count = (2_i32*plan%periodic_ewald%img_outer + 1_i32)**2
    plan%periodic_ewald%inner_count = (2_i32*plan%periodic_ewald%nimg + 1_i32)**2
    plan%periodic_ewald%k_count = (2_i32*plan%periodic_ewald%kmax + 1_i32)**2 - 1_i32

    allocate ( &
      plan%periodic_ewald%screen_shift1(plan%periodic_ewald%screen_count), &
      plan%periodic_ewald%screen_shift2(plan%periodic_ewald%screen_count), &
      plan%periodic_ewald%inner_shift1(plan%periodic_ewald%inner_count), &
      plan%periodic_ewald%inner_shift2(plan%periodic_ewald%inner_count) &
      )
    allocate ( &
      plan%periodic_ewald%k1(plan%periodic_ewald%k_count), plan%periodic_ewald%k2(plan%periodic_ewald%k_count), &
      plan%periodic_ewald%kmag(plan%periodic_ewald%k_count), plan%periodic_ewald%karg0(plan%periodic_ewald%k_count), &
      plan%periodic_ewald%kpref1(plan%periodic_ewald%k_count), plan%periodic_ewald%kpref2(plan%periodic_ewald%k_count), &
      plan%periodic_ewald%kprefz(plan%periodic_ewald%k_count) &
      )

    plan%periodic_ewald%screen_count = 0_i32
    do img_i = -plan%periodic_ewald%img_outer, plan%periodic_ewald%img_outer
      do img_j = -plan%periodic_ewald%img_outer, plan%periodic_ewald%img_outer
        plan%periodic_ewald%screen_count = plan%periodic_ewald%screen_count + 1_i32
        plan%periodic_ewald%screen_shift1(plan%periodic_ewald%screen_count) = real(img_i, dp)*plan%options%periodic_len(1)
        plan%periodic_ewald%screen_shift2(plan%periodic_ewald%screen_count) = real(img_j, dp)*plan%options%periodic_len(2)
      end do
    end do

    plan%periodic_ewald%inner_count = 0_i32
    do img_i = -plan%periodic_ewald%nimg, plan%periodic_ewald%nimg
      do img_j = -plan%periodic_ewald%nimg, plan%periodic_ewald%nimg
        plan%periodic_ewald%inner_count = plan%periodic_ewald%inner_count + 1_i32
        plan%periodic_ewald%inner_shift1(plan%periodic_ewald%inner_count) = real(img_i, dp)*plan%options%periodic_len(1)
        plan%periodic_ewald%inner_shift2(plan%periodic_ewald%inner_count) = real(img_j, dp)*plan%options%periodic_len(2)
      end do
    end do

    k_idx = 0_i32
    do h1 = -plan%periodic_ewald%kmax, plan%periodic_ewald%kmax
      k1 = two_pi_dp*real(h1, dp)/plan%options%periodic_len(1)
      do h2 = -plan%periodic_ewald%kmax, plan%periodic_ewald%kmax
        if (h1 == 0_i32 .and. h2 == 0_i32) cycle
        k2 = two_pi_dp*real(h2, dp)/plan%options%periodic_len(2)
        kmag = sqrt(k1*k1 + k2*k2)
        if (kmag <= tiny(1.0d0)) cycle
        k_idx = k_idx + 1_i32
        plan%periodic_ewald%k1(k_idx) = k1
        plan%periodic_ewald%k2(k_idx) = k2
        plan%periodic_ewald%kmag(k_idx) = kmag
        plan%periodic_ewald%karg0(k_idx) = 0.5d0*kmag/alpha
        plan%periodic_ewald%kpref1(k_idx) = pi_dp*k1/(plan%periodic_ewald%cell_area*kmag)
        plan%periodic_ewald%kpref2(k_idx) = pi_dp*k2/(plan%periodic_ewald%cell_area*kmag)
        plan%periodic_ewald%kprefz(k_idx) = pi_dp/plan%periodic_ewald%cell_area
      end do
    end do
    plan%periodic_ewald%k_count = k_idx
    plan%periodic_ewald%ready = .true.
    plan%options%periodic_ewald_alpha = alpha
  end subroutine precompute_periodic2_ewald_data

  !> 全ソース分の periodic2 Ewald 補正を加算する。
  !! @param[in] plan FMM 計画。
  !! @param[in] state ソース電荷を含む state。
  !! @param[in] target 評価点。
  !! @param[inout] e 電場。
  subroutine add_periodic2_exact_ewald_correction_all_sources(plan, state, target, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: e(3)
    integer(i32) :: idx

    if (.not. plan%periodic_ewald%ready) return
    do idx = 1_i32, plan%nsrc
      call add_periodic2_exact_ewald_correction_single_source(plan, state%src_q(idx), plan%src_pos(:, idx), target, e)
    end do
  end subroutine add_periodic2_exact_ewald_correction_all_sources

  !> 1 粒子分の periodic2 Ewald 補正を加算する。
  !! @param[in] plan FMM 計画。
  !! @param[in] q 電荷量。
  !! @param[in] src ソース位置。
  !! @param[in] target 評価点。
  !! @param[inout] e 電場。
  subroutine add_periodic2_exact_ewald_correction_single_source(plan, q, src, target, e)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: q
    real(dp), intent(in) :: src(3), target(3)
    real(dp), intent(inout) :: e(3)

    if (.not. plan%periodic_ewald%ready) return
    if (abs(q) <= tiny(1.0d0)) return

    call add_exact_periodic2_real_space_correction(plan, q, src, target, e)
    call add_exact_periodic2_reciprocal_space_correction(plan, q, src, target, e)
    call add_exact_periodic2_k0_correction(plan, q, src, target, e)
  end subroutine add_periodic2_exact_ewald_correction_single_source

  subroutine add_exact_periodic2_real_space_correction(plan, q, src, target, e)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: q, src(3), target(3)
    real(dp), intent(inout) :: e(3)
    integer(i32) :: shift_idx
    real(dp) :: shifted(3)

    do shift_idx = 1_i32, plan%periodic_ewald%screen_count
      shifted = src
      shifted(plan%periodic_ewald%axis1) = shifted(plan%periodic_ewald%axis1) + plan%periodic_ewald%screen_shift1(shift_idx)
      shifted(plan%periodic_ewald%axis2) = shifted(plan%periodic_ewald%axis2) + plan%periodic_ewald%screen_shift2(shift_idx)
      call add_screened_point_charge(q, target, shifted, plan%periodic_ewald%alpha, e)
    end do

    do shift_idx = 1_i32, plan%periodic_ewald%inner_count
      shifted = src
      shifted(plan%periodic_ewald%axis1) = shifted(plan%periodic_ewald%axis1) + plan%periodic_ewald%inner_shift1(shift_idx)
      shifted(plan%periodic_ewald%axis2) = shifted(plan%periodic_ewald%axis2) + plan%periodic_ewald%inner_shift2(shift_idx)
      call add_softened_point_charge(-q, target, shifted, plan%periodic_ewald%soft2, e)
    end do
  end subroutine add_exact_periodic2_real_space_correction

  subroutine add_exact_periodic2_reciprocal_space_correction(plan, q, src, target, e)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: q, src(3), target(3)
    real(dp), intent(inout) :: e(3)
    integer(i32) :: k_idx
    real(dp) :: theta, dz, phase_sin, phase_cos, term_p, term_m, pair_sum

    dz = target(plan%periodic_ewald%axis_free) - src(plan%periodic_ewald%axis_free)
    do k_idx = 1_i32, plan%periodic_ewald%k_count
      theta = plan%periodic_ewald%k1(k_idx)*(target(plan%periodic_ewald%axis1) - src(plan%periodic_ewald%axis1)) + &
              plan%periodic_ewald%k2(k_idx)*(target(plan%periodic_ewald%axis2) - src(plan%periodic_ewald%axis2))
      phase_sin = sin(theta)
      phase_cos = cos(theta)
      term_p = exp(plan%periodic_ewald%kmag(k_idx)*dz)*erfc(plan%periodic_ewald%karg0(k_idx) + plan%periodic_ewald%alpha*dz)
      term_m = exp(-plan%periodic_ewald%kmag(k_idx)*dz)*erfc(plan%periodic_ewald%karg0(k_idx) - plan%periodic_ewald%alpha*dz)
      pair_sum = term_p + term_m
      e(plan%periodic_ewald%axis1) = e(plan%periodic_ewald%axis1) + q*plan%periodic_ewald%kpref1(k_idx)*phase_sin*pair_sum
      e(plan%periodic_ewald%axis2) = e(plan%periodic_ewald%axis2) + q*plan%periodic_ewald%kpref2(k_idx)*phase_sin*pair_sum
      e(plan%periodic_ewald%axis_free) = e(plan%periodic_ewald%axis_free) &
                                         + q*plan%periodic_ewald%kprefz(k_idx)*phase_cos*(term_m - term_p)
    end do
  end subroutine add_exact_periodic2_reciprocal_space_correction

  subroutine add_exact_periodic2_k0_correction(plan, q, src, target, e)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: q, src(3), target(3)
    real(dp), intent(inout) :: e(3)
    real(dp) :: dz

    dz = target(plan%periodic_ewald%axis_free) - src(plan%periodic_ewald%axis_free)
    e(plan%periodic_ewald%axis_free) = e(plan%periodic_ewald%axis_free) &
                                       + plan%periodic_ewald%k0_pref*q*erf(plan%periodic_ewald%alpha*dz)
  end subroutine add_exact_periodic2_k0_correction

  subroutine add_screened_point_charge(q, target, source, alpha, e)
    real(dp), intent(in) :: q
    real(dp), intent(in) :: target(3), source(3)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: e(3)
    real(dp) :: dx(3), r2, rmag, inv_r, inv_r2, inv_r3
    real(dp) :: ar, screen, gaussian, pref

    dx = target - source
    r2 = sum(dx*dx)
    if (r2 <= tiny(1.0d0)) return
    rmag = sqrt(r2)
    inv_r = 1.0d0/rmag
    inv_r2 = inv_r*inv_r
    inv_r3 = inv_r2*inv_r
    ar = alpha*rmag
    screen = erfc(ar)
    gaussian = exp(-(ar*ar))
    pref = q*(screen*inv_r3 + 2.0d0*alpha*inv_sqrt_pi*gaussian*inv_r2)
    e = e + pref*dx
  end subroutine add_screened_point_charge

  subroutine add_softened_point_charge(q, target, source, soft2, e)
    real(dp), intent(in) :: q
    real(dp), intent(in) :: target(3), source(3)
    real(dp), intent(in) :: soft2
    real(dp), intent(inout) :: e(3)
    real(dp) :: dx(3), r2, inv_r3

    dx = target - source
    r2 = sum(dx*dx) + soft2
    if (r2 <= tiny(1.0d0)) return
    inv_r3 = 1.0d0/(sqrt(r2)*r2)
    e = e + q*inv_r3*dx
  end subroutine add_softened_point_charge

end module bem_coulomb_fmm_periodic_ewald
