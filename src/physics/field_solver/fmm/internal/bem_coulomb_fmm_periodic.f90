!> Coulomb FMM の periodic2 境界処理。
module bem_coulomb_fmm_periodic
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: inv_sqrt_pi, fmm_options_type, fmm_plan_type
  implicit none
  private

  public :: has_valid_target_box
  public :: use_periodic2_ewald_like
  public :: build_periodic_shift_values
  public :: apply_periodic2_minimum_image
  public :: wrap_periodic2_point
  public :: distance_to_source_bbox
  public :: distance_to_source_bbox_periodic
  public :: prepare_periodic2_ewald
  public :: add_screened_shifted_node_images
  public :: add_point_charge_images_field

contains

  logical function has_valid_target_box(options)
    type(fmm_options_type), intent(in) :: options

    has_valid_target_box = all(options%target_box_max > options%target_box_min)
  end function has_valid_target_box

  logical function use_periodic2_ewald_like(plan)
    type(fmm_plan_type), intent(in) :: plan

    use_periodic2_ewald_like = plan%options%use_periodic2 .and. trim(plan%options%periodic_far_correction) == 'ewald_like'
  end function use_periodic2_ewald_like

  subroutine build_periodic_shift_values(plan, shift_axis1, shift_axis2, nshift)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(out) :: shift_axis1(:), shift_axis2(:)
    integer(i32), intent(out) :: nshift
    integer(i32) :: s, img

    nshift = 1_i32
    shift_axis1(1) = 0.0d0
    shift_axis2(1) = 0.0d0
    if (.not. plan%options%use_periodic2) return

    nshift = 2_i32 * plan%options%periodic_image_layers + 1_i32
    if (size(shift_axis1) < nshift .or. size(shift_axis2) < nshift) then
      error stop 'periodic shift buffer is too small.'
    end if
    do s = 1_i32, nshift
      img = s - plan%options%periodic_image_layers - 1_i32
      shift_axis1(s) = real(img, dp) * plan%options%periodic_len(1)
      shift_axis2(s) = real(img, dp) * plan%options%periodic_len(2)
    end do
  end subroutine build_periodic_shift_values

  pure subroutine apply_periodic2_minimum_image(plan, d)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(inout) :: d(3)
    integer(i32) :: axis, k
    real(dp) :: len_k, half_len

    if (.not. plan%options%use_periodic2) return
    do k = 1_i32, 2_i32
      axis = plan%options%periodic_axes(k)
      len_k = plan%options%periodic_len(k)
      half_len = 0.5d0 * len_k
      if (d(axis) > half_len) then
        d(axis) = d(axis) - len_k
      else if (d(axis) < -half_len) then
        d(axis) = d(axis) + len_k
      end if
    end do
  end subroutine apply_periodic2_minimum_image

  subroutine wrap_periodic2_point(plan, p)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(inout) :: p(3)
    integer(i32) :: k, axis

    if (.not. plan%options%use_periodic2) return
    do k = 1_i32, 2_i32
      axis = plan%options%periodic_axes(k)
      p(axis) = plan%options%target_box_min(axis) &
                + modulo(p(axis) - plan%options%target_box_min(axis), plan%options%periodic_len(k))
    end do
  end subroutine wrap_periodic2_point

  real(dp) function distance_to_source_bbox(p, src_center, src_half)
    real(dp), intent(in) :: p(3), src_center(3), src_half(3)
    real(dp) :: d(3), q(3)
    integer(i32) :: axis

    d = p - src_center
    do axis = 1_i32, 3_i32
      q(axis) = max(0.0d0, abs(d(axis)) - src_half(axis))
    end do
    distance_to_source_bbox = sqrt(sum(q * q))
  end function distance_to_source_bbox

  real(dp) function distance_to_source_bbox_periodic(plan, p, src_center, src_half)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: p(3), src_center(3), src_half(3)
    real(dp) :: d(3), q(3)
    integer(i32) :: axis

    d = p - src_center
    call apply_periodic2_minimum_image(plan, d)
    do axis = 1_i32, 3_i32
      q(axis) = max(0.0d0, abs(d(axis)) - src_half(axis))
    end do
    distance_to_source_bbox_periodic = sqrt(sum(q * q))
  end function distance_to_source_bbox_periodic

  logical function prepare_periodic2_ewald(plan, alpha, nimg, img_outer, axis1, axis2)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(out) :: alpha
    integer(i32), intent(out) :: nimg, img_outer, axis1, axis2

    prepare_periodic2_ewald = .false.
    alpha = 0.0d0
    nimg = 0_i32
    img_outer = 0_i32
    axis1 = 0_i32
    axis2 = 0_i32

    if (.not. use_periodic2_ewald_like(plan)) return
    if (plan%options%periodic_ewald_layers <= 0_i32) return

    alpha = plan%options%periodic_ewald_alpha
    if (alpha <= 0.0d0) return

    nimg = plan%options%periodic_image_layers
    img_outer = nimg + plan%options%periodic_ewald_layers
    if (img_outer <= nimg) return

    axis1 = plan%options%periodic_axes(1)
    axis2 = plan%options%periodic_axes(2)
    if (axis1 <= 0_i32 .or. axis2 <= 0_i32) return

    prepare_periodic2_ewald = .true.
  end function prepare_periodic2_ewald

  subroutine add_screened_shifted_node_images(q, src, target, axis1, axis2, periodic_len, near_img, outer_img, alpha, e)
    real(dp), intent(in) :: q
    real(dp), intent(in) :: src(3)
    real(dp), intent(in) :: target(3)
    integer(i32), intent(in) :: axis1, axis2
    real(dp), intent(in) :: periodic_len(2)
    integer(i32), intent(in) :: near_img, outer_img
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: e(3)
    integer(i32) :: img_i, img_j
    real(dp) :: shifted(3)

    if (abs(q) <= tiny(1.0d0)) return
    do img_i = -outer_img, outer_img
      do img_j = -outer_img, outer_img
        if (abs(img_i) <= near_img .and. abs(img_j) <= near_img) cycle
        shifted = src
        shifted(axis1) = shifted(axis1) + real(img_i, dp) * periodic_len(1)
        shifted(axis2) = shifted(axis2) + real(img_j, dp) * periodic_len(2)
        call add_screened_point_charge(q, target, shifted, alpha, e)
      end do
    end do
  end subroutine add_screened_shifted_node_images

  subroutine add_point_charge_images_field(q, src, target, soft2, axis1, axis2, shift_axis1, shift_axis2, nshift, e)
    real(dp), intent(in) :: q, src(3), target(3)
    real(dp), intent(in) :: soft2
    integer(i32), intent(in) :: axis1, axis2, nshift
    real(dp), intent(in) :: shift_axis1(:), shift_axis2(:)
    real(dp), intent(inout) :: e(3)
    integer(i32) :: img_i, img_j
    real(dp) :: shifted(3), dx(3), r2, inv_r3

    if (abs(q) <= tiny(1.0d0)) return
    do img_i = 1_i32, nshift
      do img_j = 1_i32, nshift
        shifted = src
        if (axis1 > 0_i32) shifted(axis1) = shifted(axis1) + shift_axis1(img_i)
        if (axis2 > 0_i32) shifted(axis2) = shifted(axis2) + shift_axis2(img_j)
        dx = target - shifted
        r2 = sum(dx * dx) + soft2
        if (r2 <= tiny(1.0d0)) cycle
        inv_r3 = 1.0d0 / (sqrt(r2) * r2)
        e = e + q * inv_r3 * dx
      end do
    end do
  end subroutine add_point_charge_images_field

  subroutine add_screened_point_charge(q, target, source, alpha, e)
    real(dp), intent(in) :: q
    real(dp), intent(in) :: target(3), source(3)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: e(3)
    real(dp) :: dx(3), r2, rmag, inv_r, inv_r2, inv_r3
    real(dp) :: ar, screen, gaussian, pref

    if (abs(q) <= tiny(1.0d0)) return
    dx = target - source
    r2 = sum(dx * dx)
    if (r2 <= tiny(1.0d0)) return

    rmag = sqrt(r2)
    inv_r = 1.0d0 / rmag
    inv_r2 = inv_r * inv_r
    inv_r3 = inv_r2 * inv_r
    ar = alpha * rmag
    screen = erfc(ar)
    gaussian = exp(-(ar * ar))
    pref = q * (screen * inv_r3 + 2.0d0 * alpha * inv_sqrt_pi * gaussian * inv_r2)
    e = e + pref * dx
  end subroutine add_screened_point_charge

end module bem_coulomb_fmm_periodic
