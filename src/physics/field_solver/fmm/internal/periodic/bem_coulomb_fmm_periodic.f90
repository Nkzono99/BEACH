!> Coulomb FMM の periodic2 境界処理。
module bem_coulomb_fmm_periodic
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_options_type, fmm_plan_type
  implicit none
  private

  public :: has_valid_target_box
  public :: use_periodic2_m2l_root_trunc
  public :: use_periodic2_m2l_root_oracle
  public :: use_periodic2_root_operator
  public :: build_periodic_shift_values
  public :: apply_periodic2_minimum_image
  public :: wrap_periodic2_point
  public :: distance_to_source_bbox
  public :: distance_to_source_bbox_periodic
  public :: add_point_charge_images_field

contains

  !> periodic2 の target box が有効かを判定する。
  !! @param[in] options FMM 設定。
  !! @return target_box_min < target_box_max を全軸で満たすなら `.true.`。
  logical function has_valid_target_box(options)
    type(fmm_options_type), intent(in) :: options

    has_valid_target_box = all(options%target_box_max > options%target_box_min)
  end function has_valid_target_box

  !> periodic2 の far correction に trunc 版を使うか判定する。
  !! @param[in] plan FMM 計画。
  !! @return trunc 版を使うなら `.true.`。
  logical function use_periodic2_m2l_root_trunc(plan)
    type(fmm_plan_type), intent(in) :: plan

    use_periodic2_m2l_root_trunc = plan%options%use_periodic2 .and. &
                                   trim(plan%options%periodic_far_correction) == 'm2l_root_trunc'
  end function use_periodic2_m2l_root_trunc

  !> periodic2 の far correction に oracle 版を使うか判定する。
  !! @param[in] plan FMM 計画。
  !! @return oracle 版を使うなら `.true.`。
  logical function use_periodic2_m2l_root_oracle(plan)
    type(fmm_plan_type), intent(in) :: plan

    use_periodic2_m2l_root_oracle = plan%options%use_periodic2 .and. &
                                    trim(plan%options%periodic_far_correction) == 'm2l_root_oracle'
  end function use_periodic2_m2l_root_oracle

  !> periodic2 の root operator を使うか判定する。
  !! @param[in] plan FMM 計画。
  !! @return root operator を使うなら `.true.`。
  logical function use_periodic2_root_operator(plan)
    type(fmm_plan_type), intent(in) :: plan

    use_periodic2_root_operator = use_periodic2_m2l_root_trunc(plan) .or. use_periodic2_m2l_root_oracle(plan)
  end function use_periodic2_root_operator

  !> periodic2 の画像シフト値を作成する。
  !! @param[in] plan FMM 計画。
  !! @param[out] shift_axis1 軸 1 のシフト値。
  !! @param[out] shift_axis2 軸 2 のシフト値。
  !! @param[out] nshift シフト数。
  subroutine build_periodic_shift_values(plan, shift_axis1, shift_axis2, nshift)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(out) :: shift_axis1(:), shift_axis2(:)
    integer(i32), intent(out) :: nshift
    integer(i32) :: s, img

    nshift = 1_i32
    shift_axis1(1) = 0.0d0
    shift_axis2(1) = 0.0d0
    if (.not. plan%options%use_periodic2) return

    nshift = 2_i32*plan%options%periodic_image_layers + 1_i32
    if (size(shift_axis1) < nshift .or. size(shift_axis2) < nshift) then
      error stop 'periodic shift buffer is too small.'
    end if
    do s = 1_i32, nshift
      img = s - plan%options%periodic_image_layers - 1_i32
      shift_axis1(s) = real(img, dp)*plan%options%periodic_len(1)
      shift_axis2(s) = real(img, dp)*plan%options%periodic_len(2)
    end do
  end subroutine build_periodic_shift_values

  !> periodic2 の minimum image を差分ベクトルへ適用する。
  !! @param[in] plan FMM 計画。
  !! @param[inout] d 差分ベクトル。
  pure subroutine apply_periodic2_minimum_image(plan, d)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(inout) :: d(3)
    integer(i32) :: axis, k
    real(dp) :: len_k, half_len

    if (.not. plan%options%use_periodic2) return
    do k = 1_i32, 2_i32
      axis = plan%options%periodic_axes(k)
      len_k = plan%options%periodic_len(k)
      half_len = 0.5d0*len_k
      if (d(axis) > half_len) then
        d(axis) = d(axis) - len_k
      else if (d(axis) < -half_len) then
        d(axis) = d(axis) + len_k
      end if
    end do
  end subroutine apply_periodic2_minimum_image

  !> periodic2 領域内へ点座標を折り返す。
  !! @param[in] plan FMM 計画。
  !! @param[inout] p 点座標。
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

  !> 点と source BBox の距離を返す。
  !! @param[in] p 評価点。
  !! @param[in] src_center source BBox 中心。
  !! @param[in] src_half source BBox 半サイズ。
  !! @return 距離 [m]。
  real(dp) function distance_to_source_bbox(p, src_center, src_half)
    real(dp), intent(in) :: p(3), src_center(3), src_half(3)
    real(dp) :: d(3), q(3)
    integer(i32) :: axis

    d = p - src_center
    do axis = 1_i32, 3_i32
      q(axis) = max(0.0d0, abs(d(axis)) - src_half(axis))
    end do
    distance_to_source_bbox = sqrt(sum(q*q))
  end function distance_to_source_bbox

  !> periodic2 の minimum image を考慮した source BBox 距離を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] p 評価点。
  !! @param[in] src_center source BBox 中心。
  !! @param[in] src_half source BBox 半サイズ。
  !! @return 距離 [m]。
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
    distance_to_source_bbox_periodic = sqrt(sum(q*q))
  end function distance_to_source_bbox_periodic

  !> 画像電荷を足し合わせて点電荷の電場を加算する。
  !! @param[in] q 電荷量。
  !! @param[in] src 元の電荷位置。
  !! @param[in] target 評価点。
  !! @param[in] soft2 ソフトニング二乗。
  !! @param[in] axis1 画像シフト軸 1。
  !! @param[in] axis2 画像シフト軸 2。
  !! @param[in] shift_axis1 軸 1 のシフト値。
  !! @param[in] shift_axis2 軸 2 のシフト値。
  !! @param[in] nshift シフト数。
  !! @param[inout] e 電場。
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
        r2 = sum(dx*dx) + soft2
        if (r2 <= tiny(1.0d0)) cycle
        inv_r3 = 1.0d0/(sqrt(r2)*r2)
        e = e + q*inv_r3*dx
      end do
    end do
  end subroutine add_point_charge_images_field

end module bem_coulomb_fmm_periodic
