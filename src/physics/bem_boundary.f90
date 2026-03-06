!> シミュレーションボックス境界（流出/反射/周期）を適用するモジュール。
module bem_boundary
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_open, bc_reflect, bc_periodic
  implicit none
  private

  public :: apply_box_boundary
contains

  !> 1ステップの更新候補位置にボックス境界条件を適用し、生存/流出状態と位置速度を更新する。
  !! @param[in] cfg 境界条件とボックス範囲を含むシミュレーション設定。
  !! @param[inout] x 境界適用対象の粒子位置 `(x,y,z)` [m]。
  !! @param[inout] v 境界適用対象の粒子速度 `(vx,vy,vz)` [m/s]。
  !! @param[inout] alive 粒子生存フラグ（流出時は `.false.` へ更新）。
  !! @param[out] escaped_boundary この呼び出しで境界流出が発生したか。
  subroutine apply_box_boundary(cfg, x, v, alive, escaped_boundary)
    type(sim_config), intent(in) :: cfg
    real(dp), intent(inout) :: x(3)
    real(dp), intent(inout) :: v(3)
    logical, intent(inout) :: alive
    logical, intent(out) :: escaped_boundary

    integer(i32) :: axis
    integer(i32) :: bc
    real(dp) :: span, eps

    escaped_boundary = .false.
    if (.not. cfg%use_box) return
    if (.not. alive) return

    eps = 1.0d-12
    do axis = 1, 3
      span = cfg%box_max(axis) - cfg%box_min(axis)
      if (span <= 0.0_dp) then
        alive = .false.
        escaped_boundary = .true.
        return
      end if

      do while (x(axis) < cfg%box_min(axis))
        bc = cfg%bc_low(axis)
        call apply_one_side_boundary(axis, bc, cfg%box_min(axis), cfg%box_max(axis), span, eps, x, v, alive, escaped_boundary)
        if ((.not. alive) .or. (.not. escaped_boundary .and. x(axis) >= cfg%box_min(axis))) exit
      end do
      if (.not. alive) return

      do while (x(axis) > cfg%box_max(axis))
        bc = cfg%bc_high(axis)
        call apply_one_side_boundary(axis, bc, cfg%box_min(axis), cfg%box_max(axis), span, eps, x, v, alive, escaped_boundary)
        if ((.not. alive) .or. (.not. escaped_boundary .and. x(axis) <= cfg%box_max(axis))) exit
      end do
      if (.not. alive) return
    end do
  end subroutine apply_box_boundary

  !> 単一軸・単一側面の境界条件を適用する内部ヘルパ。
  !! @param[in] axis 境界判定軸（1:x, 2:y, 3:z）。
  !! @param[in] bc 適用する境界条件種別（open/reflect/periodic）。
  !! @param[in] box_min 判定軸の下限座標 [m]。
  !! @param[in] box_max 判定軸の上限座標 [m]。
  !! @param[in] span 判定軸のボックス幅 `box_max - box_min` [m]。
  !! @param[in] eps 反射/周期後に境界内へ押し戻す微小量 [m]。
  !! @param[inout] x 粒子位置 `(x,y,z)` [m]。
  !! @param[inout] v 粒子速度 `(vx,vy,vz)` [m/s]。
  !! @param[inout] alive 粒子生存フラグ。
  !! @param[inout] escaped_boundary 境界流出フラグ。
  subroutine apply_one_side_boundary(axis, bc, box_min, box_max, span, eps, x, v, alive, escaped_boundary)
    integer(i32), intent(in) :: axis, bc
    real(dp), intent(in) :: box_min, box_max, span, eps
    real(dp), intent(inout) :: x(3), v(3)
    logical, intent(inout) :: alive
    logical, intent(inout) :: escaped_boundary

    select case (bc)
    case (bc_open)
      alive = .false.
      escaped_boundary = .true.
    case (bc_reflect)
      if (x(axis) < box_min) then
        x(axis) = box_min + (box_min - x(axis))
      else if (x(axis) > box_max) then
        x(axis) = box_max - (x(axis) - box_max)
      end if
      x(axis) = min(box_max - eps, max(box_min + eps, x(axis)))
      v(axis) = -v(axis)
      escaped_boundary = .false.
    case (bc_periodic)
      x(axis) = modulo(x(axis) - box_min, span) + box_min
      x(axis) = min(box_max - eps, max(box_min + eps, x(axis)))
      escaped_boundary = .false.
    case default
      alive = .false.
      escaped_boundary = .true.
    end select
  end subroutine apply_one_side_boundary

end module bem_boundary
