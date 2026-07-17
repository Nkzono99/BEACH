!> シミュレーションボックス境界（流出/反射/周期）を適用するモジュール。
module bem_boundary
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_open, bc_reflect, bc_periodic
  implicit none
  private

  integer(i32), parameter, public :: boundary_event_ok = 0_i32
  integer(i32), parameter, public :: boundary_event_invalid_geometry = 1_i32
  integer(i32), parameter :: boundary_face_all = 63_i32

  type, public :: boundary_event_type
    logical :: has_event = .false.
    real(dp) :: fraction = 0.0_dp
    integer(i32) :: face_mask = 0_i32
    integer(i32) :: face_bc(6) = -1_i32
  end type boundary_event_type

  public :: apply_box_boundary
  public :: find_first_boundary_event
  public :: apply_escape_reflect_periodic_event
contains

  !> 閉じたbox内の始点から候補終点へ向かう線分について、最初のbox面交差を返す。
  subroutine find_first_boundary_event(cfg, x0, x1, event, status)
    type(sim_config), intent(in) :: cfg
    real(dp), intent(in) :: x0(3), x1(3)
    type(boundary_event_type), intent(out) :: event
    integer(i32), intent(out) :: status

    real(dp) :: delta(3), candidate_fraction(3), first_fraction, tie_tolerance
    integer(i32) :: axis, candidate_face(3), n_candidate

    event = boundary_event_type()
    status = boundary_event_ok
    if (.not. all(ieee_is_finite(x0)) .or. .not. all(ieee_is_finite(x1))) then
      status = boundary_event_invalid_geometry
      return
    end if
    if (.not. cfg%use_box) return
    if (.not. valid_event_box(cfg)) then
      status = boundary_event_invalid_geometry
      return
    end if
    if (any(x0 < cfg%box_min) .or. any(x0 > cfg%box_max)) then
      status = boundary_event_invalid_geometry
      return
    end if

    delta = x1 - x0
    if (.not. all(ieee_is_finite(delta))) then
      status = boundary_event_invalid_geometry
      return
    end if

    event%face_bc = boundary_face_conditions(cfg)
    candidate_fraction = 0.0_dp
    candidate_face = 0_i32
    n_candidate = 0_i32
    do axis = 1_i32, 3_i32
      if (delta(axis) > 0.0_dp .and. x1(axis) >= cfg%box_max(axis)) then
        n_candidate = n_candidate + 1_i32
        candidate_fraction(n_candidate) = (cfg%box_max(axis) - x0(axis))/delta(axis)
        candidate_face(n_candidate) = 2_i32*axis
      else if (delta(axis) < 0.0_dp .and. x1(axis) <= cfg%box_min(axis)) then
        n_candidate = n_candidate + 1_i32
        candidate_fraction(n_candidate) = (cfg%box_min(axis) - x0(axis))/delta(axis)
        candidate_face(n_candidate) = 2_i32*axis - 1_i32
      end if
    end do
    if (n_candidate == 0_i32) return
    if (.not. all(ieee_is_finite(candidate_fraction(1:n_candidate))) .or. &
        any(candidate_fraction(1:n_candidate) < 0.0_dp) .or. &
        any(candidate_fraction(1:n_candidate) > 1.0_dp)) then
      event = boundary_event_type()
      status = boundary_event_invalid_geometry
      return
    end if

    first_fraction = minval(candidate_fraction(1:n_candidate))
    if (first_fraction == 0.0_dp) first_fraction = 0.0_dp
    tie_tolerance = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(first_fraction))
    event%has_event = .true.
    event%fraction = first_fraction
    do axis = 1_i32, n_candidate
      if (abs(candidate_fraction(axis) - first_fraction) <= tie_tolerance) then
        event%face_mask = ior(event%face_mask, shiftl(1_i32, candidate_face(axis) - 1_i32))
      end if
    end do
  end subroutine find_first_boundary_event

  !> 最初のbox eventへescape/reflect/periodicを軸順序に依存せず適用する。
  subroutine apply_escape_reflect_periodic_event(cfg, event, x, v, alive, escaped, status)
    type(sim_config), intent(in) :: cfg
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(inout) :: x(3), v(3)
    logical, intent(inout) :: alive
    logical, intent(out) :: escaped
    integer(i32), intent(out) :: status

    real(dp) :: x_work(3), v_work(3)
    integer(i32) :: axis, low_face, high_face, low_bc, high_bc
    logical :: alive_work, has_open

    escaped = .false.
    status = boundary_event_ok
    if (.not. cfg%use_box .or. .not. valid_event_box(cfg)) then
      status = boundary_event_invalid_geometry
      return
    end if
    if (trim(cfg%open_boundary_model) /= 'escape') then
      status = boundary_event_invalid_geometry
      return
    end if
    if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(v))) then
      status = boundary_event_invalid_geometry
      return
    end if
    if (.not. event%has_event) then
      if (event%face_mask /= 0_i32 .or. event%fraction /= 0.0_dp) status = boundary_event_invalid_geometry
      return
    end if
    if (.not. ieee_is_finite(event%fraction) .or. event%fraction < 0.0_dp .or. event%fraction > 1.0_dp) then
      status = boundary_event_invalid_geometry
      return
    end if
    if (event%face_mask <= 0_i32 .or. iand(event%face_mask, boundary_face_all) /= event%face_mask) then
      status = boundary_event_invalid_geometry
      return
    end if
    if (any(event%face_bc /= boundary_face_conditions(cfg))) then
      status = boundary_event_invalid_geometry
      return
    end if
    do axis = 1_i32, 3_i32
      low_face = 2_i32*axis - 1_i32
      high_face = 2_i32*axis
      if (btest(event%face_mask, low_face - 1_i32) .and. btest(event%face_mask, high_face - 1_i32)) then
        status = boundary_event_invalid_geometry
        return
      end if
    end do

    x_work = x
    v_work = v
    alive_work = alive
    has_open = .false.
    do axis = 1_i32, 3_i32
      low_face = 2_i32*axis - 1_i32
      high_face = 2_i32*axis
      low_bc = event%face_bc(low_face)
      high_bc = event%face_bc(high_face)
      if (btest(event%face_mask, low_face - 1_i32)) has_open = has_open .or. low_bc == bc_open
      if (btest(event%face_mask, high_face - 1_i32)) has_open = has_open .or. high_bc == bc_open
    end do

    if (has_open) then
      alive_work = .false.
      escaped = .true.
    else
      do axis = 1_i32, 3_i32
        low_face = 2_i32*axis - 1_i32
        high_face = 2_i32*axis
        low_bc = event%face_bc(low_face)
        high_bc = event%face_bc(high_face)
        if (btest(event%face_mask, low_face - 1_i32)) then
          call apply_surviving_event_axis( &
            axis, low_bc, .false., cfg%box_min(axis), cfg%box_max(axis), x_work, v_work, status &
            )
        else if (btest(event%face_mask, high_face - 1_i32)) then
          call apply_surviving_event_axis( &
            axis, high_bc, .true., cfg%box_min(axis), cfg%box_max(axis), x_work, v_work, status &
            )
        end if
        if (status /= boundary_event_ok) return
      end do
    end if
    if (.not. all(ieee_is_finite(x_work)) .or. .not. all(ieee_is_finite(v_work))) then
      escaped = .false.
      status = boundary_event_invalid_geometry
      return
    end if

    x = x_work
    v = v_work
    alive = alive_work
  end subroutine apply_escape_reflect_periodic_event

  !> 反射または周期境界の単一軸作用をwork stateへ適用する。
  subroutine apply_surviving_event_axis(axis, bc, high_side, box_min, box_max, x, v, status)
    integer(i32), intent(in) :: axis, bc
    logical, intent(in) :: high_side
    real(dp), intent(in) :: box_min, box_max
    real(dp), intent(inout) :: x(3), v(3)
    integer(i32), intent(inout) :: status
    real(dp) :: inset

    inset = boundary_inset(box_min, box_max)

    select case (bc)
    case (bc_reflect)
      if (high_side) then
        x(axis) = box_max - inset
      else
        x(axis) = box_min + inset
      end if
      v(axis) = -v(axis)
    case (bc_periodic)
      if (high_side) then
        x(axis) = box_min + inset
      else
        x(axis) = box_max - inset
      end if
    case default
      status = boundary_event_invalid_geometry
    end select
  end subroutine apply_surviving_event_axis

  !> Event APIで扱える有限かつ正幅のbox設定かを返す。
  logical function valid_event_box(cfg) result(valid)
    type(sim_config), intent(in) :: cfg
    real(dp) :: inset
    integer(i32) :: axis

    valid = .false.
    if (.not. all(ieee_is_finite(cfg%box_min)) .or. .not. all(ieee_is_finite(cfg%box_max))) return
    if (any(cfg%box_max <= cfg%box_min)) return
    if (.not. all(valid_boundary_condition(cfg%bc_low)) .or. &
        .not. all(valid_boundary_condition(cfg%bc_high))) return
    do axis = 1_i32, 3_i32
      inset = boundary_inset(cfg%box_min(axis), cfg%box_max(axis))
      if (.not. ieee_is_finite(inset) .or. inset <= 0.0_dp) return
      if (cfg%box_min(axis) + inset >= cfg%box_max(axis) - inset) return
    end do
    valid = .true.
  end function valid_event_box

  !> 境界作用後の座標をnormal範囲かつevent toleranceより十分内側へ置く。
  pure real(dp) function boundary_inset(box_min, box_max) result(inset)
    real(dp), intent(in) :: box_min, box_max
    real(dp) :: span, scale

    span = box_max - box_min
    scale = max(abs(box_min), abs(box_max), span, tiny(1.0_dp))
    inset = max(64.0_dp*epsilon(1.0_dp)*scale, spacing(scale))
    inset = min(0.25_dp*span, inset)
  end function boundary_inset

  !> 境界条件enum配列の各値が既知かを返す。
  elemental logical function valid_boundary_condition(bc) result(valid)
    integer(i32), intent(in) :: bc

    valid = bc == bc_open .or. bc == bc_reflect .or. bc == bc_periodic
  end function valid_boundary_condition

  !> face bit順にbox設定の境界条件を並べる。
  pure function boundary_face_conditions(cfg) result(face_bc)
    type(sim_config), intent(in) :: cfg
    integer(i32) :: face_bc(6)
    integer(i32) :: axis

    do axis = 1_i32, 3_i32
      face_bc(2_i32*axis - 1_i32) = cfg%bc_low(axis)
      face_bc(2_i32*axis) = cfg%bc_high(axis)
    end do
  end function boundary_face_conditions

  !> 1ステップの更新候補位置にボックス境界条件を適用し、生存/流出状態と位置速度を更新する。
  !! @param[in] cfg 境界条件とボックス範囲を含むシミュレーション設定。
  !! @param[inout] x 境界適用対象の粒子位置 `(x,y,z)` [m]。
  !! @param[inout] v 境界適用対象の粒子速度 `(vx,vy,vz)` [m/s]。
  !! @param[inout] alive 粒子生存フラグ（流出時は `.false.` へ更新）。
  !! @param[out] escaped_boundary この呼び出しで境界流出が発生したか。
  subroutine apply_box_boundary( &
    cfg, x, v, alive, escaped_boundary, q_particle, m_particle, phi_boundary, potential_axis, potential_high_side &
    )
    type(sim_config), intent(in) :: cfg
    real(dp), intent(inout) :: x(3)
    real(dp), intent(inout) :: v(3)
    logical, intent(inout) :: alive
    logical, intent(out) :: escaped_boundary
    real(dp), intent(in), optional :: q_particle
    real(dp), intent(in), optional :: m_particle
    real(dp), intent(in), optional :: phi_boundary
    integer(i32), intent(in), optional :: potential_axis
    logical, intent(in), optional :: potential_high_side

    integer(i32) :: axis
    integer(i32) :: bc
    real(dp) :: span, eps

    escaped_boundary = .false.
    if (.not. cfg%use_box) return
    if (.not. alive) return

    eps = 1.0d-12
    if (present(potential_axis) .and. present(potential_high_side)) then
      axis = potential_axis
      if (axis >= 1_i32 .and. axis <= 3_i32) then
        span = cfg%box_max(axis) - cfg%box_min(axis)
        if (span <= 0.0_dp) then
          alive = .false.
          escaped_boundary = .true.
          return
        end if
        if ((.not. potential_high_side) .and. x(axis) < cfg%box_min(axis)) then
          call apply_one_side_boundary( &
            cfg, axis, cfg%bc_low(axis), cfg%box_min(axis), cfg%box_max(axis), span, eps, x, v, alive, &
            escaped_boundary, q_particle, m_particle, phi_boundary, potential_axis, potential_high_side &
            )
        else if (potential_high_side .and. x(axis) > cfg%box_max(axis)) then
          call apply_one_side_boundary( &
            cfg, axis, cfg%bc_high(axis), cfg%box_min(axis), cfg%box_max(axis), span, eps, x, v, alive, &
            escaped_boundary, q_particle, m_particle, phi_boundary, potential_axis, potential_high_side &
            )
        end if
        if (.not. alive) return
      end if
    end if

    do axis = 1, 3
      span = cfg%box_max(axis) - cfg%box_min(axis)
      if (span <= 0.0_dp) then
        alive = .false.
        escaped_boundary = .true.
        return
      end if

      do while (x(axis) < cfg%box_min(axis))
        bc = cfg%bc_low(axis)
        call apply_one_side_boundary( &
          cfg, axis, bc, cfg%box_min(axis), cfg%box_max(axis), span, eps, x, v, alive, escaped_boundary, &
          q_particle, m_particle, phi_boundary, potential_axis, potential_high_side &
          )
        if ((.not. alive) .or. (.not. escaped_boundary .and. x(axis) >= cfg%box_min(axis))) exit
      end do
      if (.not. alive) return

      do while (x(axis) > cfg%box_max(axis))
        bc = cfg%bc_high(axis)
        call apply_one_side_boundary( &
          cfg, axis, bc, cfg%box_min(axis), cfg%box_max(axis), span, eps, x, v, alive, escaped_boundary, &
          q_particle, m_particle, phi_boundary, potential_axis, potential_high_side &
          )
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
  subroutine apply_one_side_boundary( &
    cfg, axis, bc, box_min, box_max, span, eps, x, v, alive, escaped_boundary, q_particle, m_particle, phi_boundary, &
    potential_axis, potential_high_side &
    )
    type(sim_config), intent(in) :: cfg
    integer(i32), intent(in) :: axis, bc
    real(dp), intent(in) :: box_min, box_max, span, eps
    real(dp), intent(inout) :: x(3), v(3)
    logical, intent(inout) :: alive
    logical, intent(inout) :: escaped_boundary
    real(dp), intent(in), optional :: q_particle
    real(dp), intent(in), optional :: m_particle
    real(dp), intent(in), optional :: phi_boundary
    integer(i32), intent(in), optional :: potential_axis
    logical, intent(in), optional :: potential_high_side

    select case (bc)
    case (bc_open)
      if (open_boundary_should_reflect( &
          cfg, axis, x(axis) > box_max, v, q_particle, m_particle, phi_boundary, potential_axis, potential_high_side &
          )) then
        call reflect_one_side(axis, box_min, box_max, eps, x, v)
        escaped_boundary = .false.
      else
        alive = .false.
        escaped_boundary = .true.
      end if
    case (bc_reflect)
      call reflect_one_side(axis, box_min, box_max, eps, x, v)
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

  !> 開境界の電位障壁モデルで反射すべきかを判定する。
  logical function open_boundary_should_reflect( &
    cfg, axis, is_high_side, v, q_particle, m_particle, phi_boundary, potential_axis, potential_high_side &
    ) result(reflect)
    type(sim_config), intent(in) :: cfg
    integer(i32), intent(in) :: axis
    logical, intent(in) :: is_high_side
    real(dp), intent(in) :: v(3)
    real(dp), intent(in), optional :: q_particle
    real(dp), intent(in), optional :: m_particle
    real(dp), intent(in), optional :: phi_boundary
    integer(i32), intent(in), optional :: potential_axis
    logical, intent(in), optional :: potential_high_side

    real(dp) :: outward_v, kinetic_normal, potential_barrier

    reflect = .false.
    if (trim(cfg%open_boundary_model) /= 'potential_barrier') return
    if (.not. present(q_particle) .or. .not. present(m_particle) .or. .not. present(phi_boundary)) return
    if (present(potential_axis)) then
      if (potential_axis /= axis) return
    end if
    if (present(potential_high_side)) then
      if (potential_high_side .neqv. is_high_side) return
    end if
    if (m_particle <= 0.0d0) return

    if (is_high_side) then
      outward_v = v(axis)
    else
      outward_v = -v(axis)
    end if
    if (outward_v <= 0.0d0) return

    kinetic_normal = 0.5d0*m_particle*outward_v*outward_v
    potential_barrier = q_particle*(cfg%phi_infty - phi_boundary)
    reflect = potential_barrier > 0.0d0 .and. kinetic_normal < potential_barrier
  end function open_boundary_should_reflect

  !> 反射境界と同じ鏡像位置・法線速度反転を適用する。
  subroutine reflect_one_side(axis, box_min, box_max, eps, x, v)
    integer(i32), intent(in) :: axis
    real(dp), intent(in) :: box_min, box_max, eps
    real(dp), intent(inout) :: x(3), v(3)

    if (x(axis) < box_min) then
      x(axis) = box_min + (box_min - x(axis))
    else if (x(axis) > box_max) then
      x(axis) = box_max - (x(axis) - box_max)
    end if
    x(axis) = min(box_max - eps, max(box_min + eps, x(axis)))
    v(axis) = -v(axis)
  end subroutine reflect_one_side

end module bem_boundary
