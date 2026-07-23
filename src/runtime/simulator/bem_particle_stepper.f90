!> 同一時刻の粒子状態から、空間電場を中点評価した1ステップ候補を構築する。
module bem_particle_stepper
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, hit_info, bc_open, bc_reflect, bc_periodic
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_pusher, only: boris_push
  use bem_collision, only: collision_query_ok, find_first_hit
  use bem_boundary, only: boundary_event_type, boundary_event_ok, boundary_event_invalid_geometry, find_first_boundary_event, &
                          apply_escape_reflect_periodic_event
  use bem_interface_types, only: interface_crossing_type
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_open_escape, external_open_potential_barrier, &
    external_transport_none, external_boundary_owns_event
  implicit none
  private

  integer(i32), parameter, public :: particle_step_ok = collision_query_ok
  integer(i32), parameter, public :: particle_step_invalid_boundary = 1001_i32
  integer(i32), parameter, public :: particle_step_multiple_box_events = 1002_i32
  integer(i32), parameter, public :: particle_step_ambiguous_open_corner = 1003_i32
  integer(i32), parameter, public :: particle_step_unsupported_barrier_corner = particle_step_ambiguous_open_corner
  integer(i32), parameter, public :: particle_step_multiple_external_events = 1004_i32

  type, public :: particle_step_result
    real(dp) :: x(3) = 0.0_dp
    real(dp) :: v(3) = 0.0_dp
    logical :: absorbed = .false.
    logical :: escaped_boundary = .false.
    integer(i32) :: elem_idx = -1_i32
    integer(i32) :: status = particle_step_ok
    integer(i32) :: field_eval_count = 0_i32
    integer(i32) :: collision_query_count = 0_i32
    type(interface_crossing_type) :: interface_crossing
  end type particle_step_result

  public :: build_particle_step_candidate
  public :: advance_particle_step
  public :: resolve_particle_boundary_candidate

contains

  !> 予測中点の電場と一様磁場を使い、次時刻の位置・速度候補を返す。
  subroutine build_particle_step_candidate( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x1, v1 &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    real(dp), intent(out) :: x1(3), v1(3)
    real(dp) :: x_mid(3), e_mid(3)

    x_mid = x0 + 0.5d0*v0*dt
    call project_field_sample_to_box(sim, x_mid)
    call snapshot%eval_local_e(mesh, x_mid, e_mid)
    call boris_push(x0, v0, q, m, dt, e_mid, bfield, x1, v1)
  end subroutine build_particle_step_candidate

  !> 境界を越えるcandidateでも、場評価点はsolverのprimitive target box内に保つ。
  pure subroutine project_field_sample_to_box(sim, position)
    type(sim_config), intent(in) :: sim
    real(dp), intent(inout) :: position(3)
    integer(i32) :: axis
    real(dp) :: span

    if (.not. sim%use_box) return
    do axis = 1_i32, 3_i32
      if (.not. ieee_is_finite(sim%box_min(axis)) .or. .not. ieee_is_finite(sim%box_max(axis))) cycle
      span = sim%box_max(axis) - sim%box_min(axis)
      if (.not. ieee_is_finite(span) .or. span <= 0.0_dp) cycle
      if (sim%bc_low(axis) == bc_periodic .and. sim%bc_high(axis) == bc_periodic) then
        position(axis) = sim%box_min(axis) + modulo(position(axis) - sim%box_min(axis), span)
      else
        position(axis) = min(max(position(axis), sim%box_min(axis)), sim%box_max(axis))
      end if
    end do
  end subroutine project_field_sample_to_box

  !> 一つのouter stepについてmesh/boxの最早eventを順序付け、最大八度だけremainderを再積分する。
  subroutine advance_particle_step(mesh, sim, snapshot, bfield, x0, v0, q, m, dt, result, boundary_contract)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    type(particle_step_result), intent(out) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract

    type(hit_info) :: hit
    real(dp) :: x_candidate(3), v_candidate(3)
    integer(i32) :: query_status

    result = particle_step_result()
    result%x = x0
    result%v = v0
    if (.not. valid_particle_step_input(x0, v0, bfield, q, m, dt)) then
      result%status = particle_step_invalid_boundary
      return
    end if

    call build_particle_step_candidate(mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate)
    result%field_eval_count = 1_i32
    if (.not. all(ieee_is_finite(x_candidate)) .or. .not. all(ieee_is_finite(v_candidate))) then
      result%status = particle_step_invalid_boundary
      return
    end if

    call find_first_hit(mesh, x0, x_candidate, hit, sim=sim, status=query_status)
    result%collision_query_count = 1_i32
    if (query_status /= collision_query_ok) then
      if (sim%use_box .and. .not. point_strictly_inside_box(sim, x_candidate)) then
        call advance_particle_boundary_crossing( &
          mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, result=result, &
          boundary_contract=boundary_contract &
          )
        return
      end if
      result%status = query_status
      return
    end if
    if (.not. sim%use_box .or. point_strictly_inside_box(sim, x_candidate)) then
      if (hit%has_hit) then
        call accept_particle_hit(v0, x_candidate, v_candidate, hit, result)
        return
      end if
      result%x = x_candidate
      result%v = v_candidate
      return
    end if

    call advance_particle_boundary_crossing( &
      mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, hit, result, &
      boundary_contract=boundary_contract &
      )
  end subroutine advance_particle_step

  !> 構築済みcandidateがbox crossing候補のとき、fieldを再評価せずeventを解決する。
  subroutine resolve_particle_boundary_candidate( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, hit, result, boundary_contract &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt, x_candidate(3), v_candidate(3)
    type(hit_info), intent(in), optional :: hit
    type(particle_step_result), intent(out) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract

    result = particle_step_result()
    result%x = x0
    result%v = v0
    result%field_eval_count = 1_i32
    result%collision_query_count = merge(1_i32, 0_i32, present(hit))
    if (.not. valid_particle_step_input(x0, v0, bfield, q, m, dt) .or. &
        .not. all(ieee_is_finite(x_candidate)) .or. .not. all(ieee_is_finite(v_candidate))) then
      result%status = particle_step_invalid_boundary
      return
    end if
    call advance_particle_boundary_crossing( &
      mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, hit, result, boundary_contract &
      )
  end subroutine resolve_particle_boundary_candidate

  !> box crossing時だけevent用stateを確保し、最大八つのeventとremainderを処理する。
  subroutine advance_particle_boundary_crossing( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, hit, result, boundary_contract &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt, x_candidate(3), v_candidate(3)
    type(hit_info), intent(in), optional :: hit
    type(particle_step_result), intent(inout) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract

    integer(i32), parameter :: max_boundary_events = 8_i32
    type(boundary_event_type) :: event
    type(hit_info) :: remainder_hit
    real(dp) :: x_start(3), v_start(3), x_trial(3), v_trial(3), x_event(3), v_event(3)
    real(dp) :: dt_segment, dt_remaining, elapsed_dt, global_fraction, refined_fraction, event_order_tolerance
    integer(i32) :: query_status, boundary_status, event_count
    logical :: alive, escaped, external_event_owned, z_high_chord_candidate
    logical :: first_segment
    type(external_boundary_contract_type) :: active_boundary_contract

    active_boundary_contract = external_boundary_contract_type()
    select case (trim(sim%open_boundary_model))
    case ('escape')
      active_boundary_contract%ordinary_open_model = external_open_escape
    case ('potential_barrier')
      active_boundary_contract%ordinary_open_model = external_open_potential_barrier
    case default
      result%status = particle_step_invalid_boundary
      return
    end select
    if (present(boundary_contract)) active_boundary_contract = boundary_contract
    x_start = x0
    v_start = v0
    x_trial = x_candidate
    v_trial = v_candidate
    dt_segment = dt
    elapsed_dt = 0.0_dp
    event_count = 0_i32
    first_segment = .true.

    do
      if (point_strictly_inside_box(sim, x_trial)) then
        if (first_segment .and. present(hit)) then
          if (hit%has_hit) then
            call accept_particle_hit(v_start, x_trial, v_trial, hit, result)
            return
          end if
        else
          call query_particle_chord(mesh, sim, x_start, v_start, x_trial, v_trial, result, remainder_hit, query_status)
          if (query_status /= collision_query_ok .or. result%absorbed) return
        end if
        result%x = x_trial
        result%v = v_trial
        return
      end if

      call find_first_boundary_event(sim, x_start, x_trial, event, boundary_status)
      if (boundary_status /= boundary_event_ok) then
        result%status = particle_step_invalid_boundary
        return
      end if
      if (.not. event%has_event) then
        if (first_segment .and. present(hit)) then
          if (hit%has_hit) then
            call accept_particle_hit(v_start, x_trial, v_trial, hit, result)
            return
          end if
        else
          call query_particle_chord(mesh, sim, x_start, v_start, x_trial, v_trial, result, remainder_hit, query_status)
          if (query_status /= collision_query_ok .or. result%absorbed) return
        end if
        result%x = x_trial
        result%v = v_trial
        return
      end if

      external_event_owned = external_boundary_owns_event( &
                             active_boundary_contract, event%face_mask, event%face_bc &
                             )
      z_high_chord_candidate = external_event_owned .or. &
                               (active_boundary_contract%interface_transport /= external_transport_none .and. &
                                event%face_bc(6) == bc_open .and. x_start(3) < sim%box_max(3) .and. &
                                x_trial(3) >= sim%box_max(3))
      if (z_high_chord_candidate) then
        call refine_z_high_crossing_fraction( &
          sim, x_start, v_start, x_trial, v_trial, dt_segment, refined_fraction, boundary_status &
          )
        if (external_event_owned) then
          if (boundary_status /= boundary_event_ok) then
            result%status = particle_step_invalid_boundary
            return
          end if
          if (.not. refined_z_high_preserves_lateral_event_order( &
              sim, event, x_start, x_trial, refined_fraction &
              )) then
            result%status = particle_step_invalid_boundary
            return
          end if
          event%fraction = refined_fraction
        else if (boundary_status == boundary_event_ok) then
          event_order_tolerance = 64.0_dp*epsilon(1.0_dp)* &
                                  max(1.0_dp, abs(event%fraction), abs(refined_fraction))
          if (refined_fraction <= event%fraction + event_order_tolerance) then
            result%status = particle_step_invalid_boundary
            return
          end if
        end if
      end if

      if (first_segment .and. present(hit)) then
        if (hit%has_hit) then
          if (hit%t <= event%fraction + 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(event%fraction))) then
            call accept_particle_hit(v_start, x_trial, v_trial, hit, result)
            return
          end if
        end if
      end if

      call interpolate_boundary_state(sim, event, x_start, v_start, x_trial, v_trial, x_event, v_event)
      if (.not. (first_segment .and. present(hit))) then
        call query_particle_chord(mesh, sim, x_start, v_start, x_event, v_event, result, remainder_hit, query_status)
        if (query_status /= collision_query_ok .or. result%absorbed) return
      end if

      if (event_count >= max_boundary_events) then
        result%status = particle_step_multiple_box_events
        return
      end if
      event_count = event_count + 1_i32
      dt_remaining = (1.0_dp - event%fraction)*dt_segment
      global_fraction = (elapsed_dt + event%fraction*dt_segment)/dt

      if (external_event_owned) then
        call apply_external_companion_faces( &
          sim, event, x_event, v_event, boundary_status &
          )
        if (boundary_status /= boundary_event_ok) then
          result%status = boundary_status
          return
        end if
        result%x = x_event
        result%v = v_event
        result%interface_crossing%has_crossing = .true.
        result%interface_crossing%face_index = 6_i32
        result%interface_crossing%fraction = global_fraction
        result%interface_crossing%position = x_event
        result%interface_crossing%velocity = v_event
        result%interface_crossing%dt_remaining = dt_remaining
        return
      end if

      alive = .true.
      escaped = .false.
      if (active_boundary_contract%ordinary_open_model == external_open_potential_barrier) then
        call apply_potential_barrier_event( &
          mesh, sim, snapshot, event, q, m, x_event, v_event, alive, escaped, boundary_status &
          )
      else
        call apply_escape_reflect_periodic_event(sim, event, x_event, v_event, alive, escaped, boundary_status)
      end if
      if (boundary_status /= boundary_event_ok) then
        if (boundary_status == boundary_event_invalid_geometry) then
          result%status = particle_step_invalid_boundary
        else
          result%status = boundary_status
        end if
        return
      end if
      if (.not. alive) then
        result%x = x_event
        result%v = v_event
        result%escaped_boundary = escaped
        return
      end if
      if (dt_remaining <= 0.0_dp) then
        result%x = x_event
        result%v = v_event
        return
      end if

      elapsed_dt = elapsed_dt + event%fraction*dt_segment
      x_start = x_event
      v_start = v_event
      dt_segment = dt_remaining
      call build_particle_step_candidate( &
        mesh, sim, snapshot, bfield, x_start, v_start, q, m, dt_segment, x_trial, v_trial &
        )
      result%field_eval_count = result%field_eval_count + 1_i32
      if (.not. all(ieee_is_finite(x_trial)) .or. .not. all(ieee_is_finite(v_trial))) then
        result%status = particle_step_invalid_boundary
        return
      end if
      first_segment = .false.
    end do
  end subroutine advance_particle_boundary_crossing

  !> Boris端点と整合する二次軌道から、外向きz-high交差の時刻率を求める。
  subroutine refine_z_high_crossing_fraction(sim, x0, v0, x1, v1, dt, fraction, status)
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: x0(3), v0(3), x1(3), v1(3), dt
    real(dp), intent(out) :: fraction
    integer(i32), intent(out) :: status
    real(dp) :: a, b, c, discriminant, q_root, roots(2), velocity, coefficient_tolerance, root_tolerance
    integer(i32) :: root

    fraction = 0.0_dp
    status = boundary_event_invalid_geometry
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. x0(3) >= sim%box_max(3) .or. &
        x1(3) < sim%box_max(3)) return

    a = 0.5_dp*dt*(v1(3) - v0(3))
    b = dt*v0(3)
    c = x0(3) - sim%box_max(3)
    coefficient_tolerance = 256.0_dp*epsilon(1.0_dp)*max(tiny(1.0_dp), abs(a), abs(b), abs(c))
    root_tolerance = 256.0_dp*epsilon(1.0_dp)
    roots = huge(1.0_dp)
    if (abs(a) <= coefficient_tolerance) then
      if (b <= 0.0_dp) return
      roots(1) = -c/b
    else
      discriminant = b*b - 4.0_dp*a*c
      if (.not. ieee_is_finite(discriminant) .or. discriminant < 0.0_dp) return
      q_root = -0.5_dp*(b + sign(sqrt(discriminant), b))
      if (q_root == 0.0_dp) return
      roots = [q_root/a, c/q_root]
    end if
    do root = 1_i32, 2_i32
      velocity = v0(3) + roots(root)*(v1(3) - v0(3))
      if (roots(root) >= -root_tolerance .and. roots(root) <= 1.0_dp + root_tolerance .and. velocity > 0.0_dp) then
        fraction = min(max(roots(root), 0.0_dp), 1.0_dp)
        status = boundary_event_ok
        return
      end if
    end do
  end subroutine refine_z_high_crossing_fraction

  !> z-high時刻補正後も、既存chordが示すx/y面の先後・同時関係が保たれるか検査する。
  pure logical function refined_z_high_preserves_lateral_event_order( &
    sim, event, x0, x1, refined_fraction &
    ) result(preserves_order)
    type(sim_config), intent(in) :: sim
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(in) :: x0(3), x1(3), refined_fraction

    real(dp) :: delta, candidate_fraction, tie_tolerance
    integer(i32) :: axis, candidate_face, axis_mask
    logical :: has_original_companion

    preserves_order = .false.
    if (.not. ieee_is_finite(refined_fraction) .or. refined_fraction < 0.0_dp .or. &
        refined_fraction > 1.0_dp) return
    do axis = 1_i32, 2_i32
      axis_mask = ior(shiftl(1_i32, 2_i32*axis - 2_i32), shiftl(1_i32, 2_i32*axis - 1_i32))
      has_original_companion = iand(event%face_mask, axis_mask) /= 0_i32
      delta = x1(axis) - x0(axis)
      candidate_face = 0_i32
      candidate_fraction = huge(1.0_dp)
      if (delta > 0.0_dp .and. x1(axis) >= sim%box_max(axis)) then
        candidate_face = 2_i32*axis
        candidate_fraction = (sim%box_max(axis) - x0(axis))/delta
      else if (delta < 0.0_dp .and. x1(axis) <= sim%box_min(axis)) then
        candidate_face = 2_i32*axis - 1_i32
        candidate_fraction = (sim%box_min(axis) - x0(axis))/delta
      end if
      if (candidate_face == 0_i32) then
        if (has_original_companion) return
        cycle
      end if
      if (.not. ieee_is_finite(candidate_fraction)) return
      tie_tolerance = 64.0_dp*epsilon(1.0_dp)* &
                      max(1.0_dp, abs(candidate_fraction), abs(refined_fraction))
      if (has_original_companion) then
        if (.not. btest(event%face_mask, candidate_face - 1_i32) .or. &
            abs(candidate_fraction - refined_fraction) > tie_tolerance) return
      else if (candidate_fraction <= refined_fraction + tie_tolerance) then
        return
      end if
    end do
    preserves_order = .true.
  end function refined_z_high_preserves_lateral_event_order

  !> 既に選択済みのmesh hitをresultへ反映する。
  subroutine accept_particle_hit(va, xb, vb, hit, result)
    real(dp), intent(in) :: va(3), xb(3), vb(3)
    type(hit_info), intent(in) :: hit
    type(particle_step_result), intent(inout) :: result

    result%x = hit%pos
    result%v = va + hit%t*(vb - va)
    result%absorbed = .true.
    result%elem_idx = hit%elem_idx
  end subroutine accept_particle_hit

  !> 線分衝突を一度照会し、hit/statusを共通resultへ反映する。
  subroutine query_particle_chord(mesh, sim, xa, va, xb, vb, result, hit, query_status)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: xa(3), va(3), xb(3), vb(3)
    type(particle_step_result), intent(inout) :: result
    type(hit_info), intent(out) :: hit
    integer(i32), intent(out) :: query_status

    call find_first_hit(mesh, xa, xb, hit, sim=sim, status=query_status)
    result%collision_query_count = result%collision_query_count + 1_i32
    if (query_status /= collision_query_ok) then
      result%status = query_status
      return
    end if
    if (.not. hit%has_hit) return

    result%x = hit%pos
    result%v = va + hit%t*(vb - va)
    result%absorbed = .true.
    result%elem_idx = hit%elem_idx
  end subroutine query_particle_chord

  !> Event fractionのdense stateを作り、setされたface座標をbox値へ正確に揃える。
  subroutine interpolate_boundary_state(sim, event, x0, v0, x1, v1, x_event, v_event)
    type(sim_config), intent(in) :: sim
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(in) :: x0(3), v0(3), x1(3), v1(3)
    real(dp), intent(out) :: x_event(3), v_event(3)
    integer(i32) :: axis

    x_event = x0 + event%fraction*(x1 - x0)
    v_event = v0 + event%fraction*(v1 - v0)
    do axis = 1_i32, 3_i32
      if (btest(event%face_mask, 2_i32*axis - 2_i32)) x_event(axis) = sim%box_min(axis)
      if (btest(event%face_mask, 2_i32*axis - 1_i32)) x_event(axis) = sim%box_max(axis)
    end do
  end subroutine interpolate_boundary_state

  !> interface面と同時に交差した周期・反射面を先に合成し、別open面との曖昧なcornerは拒否する。
  subroutine apply_external_companion_faces(sim, event, x, v, status)
    type(sim_config), intent(in) :: sim
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(inout) :: x(3), v(3)
    integer(i32), intent(out) :: status

    type(sim_config) :: action_sim
    type(boundary_event_type) :: companion_event
    integer(i32) :: face, action_status
    logical :: alive, escaped

    status = boundary_event_ok
    companion_event = event
    companion_event%face_mask = ibclr(companion_event%face_mask, 5_i32)
    if (companion_event%face_mask == 0_i32) return
    do face = 1_i32, 6_i32
      if (btest(companion_event%face_mask, face - 1_i32) .and. companion_event%face_bc(face) == bc_open) then
        status = particle_step_ambiguous_open_corner
        return
      end if
    end do

    action_sim = sim
    action_sim%open_boundary_model = 'escape'
    alive = .true.
    escaped = .false.
    call apply_escape_reflect_periodic_event( &
      action_sim, companion_event, x, v, alive, escaped, action_status &
      )
    if (action_status /= boundary_event_ok .or. .not. alive .or. escaped) then
      status = particle_step_invalid_boundary
    end if
  end subroutine apply_external_companion_faces

  !> 単一open面のpotential-barrier式をevent位置と補間速度で評価する。
  subroutine apply_potential_barrier_event( &
    mesh, sim, snapshot, event, q, m, x, v, alive, escaped, status &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(in) :: q, m
    real(dp), intent(inout) :: x(3), v(3)
    logical, intent(inout) :: alive
    logical, intent(out) :: escaped
    integer(i32), intent(out) :: status

    type(sim_config) :: action_sim
    type(boundary_event_type) :: action_event
    real(dp) :: phi_boundary, outward_v, kinetic_normal, potential_barrier
    integer(i32) :: axis, face_index, open_count
    logical :: high_side, reflect_open

    status = boundary_event_ok
    escaped = .false.
    open_count = 0_i32
    face_index = 0_i32
    high_side = .false.
    do axis = 1_i32, 3_i32
      if (btest(event%face_mask, 2_i32*axis - 2_i32) .and. event%face_bc(2_i32*axis - 1_i32) == bc_open) then
        open_count = open_count + 1_i32
        face_index = 2_i32*axis - 1_i32
        high_side = .false.
      end if
      if (btest(event%face_mask, 2_i32*axis - 1_i32) .and. event%face_bc(2_i32*axis) == bc_open) then
        open_count = open_count + 1_i32
        face_index = 2_i32*axis
        high_side = .true.
      end if
    end do
    if (open_count > 1_i32) then
      status = particle_step_ambiguous_open_corner
      return
    end if

    action_sim = sim
    action_sim%open_boundary_model = 'escape'
    action_event = event
    if (open_count == 0_i32) then
      call apply_escape_reflect_periodic_event(action_sim, action_event, x, v, alive, escaped, status)
      return
    end if
    if (.not. ieee_is_finite(q) .or. .not. ieee_is_finite(m) .or. m <= 0.0_dp) then
      status = particle_step_invalid_boundary
      return
    end if

    axis = (face_index + 1_i32)/2_i32
    call snapshot%eval_local_phi(mesh, sim, x, phi_boundary)
    if (.not. ieee_is_finite(phi_boundary)) then
      status = particle_step_invalid_boundary
      return
    end if
    if (high_side) then
      outward_v = v(axis)
    else
      outward_v = -v(axis)
    end if
    kinetic_normal = 0.5_dp*m*outward_v*outward_v
    potential_barrier = q*(sim%phi_infty - phi_boundary)
    if (.not. ieee_is_finite(kinetic_normal) .or. .not. ieee_is_finite(potential_barrier)) then
      status = particle_step_invalid_boundary
      return
    end if
    reflect_open = outward_v > 0.0_dp .and. potential_barrier > 0.0_dp .and. kinetic_normal < potential_barrier
    if (.not. reflect_open) then
      alive = .false.
      escaped = .true.
      return
    end if

    action_event%face_bc(face_index) = bc_reflect
    if (high_side) then
      action_sim%bc_high(axis) = bc_reflect
    else
      action_sim%bc_low(axis) = bc_reflect
    end if
    call apply_escape_reflect_periodic_event(action_sim, action_event, x, v, alive, escaped, status)
  end subroutine apply_potential_barrier_event

  !> Candidate endpointがboxの全faceからstrictly interiorかを返すfast-path判定。
  pure logical function point_strictly_inside_box(sim, x) result(inside)
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: x(3)

    inside = all(x > sim%box_min) .and. all(x < sim%box_max)
  end function point_strictly_inside_box

  !> outer particle stepの有限値・質量・時間刻み入力を検証する。
  pure logical function valid_particle_step_input(x, v, b, q, m, dt) result(valid)
    real(dp), intent(in) :: x(3), v(3), b(3), q, m, dt

    valid = all(ieee_is_finite(x)) .and. all(ieee_is_finite(v)) .and. all(ieee_is_finite(b)) .and. &
            ieee_is_finite(q) .and. ieee_is_finite(m) .and. m > 0.0_dp .and. &
            ieee_is_finite(dt) .and. dt >= 0.0_dp
  end function valid_particle_step_input

end module bem_particle_stepper
