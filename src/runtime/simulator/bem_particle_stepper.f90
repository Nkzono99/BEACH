!> 同一時刻の粒子状態から、空間電場を中点評価した1ステップ候補を構築する。
module bem_particle_stepper
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: mesh_type, sim_config, hit_info, bc_open, bc_reflect, bc_periodic, bc_redistributed_reflect
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_pusher, only: boris_push
  use bem_collision, only: collision_query_ok, find_first_hit
  use bem_boundary, only: boundary_event_type, boundary_event_ok, boundary_event_invalid_geometry, find_first_boundary_event, &
                          apply_escape_reflect_periodic_event
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_open_escape, external_open_potential_barrier
  implicit none
  private

  integer(i32), parameter, public :: particle_step_ok = collision_query_ok
  integer(i32), parameter, public :: particle_step_invalid_boundary = 1001_i32
  integer(i32), parameter, public :: particle_step_multiple_box_events = 1002_i32
  integer(i32), parameter, public :: particle_step_ambiguous_open_corner = 1003_i32
  integer(i32), parameter, public :: particle_step_unsupported_barrier_corner = particle_step_ambiguous_open_corner

  type, public :: particle_step_result
    real(dp) :: x(3) = 0.0_dp
    real(dp) :: v(3) = 0.0_dp
    logical :: absorbed = .false.
    logical :: escaped_boundary = .false.
    integer(i32) :: elem_idx = -1_i32
    integer(i32) :: status = particle_step_ok
    integer(i32) :: field_eval_count = 0_i32
    integer(i32) :: collision_query_count = 0_i32
  end type particle_step_result

  public :: build_particle_step_candidate
  public :: advance_particle_step
  public :: resolve_particle_boundary_candidate

contains

  !> 予測中点の電場と一様磁場を使い、次時刻の位置・速度候補を返す。
  subroutine build_particle_step_candidate( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x1, v1, sampled_electric_field &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    real(dp), intent(out) :: x1(3), v1(3)
    real(dp), intent(out), optional :: sampled_electric_field(3)
    real(dp) :: x_mid(3), e_mid(3)

    x_mid = x0 + 0.5d0*v0*dt
    call project_field_sample_to_box(sim, x_mid)
    call snapshot%eval_local_e(mesh, x_mid, e_mid)
    call boris_push(x0, v0, q, m, dt, e_mid, bfield, x1, v1)
    if (present(sampled_electric_field)) sampled_electric_field = e_mid
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

  !> 一つの粒子stepについてmesh/boxの最早eventを順序付ける。
  subroutine advance_particle_step( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, result, boundary_contract, boundary_rng_counter &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    type(particle_step_result), intent(out) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in), optional :: boundary_rng_counter(4)

    call advance_particle_step_impl( &
      mesh, sim, snapshot, bfield, x0, v0, q, m, dt, result, boundary_contract, boundary_rng_counter, 0_i32 &
      )
  end subroutine advance_particle_step

  !> periodic event の適応分割深さを内部だけで引き回す1 step実装。
  recursive subroutine advance_particle_step_impl( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, result, boundary_contract, boundary_rng_counter, adaptive_depth &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    type(particle_step_result), intent(out) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in), optional :: boundary_rng_counter(4)
    integer(i32), intent(in) :: adaptive_depth

    type(hit_info) :: hit
    real(dp) :: x_candidate(3), v_candidate(3), candidate_electric_field(3)
    integer(i32) :: query_status
    integer(i64) :: active_boundary_rng_counter(4)
    logical :: has_boundary_rng_counter

    result = particle_step_result()
    result%x = x0
    result%v = v0
    active_boundary_rng_counter = 0_i64
    has_boundary_rng_counter = present(boundary_rng_counter)
    if (has_boundary_rng_counter) active_boundary_rng_counter = boundary_rng_counter
    if (.not. valid_particle_step_input(x0, v0, bfield, q, m, dt)) then
      result%status = particle_step_invalid_boundary
      return
    end if

    call build_particle_step_candidate( &
      mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, candidate_electric_field &
      )
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
          mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, candidate_electric_field, result=result, &
          boundary_contract=boundary_contract, boundary_rng_counter=active_boundary_rng_counter, &
          has_boundary_rng_counter=has_boundary_rng_counter, adaptive_depth=adaptive_depth &
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
      mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, candidate_electric_field, hit, result, &
      boundary_contract=boundary_contract, boundary_rng_counter=active_boundary_rng_counter, &
      has_boundary_rng_counter=has_boundary_rng_counter, adaptive_depth=adaptive_depth &
      )
  end subroutine advance_particle_step_impl

  !> 構築済みcandidateとその場評価値からbox eventを解決する。
  subroutine resolve_particle_boundary_candidate( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, hit, result, boundary_contract, &
    boundary_rng_counter, sampled_electric_field &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt, x_candidate(3), v_candidate(3)
    type(hit_info), intent(in), optional :: hit
    type(particle_step_result), intent(out) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in), optional :: boundary_rng_counter(4)
    real(dp), intent(in), optional :: sampled_electric_field(3)
    integer(i64) :: active_boundary_rng_counter(4)
    real(dp) :: active_electric_field(3), x_mid(3)
    logical :: has_boundary_rng_counter

    result = particle_step_result()
    result%x = x0
    result%v = v0
    active_boundary_rng_counter = 0_i64
    has_boundary_rng_counter = present(boundary_rng_counter)
    if (has_boundary_rng_counter) active_boundary_rng_counter = boundary_rng_counter
    result%field_eval_count = 1_i32
    result%collision_query_count = merge(1_i32, 0_i32, present(hit))
    if (.not. valid_particle_step_input(x0, v0, bfield, q, m, dt) .or. &
        .not. all(ieee_is_finite(x_candidate)) .or. .not. all(ieee_is_finite(v_candidate))) then
      result%status = particle_step_invalid_boundary
      return
    end if
    if (present(sampled_electric_field)) then
      active_electric_field = sampled_electric_field
    else
      x_mid = x0 + 0.5_dp*v0*dt
      call project_field_sample_to_box(sim, x_mid)
      call snapshot%eval_local_e(mesh, x_mid, active_electric_field)
      result%field_eval_count = result%field_eval_count + 1_i32
    end if
    if (.not. all(ieee_is_finite(active_electric_field))) then
      result%status = particle_step_invalid_boundary
      return
    end if
    call advance_particle_boundary_crossing( &
      mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, active_electric_field, hit, result, &
      boundary_contract, active_boundary_rng_counter, has_boundary_rng_counter, 0_i32 &
      )
  end subroutine resolve_particle_boundary_candidate

  !> box crossing時だけevent用stateを確保し、periodic過多時はremainderを適応分割する。
  subroutine advance_particle_boundary_crossing( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, candidate_electric_field, hit, result, &
    boundary_contract, boundary_rng_counter, has_boundary_rng_counter, adaptive_depth &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt, x_candidate(3), v_candidate(3)
    real(dp), intent(in) :: candidate_electric_field(3)
    type(hit_info), intent(in), optional :: hit
    type(particle_step_result), intent(inout) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in) :: boundary_rng_counter(4)
    logical, intent(in) :: has_boundary_rng_counter
    integer(i32), intent(in) :: adaptive_depth

    integer(i32), parameter :: max_boundary_events = 8_i32
    type(boundary_event_type) :: event
    type(hit_info) :: remainder_hit
    real(dp) :: x_start(3), v_start(3), x_trial(3), v_trial(3), x_event(3), v_event(3), segment_electric_field(3)
    real(dp) :: dt_segment, dt_remaining
    real(dp) :: redistribution_uniform(3)
    integer(i32) :: query_status, boundary_status, event_count
    logical :: alive, escaped
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
    segment_electric_field = candidate_electric_field
    dt_segment = dt
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

      if (first_segment .and. present(hit)) then
        if (hit%has_hit) then
          if (hit%t <= event%fraction + 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(event%fraction))) then
            call accept_particle_hit(v_start, x_trial, v_trial, hit, result)
            return
          end if
        end if
      end if

      call evaluate_boundary_state( &
        sim, event, x_start, v_start, x_trial, v_trial, bfield, q, m, dt_segment, segment_electric_field, x_event, v_event, &
        boundary_status &
        )
      if (boundary_status /= boundary_event_ok) then
        result%status = particle_step_invalid_boundary
        return
      end if
      if (.not. (first_segment .and. present(hit))) then
        call query_particle_chord(mesh, sim, x_start, v_start, x_event, v_event, result, remainder_hit, query_status)
        if (query_status /= collision_query_ok .or. result%absorbed) return
      end if

      if (event_count >= max_boundary_events) then
        if (periodic_event_can_subdivide(sim, event) .and. adaptive_depth < 12_i32) then
          call advance_periodic_substeps( &
            mesh, sim, snapshot, bfield, x_start, v_start, q, m, dt_segment, result, active_boundary_contract, &
            boundary_rng_counter, has_boundary_rng_counter, adaptive_depth &
            )
          return
        end if
        result%status = particle_step_multiple_box_events
        return
      end if
      event_count = event_count + 1_i32
      if (event_requires_redistribution_counter(event, active_boundary_contract) .and. &
          .not. has_boundary_rng_counter) then
        result%status = particle_step_invalid_boundary
        return
      end if
      call generate_redistribution_uniform( &
        sim%rng_seed, boundary_rng_counter, event_count, redistribution_uniform &
        )
      dt_remaining = (1.0_dp - event%fraction)*dt_segment
      alive = .true.
      escaped = .false.
      if (event_uses_potential_barrier(event, active_boundary_contract)) then
        call apply_potential_barrier_event( &
          mesh, sim, snapshot, event, active_boundary_contract, q, m, x_event, v_event, alive, escaped, boundary_status, &
          redistribution_uniform &
          )
      else
        call apply_escape_reflect_periodic_event( &
          sim, event, x_event, v_event, alive, escaped, boundary_status, redistribution_uniform &
          )
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

      x_start = x_event
      v_start = v_event
      dt_segment = dt_remaining
      call build_particle_step_candidate( &
        mesh, sim, snapshot, bfield, x_start, v_start, q, m, dt_segment, x_trial, v_trial, segment_electric_field &
        )
      result%field_eval_count = result%field_eval_count + 1_i32
      if (.not. all(ieee_is_finite(x_trial)) .or. .not. all(ieee_is_finite(v_trial))) then
        result%status = particle_step_invalid_boundary
        return
      end if
      first_segment = .false.
    end do
  end subroutine advance_particle_boundary_crossing

  !> 一つのperiodic remainderを二分し、各半stepのbox/collision eventを順番に解く。
  subroutine advance_periodic_substeps( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, result, boundary_contract, boundary_rng_counter, &
    has_boundary_rng_counter, adaptive_depth &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    type(particle_step_result), intent(inout) :: result
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i64), intent(in) :: boundary_rng_counter(4)
    logical, intent(in) :: has_boundary_rng_counter
    integer(i32), intent(in) :: adaptive_depth

    type(particle_step_result) :: first_half, second_half
    integer(i32) :: prior_field_evals, prior_collision_queries

    prior_field_evals = result%field_eval_count
    prior_collision_queries = result%collision_query_count
    if (has_boundary_rng_counter) then
      call advance_particle_step_impl( &
        mesh, sim, snapshot, bfield, x0, v0, q, m, 0.5_dp*dt, first_half, boundary_contract, &
        boundary_rng_counter, adaptive_depth + 1_i32 &
        )
    else
      call advance_particle_step_impl( &
        mesh, sim, snapshot, bfield, x0, v0, q, m, 0.5_dp*dt, first_half, boundary_contract=boundary_contract, &
        adaptive_depth=adaptive_depth + 1_i32 &
        )
    end if
    if (first_half%status /= particle_step_ok .or. first_half%absorbed .or. first_half%escaped_boundary) then
      result = first_half
      result%field_eval_count = prior_field_evals + first_half%field_eval_count
      result%collision_query_count = prior_collision_queries + first_half%collision_query_count
      return
    end if

    if (has_boundary_rng_counter) then
      call advance_particle_step_impl( &
        mesh, sim, snapshot, bfield, first_half%x, first_half%v, q, m, 0.5_dp*dt, second_half, boundary_contract, &
        boundary_rng_counter, adaptive_depth + 1_i32 &
        )
    else
      call advance_particle_step_impl( &
        mesh, sim, snapshot, bfield, first_half%x, first_half%v, q, m, 0.5_dp*dt, second_half, &
        boundary_contract=boundary_contract, adaptive_depth=adaptive_depth + 1_i32 &
        )
    end if
    result = second_half
    result%field_eval_count = prior_field_evals + first_half%field_eval_count + second_half%field_eval_count
    result%collision_query_count = prior_collision_queries + first_half%collision_query_count + &
                                   second_half%collision_query_count
  end subroutine advance_periodic_substeps

  !> eventの全faceがperiodicで、適応分割中に乱数境界へ入らない場合だけ再試行する。
  pure logical function periodic_event_can_subdivide(sim, event) result(can_subdivide)
    type(sim_config), intent(in) :: sim
    type(boundary_event_type), intent(in) :: event
    integer(i32) :: axis

    can_subdivide = event%has_event
    if (any(sim%bc_low == bc_redistributed_reflect) .or. any(sim%bc_high == bc_redistributed_reflect)) then
      can_subdivide = .false.
      return
    end if
    do axis = 1_i32, 3_i32
      if (btest(event%face_mask, 2_i32*axis - 2_i32)) then
        if (sim%bc_low(axis) /= bc_periodic) can_subdivide = .false.
      end if
      if (btest(event%face_mask, 2_i32*axis - 1_i32)) then
        if (sim%bc_high(axis) /= bc_periodic) can_subdivide = .false.
      end if
    end do
  end function periodic_event_can_subdivide

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

  !> Chord上のevent位置と、その接線方向・離散workに整合する速度を返す。
  subroutine evaluate_boundary_state( &
    sim, event, x0, v0, x1, v1, bfield, q, m, dt, electric_field, x_event, v_event, status &
    )
    type(sim_config), intent(in) :: sim
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(in) :: x0(3), v0(3), x1(3), v1(3), bfield(3), q, m, dt, electric_field(3)
    real(dp), intent(out) :: x_event(3), v_event(3)
    integer(i32), intent(out) :: status
    integer(i32) :: axis
    real(dp) :: chord_velocity(3), chord_speed2, event_speed2, speed_scale, speed_tolerance
    logical :: force_free

    status = boundary_event_ok
    x_event = x0 + event%fraction*(x1 - x0)
    do axis = 1_i32, 3_i32
      if (btest(event%face_mask, 2_i32*axis - 2_i32)) x_event(axis) = sim%box_min(axis)
      if (btest(event%face_mask, 2_i32*axis - 1_i32)) x_event(axis) = sim%box_max(axis)
    end do
    chord_velocity = (x1 - x0)/dt
    chord_speed2 = sum(chord_velocity*chord_velocity)
    event_speed2 = sum(v0*v0) + 2.0_dp*(q/m)*dot_product(electric_field, x_event - x0)
    speed_tolerance = 256.0_dp*epsilon(1.0_dp)*max(tiny(1.0_dp), sum(v0*v0), chord_speed2, abs(event_speed2))
    if (.not. ieee_is_finite(chord_speed2) .or. .not. ieee_is_finite(event_speed2) .or. &
        chord_speed2 <= tiny(1.0_dp) .or. event_speed2 <= speed_tolerance) then
      v_event = 0.0_dp
      status = boundary_event_invalid_geometry
      return
    end if
    force_free = q == 0.0_dp .or. (all(electric_field == 0.0_dp) .and. all(bfield == 0.0_dp))
    if (force_free .and. all(v1 == v0)) then
      v_event = v0
      if (.not. all(ieee_is_finite(x_event))) status = boundary_event_invalid_geometry
      return
    end if
    speed_scale = sqrt(event_speed2/chord_speed2)
    v_event = speed_scale*chord_velocity
    if (.not. all(ieee_is_finite(x_event)) .or. .not. all(ieee_is_finite(v_event))) then
      status = boundary_event_invalid_geometry
    end if
  end subroutine evaluate_boundary_state

  !> 単一barrier-open面のpotential-barrier式をevent位置とevent時速度で評価する。
  subroutine apply_potential_barrier_event( &
    mesh, sim, snapshot, event, boundary_contract, q, m, x, v, alive, escaped, status, redistribution_uniform &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(boundary_event_type), intent(in) :: event
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    real(dp), intent(in) :: q, m
    real(dp), intent(inout) :: x(3), v(3)
    logical, intent(inout) :: alive
    logical, intent(out) :: escaped
    integer(i32), intent(out) :: status
    real(dp), intent(in) :: redistribution_uniform(3)

    type(sim_config) :: action_sim
    type(boundary_event_type) :: action_event
    real(dp) :: phi_boundary, barrier_potential_v, outward_v, kinetic_normal, potential_barrier
    integer(i32) :: axis, face_index, barrier_open_count, ordinary_open_count
    logical :: high_side, reflect_open

    status = boundary_event_ok
    escaped = .false.
    barrier_open_count = 0_i32
    ordinary_open_count = 0_i32
    face_index = 0_i32
    high_side = .false.
    do axis = 1_i32, 3_i32
      if (btest(event%face_mask, 2_i32*axis - 2_i32) .and. event%face_bc(2_i32*axis - 1_i32) == bc_open) then
        if (open_face_uses_potential_barrier(2_i32*axis - 1_i32, boundary_contract)) then
          barrier_open_count = barrier_open_count + 1_i32
          face_index = 2_i32*axis - 1_i32
          high_side = .false.
        else
          ordinary_open_count = ordinary_open_count + 1_i32
        end if
      end if
      if (btest(event%face_mask, 2_i32*axis - 1_i32) .and. event%face_bc(2_i32*axis) == bc_open) then
        if (open_face_uses_potential_barrier(2_i32*axis, boundary_contract)) then
          barrier_open_count = barrier_open_count + 1_i32
          face_index = 2_i32*axis
          high_side = .true.
        else
          ordinary_open_count = ordinary_open_count + 1_i32
        end if
      end if
    end do
    if (ordinary_open_count > 0_i32) then
      alive = .false.
      escaped = .true.
      return
    end if
    if (barrier_open_count > 1_i32) then
      status = particle_step_ambiguous_open_corner
      return
    end if

    action_sim = sim
    action_sim%open_boundary_model = 'escape'
    action_event = event
    if (barrier_open_count == 0_i32) then
      call apply_escape_reflect_periodic_event( &
        action_sim, action_event, x, v, alive, escaped, status, redistribution_uniform &
        )
      return
    end if
    if (.not. ieee_is_finite(q) .or. .not. ieee_is_finite(m) .or. m <= 0.0_dp) then
      status = particle_step_invalid_boundary
      return
    end if

    axis = (face_index + 1_i32)/2_i32
    barrier_potential_v = sim%phi_infty
    if (high_side) then
      if (boundary_contract%barrier_override_high(axis)) then
        barrier_potential_v = boundary_contract%barrier_potential_high_v(axis)
      end if
    else
      if (boundary_contract%barrier_override_low(axis)) then
        barrier_potential_v = boundary_contract%barrier_potential_low_v(axis)
      end if
    end if
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
    potential_barrier = q*(barrier_potential_v - phi_boundary)
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
    call apply_escape_reflect_periodic_event( &
      action_sim, action_event, x, v, alive, escaped, status, redistribution_uniform &
      )
  end subroutine apply_potential_barrier_event

  !> eventに含まれるopen面へglobalまたはspecies別のpotential barrierを適用するか返す。
  pure logical function event_uses_potential_barrier(event, boundary_contract) result(uses_barrier)
    type(boundary_event_type), intent(in) :: event
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i32) :: face

    uses_barrier = .false.
    do face = 1_i32, 6_i32
      if (.not. btest(event%face_mask, face - 1_i32) .or. event%face_bc(face) /= bc_open) cycle
      uses_barrier = open_face_uses_potential_barrier(face, boundary_contract)
      if (uses_barrier) return
    end do
  end function event_uses_potential_barrier

  !> face bit順のopen面がglobalまたはspecies別barrierの対象かを返す。
  pure logical function open_face_uses_potential_barrier(face, boundary_contract) result(uses_barrier)
    integer(i32), intent(in) :: face
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i32) :: axis

    axis = (face + 1_i32)/2_i32
    uses_barrier = boundary_contract%ordinary_open_model == external_open_potential_barrier
    if (mod(face, 2_i32) == 1_i32) then
      uses_barrier = uses_barrier .or. boundary_contract%barrier_override_low(axis)
    else
      uses_barrier = uses_barrier .or. boundary_contract%barrier_override_high(axis)
    end if
  end function open_face_uses_potential_barrier

  !> event作用が面内再配置用の一意なcounterを必要とするかを返す。
  pure logical function event_requires_redistribution_counter(event, boundary_contract) result(requires)
    type(boundary_event_type), intent(in) :: event
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i32) :: face
    logical :: has_open, has_ordinary_open, has_redistributed_reflect

    has_open = .false.
    has_ordinary_open = .false.
    has_redistributed_reflect = .false.
    do face = 1_i32, 6_i32
      if (.not. btest(event%face_mask, face - 1_i32)) cycle
      has_open = has_open .or. event%face_bc(face) == bc_open
      if (event%face_bc(face) == bc_open) then
        has_ordinary_open = has_ordinary_open .or. .not. open_face_uses_potential_barrier(face, boundary_contract)
      end if
      has_redistributed_reflect = has_redistributed_reflect .or. &
                                  event%face_bc(face) == bc_redistributed_reflect
    end do
    requires = has_redistributed_reflect .and. .not. has_ordinary_open .and. &
               (.not. has_open .or. event_uses_potential_barrier(event, boundary_contract))
  end function event_requires_redistribution_counter

  !> 粒子event識別子からOpenMP実行順序に依存しない一様乱数を構築する。
  pure subroutine generate_redistribution_uniform(seed, counter, event_index, uniform)
    integer(i32), intent(in) :: seed, event_index
    integer(i64), intent(in) :: counter(4)
    real(dp), intent(out) :: uniform(3)
    integer(i32) :: axis, item
    integer(i64) :: state
    integer(i64), parameter :: hash_modulus = 2147483647_i64

    do axis = 1_i32, 3_i32
      state = counter_hash_mix(104729_i64, int(seed, i64), 37_i64*int(axis, i64))
      state = counter_hash_mix(state, int(axis, i64), 1009_i64)
      do item = 1_i32, 4_i32
        state = counter_hash_mix(state, counter(item), 7919_i64*int(item, i64))
      end do
      state = counter_hash_mix(state, int(event_index, i64), 104729_i64)
      state = counter_hash_mix(state, counter(3), 13007_i64*int(axis, i64))
      state = counter_hash_mix(state, counter(1), 433494437_i64)
      uniform(axis) = (real(state, dp) - 0.5_dp)/real(hash_modulus, dp)
    end do
  end subroutine generate_redistribution_uniform

  !> 31-bit素数体上の非線形写像でtuple成分をcounter hashへ混合する。
  pure integer(i64) function counter_hash_mix(state, value, salt) result(mixed)
    integer(i64), intent(in) :: state, value, salt
    integer(i64), parameter :: hash_modulus = 2147483647_i64
    integer(i64), parameter :: hash_multiplier = 48271_i64

    mixed = modulo(state*state, hash_modulus)
    mixed = modulo(mixed + modulo(hash_multiplier*state, hash_modulus), hash_modulus)
    mixed = modulo(mixed + modulo(value, hash_modulus), hash_modulus)
    mixed = modulo(mixed + modulo(salt, hash_modulus) + 1_i64, hash_modulus)
    if (mixed == 0_i64) mixed = 1_i64
  end function counter_hash_mix

  !> Candidate endpointがboxの全faceからstrictly interiorかを返すfast-path判定。
  pure logical function point_strictly_inside_box(sim, x) result(inside)
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: x(3)

    inside = all(x > sim%box_min) .and. all(x < sim%box_max)
  end function point_strictly_inside_box

  !> 粒子stepの有限値・質量・時間刻み入力を検証する。
  pure logical function valid_particle_step_input(x, v, b, q, m, dt) result(valid)
    real(dp), intent(in) :: x(3), v(3), b(3), q, m, dt

    valid = all(ieee_is_finite(x)) .and. all(ieee_is_finite(v)) .and. all(ieee_is_finite(b)) .and. &
            ieee_is_finite(q) .and. ieee_is_finite(m) .and. m > 0.0_dp .and. &
            ieee_is_finite(dt) .and. dt >= 0.0_dp
  end function valid_particle_step_input

end module bem_particle_stepper
