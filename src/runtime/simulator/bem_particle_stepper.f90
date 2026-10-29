!> 同一時刻の粒子状態から、空間電場を中点評価した1ステップ候補を構築する。
module bem_particle_stepper
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: mesh_type, particle_domain_config, hit_info, bc_open, bc_reflect, bc_periodic, bc_redistributed_reflect
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_pusher, only: boris_push
  use bem_collision, only: find_first_hit
  use bem_status_codes, only: collision_query_ok, boundary_event_ok, boundary_event_invalid_geometry, particle_step_ok, &
                              particle_step_invalid_boundary, particle_step_multiple_box_events, &
                              particle_step_ambiguous_open_corner
  use bem_boundary, only: boundary_event_type, find_first_boundary_event, &
                          apply_escape_reflect_periodic_event
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type
  use bem_irradiation_cell_contract, only: irradiation_cell_contract_type
  use bem_counter_rng, only: counter_uniforms
  implicit none
  private

  !> 適応サブステップが step を半分にする最小の仕事 [J]。電場変化が相対的に大きくても、
  !! step の仕事に換算してこれより小さければ粒子のエネルギーには効かないので分けない。
  real(dp), parameter :: particle_substep_energy_floor_j = 1.0e-4_dp*1.602176634e-19_dp

  type, public :: particle_step_result
    real(dp) :: x(3) = 0.0_dp
    real(dp) :: v(3) = 0.0_dp
    logical :: absorbed = .false.
    logical :: escaped_boundary = .false.
    integer(i32) :: elem_idx = -1_i32
    integer(i32) :: status = particle_step_ok
    integer(i32) :: wall_outward_event_count = 0_i32
    integer(i32) :: outer_barrier_return_count = 0_i32
    integer(i32) :: outer_barrier_escape_count = 0_i32
    real(dp) :: wall_outward_normal_kinetic_energy_j_sum = 0.0_dp
  end type particle_step_result

  !> 1 粒子 step の間変わらない文脈。一様磁場・粒子物性・境界 contract・面内再配置の RNG counter・照射セルと、
  !! 上部真空 Fourier 場での再試行か、periodic remainder の適応分割の深さを持つ。mesh・snapshot と解決済みの
  !! 領域設定は引数で参照し、ここには公開 routine の入口で一度だけ組み立てる小さな値だけを置く。
  type :: particle_step_context
    real(dp) :: bfield(3) = 0.0_dp
    real(dp) :: q = 0.0_dp
    real(dp) :: m = 0.0_dp
    type(external_boundary_contract_type) :: boundary_contract
    integer(i64) :: boundary_rng_counter(4) = 0_i64
    !> 呼び出し側が RNG counter を渡したか。面内再配置が要る event は counter なしでは解けない。
    logical :: has_boundary_rng_counter = .false.
    logical :: use_upper_panel_fourier = .false.
    integer(i32) :: adaptive_depth = 0_i32
    type(irradiation_cell_contract_type) :: cell
  end type particle_step_context

  public :: build_particle_step_candidate
  public :: particle_substep_needs_halving
  public :: particle_substep_energy_floor_j
  public :: advance_particle_step
  public :: advance_particle_step_upper_panel_fourier
  public :: resolve_particle_boundary_candidate

contains

  !> step 候補の中点電場 e_mid が参照電場 e_ref から change_max の割合より大きく変わり、
  !! その差が変位 x1-x0 にわたる仕事に換算して床値を超えるとき真。
  !! @param[in] energy_floor 仕事の床値 [J]（省略時は `particle_substep_energy_floor_j`）。
  pure logical function particle_substep_needs_halving(e_mid, e_ref, q, x0, x1, change_max, energy_floor) result(needs)
    real(dp), intent(in) :: e_mid(3), e_ref(3), q, x0(3), x1(3), change_max
    real(dp), intent(in), optional :: energy_floor
    real(dp) :: change, floor_j

    floor_j = particle_substep_energy_floor_j
    if (present(energy_floor)) floor_j = max(floor_j, energy_floor)
    change = norm2(e_mid - e_ref)
    needs = change > change_max*max(norm2(e_mid), norm2(e_ref)) .and. &
            abs(q)*change*norm2(x1 - x0) > floor_j
  end function particle_substep_needs_halving

  !> 予測中点の電場と一様磁場を使い、次時刻の位置・速度候補を返す。
  subroutine build_particle_step_candidate( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x1, v1, sampled_electric_field, &
    use_upper_panel_fourier, field_available &
    )
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    real(dp), intent(out) :: x1(3), v1(3)
    real(dp), intent(out), optional :: sampled_electric_field(3)
    logical, intent(in), optional :: use_upper_panel_fourier
    logical, intent(out), optional :: field_available
    real(dp) :: x_mid(3), e_mid(3)
    logical :: available, use_upper

    x_mid = x0 + 0.5d0*v0*dt
    call project_field_sample_to_box(sim, x_mid)
    use_upper = .false.
    if (present(use_upper_panel_fourier)) use_upper = use_upper_panel_fourier
    available = .true.
    if (use_upper) then
      call snapshot%eval_upper_panel_fourier_e(mesh, x_mid, e_mid, available)
    else
      call snapshot%eval_local_e(mesh, x_mid, e_mid)
    end if
    if (present(field_available)) field_available = available
    if (.not. available) then
      x1 = x0
      v1 = v0
      if (present(sampled_electric_field)) sampled_electric_field = 0.0_dp
      return
    end if
    call boris_push(x0, v0, q, m, dt, e_mid, bfield, x1, v1)
    if (present(sampled_electric_field)) sampled_electric_field = e_mid
  end subroutine build_particle_step_candidate

  !> 境界を越えるcandidateでも、場評価点はsolverのprimitive target box内に保つ。
  pure subroutine project_field_sample_to_box(sim, position)
    class(particle_domain_config), intent(in) :: sim
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
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, result, boundary_contract, boundary_rng_counter, &
    cell &
    )
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    type(particle_step_result), intent(out) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in), optional :: boundary_rng_counter(4)
    logical :: field_available

    field_available = .true.
    call advance_particle_step_impl( &
      mesh, sim, snapshot, &
      make_step_context(bfield, q, m, cell, .false., boundary_contract, boundary_rng_counter), &
      x0, v0, dt, result, field_available &
      )
  end subroutine advance_particle_step

  !> 上部真空域のfactorized P0-panel Fourier場だけで、失敗した1 stepを最初から再試行する。
  subroutine advance_particle_step_upper_panel_fourier( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, result, field_available, &
    boundary_contract, boundary_rng_counter, &
    cell &
    )
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    type(particle_step_result), intent(out) :: result
    logical, intent(out) :: field_available
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in), optional :: boundary_rng_counter(4)

    field_available = .true.
    call advance_particle_step_impl( &
      mesh, sim, snapshot, &
      make_step_context(bfield, q, m, cell, .true., boundary_contract, boundary_rng_counter), &
      x0, v0, dt, result, field_available &
      )
  end subroutine advance_particle_step_upper_panel_fourier

  !> 公開 routine の引数から step 文脈を組み立てる。contract を省略すると override なし、counter を省略すると
  !! 面内再配置の event は解けない（has_boundary_rng_counter=.false.）。
  function make_step_context( &
    bfield, q, m, cell, use_upper_panel_fourier, boundary_contract, boundary_rng_counter &
    ) result(ctx)
    real(dp), intent(in) :: bfield(3), q, m
    type(irradiation_cell_contract_type), intent(in) :: cell
    logical, intent(in) :: use_upper_panel_fourier
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in), optional :: boundary_rng_counter(4)
    type(particle_step_context) :: ctx

    ctx%bfield = bfield
    ctx%q = q
    ctx%m = m
    ctx%cell = cell
    ctx%use_upper_panel_fourier = use_upper_panel_fourier
    if (present(boundary_contract)) ctx%boundary_contract = boundary_contract
    ctx%has_boundary_rng_counter = present(boundary_rng_counter)
    if (ctx%has_boundary_rng_counter) ctx%boundary_rng_counter = boundary_rng_counter
  end function make_step_context

  !> periodic event の適応分割深さ（ctx%adaptive_depth）を内部だけで引き回す1 step実装。
  recursive subroutine advance_particle_step_impl(mesh, sim, snapshot, ctx, x0, v0, dt, result, field_available)
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(particle_step_context), intent(in) :: ctx
    real(dp), intent(in) :: x0(3), v0(3), dt
    type(particle_step_result), intent(out) :: result
    logical, intent(inout) :: field_available

    type(hit_info) :: hit
    real(dp) :: x_candidate(3), v_candidate(3), candidate_electric_field(3)
    integer(i32) :: query_status

    result = particle_step_result()
    result%x = x0
    result%v = v0
    if (.not. valid_particle_step_input(x0, v0, ctx%bfield, ctx%q, ctx%m, dt)) then
      result%status = particle_step_invalid_boundary
      return
    end if

    call build_particle_step_candidate( &
      mesh, sim, snapshot, ctx%bfield, x0, v0, ctx%q, ctx%m, dt, x_candidate, v_candidate, candidate_electric_field, &
      ctx%use_upper_panel_fourier, field_available &
      )
    if (.not. field_available) then
      result%status = particle_step_invalid_boundary
      return
    end if
    if (.not. all(ieee_is_finite(x_candidate)) .or. .not. all(ieee_is_finite(v_candidate))) then
      result%status = particle_step_invalid_boundary
      return
    end if

    call find_first_hit(mesh, x0, x_candidate, hit, sim=sim, status=query_status)
    if (query_status /= collision_query_ok) then
      if (sim%use_box .and. .not. point_strictly_inside_box(sim, x_candidate)) then
        call advance_particle_boundary_crossing( &
          mesh, sim, snapshot, ctx, x0, v0, dt, x_candidate, v_candidate, candidate_electric_field, result=result, &
          field_available=field_available &
          )
        return
      end if
      result%status = query_status
      return
    end if
    if (.not. sim%use_box .or. point_strictly_inside_box(sim, x_candidate)) then
      if (hit%has_hit) then
        call accept_particle_hit(v0, v_candidate, hit, result)
        return
      end if
      result%x = x_candidate
      result%v = v_candidate
      return
    end if

    call advance_particle_boundary_crossing( &
      mesh, sim, snapshot, ctx, x0, v0, dt, x_candidate, v_candidate, candidate_electric_field, hit, result, &
      field_available &
      )
  end subroutine advance_particle_step_impl

  !> 構築済みcandidateとその場評価値からbox eventを解決する。
  subroutine resolve_particle_boundary_candidate( &
    mesh, sim, snapshot, bfield, x0, v0, q, m, dt, x_candidate, v_candidate, hit, result, boundary_contract, &
    boundary_rng_counter, sampled_electric_field, &
    cell &
    )
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt, x_candidate(3), v_candidate(3)
    type(hit_info), intent(in), optional :: hit
    type(particle_step_result), intent(out) :: result
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract
    integer(i64), intent(in), optional :: boundary_rng_counter(4)
    real(dp), intent(in), optional :: sampled_electric_field(3)
    real(dp) :: active_electric_field(3), x_mid(3)
    logical :: field_available

    result = particle_step_result()
    result%x = x0
    result%v = v0
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
    end if
    if (.not. all(ieee_is_finite(active_electric_field))) then
      result%status = particle_step_invalid_boundary
      return
    end if
    field_available = .true.
    call advance_particle_boundary_crossing( &
      mesh, sim, snapshot, &
      make_step_context(bfield, q, m, cell, .false., boundary_contract, boundary_rng_counter), &
      x0, v0, dt, x_candidate, v_candidate, active_electric_field, hit, result, field_available &
      )
  end subroutine resolve_particle_boundary_candidate

  !> box crossing時だけevent用stateを確保し、periodic過多時はremainderを適応分割する。
  subroutine advance_particle_boundary_crossing( &
    mesh, sim, snapshot, ctx, x0, v0, dt, x_candidate, v_candidate, candidate_electric_field, hit, result, &
    field_available &
    )
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(particle_step_context), intent(in) :: ctx
    real(dp), intent(in) :: x0(3), v0(3), dt, x_candidate(3), v_candidate(3)
    real(dp), intent(in) :: candidate_electric_field(3)
    type(hit_info), intent(in), optional :: hit
    type(particle_step_result), intent(inout) :: result
    logical, intent(inout) :: field_available

    integer(i32), parameter :: max_boundary_events = 8_i32
    type(boundary_event_type) :: event
    type(hit_info) :: remainder_hit
    real(dp) :: wall_normal_velocity
    real(dp) :: x_start(3), v_start(3), x_trial(3), v_trial(3), x_event(3), v_event(3), segment_electric_field(3)
    real(dp) :: dt_segment, dt_remaining
    real(dp) :: redistribution_uniform(3)
    integer(i32) :: query_status, boundary_status, event_count
    logical :: alive, escaped, wall_outward_event, wall_barrier_event
    logical :: first_segment, requires_redistribution

    x_start = x0
    v_start = v0
    x_trial = x_candidate
    v_trial = v_candidate
    segment_electric_field = candidate_electric_field
    dt_segment = dt
    event_count = 0_i32
    first_segment = .true.

    do
      event = boundary_event_type()
      if (.not. point_strictly_inside_box(sim, x_trial)) then
        call find_first_boundary_event(sim, x_start, x_trial, event, boundary_status)
        if (boundary_status /= boundary_event_ok) then
          result%status = particle_step_invalid_boundary
          return
        end if
      end if
      if (.not. event%has_event) then
        if (first_segment .and. present(hit)) then
          if (hit%has_hit) then
            call accept_particle_hit(v_start, v_trial, hit, result)
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
            call accept_particle_hit(v_start, v_trial, hit, result)
            return
          end if
        end if
      end if

      call evaluate_boundary_state( &
        sim, ctx, event, x_start, v_start, x_trial, v_trial, dt_segment, segment_electric_field, x_event, v_event, &
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
        if (periodic_event_can_subdivide(sim, event) .and. ctx%adaptive_depth < 12_i32) then
          call advance_periodic_substeps(mesh, sim, snapshot, ctx, x_start, v_start, dt_segment, result, field_available)
          return
        end if
        result%status = particle_step_multiple_box_events
        return
      end if
      event_count = event_count + 1_i32
      wall_outward_event = ctx%cell%crosses_wall_outward(event%face_mask, event%face_bc, v_event)
      wall_barrier_event = wall_outward_event .and. &
                           open_face_uses_kinetic_barrier(ctx%cell%wall_face, ctx%boundary_contract)
      if (wall_outward_event) then
        wall_normal_velocity = ctx%cell%outward_normal_velocity(v_event)
        result%wall_outward_event_count = result%wall_outward_event_count + 1_i32
        result%wall_outward_normal_kinetic_energy_j_sum = &
          result%wall_outward_normal_kinetic_energy_j_sum + 0.5_dp*ctx%m*wall_normal_velocity*wall_normal_velocity
      end if
      requires_redistribution = event_requires_redistribution_counter(event, ctx%boundary_contract)
      redistribution_uniform = 0.0_dp
      if (requires_redistribution) then
        if (.not. ctx%has_boundary_rng_counter) then
          result%status = particle_step_invalid_boundary
          return
        end if
        call generate_redistribution_uniform( &
          sim%rng_seed, ctx%boundary_rng_counter, event_count, redistribution_uniform &
          )
      end if
      dt_remaining = (1.0_dp - event%fraction)*dt_segment
      alive = .true.
      escaped = .false.
      if (event_uses_kinetic_barrier(event, ctx%boundary_contract)) then
        call apply_kinetic_barrier_event( &
          mesh, sim, snapshot, ctx, event, x_event, v_event, alive, escaped, boundary_status, redistribution_uniform, &
          field_available &
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
      if (wall_barrier_event) then
        if (alive) then
          result%outer_barrier_return_count = result%outer_barrier_return_count + 1_i32
        else if (escaped) then
          result%outer_barrier_escape_count = result%outer_barrier_escape_count + 1_i32
        end if
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
        mesh, sim, snapshot, ctx%bfield, x_start, v_start, ctx%q, ctx%m, dt_segment, x_trial, v_trial, &
        segment_electric_field, ctx%use_upper_panel_fourier, field_available &
        )
      if (.not. field_available) then
        result%status = particle_step_invalid_boundary
        return
      end if
      if (.not. all(ieee_is_finite(x_trial)) .or. .not. all(ieee_is_finite(v_trial))) then
        result%status = particle_step_invalid_boundary
        return
      end if
      first_segment = .false.
    end do
  end subroutine advance_particle_boundary_crossing

  !> 一つのperiodic remainderを二分し、各半stepのbox/collision eventを順番に解く。
  subroutine advance_periodic_substeps(mesh, sim, snapshot, ctx, x0, v0, dt, result, field_available)
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(particle_step_context), intent(in) :: ctx
    real(dp), intent(in) :: x0(3), v0(3), dt
    type(particle_step_result), intent(inout) :: result
    logical, intent(inout) :: field_available

    type(particle_step_context) :: half_ctx
    type(particle_step_result) :: first_half, second_half
    integer(i32) :: prior_wall_events, prior_outer_returns, prior_outer_escapes
    real(dp) :: prior_wall_energy

    prior_wall_events = result%wall_outward_event_count
    prior_outer_returns = result%outer_barrier_return_count
    prior_outer_escapes = result%outer_barrier_escape_count
    prior_wall_energy = result%wall_outward_normal_kinetic_energy_j_sum
    half_ctx = ctx
    half_ctx%adaptive_depth = ctx%adaptive_depth + 1_i32
    call advance_particle_step_impl(mesh, sim, snapshot, half_ctx, x0, v0, 0.5_dp*dt, first_half, field_available)
    if (.not. field_available) then
      result = first_half
      return
    end if
    if (first_half%status /= particle_step_ok .or. first_half%absorbed .or. first_half%escaped_boundary) then
      result = first_half
      result%wall_outward_event_count = prior_wall_events + first_half%wall_outward_event_count
      result%outer_barrier_return_count = prior_outer_returns + first_half%outer_barrier_return_count
      result%outer_barrier_escape_count = prior_outer_escapes + first_half%outer_barrier_escape_count
      result%wall_outward_normal_kinetic_energy_j_sum = prior_wall_energy + &
                                                        first_half%wall_outward_normal_kinetic_energy_j_sum
      return
    end if

    call advance_particle_step_impl( &
      mesh, sim, snapshot, half_ctx, first_half%x, first_half%v, 0.5_dp*dt, second_half, field_available &
      )
    result = second_half
    result%wall_outward_event_count = prior_wall_events + first_half%wall_outward_event_count + &
                                      second_half%wall_outward_event_count
    result%outer_barrier_return_count = prior_outer_returns + first_half%outer_barrier_return_count + &
                                        second_half%outer_barrier_return_count
    result%outer_barrier_escape_count = prior_outer_escapes + first_half%outer_barrier_escape_count + &
                                        second_half%outer_barrier_escape_count
    result%wall_outward_normal_kinetic_energy_j_sum = prior_wall_energy + &
                                                      first_half%wall_outward_normal_kinetic_energy_j_sum + &
                                                      second_half%wall_outward_normal_kinetic_energy_j_sum
  end subroutine advance_periodic_substeps

  !> eventの全faceがperiodicで、適応分割中に乱数境界へ入らない場合だけ再試行する。
  pure logical function periodic_event_can_subdivide(sim, event) result(can_subdivide)
    class(particle_domain_config), intent(in) :: sim
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
  subroutine accept_particle_hit(va, vb, hit, result)
    real(dp), intent(in) :: va(3), vb(3)
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
    class(particle_domain_config), intent(in) :: sim
    real(dp), intent(in) :: xa(3), va(3), xb(3), vb(3)
    type(particle_step_result), intent(inout) :: result
    type(hit_info), intent(out) :: hit
    integer(i32), intent(out) :: query_status

    call find_first_hit(mesh, xa, xb, hit, sim=sim, status=query_status)
    if (query_status /= collision_query_ok) then
      result%status = query_status
      return
    end if
    if (.not. hit%has_hit) return

    call accept_particle_hit(va, vb, hit, result)
  end subroutine query_particle_chord

  !> Chord上のevent位置と、その接線方向・離散workに整合する速度を返す。
  subroutine evaluate_boundary_state( &
    sim, ctx, event, x0, v0, x1, v1, dt, electric_field, x_event, v_event, status &
    )
    class(particle_domain_config), intent(in) :: sim
    type(particle_step_context), intent(in) :: ctx
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(in) :: x0(3), v0(3), x1(3), v1(3), dt, electric_field(3)
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
    event_speed2 = sum(v0*v0) + 2.0_dp*(ctx%q/ctx%m)*dot_product(electric_field, x_event - x0)
    speed_tolerance = 256.0_dp*epsilon(1.0_dp)*max(tiny(1.0_dp), sum(v0*v0), chord_speed2, abs(event_speed2))
    if (.not. ieee_is_finite(chord_speed2) .or. .not. ieee_is_finite(event_speed2) .or. &
        chord_speed2 <= tiny(1.0_dp) .or. event_speed2 <= speed_tolerance) then
      v_event = 0.0_dp
      status = boundary_event_invalid_geometry
      return
    end if
    force_free = ctx%q == 0.0_dp .or. (all(electric_field == 0.0_dp) .and. all(ctx%bfield == 0.0_dp))
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

  !> 種別 kinetic 障壁を event 位置と event 時速度で評価する。
  subroutine apply_kinetic_barrier_event( &
    mesh, sim, snapshot, ctx, event, x, v, alive, escaped, status, redistribution_uniform, field_available &
    )
    type(mesh_type), intent(in) :: mesh
    class(particle_domain_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(particle_step_context), intent(in) :: ctx
    type(boundary_event_type), intent(in) :: event
    real(dp), intent(inout) :: x(3), v(3)
    logical, intent(inout) :: alive
    logical, intent(out) :: escaped
    integer(i32), intent(out) :: status
    real(dp), intent(in) :: redistribution_uniform(3)
    logical, intent(inout) :: field_available

    type(boundary_event_type) :: action_event
    real(dp) :: phi_boundary, barrier_potential_v, outward_v, kinetic_normal, barrier_energy_j
    real(dp) :: x_crossing(3), v_crossing(3), phi_return, kinetic_return
    integer(i32) :: axis, face_index, barrier_open_count, ordinary_open_count, reflect_bc, action_face_bc(6)
    logical :: high_side, reflect_open

    status = boundary_event_ok
    escaped = .false.
    barrier_open_count = 0_i32
    ordinary_open_count = 0_i32
    face_index = 0_i32
    high_side = .false.
    do axis = 1_i32, 3_i32
      if (btest(event%face_mask, 2_i32*axis - 2_i32) .and. event%face_bc(2_i32*axis - 1_i32) == bc_open) then
        if (open_face_uses_kinetic_barrier(2_i32*axis - 1_i32, ctx%boundary_contract)) then
          barrier_open_count = barrier_open_count + 1_i32
          face_index = 2_i32*axis - 1_i32
          high_side = .false.
        else
          ordinary_open_count = ordinary_open_count + 1_i32
        end if
      end if
      if (btest(event%face_mask, 2_i32*axis - 1_i32) .and. event%face_bc(2_i32*axis) == bc_open) then
        if (open_face_uses_kinetic_barrier(2_i32*axis, ctx%boundary_contract)) then
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

    ! 障壁で戻す面だけを反射に置き換えた face 別境界条件。他の面は box 設定のまま event と照合する。
    do axis = 1_i32, 3_i32
      action_face_bc(2_i32*axis - 1_i32) = sim%bc_low(axis)
      action_face_bc(2_i32*axis) = sim%bc_high(axis)
    end do
    action_event = event
    if (.not. ieee_is_finite(ctx%q) .or. .not. ieee_is_finite(ctx%m) .or. ctx%m <= 0.0_dp) then
      status = particle_step_invalid_boundary
      return
    end if

    axis = (face_index + 1_i32)/2_i32
    if (high_side) then
      barrier_potential_v = ctx%boundary_contract%barrier_potential_high_v(axis)
    else
      barrier_potential_v = ctx%boundary_contract%barrier_potential_low_v(axis)
    end if
    call evaluate_face_potential(x, phi_boundary)
    if (status /= boundary_event_ok) return
    if (high_side) then
      outward_v = v(axis)
    else
      outward_v = -v(axis)
    end if
    kinetic_normal = 0.5_dp*ctx%m*outward_v*outward_v
    barrier_energy_j = ctx%q*(barrier_potential_v - phi_boundary)
    if (.not. ieee_is_finite(kinetic_normal) .or. .not. ieee_is_finite(barrier_energy_j)) then
      status = particle_step_invalid_boundary
      return
    end if
    reflect_open = outward_v > 0.0_dp .and. barrier_energy_j > 0.0_dp .and. kinetic_normal < barrier_energy_j
    if (.not. reflect_open) then
      alive = .false.
      escaped = .true.
      return
    end if

    reflect_bc = bc_reflect
    if (open_face_redistributes_barrier_return(face_index, ctx%boundary_contract)) reflect_bc = bc_redistributed_reflect
    x_crossing = x
    v_crossing = v
    call reflect_on_barrier_face(reflect_bc)
    if (reflect_bc /= bc_redistributed_reflect .or. status /= boundary_event_ok .or. .not. alive) return

    ! 外部の静電場では全エネルギーが保存する。接線速度は保ったまま、戻り位置の局所電位で法線速度を合わせる。
    call evaluate_face_potential(x, phi_return)
    if (status /= boundary_event_ok) return
    kinetic_return = kinetic_normal + ctx%q*(phi_boundary - phi_return)
    if (.not. ieee_is_finite(kinetic_return)) then
      status = particle_step_invalid_boundary
      return
    end if
    if (kinetic_return > 0.0_dp) then
      v(axis) = sign(sqrt(2.0_dp*kinetic_return/ctx%m), v(axis))
    else
      ! 戻り位置へ届かない遅い粒子は外部での滞空が短く横へほとんど動かないので、横切り位置へ戻す。
      x = x_crossing
      v = v_crossing
      call reflect_on_barrier_face(bc_reflect)
    end if

  contains

    !> 障壁を判定する面上の局所電位。面平均ではなく、その位置の値を使う。
    subroutine evaluate_face_potential(position, phi)
      real(dp), intent(in) :: position(3)
      real(dp), intent(out) :: phi

      phi = 0.0_dp
      if (ctx%use_upper_panel_fourier) then
        call snapshot%eval_upper_panel_fourier_phi(mesh, position, phi, field_available)
        if (.not. field_available) then
          status = particle_step_invalid_boundary
          return
        end if
      else
        call snapshot%eval_local_phi(mesh, sim, position, phi)
      end if
      if (.not. ieee_is_finite(phi)) status = particle_step_invalid_boundary
    end subroutine evaluate_face_potential

    subroutine reflect_on_barrier_face(bc)
      integer(i32), intent(in) :: bc

      action_event%face_bc(face_index) = bc
      action_face_bc(face_index) = bc
      call apply_escape_reflect_periodic_event( &
        sim, action_event, x, v, alive, escaped, status, redistribution_uniform, face_bc=action_face_bc &
        )
    end subroutine reflect_on_barrier_face
  end subroutine apply_kinetic_barrier_event

  !> eventに含まれるopen面へspecies別のkinetic barrierを適用するか返す。
  pure logical function event_uses_kinetic_barrier(event, boundary_contract) result(uses_barrier)
    type(boundary_event_type), intent(in) :: event
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i32) :: face

    uses_barrier = .false.
    do face = 1_i32, 6_i32
      if (.not. btest(event%face_mask, face - 1_i32) .or. event%face_bc(face) /= bc_open) cycle
      uses_barrier = open_face_uses_kinetic_barrier(face, boundary_contract)
      if (uses_barrier) return
    end do
  end function event_uses_kinetic_barrier

  !> face bit順のopen面がspecies別kinetic barrierの対象かを返す。
  pure logical function open_face_uses_kinetic_barrier(face, boundary_contract) result(uses_barrier)
    integer(i32), intent(in) :: face
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i32) :: axis

    axis = (face + 1_i32)/2_i32
    if (mod(face, 2_i32) == 1_i32) then
      uses_barrier = boundary_contract%barrier_override_low(axis)
    else
      uses_barrier = boundary_contract%barrier_override_high(axis)
    end if
  end function open_face_uses_kinetic_barrier

  !> 外部障壁で戻る粒子を、障壁面に平行な座標について面内一様へ戻すかを返す。
  pure logical function open_face_redistributes_barrier_return(face, boundary_contract) result(redistributes)
    integer(i32), intent(in) :: face
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i32) :: axis

    axis = (face + 1_i32)/2_i32
    if (mod(face, 2_i32) == 1_i32) then
      redistributes = boundary_contract%barrier_override_low(axis) .and. boundary_contract%barrier_redistribute_low(axis)
    else
      redistributes = boundary_contract%barrier_override_high(axis) .and. &
                      boundary_contract%barrier_redistribute_high(axis)
    end if
  end function open_face_redistributes_barrier_return

  !> event作用が面内再配置用の一意なcounterを必要とするかを返す。
  pure logical function event_requires_redistribution_counter(event, boundary_contract) result(requires)
    type(boundary_event_type), intent(in) :: event
    type(external_boundary_contract_type), intent(in) :: boundary_contract
    integer(i32) :: face
    logical :: has_open, has_ordinary_open, has_redistributed_reflect, has_redistributed_barrier

    has_open = .false.
    has_ordinary_open = .false.
    has_redistributed_reflect = .false.
    has_redistributed_barrier = .false.
    do face = 1_i32, 6_i32
      if (.not. btest(event%face_mask, face - 1_i32)) cycle
      has_open = has_open .or. event%face_bc(face) == bc_open
      if (event%face_bc(face) == bc_open) then
        has_ordinary_open = has_ordinary_open .or. .not. open_face_uses_kinetic_barrier(face, boundary_contract)
        has_redistributed_barrier = has_redistributed_barrier .or. &
                                    open_face_redistributes_barrier_return(face, boundary_contract)
      end if
      has_redistributed_reflect = has_redistributed_reflect .or. &
                                  event%face_bc(face) == bc_redistributed_reflect
    end do
    requires = (has_redistributed_reflect .or. has_redistributed_barrier) .and. .not. has_ordinary_open .and. &
               (.not. has_open .or. event_uses_kinetic_barrier(event, boundary_contract))
  end function event_requires_redistribution_counter

  !> 粒子event識別子からOpenMP実行順序に依存しない一様乱数を構築する（PCG64 counter-based）。
  subroutine generate_redistribution_uniform(seed, counter, event_index, uniform)
    integer(i32), intent(in) :: seed, event_index
    integer(i64), intent(in) :: counter(4)
    real(dp), intent(out) :: uniform(3)

    call counter_uniforms(seed, [counter, int(event_index, i64), 1009_i64], 1_i64, uniform)
  end subroutine generate_redistribution_uniform

  !> Candidate endpointがboxの全faceからstrictly interiorかを返すfast-path判定。
  pure logical function point_strictly_inside_box(sim, x) result(inside)
    class(particle_domain_config), intent(in) :: sim
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
