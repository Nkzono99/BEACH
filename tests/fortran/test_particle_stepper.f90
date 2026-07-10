!> 粒子ステップ候補の外部電場加算と空間電場の時間精度を検証する。
program test_particle_stepper
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, bc_open, bc_reflect, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_field_solver, only: field_solver_type
  use bem_particle_stepper, only: build_particle_step_candidate, advance_particle_step, particle_step_result, &
                                  particle_step_ok, particle_step_invalid_boundary, particle_step_multiple_box_events, &
                                  particle_step_unsupported_barrier_corner
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_allclose_1d
  implicit none

  call test_init(11)

  call test_begin('uniform_e0_included_once')
  call test_uniform_e0_included_once()
  call test_end()

  call test_begin('charged_mesh_second_order_convergence')
  call test_charged_mesh_second_order_convergence()
  call test_end()

  call test_begin('advance_no_crossing_fast_path')
  call test_advance_no_crossing_fast_path()
  call test_end()

  call test_begin('advance_mesh_before_box_absorbs')
  call test_advance_mesh_before_box_absorbs()
  call test_end()

  call test_begin('advance_box_before_mesh_escapes')
  call test_advance_box_before_mesh_escapes()
  call test_end()

  call test_begin('advance_reflected_remainder_absorbs')
  call test_advance_reflected_remainder_absorbs()
  call test_end()

  call test_begin('advance_periodic_remainder')
  call test_advance_periodic_remainder()
  call test_end()

  call test_begin('advance_corner_reflection')
  call test_advance_corner_reflection()
  call test_end()

  call test_begin('advance_second_box_event_fails')
  call test_advance_second_box_event_fails()
  call test_end()

  call test_begin('advance_legacy_barrier_single_face_only')
  call test_advance_legacy_barrier_single_face_only()
  call test_end()

  call test_begin('advance_invalid_boundary_status_namespace')
  call test_advance_invalid_boundary_status_namespace()
  call test_end()

  call test_summary()

contains

  subroutine test_uniform_e0_included_once()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver = field_solver_type()
    real(dp) :: position_new(3), velocity_new(3)

    call init_single_element_mesh(mesh, 0.0d0)
    sim = sim_config()
    sim%field_solver = 'direct'
    sim%field_normalization = 'si'
    sim%softening = 0.5d0
    sim%e0 = [1.0d0, 0.0d0, 0.0d0]
    call field_solver%init(mesh, sim)
    call field_solver%refresh(mesh)

    call build_particle_step_candidate( &
      mesh, sim, field_solver, [0.0d0, 0.0d0, 0.0d0], &
      [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], &
      2.0d0, 1.0d0, 0.5d0, position_new, velocity_new &
      )

    call assert_allclose_1d( &
      velocity_new, [1.0d0, 0.0d0, 0.0d0], 0.0d0, &
      'uniform e0 should accelerate the candidate exactly once' &
      )
    call assert_allclose_1d( &
      position_new, [0.25d0, 0.0d0, 0.0d0], 0.0d0, &
      'uniform e0 candidate displacement should use same-time states' &
      )
  end subroutine test_uniform_e0_included_once

  subroutine test_charged_mesh_second_order_convergence()
    integer(i32), parameter :: nstep_values(3) = [4_i32, 8_i32, 16_i32]
    integer(i32), parameter :: reference_nstep = 256_i32
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver = field_solver_type()
    integer(i32) :: refinement
    real(dp) :: position(3), velocity(3), reference_position(3), reference_velocity(3)
    real(dp) :: position_error(3), error_ratio(2)

    call init_single_element_mesh(mesh, 1.0d-10)
    sim = sim_config()
    sim%field_solver = 'direct'
    sim%field_normalization = 'si'
    sim%softening = 0.5d0
    sim%e0 = 0.0d0
    call field_solver%init(mesh, sim)
    call field_solver%refresh(mesh)

    call integrate_candidate( &
      mesh, sim, field_solver, reference_nstep, reference_position, reference_velocity &
      )
    do refinement = 1_i32, 3_i32
      call integrate_candidate( &
        mesh, sim, field_solver, nstep_values(refinement), position, velocity &
        )
      position_error(refinement) = sqrt(sum((position - reference_position)**2))
    end do
    error_ratio = position_error(:2)/position_error(2:)

    call assert_true(all(ieee_is_finite(position_error)), 'charged-mesh position errors must be finite')
    call assert_true(all(position_error > 0.0d0), 'charged-mesh position errors must be positive')
    call assert_true( &
      position_error(1) > position_error(2) .and. position_error(2) > position_error(3), &
      'charged-mesh position errors must decrease under refinement' &
      )
    call assert_true(all(ieee_is_finite(error_ratio)), 'charged-mesh position error ratios must be finite')
    call assert_true(error_ratio(1) >= 3.2d0, 'dt to dt/2 charged-mesh position error ratio must be at least 3.2')
    call assert_true(error_ratio(2) >= 3.2d0, 'dt/2 to dt/4 charged-mesh position error ratio must be at least 3.2')
  end subroutine test_charged_mesh_second_order_convergence

  subroutine integrate_candidate(mesh, sim, field_solver, nstep, position, velocity)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_solver_type), intent(inout) :: field_solver
    integer(i32), intent(in) :: nstep
    real(dp), intent(out) :: position(3), velocity(3)
    integer(i32) :: step
    real(dp) :: dt, position_new(3), velocity_new(3)

    dt = 1.0d0/real(nstep, dp)
    position = [1.0d0, 0.25d0, 0.1d0]
    velocity = [0.1d0, -0.2d0, 0.05d0]
    do step = 1_i32, nstep
      call build_particle_step_candidate( &
        mesh, sim, field_solver, [0.0d0, 0.0d0, 0.0d0], &
        position, velocity, 1.0d0, 1.0d0, dt, position_new, velocity_new &
        )
      position = position_new
      velocity = velocity_new
    end do
  end subroutine integrate_candidate

  subroutine init_single_element_mesh(mesh, charge)
    type(mesh_type), intent(out) :: mesh
    real(dp), intent(in) :: charge
    real(dp) :: vertex0(3, 1), vertex1(3, 1), vertex2(3, 1), element_charge(1)

    vertex0(:, 1) = [-1.0d0, -1.0d0, 0.0d0]
    vertex1(:, 1) = [2.0d0, -1.0d0, 0.0d0]
    vertex2(:, 1) = [-1.0d0, 2.0d0, 0.0d0]
    element_charge(1) = charge
    call init_mesh(mesh, vertex0, vertex1, vertex2, q0=element_charge)
  end subroutine init_single_element_mesh

  subroutine init_x_plane_mesh(mesh, x_plane)
    type(mesh_type), intent(out) :: mesh
    real(dp), intent(in) :: x_plane
    real(dp) :: vertex0(3, 1), vertex1(3, 1), vertex2(3, 1), element_charge(1)

    vertex0(:, 1) = [x_plane, -10.0_dp, -10.0_dp]
    vertex1(:, 1) = [x_plane, 10.0_dp, -10.0_dp]
    vertex2(:, 1) = [x_plane, 0.0_dp, 10.0_dp]
    element_charge = 0.0_dp
    call init_mesh(mesh, vertex0, vertex1, vertex2, q0=element_charge)
  end subroutine init_x_plane_mesh

  subroutine init_box_stepper(mesh, sim, field_solver, mesh_x)
    type(mesh_type), intent(out) :: mesh
    type(sim_config), intent(out) :: sim
    type(field_solver_type), intent(out) :: field_solver
    real(dp), intent(in) :: mesh_x

    call init_x_plane_mesh(mesh, mesh_x)
    sim = sim_config()
    sim%field_solver = 'direct'
    sim%field_normalization = 'si'
    sim%softening = 0.5_dp
    sim%use_box = .true.
    sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    call field_solver%init(mesh, sim)
    call field_solver%refresh(mesh)
  end subroutine init_box_stepper

  subroutine test_advance_no_crossing_fast_path()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 10.0_dp)
    call advance_particle_step( &
      mesh, sim, field_solver, 0.0_dp*[1.0_dp, 1.0_dp, 1.0_dp], &
      [0.2_dp, 0.2_dp, 0.2_dp], [0.1_dp, 0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, 1.0_dp, result &
      )

    call assert_true(result%status == particle_step_ok, 'no-crossing step status mismatch')
    call assert_allclose_1d(result%x, [0.3_dp, 0.2_dp, 0.2_dp], 1.0e-14_dp, 'no-crossing position mismatch')
    call assert_allclose_1d(result%v, [0.1_dp, 0.0_dp, 0.0_dp], 0.0_dp, 'no-crossing velocity mismatch')
    call assert_true(.not. result%absorbed .and. .not. result%escaped_boundary, 'no-crossing terminal flags mismatch')
    call assert_true(result%field_eval_count == 1_i32, 'no-crossing path must evaluate E once')
    call assert_true(result%collision_query_count == 1_i32, 'no-crossing path must query collision once')
  end subroutine test_advance_no_crossing_fast_path

  subroutine test_advance_mesh_before_box_absorbs()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 0.8_dp)
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.5_dp, 0.2_dp, 0.2_dp], [1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, 1.0_dp, result &
      )

    call assert_true(result%status == particle_step_ok, 'mesh-before-box status mismatch')
    call assert_true(result%absorbed .and. .not. result%escaped_boundary, 'mesh-before-box should absorb')
    call assert_true(result%elem_idx == 1_i32, 'mesh-before-box element mismatch')
    call assert_close_dp(result%x(1), 0.8_dp, 1.0e-12_dp, 'mesh-before-box hit position mismatch')
  end subroutine test_advance_mesh_before_box_absorbs

  subroutine test_advance_box_before_mesh_escapes()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 1.2_dp)
    sim%bc_high(1) = bc_open
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.5_dp, 0.2_dp, 0.2_dp], [1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, 1.0_dp, result &
      )

    call assert_true(result%status == particle_step_ok, 'box-before-mesh status mismatch')
    call assert_true(.not. result%absorbed .and. result%escaped_boundary, 'box-before-mesh should escape')
    call assert_close_dp(result%x(1), 1.0_dp, 0.0_dp, 'open escape should stop at boundary event')
    call assert_true(result%field_eval_count == 1_i32, 'open crossing should not rebuild a remainder')
    call assert_true(result%collision_query_count == 1_i32, 'open crossing should query only to the boundary')
  end subroutine test_advance_box_before_mesh_escapes

  subroutine test_advance_reflected_remainder_absorbs()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 0.5_dp)
    sim%bc_high(1) = bc_reflect
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.8_dp, 0.2_dp, 0.2_dp], [1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, 1.0_dp, result &
      )

    call assert_true(result%status == particle_step_ok, 'reflected remainder status mismatch')
    call assert_true(result%absorbed, 'reflected remainder should hit the mesh')
    call assert_close_dp(result%x(1), 0.5_dp, 1.0e-12_dp, 'reflected remainder hit position mismatch')
    call assert_true(result%field_eval_count == 2_i32, 'reflected remainder should evaluate E twice')
    call assert_true(result%collision_query_count == 2_i32, 'reflected remainder should query both chords')
  end subroutine test_advance_reflected_remainder_absorbs

  subroutine test_advance_periodic_remainder()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 10.0_dp)
    sim%bc_low(1) = bc_periodic
    sim%bc_high(1) = bc_periodic
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.8_dp, 0.2_dp, 0.2_dp], [1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, 1.0_dp, result &
      )

    call assert_true(result%status == particle_step_ok, 'periodic remainder status mismatch')
    call assert_true(.not. result%absorbed .and. .not. result%escaped_boundary, 'periodic remainder should survive')
    call assert_close_dp(result%x(1), 0.8_dp, 1.0e-12_dp, 'periodic remainder position mismatch')
    call assert_close_dp(result%v(1), 1.0_dp, 0.0_dp, 'periodic remainder velocity mismatch')
  end subroutine test_advance_periodic_remainder

  subroutine test_advance_corner_reflection()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 10.0_dp)
    sim%bc_high(1:2) = bc_reflect
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.8_dp, 0.8_dp, 0.2_dp], [1.0_dp, 1.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, 1.0_dp, result &
      )

    call assert_true(result%status == particle_step_ok, 'corner reflection status mismatch')
    call assert_allclose_1d(result%x, [0.2_dp, 0.2_dp, 0.2_dp], 1.0e-12_dp, 'corner reflection position mismatch')
    call assert_allclose_1d(result%v, [-1.0_dp, -1.0_dp, 0.0_dp], 0.0_dp, 'corner reflection velocity mismatch')
  end subroutine test_advance_corner_reflection

  subroutine test_advance_second_box_event_fails()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result
    real(dp) :: x0(3), v0(3)

    call init_box_stepper(mesh, sim, field_solver, 10.0_dp)
    sim%bc_low(1) = bc_reflect
    sim%bc_high(1) = bc_reflect
    x0 = [0.9_dp, 0.2_dp, 0.2_dp]
    v0 = [3.0_dp, 0.0_dp, 0.0_dp]
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], x0, v0, 0.0_dp, 1.0_dp, 1.0_dp, result &
      )

    call assert_true(result%status == particle_step_multiple_box_events, 'second box event must fail closed')
    call assert_allclose_1d(result%x, x0, 0.0_dp, 'failed step must preserve initial position')
    call assert_allclose_1d(result%v, v0, 0.0_dp, 'failed step must preserve initial velocity')
    call assert_true(result%collision_query_count == 2_i32, 'second event path should query mesh before failing')
  end subroutine test_advance_second_box_event_fails

  subroutine test_advance_legacy_barrier_single_face_only()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 10.0_dp)
    sim%open_boundary_model = 'potential_barrier'
    sim%phi_infty = 10.0_dp
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.8_dp, 0.2_dp, 0.2_dp], [1.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, result &
      )
    call assert_true(result%status == particle_step_ok, 'single-face legacy barrier status mismatch')
    call assert_true(.not. result%escaped_boundary, 'large single-face barrier should reflect')
    call assert_close_dp(result%x(1), 0.2_dp, 1.0e-12_dp, 'legacy barrier reflected remainder mismatch')

    sim%bc_high(2) = bc_open
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.8_dp, 0.8_dp, 0.2_dp], [1.0_dp, 1.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, result &
      )
    call assert_true( &
      result%status == particle_step_unsupported_barrier_corner, &
      'simultaneous multi-open legacy barrier should fail closed' &
      )
  end subroutine test_advance_legacy_barrier_single_face_only

  subroutine test_advance_invalid_boundary_status_namespace()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_solver_type) :: field_solver
    type(particle_step_result) :: result

    call init_box_stepper(mesh, sim, field_solver, 10.0_dp)
    sim%open_boundary_model = 'unknown'
    call advance_particle_step( &
      mesh, sim, field_solver, [0.0_dp, 0.0_dp, 0.0_dp], &
      [0.8_dp, 0.2_dp, 0.2_dp], [1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, 1.0_dp, result &
      )
    call assert_true( &
      result%status == particle_step_invalid_boundary, &
      'boundary geometry status must map into the particle-step namespace' &
      )
  end subroutine test_advance_invalid_boundary_status_namespace

end program test_particle_stepper
