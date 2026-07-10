!> 粒子ステップ候補の外部電場加算と空間電場の時間精度を検証する。
program test_particle_stepper
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config
  use bem_mesh, only: init_mesh
  use bem_field_solver, only: field_solver_type
  use bem_particle_stepper, only: build_particle_step_candidate
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_allclose_1d
  implicit none

  call test_init(2)

  call test_begin('uniform_e0_included_once')
  call test_uniform_e0_included_once()
  call test_end()

  call test_begin('charged_mesh_second_order_convergence')
  call test_charged_mesh_second_order_convergence()
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

end program test_particle_stepper
