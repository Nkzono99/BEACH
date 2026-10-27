!> free-space の導体再配分と絶縁体電荷の保存を検証する。
program test_surface_models
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_mesh, only: init_mesh
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_surface_models, only: floating_conductor_solver_type
  use bem_types, only: mesh_type, sim_config, surface_model_conductor, surface_model_insulator
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  call test_init(5)

  call test_begin('all_insulator_mesh_leaves_the_solver_inactive')
  call test_all_insulator_mesh_leaves_the_solver_inactive()
  call test_end()

  call test_begin('fixed_charge_induces_conductor_dipole')
  call test_fixed_charge_induces_conductor_dipole()
  call test_end()

  call test_begin('uniform_field_induces_conductor_dipole')
  call test_uniform_field_induces_conductor_dipole()
  call test_end()

  call test_begin('interleaved_conductor_groups_conserve_charge_and_equalize')
  call test_interleaved_conductor_groups_conserve_charge_and_equalize()
  call test_end()

  call test_begin('large_charged_conductor_equalizes_after_pivoting')
  call test_large_conductor_pivot()
  call test_end()

  call test_summary()

contains

  !> 異なる面積の大きな要素と非ゼロ総電荷で、LU の途中の行交換を検証する。
  subroutine test_large_conductor_pivot()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(electrostatic_snapshot_type) :: snapshot
    type(floating_conductor_solver_type) :: solver
    real(dp) :: centers(3, 3), q0(3), phi1, phi2, v0(3, 3), v1(3, 3), v2(3, 3), half_size
    integer(i32) :: i
    integer(i32) :: mesh_ids(3), surface_models(3)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    centers(:, 3) = [-1.0d0, 0.0d0, 0.0d0]
    q0 = [1.0d-10, 0.0d0, 1.0d-10]
    mesh_ids = [1_i32, 1_i32, 2_i32]
    surface_models = [surface_model_conductor, surface_model_conductor, surface_model_insulator]
    centers = 100.0_dp*centers
    do i = 1, 3
      half_size = 5.0_dp*real(i, dp)
      v0(:, i) = centers(:, i) + [-half_size, -half_size, 0.0_dp]
      v1(:, i) = centers(:, i) + [half_size, 0.0_dp, 0.0_dp]
      v2(:, i) = centers(:, i) + [0.0_dp, half_size, 0.0_dp]
    end do
    call init_mesh(mesh, v0, v1, v2, q0=q0, elem_mesh_id0=mesh_ids, elem_surface_model0=surface_models)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals
    call init_free_space_snapshot(mesh, [0.0d0, 0.0d0, 0.0d0], sim, snapshot)

    call solver%init(mesh, sim, snapshot)
    call assert_true(solver%active(), 'conductor mesh must activate the solver')
    call solver%relax(mesh, sim, snapshot)

    phi1 = scaled_potential_at_elem(mesh, 1_i32, [0.0d0, 0.0d0, 0.0d0])
    phi2 = scaled_potential_at_elem(mesh, 2_i32, [0.0d0, 0.0d0, 0.0d0])
    call assert_close_dp(mesh%q_elem(1) + mesh%q_elem(2), 1.0d-10, 1.0d-22, 'induced conductor net charge mismatch')
    call assert_close_dp(phi1, phi2, 1.0d-20, 'conductor elements should be equipotential')
    call assert_close_dp(mesh%q_elem(3), q0(3), 0.0d0, 'fixed insulator charge must remain unchanged')
  end subroutine test_large_conductor_pivot

  subroutine test_all_insulator_mesh_leaves_the_solver_inactive()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(electrostatic_snapshot_type) :: snapshot
    type(floating_conductor_solver_type) :: solver
    real(dp) :: centers(3, 2), q0(2)
    integer(i32) :: mesh_ids(2), surface_models(2)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [1.0d0, 0.0d0, 0.0d0]
    q0 = [2.0d-12, -1.0d-12]
    mesh_ids = [1_i32, 1_i32]
    surface_models = [surface_model_insulator, surface_model_insulator]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)
    call init_free_space_snapshot(mesh, [1.0d0, 0.0d0, 0.0d0], sim, snapshot)

    call solver%init(mesh, sim, snapshot)
    call assert_true(.not. solver%active(), 'all-insulator mesh must not activate the conductor solver')
    call solver%relax(mesh, sim, snapshot)

    call assert_close_dp(mesh%q_elem(1), q0(1), 0.0d0, 'first insulator charge must remain unchanged')
    call assert_close_dp(mesh%q_elem(2), q0(2), 0.0d0, 'second insulator charge must remain unchanged')
  end subroutine test_all_insulator_mesh_leaves_the_solver_inactive

  subroutine test_fixed_charge_induces_conductor_dipole()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(electrostatic_snapshot_type) :: snapshot
    type(floating_conductor_solver_type) :: solver
    real(dp) :: centers(3, 3), q0(3), phi1, phi2
    integer(i32) :: mesh_ids(3), surface_models(3)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    centers(:, 3) = [-1.0d0, 0.0d0, 0.0d0]
    q0 = [0.0d0, 0.0d0, 1.0d-10]
    mesh_ids = [1_i32, 1_i32, 2_i32]
    surface_models = [surface_model_conductor, surface_model_conductor, surface_model_insulator]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)
    call init_free_space_snapshot(mesh, [0.0d0, 0.0d0, 0.0d0], sim, snapshot)

    call solver%init(mesh, sim, snapshot)
    call assert_true(solver%active(), 'conductor mesh must activate the solver')
    call solver%relax(mesh, sim, snapshot)

    phi1 = scaled_potential_at_elem(mesh, 1_i32, [0.0d0, 0.0d0, 0.0d0])
    phi2 = scaled_potential_at_elem(mesh, 2_i32, [0.0d0, 0.0d0, 0.0d0])
    call assert_close_dp(mesh%q_elem(1) + mesh%q_elem(2), 0.0d0, 1.0d-22, 'induced conductor net charge mismatch')
    call assert_close_dp(phi1, phi2, 1.0d-20, 'conductor elements should be equipotential')
    call assert_true(mesh%q_elem(1) < 0.0d0, 'near fixed positive charge should induce negative conductor charge')
    call assert_true(mesh%q_elem(2) > 0.0d0, 'far conductor side should induce positive charge')
    call assert_close_dp(mesh%q_elem(3), q0(3), 0.0d0, 'fixed insulator charge must remain unchanged')
  end subroutine test_fixed_charge_induces_conductor_dipole

  subroutine test_uniform_field_induces_conductor_dipole()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(electrostatic_snapshot_type) :: snapshot
    type(floating_conductor_solver_type) :: solver
    real(dp) :: centers(3, 2), q0(2), external_e(3), phi1, phi2
    integer(i32) :: mesh_ids(2), surface_models(2)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    q0 = [0.0d0, 0.0d0]
    mesh_ids = [1_i32, 1_i32]
    surface_models = [surface_model_conductor, surface_model_conductor]
    external_e = [1.0d0, 0.0d0, 0.0d0]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)
    call init_free_space_snapshot(mesh, external_e, sim, snapshot)

    call solver%init(mesh, sim, snapshot)
    call solver%relax(mesh, sim, snapshot)

    phi1 = scaled_potential_at_elem(mesh, 1_i32, external_e)
    phi2 = scaled_potential_at_elem(mesh, 2_i32, external_e)
    call assert_close_dp(mesh%q_elem(1) + mesh%q_elem(2), 0.0d0, 1.0d-24, 'uniform-field net charge mismatch')
    call assert_close_dp(phi1, phi2, 1.0d-20, 'uniform-field conductor should be equipotential')
    call assert_true(mesh%q_elem(1) < 0.0d0, 'left side should become negative for +x external field')
    call assert_true(mesh%q_elem(2) > 0.0d0, 'right side should become positive for +x external field')
  end subroutine test_uniform_field_induces_conductor_dipole

  subroutine test_interleaved_conductor_groups_conserve_charge_and_equalize()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(electrostatic_snapshot_type) :: snapshot
    type(floating_conductor_solver_type) :: solver
    real(dp) :: centers(3, 4), q0(4), phi1, phi2, phi3, phi4
    integer(i32) :: mesh_ids(4), surface_models(4)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [4.0d0, 0.0d0, 0.0d0]
    centers(:, 3) = [1.0d0, 0.0d0, 0.0d0]
    centers(:, 4) = [5.0d0, 0.0d0, 0.0d0]
    q0 = [2.0d-12, 4.0d-12, 0.0d0, 0.0d0]
    mesh_ids = [1_i32, 2_i32, 1_i32, 2_i32]
    surface_models = [surface_model_conductor, surface_model_conductor, &
                      surface_model_conductor, surface_model_conductor]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)
    call init_free_space_snapshot(mesh, [0.0d0, 0.0d0, 0.0d0], sim, snapshot)

    call solver%init(mesh, sim, snapshot)
    call solver%relax(mesh, sim, snapshot)

    phi1 = scaled_potential_at_elem(mesh, 1_i32, [0.0d0, 0.0d0, 0.0d0])
    phi2 = scaled_potential_at_elem(mesh, 2_i32, [0.0d0, 0.0d0, 0.0d0])
    phi3 = scaled_potential_at_elem(mesh, 3_i32, [0.0d0, 0.0d0, 0.0d0])
    phi4 = scaled_potential_at_elem(mesh, 4_i32, [0.0d0, 0.0d0, 0.0d0])
    call assert_close_dp(mesh%q_elem(1) + mesh%q_elem(3), 2.0d-12, 1.0d-23, &
                         'first conductor mesh_id total charge mismatch')
    call assert_close_dp(mesh%q_elem(2) + mesh%q_elem(4), 4.0d-12, 1.0d-23, &
                         'second conductor mesh_id total charge mismatch')
    call assert_close_dp(phi1, phi3, 1.0d-20, 'first conductor mesh_id should be equipotential')
    call assert_close_dp(phi2, phi4, 1.0d-20, 'second conductor mesh_id should be equipotential')
  end subroutine test_interleaved_conductor_groups_conserve_charge_and_equalize

  subroutine init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)
    type(mesh_type), intent(out) :: mesh
    real(dp), intent(in) :: centers(:, :)
    real(dp), intent(in) :: q0(:)
    integer(i32), intent(in) :: mesh_ids(:), surface_models(:)
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), parameter :: half_size = 0.05d0
    integer(i32) :: i, n

    n = int(size(centers, 2), kind=i32)
    allocate (v0(3, n), v1(3, n), v2(3, n))
    do i = 1, n
      v0(:, i) = centers(:, i) + [-half_size, -half_size, 0.0d0]
      v1(:, i) = centers(:, i) + [half_size, 0.0d0, 0.0d0]
      v2(:, i) = centers(:, i) + [0.0d0, half_size, 0.0d0]
    end do
    call init_mesh(mesh, v0, v1, v2, q0=q0, elem_mesh_id0=mesh_ids, elem_surface_model0=surface_models)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals
  end subroutine init_test_mesh

  subroutine init_free_space_snapshot(mesh, external_e, sim, snapshot)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: external_e(3)
    type(sim_config), intent(out) :: sim
    type(electrostatic_snapshot_type), intent(out) :: snapshot
    type(field_physics_config) :: field_config
    type(periodic2_physics_config) :: periodic_config
    type(panel_kernel_config) :: panel_config

    sim = sim_config()
    sim%field_solver = 'direct'
    sim%e0 = external_e
    field_config = field_physics_config(backend='direct', normalization='si')
    periodic_config = periodic2_physics_config()
    panel_config = panel_kernel_config(kernel_id='triangle_p0_exact_direct')
    call snapshot%init(mesh, sim, field_config, periodic_config, panel_config)
    call snapshot%refresh(mesh)
  end subroutine init_free_space_snapshot

  real(dp) function scaled_potential_at_elem(mesh, elem_i, external_e) result(phi)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elem_i
    real(dp), intent(in) :: external_e(3)
    integer(i32) :: elem_j

    phi = -dot_product(external_e, mesh%centers(:, elem_i))/k_coulomb
    do elem_j = 1, mesh%nelem
      phi = phi + mesh%q_elem(elem_j)*potential_coeff(mesh, elem_i, elem_j)
    end do
  end function scaled_potential_at_elem

  real(dp) function potential_coeff(mesh, elem_i, elem_j) result(coeff)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elem_i, elem_j
    type(panel_geometry_type) :: geometry
    real(dp) :: potential, field(3)
    integer(i32) :: status

    call init_panel_geometry(mesh%v0(:, elem_j), mesh%v1(:, elem_j), mesh%v2(:, elem_j), geometry, status)
    if (status /= panel_geometry_ok) error stop 'test conductor panel geometry initialization failed.'
    call panel_potential_field( &
      geometry, 1.0_dp, mesh%centers(:, elem_i), panel_side_principal_value, potential, field &
      )
    coeff = potential/k_coulomb
  end function potential_coeff

end program test_surface_models
