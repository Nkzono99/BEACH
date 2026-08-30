program test_electrostatic_snapshot
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_allclose_1d
  implicit none

  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(field_physics_config) :: field_config
  type(periodic2_physics_config) :: periodic_config
  type(panel_kernel_config) :: panel_config
  type(electrostatic_snapshot_type) :: snapshot
  type(panel_geometry_type) :: geometry
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), target(3), electric_field(3), potential, mesh_potential(1)
  real(dp) :: expected_field(3), nonzero_field(3), expected_potential, nonzero_potential
  real(dp) :: potential_without_self, source_potential, source_field(3), zero_field, zero_potential
  integer(i32) :: panel_status

  call test_init(4)
  v0(:, 1) = [-0.5_dp, -0.5_dp, 0.0_dp]
  v1(:, 1) = [0.5_dp, -0.5_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 0.5_dp, 0.0_dp]
  call init_mesh(mesh, v0(:, :1), v1(:, :1), v2(:, :1), q0=[0.0_dp])
  mesh%elem_vacuum_sign = 1_i32
  mesh%vacuum_normals = mesh%normals
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%e0 = [2.0_dp, -3.0_dp, 4.0_dp]
  field_config = field_physics_config(backend='direct', normalization='si')
  periodic_config = periodic2_physics_config()
  panel_config = panel_kernel_config(kernel_id='triangle_p0_exact_direct')
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config)
  call snapshot%refresh(mesh)

  call test_begin('free_snapshot_composes_prescribed_field_and_mesh_potential')
  target = [0.25_dp, -0.5_dp, 0.75_dp]
  call snapshot%eval_local_e(mesh, target, electric_field)
  call assert_allclose_1d(electric_field, sim%e0, 1.0e-14_dp, 'snapshot electric field mismatch')
  call snapshot%eval_local_phi(mesh, sim, target, potential)
  call assert_close_dp(potential, -dot_product(sim%e0, target), 1.0e-14_dp, 'snapshot potential mismatch')
  call snapshot%compute_mesh_potential(mesh, sim, mesh_potential)
  call assert_close_dp( &
    mesh_potential(1), -dot_product(sim%e0, mesh%centers(:, 1)), 1.0e-14_dp, &
    'snapshot mesh potential mismatch' &
    )
  call test_end()

  v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
  call init_mesh(mesh, v0(:, :1), v1(:, :1), v2(:, :1), q0=[2.0e-12_dp])
  mesh%elem_vacuum_sign = 1
  mesh%vacuum_normals = mesh%normals
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%field_bc_mode = 'periodic2'
  sim%use_box = .true.
  sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  field_config = field_physics_config(backend='direct', normalization='si')
  panel_config = panel_kernel_config( &
                 kernel_id='triangle_p0_exact_direct', surface_side_policy='per_element' &
                 )
  periodic_config = periodic2_physics_config( &
                    nonzero_mode_backend='panel_spectral_reference', zero_mode_policy='exclude_k0', &
                    lower_boundary_model='e_bottom_zero', reference_mode_layers=4, panel_quadrature_order=16 &
                    )
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config)
  call snapshot%refresh(mesh)
  call test_begin('panel_spectral_snapshot_composes_kneq0_and_k0')
  target = [0.37_dp, 0.43_dp, 0.75_dp]
  call eval_periodic_nonzero_panel_reference( &
    mesh, target, 1.0_dp, 1.0_dp, 4_i32, 16_i32, nonzero_potential, nonzero_field &
    )
  zero_field = mesh%q_elem(1)/eps0
  zero_potential = -zero_field*(target(3) - 0.25_dp)
  expected_field = nonzero_field
  expected_field(3) = expected_field(3) + zero_field
  expected_potential = nonzero_potential + zero_potential
  call snapshot%eval_local_e(mesh, target, electric_field)
  call snapshot%eval_local_phi(mesh, sim, target, potential)
  call assert_allclose_1d(electric_field, expected_field, 1.0e-12_dp, 'panel spectral field composition mismatch')
  call assert_close_dp( &
    potential, expected_potential, 1.0e-12_dp, 'panel spectral potential composition mismatch' &
    )
  call test_end()

  call test_begin('matching_plane_gauge_survives_refresh_and_clears')
  target = [0.37_dp, 0.43_dp, sim%box_max(3)]
  call snapshot%set_matching_plane_gauge(mesh, sim%box_max(3), -2.5_dp)
  call eval_periodic_nonzero_panel_reference( &
    mesh, target, 1.0_dp, 1.0_dp, 4_i32, 16_i32, nonzero_potential, nonzero_field &
    )
  call snapshot%eval_local_phi(mesh, sim, target, potential)
  call assert_close_dp(potential, nonzero_potential - 2.5_dp, 1.0e-12_dp, 'matching-plane gauge mismatch')
  call assert_close_dp( &
    snapshot%get_matching_plane_displacement(), 2.0e-12_dp, 1.0e-24_dp, &
    'matching-plane displacement mismatch' &
    )
  mesh%q_elem = 3.0e-12_dp
  call snapshot%refresh(mesh)
  call eval_periodic_nonzero_panel_reference( &
    mesh, target, 1.0_dp, 1.0_dp, 4_i32, 16_i32, nonzero_potential, nonzero_field &
    )
  call snapshot%eval_local_phi(mesh, sim, target, potential)
  call assert_close_dp(potential, nonzero_potential - 2.5_dp, 1.0e-12_dp, 'matching-plane gauge was not retained')
  call assert_close_dp( &
    snapshot%get_matching_plane_displacement(), 3.0e-12_dp, 1.0e-24_dp, &
    'matching-plane displacement was not refreshed' &
    )
  call snapshot%clear_matching_plane_gauge(mesh)
  zero_potential = -(mesh%q_elem(1)/eps0)*(target(3) - 0.25_dp)
  call snapshot%eval_local_phi(mesh, sim, target, potential)
  call assert_close_dp( &
    potential, nonzero_potential + zero_potential, 1.0e-12_dp, &
    'clearing matching-plane gauge did not restore the default gauge' &
    )
  call test_end()

  call test_begin('triangle_p0_primary_self_is_excluded')
  v0(:, 1) = [-0.5_dp, -0.5_dp, 0.0_dp]
  v1(:, 1) = [0.5_dp, -0.5_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 0.5_dp, 0.0_dp]
  v0(:, 2) = [-0.5_dp, -0.5_dp, 1.0_dp]
  v1(:, 2) = [0.5_dp, -0.5_dp, 1.0_dp]
  v2(:, 2) = [0.0_dp, 0.5_dp, 1.0_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[2.0e-12_dp, -0.5e-12_dp])
  mesh%elem_vacuum_sign = 1_i32
  mesh%vacuum_normals = mesh%normals
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%e0 = [1.0_dp, -2.0_dp, 3.0_dp]
  field_config = field_physics_config(backend='direct', normalization='si')
  periodic_config = periodic2_physics_config()
  panel_config = panel_kernel_config( &
                 kernel_id='triangle_p0_exact_direct', surface_side_policy='per_element' &
                 )
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config)
  call snapshot%refresh(mesh)
  call snapshot%eval_local_phi_without_primary_self(mesh, sim, 1_i32, potential_without_self)
  call init_panel_geometry(mesh%v0(:, 2), mesh%v1(:, 2), mesh%v2(:, 2), geometry, panel_status)
  if (panel_status /= panel_geometry_ok) error stop 'test panel geometry initialization failed.'
  call panel_potential_field( &
    geometry, mesh%q_elem(2), mesh%centers(:, 1), panel_side_principal_value, source_potential, source_field &
    )
  expected_potential = source_potential - dot_product(sim%e0, mesh%centers(:, 1))
  call assert_close_dp( &
    potential_without_self, expected_potential, 1.0e-14_dp*max(1.0_dp, abs(expected_potential)), &
    'triangle_p0 primary self exclusion mismatch' &
    )
  call test_end()
  call test_summary()
end program test_electrostatic_snapshot
