program test_electrostatic_snapshot
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, k_coulomb
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      outer_plasma_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_diagnostics_type
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_allclose_1d
  implicit none

  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(field_physics_config) :: field_config
  type(periodic2_physics_config) :: periodic_config
  type(panel_kernel_config) :: panel_config
  type(outer_plasma_config) :: outer_config
  type(kinetic_outer_plasma_options_type) :: kinetic_options
  type(electrostatic_snapshot_type) :: snapshot
  type(electrostatic_diagnostics_type) :: diagnostics
  type(panel_geometry_type) :: geometry
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), target(3), electric_field(3), potential, mesh_potential(1)
  real(dp) :: potential_without_self, expected_potential, self_potential, self_field(3)
  integer(i32) :: panel_status

  call test_init(5)
  v0(:, 1) = [-0.5_dp, -0.5_dp, 0.0_dp]
  v1(:, 1) = [0.5_dp, -0.5_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 0.5_dp, 0.0_dp]
  call init_mesh(mesh, v0(:, :1), v1(:, :1), v2(:, :1), q0=[0.0_dp])
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%e0 = [2.0_dp, -3.0_dp, 4.0_dp]
  field_config = field_physics_config(backend='direct', normalization='si')
  periodic_config = periodic2_physics_config()
  panel_config = panel_kernel_config()
  outer_config = outer_plasma_config()
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config, outer_config)
  call snapshot%refresh(mesh)

  call test_begin('prescribed_field_is_composed_once')
  target = [0.25_dp, -0.5_dp, 0.75_dp]
  call snapshot%eval_local_e(mesh, target, electric_field)
  call assert_allclose_1d(electric_field, sim%e0, 1.0e-14_dp, 'snapshot electric field mismatch')
  call snapshot%eval_local_phi(mesh, sim, target, potential)
  call assert_close_dp(potential, -dot_product(sim%e0, target), 1.0e-14_dp, 'snapshot potential mismatch')
  call test_end()

  call test_begin('periodic_split_gauss_closure')
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
  call init_mesh(mesh, v0(:, :1), v1(:, :1), v2(:, :1), q0=[2.0e-12_dp])
  mesh%elem_vacuum_sign = 1
  mesh%vacuum_normals = mesh%normals
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%field_bc_mode = 'periodic2'
  sim%softening = 0.0_dp
  sim%use_box = .true.
  sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  field_config = field_physics_config(backend='direct', normalization='si')
  panel_config = panel_kernel_config( &
                 source_model='triangle_p0', kernel_id='triangle_p0_exact_direct', &
                 surface_side_policy='per_element' &
                 )
  periodic_config = periodic2_physics_config( &
                    nonzero_mode_backend='panel_spectral_reference', zero_mode_policy='exclude_k0', &
                    lower_boundary_model='e_bottom_zero', reference_mode_layers=4, panel_quadrature_order=16, &
                    interface_phi_tolerance=1.0e12_dp, interface_field_tolerance=1.0e12_dp &
                    )
  outer_config = outer_plasma_config( &
                 model='kinetic_1d', interface_z=1.0_dp, &
                 debye_length=0.2_dp, thermal_voltage=10.0_dp, max_linearity_ratio=0.5_dp &
                 )
  kinetic_options = kinetic_outer_plasma_options_type( &
                    grid_points=17_i32, domain_length=1.0_dp, tail_length=0.2_dp, &
                    electron_charge=-1.0_dp, electron_mass=1.0_dp, electron_density_infinity=0.0_dp, &
                    electron_temperature_j=1.0_dp, ion_charge=1.0_dp, ion_mass=1.0_dp, &
                    ion_density_infinity=0.0_dp, ion_temperature_j=0.0_dp, ion_drift_infinity=1.0_dp &
                    )
  call snapshot%init( &
    mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options=kinetic_options &
    )
  call snapshot%refresh(mesh)
  call snapshot%get_diagnostics(diagnostics)
  call assert_close_dp( &
    sum(mesh%q_elem) - sim%box_max(1)*sim%box_max(2)*eps0*snapshot%outer%interface_field, &
    0.0_dp, 1.0e-24_dp, 'surface plus outer charge must close' &
    )
  call assert_close_dp(snapshot%gauss_residual, 0.0_dp, 1.0e-24_dp, 'snapshot Gauss residual mismatch')
  call assert_close_dp( &
    diagnostics%interface_potential, snapshot%outer%interface_potential, 1.0e-24_dp, &
    'diagnostic interface potential mismatch' &
    )
  call snapshot%eval_local_phi(mesh, sim, mesh%centers(:, 1), potential)
  call snapshot%eval_local_phi_without_primary_self(mesh, sim, 1_i32, potential_without_self)
  call init_panel_geometry(mesh%v0(:, 1), mesh%v1(:, 1), mesh%v2(:, 1), geometry, panel_status)
  if (panel_status /= panel_geometry_ok) error stop 'periodic test panel geometry initialization failed.'
  call panel_potential_field( &
    geometry, mesh%q_elem(1), mesh%centers(:, 1), panel_side_principal_value, self_potential, self_field &
    )
  call assert_close_dp( &
    potential_without_self, potential - self_potential, 1.0e-12_dp*max(1.0_dp, abs(potential)), &
    'panel spectral primary self exclusion mismatch' &
    )
  call test_end()

  call test_begin('mesh_potential_uses_same_snapshot')
  call snapshot%compute_mesh_potential(mesh, sim, mesh_potential)
  call snapshot%eval_local_phi(mesh, sim, mesh%centers(:, 1), potential)
  call assert_close_dp( &
    mesh_potential(1), potential, 1.0e-14_dp, &
    'snapshot mesh potential mismatch' &
    )
  call test_end()

  call test_begin('point_primary_self_is_excluded')
  v0(:, 1) = [-0.5_dp, -0.5_dp, 0.0_dp]
  v1(:, 1) = [0.5_dp, -0.5_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 0.5_dp, 0.0_dp]
  v0(:, 2) = [-0.5_dp, -0.5_dp, 1.0_dp]
  v1(:, 2) = [0.5_dp, -0.5_dp, 1.0_dp]
  v2(:, 2) = [0.0_dp, 0.5_dp, 1.0_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[2.0e-12_dp, -0.5e-12_dp])
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%softening = 0.2_dp
  sim%e0 = [0.0_dp, 0.0_dp, 3.0_dp]
  field_config = field_physics_config(backend='direct', normalization='si')
  periodic_config = periodic2_physics_config()
  panel_config = panel_kernel_config()
  outer_config = outer_plasma_config()
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config, outer_config)
  call snapshot%refresh(mesh)
  call snapshot%eval_local_phi(mesh, sim, mesh%centers(:, 1), potential)
  call snapshot%eval_local_phi_without_primary_self(mesh, sim, 1_i32, potential_without_self)
  expected_potential = potential - k_coulomb*mesh%q_elem(1)/sim%softening
  call assert_close_dp( &
    potential_without_self, expected_potential, 1.0e-14_dp*max(1.0_dp, abs(expected_potential)), &
    'point primary self exclusion mismatch' &
    )
  call test_end()

  call test_begin('triangle_p0_primary_self_is_excluded')
  call init_mesh(mesh, v0, v1, v2, q0=[2.0e-12_dp, -0.5e-12_dp])
  mesh%elem_vacuum_sign = 1_i32
  mesh%vacuum_normals = mesh%normals
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%softening = 0.0_dp
  sim%e0 = [0.0_dp, 0.0_dp, 3.0_dp]
  field_config = field_physics_config(backend='direct', normalization='si')
  panel_config = panel_kernel_config( &
                 source_model='triangle_p0', kernel_id='triangle_p0_exact_direct', &
                 surface_side_policy='per_element' &
                 )
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config, outer_config)
  call snapshot%refresh(mesh)
  call snapshot%eval_local_phi(mesh, sim, mesh%centers(:, 1), potential)
  call snapshot%eval_local_phi_without_primary_self(mesh, sim, 1_i32, potential_without_self)
  call init_panel_geometry(mesh%v0(:, 1), mesh%v1(:, 1), mesh%v2(:, 1), geometry, panel_status)
  if (panel_status /= panel_geometry_ok) error stop 'test panel geometry initialization failed.'
  call panel_potential_field( &
    geometry, mesh%q_elem(1), mesh%centers(:, 1), panel_side_principal_value, self_potential, self_field &
    )
  expected_potential = potential - self_potential
  call assert_close_dp( &
    potential_without_self, expected_potential, 1.0e-14_dp*max(1.0_dp, abs(expected_potential)), &
    'triangle_p0 primary self exclusion mismatch' &
    )
  call test_end()
  call test_summary()
end program test_electrostatic_snapshot
