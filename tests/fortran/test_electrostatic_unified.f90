program test_electrostatic_unified
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, &
                                      panel_kernel_config, outer_plasma_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_diagnostics_type, &
                                        electrostatic_restart_state_type
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, &
                          assert_close_dp, assert_allclose_1d
  implicit none

  type(mesh_type) :: mesh
  type(sim_config) :: sim_low, sim_high
  type(field_physics_config) :: field
  type(periodic2_physics_config) :: periodic
  type(panel_kernel_config) :: panel
  type(outer_plasma_config) :: outer_low, outer_high
  type(electrostatic_snapshot_type) :: snapshot_low, snapshot_high, restarted_snapshot
  type(electrostatic_diagnostics_type) :: diagnostics
  type(electrostatic_restart_state_type) :: restart_state
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), target(3)
  real(dp) :: field_low(3), field_high(3), potential_low, potential_high

  call configure_fixture()
  call test_init(3)

  call snapshot_low%init(mesh, sim_low, field, periodic, panel, outer_low)
  call snapshot_high%init(mesh, sim_high, field, periodic, panel, outer_high)
  call snapshot_low%refresh(mesh)
  call snapshot_high%refresh(mesh)

  call test_begin('ownership_interface_does_not_truncate_field_domain')
  target = [0.37_dp, 0.61_dp, 0.55_dp]
  call snapshot_low%eval_local_e(mesh, target, field_low)
  call snapshot_high%eval_local_e(mesh, target, field_high)
  call snapshot_low%eval_local_phi(mesh, sim_low, target, potential_low)
  call snapshot_high%eval_local_phi(mesh, sim_high, target, potential_high)
  call assert_allclose_1d(field_low, field_high, 2.0e-11_dp, 'ownership height changed local field')
  call assert_close_dp(potential_low, potential_high, 2.0e-11_dp, 'ownership height changed local potential')
  call assert_close_dp(snapshot_low%nonzero_tail%handoff_z, snapshot_high%nonzero_tail%handoff_z, &
                       1.0e-14_dp, 'physical response plane must be geometry-derived')
  call assert_close_dp( &
    snapshot_high%outer%interface_potential, &
    snapshot_high%outer%potential(snapshot_high%outer%profile_n)* &
    exp(-(outer_high%interface_z - snapshot_high%unified_grid%z(snapshot_high%unified_grid%n))/ &
        snapshot_high%unified_options%tail_length), &
    1.0e-14_dp, 'far handoff potential must continue the Robin tail' &
    )
  call test_end()

  call test_begin('unified_profile_closes_gauss_law')
  call snapshot_low%get_diagnostics(diagnostics)
  call assert_true(snapshot_low%use_unified_outer, 'unified outer flag must be active')
  call assert_true(diagnostics%applicable, 'unified linear response must be applicable')
  call assert_true(diagnostics%accessible_fraction_min < diagnostics%accessible_fraction_max, &
                   'rough fixture must exercise accessible-area variation')
  call assert_close_dp(diagnostics%response_start_z, snapshot_low%nonzero_tail%handoff_z, &
                       0.0_dp, 'response-plane diagnostic mismatch')
  call assert_close_dp(snapshot_low%gauss_residual, 0.0_dp, 2.0e-24_dp, 'unified Gauss residual mismatch')
  call assert_close_dp(snapshot_low%outer%interface_z, outer_low%interface_z, 1.0e-14_dp, &
                       'handoff state must be sampled at ownership interface')
  call test_end()

  call test_begin('unified_profile_restart_is_compatible')
  call snapshot_low%export_restart_state(4_i32, restart_state)
  call restarted_snapshot%init(mesh, sim_low, field, periodic, panel, outer_low)
  call restarted_snapshot%restore_outer_state(restart_state)
  call assert_allclose_1d(restarted_snapshot%outer%potential, snapshot_low%outer%potential, &
                          0.0_dp, 'unified restart profile mismatch')
  call restarted_snapshot%refresh(mesh)
  call assert_close_dp(restarted_snapshot%gauss_residual, 0.0_dp, 2.0e-24_dp, &
                       'restarted unified Gauss residual mismatch')
  call test_end()
  call test_summary()

contains

  subroutine configure_fixture()
    v0(:, 1) = [0.0_dp, 0.0_dp, 0.18_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.24_dp]
    v2(:, 1) = [1.0_dp, 1.0_dp, 0.30_dp]
    v0(:, 2) = [0.0_dp, 0.0_dp, 0.18_dp]
    v1(:, 2) = [1.0_dp, 1.0_dp, 0.30_dp]
    v2(:, 2) = [0.0_dp, 1.0_dp, 0.22_dp]
    call init_mesh(mesh, v0, v1, v2, q0=[2.0e-14_dp, -0.5e-14_dp])
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals

    sim_low = sim_config()
    sim_low%field_solver = 'direct'
    sim_low%field_bc_mode = 'periodic2'
    sim_low%softening = 0.0_dp
    sim_low%use_box = .true.
    sim_low%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    sim_low%box_max = [1.0_dp, 1.0_dp, 0.70_dp]
    sim_low%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim_low%bc_high = [bc_periodic, bc_periodic, bc_open]
    sim_high = sim_low
    sim_high%box_max(3) = 3.0_dp

    field = field_physics_config(backend='direct', normalization='si')
    panel = panel_kernel_config(source_model='triangle_p0', kernel_id='triangle_p0_exact_direct', &
                                surface_side_policy='per_element')
    periodic = periodic2_physics_config( &
               nonzero_mode_backend='panel_spectral_reference', zero_mode_policy='exclude_k0', &
               lower_boundary_model='e_bottom_zero', reference_mode_layers=3_i32, &
               panel_quadrature_order=8_i32, interface_sample_n=8_i32)
    outer_low = outer_plasma_config(model='unified_linear_response', interface_z=sim_low%box_max(3), &
                                    debye_length=0.20_dp, thermal_voltage=5.0_dp, &
                                    max_linearity_ratio=0.5_dp)
    outer_high = outer_low
    outer_high%interface_z = sim_high%box_max(3)
  end subroutine configure_fixture

end program test_electrostatic_unified
