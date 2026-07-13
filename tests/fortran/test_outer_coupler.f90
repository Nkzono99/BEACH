program test_outer_coupler
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      outer_plasma_config, coupling_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_restart_state_type
  use bem_outer_coupler, only: outer_coupler_type
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(field_physics_config) :: field_config
  type(periodic2_physics_config) :: periodic_config
  type(panel_kernel_config) :: panel_config
  type(outer_plasma_config) :: outer_config
  type(coupling_config) :: coupling
  type(electrostatic_snapshot_type) :: snapshot
  type(outer_coupler_type) :: coupler
  type(electrostatic_snapshot_type) :: restarted_snapshot
  type(outer_coupler_type) :: restarted_coupler
  type(electrostatic_restart_state_type) :: restart_state
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  logical :: updated

  call test_init(3)
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[1.0e-12_dp])
  mesh%elem_vacuum_sign = 1_i32
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
                    lower_boundary_model='symmetric_vacuum', reference_mode_layers=2, panel_quadrature_order=8, &
                    interface_phi_tolerance=1.0e12_dp, interface_field_tolerance=1.0e12_dp &
                    )
  outer_config = outer_plasma_config( &
                 model='linear_debye', interface_z=1.0_dp, debye_length=0.2_dp, thermal_voltage=10.0_dp, &
                 max_linearity_ratio=0.5_dp &
                 )
  coupling = coupling_config(outer_update_stride=2_i32)
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config, outer_config)
  call coupler%init(coupling)

  call test_begin('first_batch_and_stride_ownership')
  call coupler%refresh(snapshot, mesh, 1_i32, updated)
  call assert_true(updated, 'first batch must update outer state')
  call assert_close_dp(snapshot%gauss_residual, 0.0_dp, 1.0e-24_dp, 'first batch closure mismatch')
  mesh%q_elem = 2.0e-12_dp
  call coupler%refresh(snapshot, mesh, 2_i32, updated)
  call assert_true(.not. updated, 'stride=2 must hold outer state on batch 2')
  call assert_close_dp(snapshot%gauss_residual, 0.5e-12_dp, 1.0e-24_dp, 'held-state residual mismatch')
  call test_end()

  call test_begin('restart_preserves_stride_phase')
  mesh%q_elem = 1.0e-12_dp
  call coupler%init(coupling)
  call snapshot%refresh(mesh)
  call coupler%refresh(snapshot, mesh, 1_i32, updated)
  mesh%q_elem = 2.0e-12_dp
  call coupler%refresh(snapshot, mesh, 2_i32, updated)
  call snapshot%export_restart_state(coupler%last_outer_update_batch, restart_state)
  call restarted_snapshot%init(mesh, sim, field_config, periodic_config, panel_config, outer_config)
  call restarted_snapshot%restore_outer_state(restart_state)
  call restarted_coupler%init(coupling, restart_state%last_outer_update_batch)
  call restarted_coupler%refresh(restarted_snapshot, mesh, 3_i32, updated)
  call assert_true(updated, 'restart must preserve the batch-3 update schedule')
  call assert_close_dp(restarted_snapshot%gauss_residual, 0.0_dp, 1.0e-24_dp, 'restart closure mismatch')
  call test_end()

  call test_begin('scheduled_update_restores_closure')
  call coupler%refresh(snapshot, mesh, 3_i32, updated)
  call assert_true(updated, 'batch 3 must update outer state')
  call assert_close_dp(snapshot%gauss_residual, 0.0_dp, 1.0e-24_dp, 'scheduled closure mismatch')
  call test_end()

  call test_summary()
end program test_outer_coupler
