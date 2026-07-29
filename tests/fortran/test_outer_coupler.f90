program test_outer_coupler
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_app_config, only: app_config, default_app_config, build_mesh_from_config
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      outer_plasma_config, coupling_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_restart_state_type
  use bem_outer_coupler, only: outer_coupler_type
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid
  use bem_zhao_steady_start, only: initialize_zhao_floating_steady_start
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  type(mesh_type) :: mesh
  type(mesh_type) :: steady_mesh
  type(sim_config) :: sim
  type(sim_config) :: steady_sim
  type(app_config) :: steady_app
  type(field_physics_config) :: field_config
  type(periodic2_physics_config) :: periodic_config
  type(periodic2_physics_config) :: steady_periodic
  type(panel_kernel_config) :: panel_config
  type(outer_plasma_config) :: outer_config
  type(coupling_config) :: coupling
  type(coupling_config) :: implicit_coupling
  type(coupling_config) :: steady_coupling
  type(kinetic_outer_plasma_options_type) :: kinetic_options
  type(outer_plasma_state_type) :: steady_state
  type(electrostatic_snapshot_type) :: snapshot
  type(electrostatic_snapshot_type) :: implicit_snapshot
  type(outer_coupler_type) :: coupler
  type(outer_coupler_type) :: implicit_coupler
  type(electrostatic_snapshot_type) :: restarted_snapshot
  type(outer_coupler_type) :: restarted_coupler
  type(electrostatic_restart_state_type) :: restart_state
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  real(dp) :: symmetric_charge, zero_bottom_charge, expected_charge, selected_area
  integer(i32) :: status, element
  character(len=256) :: message
  logical :: updated

  call test_init(6)
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[1.0e-12_dp])
  mesh%elem_vacuum_sign = 1_i32
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
                    lower_boundary_model='symmetric_vacuum', reference_mode_layers=2, panel_quadrature_order=8, &
                    interface_phi_tolerance=1.0e12_dp, interface_field_tolerance=1.0e12_dp &
                    )
  outer_config = outer_plasma_config( &
                 model='kinetic_1d', interface_z=1.0_dp, debye_length=0.2_dp, thermal_voltage=10.0_dp &
                 )
  kinetic_options = kinetic_outer_plasma_options_type( &
                    grid_points=17_i32, domain_length=1.0_dp, tail_length=0.2_dp, &
                    electron_charge=-1.0_dp, electron_mass=1.0_dp, electron_density_infinity=0.0_dp, &
                    electron_temperature_j=1.0_dp, ion_charge=1.0_dp, ion_mass=1.0_dp, &
                    ion_density_infinity=0.0_dp, ion_temperature_j=0.0_dp, ion_drift_infinity=1.0_dp &
                    )
  coupling = coupling_config(outer_update_stride=2_i32)
  call snapshot%init( &
    mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options=kinetic_options &
    )
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

  call test_begin('forced_same_batch_refresh_bypasses_stride')
  call coupler%refresh( &
    snapshot, mesh, 2_i32, updated, continuation_stage='post_commit', force_outer_update=.true. &
    )
  call assert_true(updated, 'post-commit corrector must force an outer update in the same batch')
  call assert_equal_i32(coupler%last_outer_update_batch, 2_i32, &
                        'forced refresh must record the actual outer update batch')
  call assert_close_dp(snapshot%gauss_residual, 0.0_dp, 1.0e-24_dp, &
                       'forced post-commit refresh did not restore closure')
  mesh%q_elem = 3.0e-12_dp
  call coupler%refresh(snapshot, mesh, 2_i32, updated)
  call assert_true(.not. updated, 'ordinary refresh must still obey the stride after a forced refresh')
  call assert_close_dp(snapshot%gauss_residual, 0.5e-12_dp, 1.0e-24_dp, &
                       'ordinary same-batch refresh unexpectedly changed the outer state')
  call test_end()

  call test_begin('restart_preserves_stride_phase')
  mesh%q_elem = 1.0e-12_dp
  call coupler%init(coupling)
  call snapshot%refresh(mesh)
  call coupler%refresh(snapshot, mesh, 1_i32, updated)
  mesh%q_elem = 2.0e-12_dp
  call coupler%refresh(snapshot, mesh, 2_i32, updated)
  call snapshot%export_restart_state(coupler%last_outer_update_batch, restart_state)
  call restarted_snapshot%init( &
    mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options=kinetic_options &
    )
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

  call test_begin('implicit_mean_reuses_post_refresh_until_mesh_changes')
  mesh%q_elem = 1.0e-12_dp
  implicit_coupling = coupling_config(update_mode='implicit_mean', outer_update_stride=1_i32)
  call implicit_snapshot%init( &
    mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options=kinetic_options &
    )
  call implicit_coupler%init(implicit_coupling)
  call implicit_coupler%refresh( &
    implicit_snapshot, mesh, 1_i32, updated, continuation_stage='pre_batch' &
    )
  call assert_true(updated, 'first implicit-mean pre-batch refresh must initialize the snapshot')
  call assert_true(implicit_coupler%snapshot_matches_mesh, &
                   'first implicit-mean refresh must mark the snapshot current')

  mesh%q_elem = 2.0e-12_dp
  call implicit_coupler%mark_mesh_changed()
  call implicit_coupler%refresh( &
    implicit_snapshot, mesh, 1_i32, updated, continuation_stage='post_implicit_mean', force_outer_update=.true. &
    )
  call assert_true(updated, 'implicit-mean post update must refresh the changed mesh')
  call assert_close_dp(implicit_snapshot%gauss_residual, 0.0_dp, 1.0e-24_dp, &
                       'implicit-mean post update did not restore closure')

  call implicit_coupler%refresh( &
    implicit_snapshot, mesh, 2_i32, updated, continuation_stage='pre_batch' &
    )
  call assert_true(.not. updated, 'current implicit-mean snapshot must be reused at the next pre-batch stage')
  call assert_equal_i32(implicit_coupler%last_outer_update_batch, 1_i32, &
                        'reused pre-batch snapshot must not record a duplicate outer update')

  mesh%q_elem = 3.0e-12_dp
  call implicit_coupler%mark_mesh_changed()
  call implicit_coupler%refresh( &
    implicit_snapshot, mesh, 2_i32, updated, continuation_stage='pre_batch' &
    )
  call assert_true(updated, 'stale implicit-mean snapshot must refresh before particle tracing')
  call assert_close_dp(implicit_snapshot%gauss_residual, 0.0_dp, 1.0e-24_dp, &
                       'stale implicit-mean pre-batch refresh did not restore closure')

  call implicit_coupler%init(implicit_coupling, last_outer_update_batch=2_i32)
  call assert_true(.not. implicit_coupler%snapshot_matches_mesh, &
                   'restart initialization must not assume that the restored snapshot field is current')
  call implicit_coupler%refresh( &
    implicit_snapshot, mesh, 3_i32, updated, continuation_stage='pre_batch' &
    )
  call assert_true(updated, 'restart must perform its first implicit-mean pre-batch refresh')

  call implicit_snapshot%export_restart_state(implicit_coupler%last_outer_update_batch, restart_state)
  call restarted_snapshot%init( &
    mesh, sim, field_config, periodic_config, panel_config, outer_config, kinetic_options=kinetic_options &
    )
  call restarted_snapshot%restore_outer_state(restart_state)
  call restarted_coupler%init(implicit_coupling, restart_state%last_outer_update_batch)
  call restarted_coupler%accept_restored_snapshot()
  call restarted_coupler%refresh( &
    restarted_snapshot, mesh, 4_i32, updated, continuation_stage='pre_batch' &
    )
  call assert_true(.not. updated, &
                   'an explicitly accepted restored implicit-mean snapshot must be reused')
  call assert_true(restarted_coupler%snapshot_matches_mesh, &
                   'accepted restored snapshot must remain marked current')
  call test_end()

  call test_begin('Zhao steady start seeds only the selected plane by panel area')
  call default_app_config(steady_app)
  steady_app%mesh_mode = 'template'
  steady_app%n_templates = 2_i32
  steady_app%templates(1)%enabled = .true.
  steady_app%templates(1)%kind = 'plane'
  steady_app%templates(1)%size_x = 1.0_dp
  steady_app%templates(1)%size_y = 1.0_dp
  steady_app%templates(1)%nx = 2_i32
  steady_app%templates(1)%ny = 1_i32
  steady_app%templates(1)%center = [0.5_dp, 0.5_dp, 0.25_dp]
  steady_app%templates(1)%surface_side_policy = 'normal_plus'
  steady_app%templates(2)%enabled = .true.
  steady_app%templates(2)%kind = 'sphere'
  steady_app%templates(2)%radius = 0.05_dp
  steady_app%templates(2)%n_lon = 8_i32
  steady_app%templates(2)%n_lat = 4_i32
  steady_app%templates(2)%center = [0.5_dp, 0.5_dp, 0.5_dp]
  steady_app%templates(2)%surface_side_policy = 'outward_closed'
  call build_mesh_from_config(steady_app, steady_mesh)

  steady_sim = sim_config()
  steady_sim%use_box = .true.
  steady_sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  steady_sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  steady_periodic = periodic2_physics_config(lower_boundary_model='symmetric_vacuum')
  steady_coupling = coupling_config(steady_start_mode='zhao_floating', steady_start_mesh_id=1_i32)
  kinetic_options = zhao_stationary_options()

  call initialize_zhao_floating_steady_start( &
    steady_mesh, steady_app%mesh_mode, steady_sim, steady_periodic, steady_coupling, kinetic_options, steady_state, &
    symmetric_charge, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'symmetric Zhao steady start failed: '//trim(message))
  expected_charge = 2.0_dp*eps0*steady_state%interface_field
  call assert_close_dp(symmetric_charge, expected_charge, 1.0e-30_dp, &
                       'symmetric-vacuum Zhao seed coefficient mismatch')
  call assert_close_dp(sum(steady_mesh%q_elem), symmetric_charge, 1.0e-30_dp, &
                       'distributed Zhao seed charge mismatch')
  call assert_true(maxval(abs(pack(steady_mesh%q_elem, steady_mesh%elem_mesh_id == 2_i32))) == 0.0_dp, &
                   'unselected sphere must remain neutral')
  selected_area = sum(steady_mesh%panel_area, mask=steady_mesh%elem_mesh_id == 1_i32)
  do element = 1_i32, steady_mesh%nelem
    if (steady_mesh%elem_mesh_id(element) /= 1_i32) cycle
    call assert_close_dp( &
      steady_mesh%q_elem(element), symmetric_charge*steady_mesh%panel_area(element)/selected_area, &
      1.0e-30_dp, 'Zhao seed was not distributed in proportion to panel area' &
      )
  end do

  steady_mesh%q_elem = 0.0_dp
  steady_periodic%lower_boundary_model = 'e_bottom_zero'
  call initialize_zhao_floating_steady_start( &
    steady_mesh, steady_app%mesh_mode, steady_sim, steady_periodic, steady_coupling, kinetic_options, steady_state, &
    zero_bottom_charge, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'zero-bottom Zhao steady start failed: '//trim(message))
  expected_charge = eps0*steady_state%interface_field
  call assert_close_dp(zero_bottom_charge, expected_charge, 1.0e-30_dp, &
                       'e_bottom_zero Zhao seed coefficient mismatch')
  call assert_close_dp(symmetric_charge, 2.0_dp*zero_bottom_charge, 1.0e-30_dp, &
                       'lower-boundary Zhao seed factors must differ by two')
  call assert_true(maxval(abs(pack(steady_mesh%q_elem, steady_mesh%elem_mesh_id == 2_i32))) == 0.0_dp, &
                   'unselected sphere must remain neutral after e_bottom_zero initialization')

  steady_mesh%q_elem = 0.0_dp
  steady_coupling%steady_start_mesh_id = 2_i32
  call initialize_zhao_floating_steady_start( &
    steady_mesh, steady_app%mesh_mode, steady_sim, steady_periodic, steady_coupling, kinetic_options, steady_state, &
    zero_bottom_charge, status, message &
    )
  call assert_equal_i32(status, outer_plasma_invalid, 'a sphere must be rejected as the steady-start plane')
  call assert_true(index(message, 'must be horizontal') > 0, &
                   'non-plane steady-start rejection lost its actionable reason')
  call assert_true(all(steady_mesh%q_elem == 0.0_dp), &
                   'failed steady-start validation must not mutate mesh charge')

  steady_coupling%steady_start_mesh_id = 1_i32
  call initialize_zhao_floating_steady_start( &
    steady_mesh, 'obj', steady_sim, steady_periodic, steady_coupling, kinetic_options, steady_state, &
    zero_bottom_charge, status, message &
    )
  call assert_equal_i32(status, outer_plasma_invalid, 'non-template mesh mode must be rejected')
  call assert_true(index(message, 'mesh.mode="template"') > 0, &
                   'non-template steady-start rejection lost its actionable reason')
  call assert_true(all(steady_mesh%q_elem == 0.0_dp), &
                   'non-template steady-start rejection must not mutate mesh charge')
  call test_end()

  call test_summary()

contains

  function zhao_stationary_options() result(value)
    type(kinetic_outer_plasma_options_type) :: value

    value%kinetic_closure = 'zhao_charge_driven'
    value%zhao_branch = 'auto'
    value%grid_points = 65_i32
    value%electron_charge = -qe
    value%electron_mass = 9.1093837139e-31_dp
    value%electron_density_infinity = 8.7e6_dp
    value%electron_temperature_j = 12.0_dp*qe
    value%electron_drift_infinity = 4.0529988897111727e5_dp
    value%ion_charge = qe
    value%ion_mass = 1.67262192369e-27_dp
    value%ion_density_infinity = 8.7e6_dp
    value%ion_temperature_j = 0.1_dp*qe
    value%ion_drift_infinity = 4.0529988897111727e5_dp
    value%photoelectron_charge = -qe
    value%photoelectron_mass = value%electron_mass
    value%photoelectron_temperature_j = 2.2_dp*qe
    value%photoelectron_reference_density = 64.0e6_dp
    value%zhao_alpha_deg = 60.0_dp
  end function zhao_stationary_options

end program test_outer_coupler
