!> Public snapshot switching, geometry-only selection and sheath Neumann contracts.
program test_upper_fourier_switch
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi
  use bem_types, only: mesh_type, sim_config, bc_periodic, bc_open
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      validate_active_physics_config, validate_upper_fourier_height, physics_config_ok
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_allclose_1d
  implicit none
  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(field_physics_config) :: field
  type(periodic2_physics_config) :: periodic
  type(panel_kernel_config) :: panel
  type(electrostatic_snapshot_type) :: switched, baseline, automatic, repeated
  real(dp) :: a(3, 2), b(3, 2), c(3, 2), target(3), e(3), expected_e(3), phi, expected_phi, phi0, ez0
  real(dp) :: height, mismatch, phi_rms, ez_rms, top_field, old_phi, bottom_field, delta
  integer(i32) :: status, point
  character(len=256) :: message

  call test_init(4)
  a(:, 1) = [0.1_dp, 0.1_dp, 0.2_dp]
  b(:, 1) = [0.4_dp, 0.1_dp, 0.2_dp]
  c(:, 1) = [0.1_dp, 0.4_dp, 0.2_dp]
  a(:, 2) = [0.6_dp, 0.6_dp, 0.25_dp]
  b(:, 2) = [0.9_dp, 0.6_dp, 0.3_dp]
  c(:, 2) = [0.6_dp, 0.9_dp, 0.35_dp]
  call init_mesh(mesh, a, b, c, q0=[1.0e-12_dp, -0.3e-12_dp], elem_mesh_id0=[1_i32, 2_i32])
  mesh%elem_vacuum_sign = 1_i32
  mesh%vacuum_normals = mesh%normals
  sim%field_solver = 'fmm'
  sim%field_bc_mode = 'periodic2'
  sim%field_periodic_far_correction = 'cached_kneq0'
  sim%field_periodic_ewald_layers = 3_i32
  sim%field_periodic_cache_dir = '.beach_cache/test_upper_fourier_switch'
  sim%use_box = .true.
  sim%box_min = 0.0_dp
  sim%box_max = 1.0_dp
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = sim%bc_low
  sim%e0 = [0.1_dp, -0.2_dp, 0.3_dp]
  field%backend = 'fmm'
  panel%kernel_id = 'triangle_p0_exact_p2m_near'
  periodic%nonzero_mode_backend = 'cached_kneq0'
  periodic%zero_mode_policy = 'exclude_k0'
  periodic%lower_boundary_model = 'e_bottom_zero'

  call test_begin('switch_dispatches_all_local_fields_and_preserves_lower_region')
  call baseline%init(mesh, sim, field, periodic, panel)
  periodic%upper_fourier_switch = 'height'
  periodic%upper_fourier_switch_height = 0.7_dp
  call switched%init(mesh, sim, field, periodic, panel)
  call assert_true(.not. baseline%diagnostics%upper_fourier_switch_active, 'default switch must be off')
  do point = 1, 3
    target = [0.31_dp, 0.53_dp, 0.7_dp + 0.15_dp*real(point - 1, dp)]
    ! Independent panel Fourier quadrature, rather than the switch plan itself.
    call eval_periodic_nonzero_panel_reference(mesh, target, 1.0_dp, 1.0_dp, &
                                               switched%upper_fourier_plan%mode_layers, 32_i32, expected_phi, expected_e)
    call eval_periodic_zero_mode(switched%zero_plan, switched%zero_state, target(3), zero_mode_trace_plus, phi0, ez0)
    expected_phi = expected_phi + phi0 - dot_product(sim%e0, target)
    expected_e = expected_e + sim%e0
    expected_e(3) = expected_e(3) + ez0
    call switched%eval_local_phi(mesh, sim, target, phi)
    call switched%eval_local_e(mesh, target, e)
    ! Existing upper-vacuum plan uses seven-point panel quadrature; bound its independent quadrature error.
    call assert_close_dp(phi, expected_phi, 2.0e-6_dp, 'upper potential independent Fourier oracle')
    call assert_allclose_1d(e, expected_e, 3.0e-5_dp, 'upper field independent Fourier oracle')
    call switched%upper_fourier_plan%eval(mesh%q_elem, target, expected_phi, expected_e, status)
    call assert_close_dp(phi, expected_phi + phi0 - dot_product(sim%e0, target), 1.0e-13_dp, 'exact dispatch phi')
    expected_e = expected_e + sim%e0
    expected_e(3) = expected_e(3) + ez0
    call assert_allclose_1d(e, expected_e, 1.0e-13_dp, 'exact dispatch E')
  end do
  target(3) = 0.69_dp
  call baseline%eval_local_phi(mesh, sim, target, expected_phi)
  call switched%eval_local_phi(mesh, sim, target, phi)
  call assert_close_dp(phi, expected_phi, 0.0_dp, 'below switch legacy phi')
  call baseline%eval_local_e(mesh, target, expected_e)
  call switched%eval_local_e(mesh, target, e)
  call assert_allclose_1d(e, expected_e, 0.0_dp, 'below switch legacy E')
  call switched%measure_switch_mismatch(mesh, sim, phi_rms, ez_rms)
  call assert_true(phi_rms >= 0.0_dp .and. ez_rms >= 0.0_dp, 'real-charge mismatch diagnostics')
  call test_end()

  call test_begin('auto_selection_is_deterministic_and_obeys_decay_order')
  periodic%upper_fourier_switch = 'auto'
  call automatic%init(mesh, sim, field, periodic, panel)
  mesh%q_elem = [-5.0e-12_dp, 9.0e-12_dp]
  call repeated%init(mesh, sim, field, periodic, panel)
  mesh%q_elem = [1.0e-12_dp, -0.3e-12_dp]
  height = automatic%diagnostics%upper_fourier_switch_height_m
  mismatch = automatic%diagnostics%upper_fourier_switch_selected_mismatch_v
  call assert_true(height >= 0.6_dp .and. height < 1.0_dp, 'auto selected range')
  call assert_close_dp(height, repeated%diagnostics%upper_fourier_switch_height_m, 0.0_dp, 'deterministic height')
  call assert_close_dp(mismatch, repeated%diagnostics%upper_fourier_switch_selected_mismatch_v, 0.0_dp, &
                       'deterministic proxy mismatch')
  call assert_true(exp(-2.0_dp*pi*real(automatic%upper_fourier_plan%mode_layers + 1, dp)*(height - 0.35_dp)) &
                   <= 1.0e-6_dp, 'automatic truncation decay')
  call test_end()

  call test_begin('invalid_switches_are_preflight_errors')
  periodic%upper_fourier_switch = 'height'
  periodic%upper_fourier_switch_height = 0.44_dp
  call validate_upper_fourier_height(sim, periodic, 0.35_dp, status, message)
  call assert_true(status /= physics_config_ok, 'insufficient clearance accepted')
  periodic%upper_fourier_switch_height = 1.0_dp
  call validate_upper_fourier_height(sim, periodic, 0.35_dp, status, message)
  call assert_true(status /= physics_config_ok, 'H accepted as switch')
  periodic%upper_fourier_switch = 'auto'
  call validate_upper_fourier_height(sim, periodic, 0.9_dp, status, message)
  call assert_true(status /= physics_config_ok, 'empty auto interval accepted')
  periodic%nonzero_mode_backend = 'legacy_finite_images'
  call validate_active_physics_config(sim, field, periodic, panel, status, message)
  call assert_true(status /= physics_config_ok, 'wrong backend accepted')
  periodic%nonzero_mode_backend = 'cached_kneq0'
  periodic%upper_fourier_mode_layers = -1_i32
  call validate_active_physics_config(sim, field, periodic, panel, status, message)
  call assert_true(status /= physics_config_ok, 'negative mode layers accepted')
  call test_end()

  call test_begin('snapshot_sheath_top_field_matches_gauss_and_layer_potential')
  periodic%upper_fourier_switch = 'off'
  periodic%upper_fourier_mode_layers = 0_i32
  periodic%lower_boundary_model = 'e_top_from_sheath'
  call switched%init(mesh, sim, field, periodic, panel)
  call switched%set_matching_plane_gauge(mesh, 1.0_dp, 2.0_dp)
  top_field = -2.7_dp
  call switched%set_matching_plane_top_field(mesh, top_field)
  call eval_periodic_zero_mode(switched%zero_plan, switched%zero_state, 1.0_dp, zero_mode_trace_plus, phi, ez0)
  call assert_close_dp(ez0, top_field, 1.0e-12_dp*abs(top_field), 'Neumann field at H')
  call baseline%set_matching_plane_gauge(mesh, 1.0_dp, 2.0_dp)
  bottom_field = top_field - sum(mesh%q_elem)/(eps0*switched%zero_plan%area_xy)
  do point = 1, 5
    target(3) = 0.15_dp*real(point, dp)
    call eval_periodic_zero_mode(switched%zero_plan, switched%zero_state, target(3), zero_mode_trace_plus, phi, ez0)
    call eval_periodic_zero_mode(baseline%zero_plan, baseline%zero_state, target(3), zero_mode_trace_plus, old_phi, ez0)
    delta = bottom_field*(1.0_dp - target(3))
    call assert_close_dp(phi - old_phi, delta, 1.0e-12_dp, 'shift from e_bottom_zero')
  end do
  call switched%refresh(mesh)
  call eval_periodic_zero_mode(switched%zero_plan, switched%zero_state, 1.0_dp, zero_mode_trace_plus, phi, ez0)
  call assert_close_dp(ez0, top_field, 1.0e-12_dp*abs(top_field), 'refresh retains top field')
  call test_end()
  call test_summary()
end program test_upper_fourier_switch
