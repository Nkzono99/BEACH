program test_periodic2_cached_snapshot
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      outer_plasma_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, &
                          assert_allclose_1d, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(field_physics_config) :: field_config
  type(periodic2_physics_config) :: periodic_config
  type(panel_kernel_config) :: panel_config
  type(outer_plasma_config) :: outer_config
  type(electrostatic_snapshot_type) :: snapshot
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), target(3)
  real(dp) :: total_field(3), expected_field(3), nonzero_field(3), zero_field, zero_potential
  real(dp) :: total_potential, expected_potential, nonzero_potential
  real(dp) :: reference_field(3), reference_potential, field_error, potential_error, charge_scale
  character(len=512) :: cache_path
  character(len=*), parameter :: cache_dir = 'test_periodic2_cached_snapshot_tmp'

  call configure_fixture(mesh, sim, field_config, periodic_config, panel_config, outer_config, v0, v1, v2)
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config, outer_config)
  cache_path = snapshot%nonzero_solver%fmm_core_plan%periodic_cache_path
  call snapshot%refresh(mesh)
  target = [0.37_dp, 0.61_dp, 0.42_dp]

  call test_init(2)
  call test_begin('cached_snapshot_composes_kneq0_and_k0_once')
  call assert_true(snapshot%use_cached_kneq0 .and. snapshot%use_zero_mode, 'cached split flags must be active')
  call assert_true( &
    trim(snapshot%diagnostics%periodic_cache_fingerprint) == &
    trim(snapshot%nonzero_solver%fmm_core_plan%periodic_cache_fingerprint), &
    'snapshot must expose the periodic cache fingerprint' &
    )
  call evaluate_components(expected_field, expected_potential, nonzero_field, nonzero_potential, zero_field, zero_potential)
  call snapshot%eval_local_e(mesh, target, total_field)
  call snapshot%eval_local_phi(mesh, sim, target, total_potential)
  call assert_true(abs(zero_field) > 1.0e-12_dp, 'fixture must exercise a nonzero k=0 field')
  call assert_allclose_1d(total_field, expected_field, 1.0e-10_dp, 'snapshot field composition mismatch')
  call assert_close_dp(total_potential, expected_potential, 1.0e-10_dp, 'snapshot potential composition mismatch')
  call eval_periodic_nonzero_panel_reference(mesh, target, 1.0_dp, 1.0_dp, 12_i32, 16_i32, &
                                             reference_potential, reference_field)
  charge_scale = sum(abs(mesh%q_elem))/eps0
  field_error = sqrt(sum((nonzero_field - reference_field)**2))/charge_scale
  potential_error = abs(nonzero_potential - reference_potential)/charge_scale
  write (*, '(a,2(es12.5,1x))') 'panel cached kneq0 errors(field,potential)=', field_error, potential_error
  write (*, '(a,8(es12.5,1x))') 'panel cached/ref values=', nonzero_field, reference_field, &
    nonzero_potential, reference_potential
  call assert_true(field_error < 1.0e-1_dp, 'panel cached kneq0 field exceeds the charge-scale error contract')
  call assert_true(potential_error < 1.0e-1_dp, 'panel cached kneq0 potential exceeds the charge-scale error contract')
  call test_end()

  call test_begin('cached_snapshot_refreshes_both_components')
  mesh%q_elem = [1.0e-12_dp, -1.0e-12_dp]
  call snapshot%refresh(mesh)
  call evaluate_components(expected_field, expected_potential, nonzero_field, nonzero_potential, zero_field, zero_potential)
  call snapshot%eval_local_e(mesh, target, total_field)
  call snapshot%eval_local_phi(mesh, sim, target, total_potential)
  call assert_allclose_1d(total_field, expected_field, 1.0e-10_dp, 'refreshed snapshot field composition mismatch')
  call assert_close_dp(total_potential, expected_potential, 1.0e-10_dp, 'refreshed snapshot potential composition mismatch')
  call eval_periodic_nonzero_panel_reference(mesh, target, 1.0_dp, 1.0_dp, 12_i32, 16_i32, &
                                             reference_potential, reference_field)
  charge_scale = sum(abs(mesh%q_elem))/eps0
  field_error = sqrt(sum((nonzero_field - reference_field)**2))/charge_scale
  potential_error = abs(nonzero_potential - reference_potential)/charge_scale
  call assert_true(field_error < 1.0e-1_dp, 'neutral panel field exceeds the charge-scale error contract')
  call assert_true(potential_error < 1.0e-1_dp, 'neutral panel potential exceeds the charge-scale error contract')
  call test_end()

  call delete_file_if_exists(cache_path)
  call delete_file_if_exists(trim(cache_path)//'.lock')
  call remove_empty_directory(cache_dir)
  call test_summary()

contains

  subroutine evaluate_components(expected_e, expected_phi, nonzero_e, nonzero_phi, zero_e, zero_phi)
    real(dp), intent(out) :: expected_e(3), expected_phi, nonzero_e(3), nonzero_phi, zero_e, zero_phi
    call snapshot%nonzero_solver%eval_e(mesh, target, nonzero_e)
    call snapshot%nonzero_solver%eval_potential(mesh, sim, target, nonzero_phi)
    call eval_periodic_zero_mode(snapshot%zero_plan, snapshot%zero_state, target(3), zero_mode_trace_plus, zero_phi, zero_e)
    expected_e = nonzero_e + sim%e0
    expected_e(3) = expected_e(3) + zero_e
    expected_phi = nonzero_phi + zero_phi - dot_product(sim%e0, target)
  end subroutine evaluate_components

  subroutine configure_fixture(mesh_out, sim_out, field_out, periodic_out, panel_out, outer_out, a, b, c)
    type(mesh_type), intent(out) :: mesh_out
    type(sim_config), intent(out) :: sim_out
    type(field_physics_config), intent(out) :: field_out
    type(periodic2_physics_config), intent(out) :: periodic_out
    type(panel_kernel_config), intent(out) :: panel_out
    type(outer_plasma_config), intent(out) :: outer_out
    real(dp), intent(out) :: a(3, 2), b(3, 2), c(3, 2)

    a(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    b(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    c(:, 1) = [1.0_dp, 1.0_dp, 0.25_dp]
    a(:, 2) = [0.0_dp, 0.0_dp, 0.25_dp]
    b(:, 2) = [1.0_dp, 1.0_dp, 0.25_dp]
    c(:, 2) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh(mesh_out, a, b, c, q0=[1.0e-12_dp, 2.0e-12_dp])
    mesh_out%elem_vacuum_sign = 1_i32
    mesh_out%vacuum_normals = mesh_out%normals

    sim_out = sim_config()
    sim_out%field_solver = 'fmm'
    sim_out%field_bc_mode = 'periodic2'
    sim_out%field_periodic_far_correction = 'cached_kneq0'
    sim_out%field_periodic_image_layers = 1_i32
    sim_out%field_periodic_ewald_layers = 3_i32
    sim_out%field_periodic_cache_dir = cache_dir
    sim_out%field_periodic_generation_tolerance = 1.0e-8_dp
    sim_out%softening = 0.0_dp
    sim_out%e0 = [0.2_dp, -0.1_dp, 0.3_dp]
    sim_out%use_box = .true.
    sim_out%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    sim_out%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    sim_out%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim_out%bc_high = [bc_periodic, bc_periodic, bc_open]
    field_out = field_physics_config(backend='fmm', normalization='si')
    periodic_out = periodic2_physics_config( &
                   nonzero_mode_backend='cached_kneq0', zero_mode_policy='exclude_k0', &
                   lower_boundary_model='e_bottom_zero' &
                   )
    panel_out = panel_kernel_config( &
                source_model='triangle_p0', kernel_id='triangle_p0_exact_p2m_near', &
                surface_side_policy='per_element' &
                )
    outer_out = outer_plasma_config(model='none')
  end subroutine configure_fixture

end program test_periodic2_cached_snapshot
