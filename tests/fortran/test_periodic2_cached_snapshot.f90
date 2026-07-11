program test_periodic2_cached_snapshot
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_loc, c_ptr
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_field_kernel_c, only: beach_kernel_build_panel, beach_kernel_create, beach_kernel_destroy, &
                                beach_kernel_eval_e, beach_kernel_eval_phi, &
                                beach_kernel_get_periodic_cache_info, beach_kernel_ok, &
                                beach_kernel_set_periodic_cache, beach_kernel_update_charges
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config, &
                                      outer_plasma_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, &
                          assert_allclose_1d, assert_equal_i32, delete_file_if_exists, remove_empty_directory
  implicit none

  interface
    integer(c_int) function beach_zero_mode_create(handle) bind(C, name='beach_zero_mode_create') result(status)
      import :: c_int, c_ptr
      type(c_ptr), intent(out) :: handle
    end function beach_zero_mode_create

    integer(c_int) function beach_zero_mode_destroy(handle) bind(C, name='beach_zero_mode_destroy') result(status)
      import :: c_int, c_ptr
      type(c_ptr), value :: handle
    end function beach_zero_mode_destroy

    integer(c_int) function beach_zero_mode_build(handle, nsrc, source_heights_ptr, area_xy) &
      bind(C, name='beach_zero_mode_build') result(status)
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: nsrc
      type(c_ptr), value :: source_heights_ptr
      real(c_double), value :: area_xy
    end function beach_zero_mode_build

    integer(c_int) function beach_zero_mode_update( &
      handle, nsrc, charge_ptr, e_bottom, z_gauge, phi_gauge &
      ) bind(C, name='beach_zero_mode_update') result(status)
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: nsrc
      type(c_ptr), value :: charge_ptr
      real(c_double), value :: e_bottom, z_gauge, phi_gauge
    end function beach_zero_mode_update

    integer(c_int) function beach_zero_mode_eval(handle, ntarget, z_ptr, trace, phi_ptr, ez_ptr) &
      bind(C, name='beach_zero_mode_eval') result(status)
      import :: c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: ntarget
      type(c_ptr), value :: z_ptr
      integer(c_int), value :: trace
      type(c_ptr), value :: phi_ptr, ez_ptr
    end function beach_zero_mode_eval
  end interface

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

  call test_init(3)
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

  call test_public_c_abi_acceptance()

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

  subroutine test_public_c_abi_acceptance()
    integer, parameter :: ntarget = 6
    type(c_ptr) :: kernel_handle, zero_handle
    integer(c_int) :: status
    integer(c_int), target :: periodic_axes(2), cache_hit, cache_build_count
    integer(c_int), target :: fingerprint_length, path_length
    real(c_double), target :: panel_v0(3, 2), panel_v1(3, 2), panel_v2(3, 2)
    real(c_double), target :: charge(2), source_heights(3, 2)
    real(c_double), target :: periodic_length(2), box_min(3), box_max(3)
    real(c_double), target :: target_points(3, ntarget), target_z(ntarget)
    real(c_double), target :: nonzero_e(3, ntarget), nonzero_phi(ntarget)
    real(c_double), target :: zero_e(ntarget), zero_phi(ntarget)
    character(kind=c_char), target :: cache_bytes(len(cache_dir))
    character(kind=c_char), target :: fingerprint_buffer(129), path_buffer(513)
    real(dp) :: snapshot_e(3), expected_e(3), snapshot_phi, expected_phi
    character(len=160) :: message
    integer :: point

    call test_begin('cached_snapshot_matches_public_c_abis_at_and_off_surface')

    panel_v0 = real(mesh%v0, c_double)
    panel_v1 = real(mesh%v1, c_double)
    panel_v2 = real(mesh%v2, c_double)
    charge = real(mesh%q_elem, c_double)
    source_heights(1, :) = panel_v0(3, :)
    source_heights(2, :) = panel_v1(3, :)
    source_heights(3, :) = panel_v2(3, :)
    periodic_axes = int(snapshot%nonzero_solver%fmm_core_options%periodic_axes, c_int)
    periodic_length = real(snapshot%nonzero_solver%fmm_core_options%periodic_len, c_double)
    box_min = real(snapshot%nonzero_solver%fmm_core_options%target_box_min, c_double)
    box_max = real(snapshot%nonzero_solver%fmm_core_options%target_box_max, c_double)
    call set_c_text(cache_bytes, cache_dir)

    status = beach_kernel_create(kernel_handle)
    call assert_c_ok(status, 'FieldKernel create status')
    status = beach_kernel_set_periodic_cache( &
             kernel_handle, c_loc(cache_bytes), int(size(cache_bytes), c_int), &
             real(sim%field_periodic_generation_tolerance, c_double) &
             )
    call assert_c_ok(status, 'FieldKernel cache status')
    status = beach_kernel_build_panel( &
             kernel_handle, int(mesh%nelem, c_int), c_loc(panel_v0), c_loc(panel_v1), c_loc(panel_v2), &
             real(snapshot%nonzero_solver%fmm_core_options%theta, c_double), &
             int(snapshot%nonzero_solver%fmm_core_options%leaf_max, c_int), &
             int(snapshot%nonzero_solver%fmm_core_options%order, c_int), &
             real(snapshot%nonzero_solver%fmm_core_options%softening, c_double), &
             1_c_int, c_loc(periodic_axes), c_loc(periodic_length), &
             int(snapshot%nonzero_solver%fmm_core_options%periodic_image_layers, c_int), 3_c_int, &
             real(snapshot%nonzero_solver%fmm_core_options%periodic_ewald_alpha, c_double), &
             int(snapshot%nonzero_solver%fmm_core_options%periodic_ewald_layers, c_int), &
             c_loc(box_min), c_loc(box_max) &
             )
    call assert_c_ok(status, 'FieldKernel cached panel build status')
    status = beach_kernel_get_periodic_cache_info( &
             kernel_handle, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), &
             int(size(fingerprint_buffer), c_int), c_loc(fingerprint_length), c_loc(path_buffer), &
             int(size(path_buffer), c_int), c_loc(path_length) &
             )
    call assert_c_ok(status, 'FieldKernel cache info status')
    call assert_equal_i32(int(cache_hit, i32), 1_i32, 'FieldKernel must reuse the snapshot cache')
    call assert_equal_i32(int(cache_build_count, i32), 0_i32, 'FieldKernel must not rebuild the snapshot cache')
    call assert_true( &
      trim(c_buffer_text(fingerprint_buffer, fingerprint_length)) == &
      trim(snapshot%diagnostics%periodic_cache_fingerprint), &
      'FieldKernel and snapshot cache fingerprints must match' &
      )
    call assert_true( &
      trim(c_buffer_text(path_buffer, path_length)) == trim(snapshot%diagnostics%periodic_cache_path), &
      'FieldKernel and snapshot cache paths must match' &
      )
    status = beach_kernel_update_charges(kernel_handle, int(mesh%nelem, c_int), c_loc(charge))
    call assert_c_ok(status, 'FieldKernel charge update status')

    status = beach_zero_mode_create(zero_handle)
    call assert_c_ok(status, 'zero-mode create status')
    status = beach_zero_mode_build( &
             zero_handle, int(mesh%nelem, c_int), c_loc(source_heights), &
             real(product(sim%box_max(1:2) - sim%box_min(1:2)), c_double) &
             )
    call assert_c_ok(status, 'zero-mode build status')
    status = beach_zero_mode_update( &
             zero_handle, int(mesh%nelem, c_int), c_loc(charge), 0.0_c_double, &
             minval(source_heights), 0.0_c_double &
             )
    call assert_c_ok(status, 'zero-mode charge update status')

    ! The first two points share a source height but lie outside the corresponding panel.
    target_points(:, 1) = [0.85_c_double, 0.15_c_double, 0.25_c_double]
    target_points(:, 2) = [0.15_c_double, 0.85_c_double, 0.65_c_double]
    target_points(:, 3) = [0.37_c_double, 0.61_c_double, 0.42_c_double]
    target_points(:, 4) = [0.19_c_double, 0.73_c_double, 0.88_c_double]
    target_points(:, 5) = [0.25_c_double, 0.25_c_double, 0.25_c_double]
    target_points(:, 6) = [0.75_c_double, 0.60_c_double, 0.65_c_double]
    target_z = target_points(3, :)
    status = beach_kernel_eval_e( &
             kernel_handle, int(ntarget, c_int), c_loc(target_points), c_loc(nonzero_e) &
             )
    call assert_c_ok(status, 'FieldKernel field evaluation status')
    status = beach_kernel_eval_phi( &
             kernel_handle, int(ntarget, c_int), c_loc(target_points), c_loc(nonzero_phi) &
             )
    call assert_c_ok(status, 'FieldKernel potential evaluation status')
    status = beach_zero_mode_eval( &
             zero_handle, int(ntarget, c_int), c_loc(target_z), int(zero_mode_trace_plus, c_int), &
             c_loc(zero_phi), c_loc(zero_e) &
             )
    call assert_c_ok(status, 'zero-mode plus-trace evaluation status')

    do point = 1, ntarget
      call snapshot%eval_local_e(mesh, real(target_points(:, point), dp), snapshot_e)
      call snapshot%eval_local_phi(mesh, sim, real(target_points(:, point), dp), snapshot_phi)
      expected_e = real(nonzero_e(:, point), dp) + sim%e0
      expected_e(3) = expected_e(3) + real(zero_e(point), dp)
      expected_phi = real(nonzero_phi(point) + zero_phi(point), dp) - &
                     dot_product(sim%e0, real(target_points(:, point), dp))
      write (message, '(a,i0)') 'snapshot/C-ABI field mismatch at point ', point
      call assert_allclose_1d(snapshot_e, expected_e, 1.0e-10_dp, trim(message))
      write (message, '(a,i0)') 'snapshot/C-ABI potential mismatch at point ', point
      call assert_close_dp(snapshot_phi, expected_phi, 1.0e-10_dp, trim(message))
    end do

    status = beach_zero_mode_destroy(zero_handle)
    call assert_c_ok(status, 'zero-mode destroy status')
    status = beach_kernel_destroy(kernel_handle)
    call assert_c_ok(status, 'FieldKernel destroy status')
    call test_end()
  end subroutine test_public_c_abi_acceptance

  subroutine assert_c_ok(actual, message)
    integer(c_int), intent(in) :: actual
    character(len=*), intent(in) :: message

    call assert_equal_i32(int(actual, i32), int(beach_kernel_ok, i32), message)
  end subroutine assert_c_ok

  subroutine set_c_text(output, value)
    character(kind=c_char), intent(out) :: output(:)
    character(len=*), intent(in) :: value
    integer :: i

    do i = 1, size(output)
      output(i) = achar(iachar(value(i:i)), kind=c_char)
    end do
  end subroutine set_c_text

  function c_buffer_text(buffer, length) result(value)
    character(kind=c_char), intent(in) :: buffer(:)
    integer(c_int), intent(in) :: length
    character(len=size(buffer)) :: value
    integer :: i

    value = ''
    do i = 1, int(length)
      value(i:i) = achar(iachar(buffer(i)))
    end do
  end function c_buffer_text

  subroutine configure_fixture(mesh_out, sim_out, field_out, periodic_out, panel_out, outer_out, a, b, c)
    type(mesh_type), intent(out) :: mesh_out
    type(sim_config), intent(out) :: sim_out
    type(field_physics_config), intent(out) :: field_out
    type(periodic2_physics_config), intent(out) :: periodic_out
    type(panel_kernel_config), intent(out) :: panel_out
    type(outer_plasma_config), intent(out) :: outer_out
    real(dp), intent(out) :: a(3, 2), b(3, 2), c(3, 2)

    a(:, 1) = [0.10_dp, 0.10_dp, 0.25_dp]
    b(:, 1) = [0.55_dp, 0.10_dp, 0.25_dp]
    c(:, 1) = [0.10_dp, 0.55_dp, 0.25_dp]
    a(:, 2) = [0.45_dp, 0.45_dp, 0.65_dp]
    b(:, 2) = [0.90_dp, 0.45_dp, 0.65_dp]
    c(:, 2) = [0.90_dp, 0.90_dp, 0.65_dp]
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
    sim_out%field_normalization = 'si'
    sim_out%tree_theta = 0.5_dp
    sim_out%has_tree_theta = .true.
    sim_out%tree_leaf_max = 8_i32
    sim_out%has_tree_leaf_max = .true.
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
