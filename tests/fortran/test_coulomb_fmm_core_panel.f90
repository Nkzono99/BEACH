program test_coulomb_fmm_core_panel
  use bem_constants, only: k_coulomb
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_panel_plan, &
                                  update_state, eval_point, eval_potential_point, destroy_plan, destroy_state
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use test_support, only: test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  integer(i32), parameter :: nx = 8_i32, ny = 8_i32, nsrc = 2_i32*nx*ny
  type(fmm_options_type) :: options
  type(fmm_plan_type) :: plan
  type(fmm_state_type) :: state
  type(panel_geometry_type) :: image_geometry
  real(dp) :: v0(3, nsrc), v1(3, nsrc), v2(3, nsrc), q(nsrc)
  real(dp) :: targets(3, 5), field(3), field_ref(3), phi, phi_ref
  real(dp) :: field_error, potential_error, max_field_error, max_potential_error
  real(dp) :: source_field(3), source_phi, image_shift(3)
  integer(i32) :: target_idx, image_x, image_y, status

  call test_begin('full_panel_fmm_matches_direct_panel_oracle')
  call make_panel_grid(v0, v1, v2, q)
  options%theta = 0.45_dp
  options%leaf_max = 4_i32
  options%order = 5_i32
  options%target_box_min = [-1.2_dp, -1.2_dp, -0.5_dp]
  options%target_box_max = [1.2_dp, 1.2_dp, 2.5_dp]
  call build_panel_plan(plan, v0, v1, v2, options)
  call update_state(plan, state, q)
  targets(:, 1) = [0.03_dp, -0.07_dp, 0.18_dp]
  targets(:, 2) = [0.42_dp, 0.31_dp, 0.45_dp]
  targets(:, 3) = [-0.55_dp, 0.22_dp, 0.8_dp]
  targets(:, 4) = [0.15_dp, -0.63_dp, 1.4_dp]
  targets(:, 5) = [-0.2_dp, 0.1_dp, 2.1_dp]
  max_field_error = 0.0_dp
  max_potential_error = 0.0_dp
  do target_idx = 1_i32, size(targets, 2)
    call eval_point(plan, state, targets(:, target_idx), field)
    call eval_potential_point(plan, state, targets(:, target_idx), phi)
    call direct_panel_oracle(v0, v1, v2, q, targets(:, target_idx), field_ref, phi_ref)
    field_error = sqrt(sum((field - field_ref)**2))/max(1.0e-14_dp, sqrt(sum(field_ref*field_ref)))
    potential_error = abs(phi - phi_ref)/max(1.0e-14_dp, abs(phi_ref))
    max_field_error = max(max_field_error, field_error)
    max_potential_error = max(max_potential_error, potential_error)
  end do
  call assert_true(max_field_error < 2.0e-3_dp, 'panel FMM field error exceeds 2e-3')
  call assert_true(max_potential_error < 2.0e-3_dp, 'panel FMM potential error exceeds 2e-3')
  call destroy_state(state)
  call destroy_plan(plan)
  call test_end()

  call test_begin('periodic_panel_near_images_match_finite_shell')
  v0(:, 1) = [0.2_dp, 0.2_dp, 0.0_dp]
  v1(:, 1) = [0.4_dp, 0.2_dp, 0.0_dp]
  v2(:, 1) = [0.2_dp, 0.45_dp, 0.0_dp]
  q(1) = 2.0_dp
  options = fmm_options_type()
  options%leaf_max = 1_i32
  options%order = 4_i32
  options%use_periodic2 = .true.
  options%periodic_axes = [1_i32, 2_i32]
  options%periodic_len = [1.0_dp, 1.0_dp]
  options%periodic_image_layers = 1_i32
  options%periodic_far_correction = 'none'
  options%target_box_min = [0.0_dp, 0.0_dp, -1.0_dp]
  options%target_box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  call build_panel_plan(plan, v0(:, :1), v1(:, :1), v2(:, :1), options)
  call assert_true(any(abs(plan%near_source_shift1) > 0.5_dp) .or. any(abs(plan%near_source_shift2) > 0.5_dp), &
                   'periodic panel plan must retain shifted near images')
  call update_state(plan, state, q(:1))
  targets(:, 1) = [0.92_dp, 0.15_dp, 0.12_dp]
  call eval_point(plan, state, targets(:, 1), field)
  call eval_potential_point(plan, state, targets(:, 1), phi)
  field_ref = 0.0_dp
  phi_ref = 0.0_dp
  do image_y = -1_i32, 1_i32
    do image_x = -1_i32, 1_i32
      image_shift = [real(image_x, dp), real(image_y, dp), 0.0_dp]
      call init_panel_geometry( &
        v0(:, 1) + image_shift, v1(:, 1) + image_shift, v2(:, 1) + image_shift, image_geometry, status &
        )
      call assert_equal_i32(status, panel_geometry_ok, 'periodic image geometry status')
      call panel_potential_field( &
        image_geometry, q(1), targets(:, 1), panel_side_principal_value, source_phi, source_field &
        )
      field_ref = field_ref + source_field/k_coulomb
      phi_ref = phi_ref + source_phi/k_coulomb
    end do
  end do
  field_error = sqrt(sum((field - field_ref)**2))/sqrt(sum(field_ref*field_ref))
  potential_error = abs(phi - phi_ref)/abs(phi_ref)
  write (*, '(a,es12.5,a,es12.5)') &
    'periodic panel finite-shell errors: field=', field_error, ', potential=', potential_error
  call assert_true(field_error < 2.0e-3_dp, 'periodic panel field error exceeds 2e-3')
  call assert_true(potential_error < 2.0e-3_dp, 'periodic panel potential error exceeds 2e-3')
  call destroy_state(state)
  call destroy_plan(plan)
  call test_end()
  call test_summary()

contains

  subroutine make_panel_grid(v0_out, v1_out, v2_out, q_out)
    real(dp), intent(out) :: v0_out(3, nsrc), v1_out(3, nsrc), v2_out(3, nsrc), q_out(nsrc)
    integer(i32) :: ix, iy, element
    real(dp) :: x0, x1, y0, y1, dx, dy

    dx = 1.6_dp/real(nx, dp)
    dy = 1.6_dp/real(ny, dp)
    element = 0_i32
    do iy = 1_i32, ny
      y0 = -0.8_dp + real(iy - 1_i32, dp)*dy
      y1 = y0 + dy
      do ix = 1_i32, nx
        x0 = -0.8_dp + real(ix - 1_i32, dp)*dx
        x1 = x0 + dx
        element = element + 1_i32
        v0_out(:, element) = [x0, y0, 0.0_dp]
        v1_out(:, element) = [x1, y0, 0.0_dp]
        v2_out(:, element) = [x1, y1, 0.0_dp]
        q_out(element) = 1.0_dp + 0.2_dp*sin(real(element, dp))
        element = element + 1_i32
        v0_out(:, element) = [x0, y0, 0.0_dp]
        v1_out(:, element) = [x1, y1, 0.0_dp]
        v2_out(:, element) = [x0, y1, 0.0_dp]
        q_out(element) = 1.0_dp + 0.2_dp*cos(real(element, dp))
      end do
    end do
  end subroutine make_panel_grid

  subroutine direct_panel_oracle(v0_in, v1_in, v2_in, q_in, target, field_out, phi_out)
    real(dp), intent(in) :: v0_in(3, nsrc), v1_in(3, nsrc), v2_in(3, nsrc), q_in(nsrc), target(3)
    real(dp), intent(out) :: field_out(3), phi_out
    type(panel_geometry_type) :: geometry
    real(dp) :: source_field(3), source_phi
    integer(i32) :: element, status

    field_out = 0.0_dp
    phi_out = 0.0_dp
    do element = 1_i32, nsrc
      call init_panel_geometry(v0_in(:, element), v1_in(:, element), v2_in(:, element), geometry, status)
      call assert_equal_i32(status, panel_geometry_ok, 'grid panel geometry status')
      call panel_potential_field(geometry, q_in(element), target, panel_side_principal_value, source_phi, source_field)
      field_out = field_out + source_field/k_coulomb
      phi_out = phi_out + source_phi/k_coulomb
    end do
  end subroutine direct_panel_oracle
end program test_coulomb_fmm_core_panel
