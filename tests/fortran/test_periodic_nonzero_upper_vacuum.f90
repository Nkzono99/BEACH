program test_periodic_nonzero_upper_vacuum
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_coulomb_fmm_periodic_nonzero_upper_vacuum, only: &
    periodic_nonzero_upper_vacuum_plan_type, upper_vacuum_eval_ok, &
    upper_vacuum_target_not_above_source
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, &
                          assert_allclose_1d, assert_equal_i32
  implicit none

  integer(i32), parameter :: mode_layers = 3_i32
  integer(i32), parameter :: high_mode_layers = 24_i32
  real(dp), parameter :: length_x = 1.3_dp, length_y = 0.9_dp
  real(dp), parameter :: charge(2) = [1.1e-12_dp, -0.4e-12_dp]
  real(dp), parameter :: updated_charge(2) = [-0.2e-12_dp, 0.9e-12_dp]
  type(mesh_type) :: mesh
  type(periodic_nonzero_upper_vacuum_plan_type) :: plan, high_mode_plan
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), target(3)
  real(dp) :: potential, direct_potential, seam_potential, field(3), direct_field(3), seam_field(3)
  real(dp) :: reference_potential, reference_field(3), source_top
  integer(i32) :: status

  call test_init(4)

  v0(:, 1) = [0.12_dp, 0.08_dp, 0.10_dp]
  v1(:, 1) = [0.58_dp, 0.11_dp, 0.24_dp]
  v2(:, 1) = [0.19_dp, 0.46_dp, 0.17_dp]
  v0(:, 2) = [0.67_dp, 0.22_dp, 0.06_dp]
  v1(:, 2) = [1.02_dp, 0.31_dp, 0.28_dp]
  v2(:, 2) = [0.74_dp, 0.72_dp, 0.35_dp]
  source_top = max(maxval(v0(3, :)), maxval(v1(3, :)), maxval(v2(3, :)))
  call init_mesh(mesh, v0, v1, v2, q0=charge)
  call plan%build(mesh, length_x, length_y, mode_layers)

  call test_begin('high_modes_near_source_top_match_duffy')
  call high_mode_plan%build(mesh, length_x, length_y, high_mode_layers)
  target = [0.37_dp, 0.41_dp, source_top + 1.0e-6_dp]
  call high_mode_plan%eval(mesh%q_elem, target, potential, field, status)
  call eval_periodic_nonzero_panel_reference( &
    mesh, target, length_x, length_y, high_mode_layers, 32_i32, reference_potential, reference_field &
    )
  call assert_equal_i32(status, upper_vacuum_eval_ok, 'near-top high-mode target was rejected')
  call assert_close_dp( &
    potential, reference_potential, 2.0e-3_dp*max(tiny(1.0_dp), abs(reference_potential)), &
    'high-mode near-top potential does not match the Duffy reference' &
    )
  call assert_allclose_1d( &
    field, reference_field, 3.0e-2_dp*max(tiny(1.0_dp), maxval(abs(reference_field))), &
    'high-mode near-top field does not match the Duffy reference' &
    )
  call test_end()

  call test_begin('factorized_plan_uses_current_charge_without_refresh')
  target = [0.37_dp, 0.41_dp, 0.83_dp]
  call plan%eval(updated_charge, target, potential, field, status)
  call assert_equal_i32(status, upper_vacuum_eval_ok, 'updated-charge target was rejected')
  call eval_unfactorized(mesh, updated_charge, target, length_x, length_y, mode_layers, direct_potential, direct_field)
  call assert_close_dp( &
    potential, direct_potential, 2.0e-13_dp*max(1.0_dp, abs(direct_potential)), &
    'updated charge potential mismatch' &
    )
  call assert_allclose_1d( &
    field, direct_field, 2.0e-12_dp*max(1.0_dp, maxval(abs(direct_field))), &
    'updated charge field mismatch' &
    )
  call test_end()

  call test_begin('periodic_faces_share_one_value')
  target = [0.0_dp, 0.29_dp, 0.71_dp]
  call plan%eval(mesh%q_elem, target, potential, field, status)
  call assert_equal_i32(status, upper_vacuum_eval_ok, 'x-face upper-vacuum target was rejected')
  call plan%eval(mesh%q_elem, target + [length_x, 0.0_dp, 0.0_dp], seam_potential, seam_field, status)
  call assert_equal_i32(status, upper_vacuum_eval_ok, 'translated upper-vacuum target was rejected')
  call assert_close_dp(seam_potential, potential, 2.0e-13_dp, 'x-periodic seam potential mismatch')
  call assert_allclose_1d(seam_field, field, 2.0e-12_dp, 'x-periodic seam field mismatch')
  target = [0.43_dp, 0.0_dp, 0.71_dp]
  call plan%eval(mesh%q_elem, target, potential, field, status)
  call assert_equal_i32(status, upper_vacuum_eval_ok, 'y-face upper-vacuum target was rejected')
  call plan%eval(mesh%q_elem, target + [0.0_dp, length_y, 0.0_dp], seam_potential, seam_field, status)
  call assert_equal_i32(status, upper_vacuum_eval_ok, 'translated y-face upper-vacuum target was rejected')
  call assert_close_dp(seam_potential, potential, 2.0e-13_dp, 'y-periodic seam potential mismatch')
  call assert_allclose_1d(seam_field, field, 2.0e-12_dp, 'y-periodic seam field mismatch')
  call test_end()

  call test_begin('target_at_or_below_source_top_is_rejected')
  target = [0.37_dp, 0.41_dp, source_top]
  call plan%eval(mesh%q_elem, target, potential, field, status)
  call assert_equal_i32(status, upper_vacuum_target_not_above_source, 'z=H target must be rejected')
  call assert_close_dp(potential, 0.0_dp, 0.0_dp, 'rejected target potential must be zero')
  call assert_allclose_1d(field, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 'rejected target field must be zero')
  target(3) = source_top - 0.01_dp
  call plan%eval(mesh%q_elem, target, potential, field, status)
  call assert_equal_i32(status, upper_vacuum_target_not_above_source, 'z<H target must be rejected')
  call test_end()

  call test_summary()

contains

  subroutine eval_unfactorized(mesh, source_charge, target, lx, ly, layers, potential, field)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: source_charge(:), target(3), lx, ly
    integer(i32), intent(in) :: layers
    real(dp), intent(out) :: potential, field(3)
    real(dp) :: area_xy, coefficient, decay, dz, kx, ky, phase, sample_charge, wave_number
    integer(i32) :: element, nx, ny, point

    area_xy = lx*ly
    potential = 0.0_dp
    field = 0.0_dp
    do element = 1, mesh%nelem
      do point = 1, 7
        sample_charge = source_charge(element)*mesh%panel_quad_weight(point, element)/mesh%panel_area(element)
        do nx = -layers, layers
          kx = 2.0_dp*pi*real(nx, dp)/lx
          do ny = -layers, layers
            if (nx == 0_i32 .and. ny == 0_i32) cycle
            ky = 2.0_dp*pi*real(ny, dp)/ly
            wave_number = sqrt(kx*kx + ky*ky)
            phase = kx*(target(1) - mesh%panel_quad_position(1, point, element)) + &
                    ky*(target(2) - mesh%panel_quad_position(2, point, element))
            dz = target(3) - mesh%panel_quad_position(3, point, element)
            decay = exp(-wave_number*abs(dz))
            coefficient = sample_charge*decay/(2.0_dp*eps0*area_xy)
            potential = potential + coefficient*cos(phase)/wave_number
            field(1) = field(1) + coefficient*kx*sin(phase)/wave_number
            field(2) = field(2) + coefficient*ky*sin(phase)/wave_number
            field(3) = field(3) + coefficient*cos(phase)
          end do
        end do
      end do
    end do
  end subroutine eval_unfactorized

end program test_periodic_nonzero_upper_vacuum
