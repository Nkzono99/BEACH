module bem_coulomb_fmm_periodic_nonzero_reference
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi
  use bem_types, only: mesh_type
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_quadrature, only: panel_quadrature_plan_type, build_panel_duffy_quadrature
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: eval_periodic_nonzero_panel_reference

contains

  subroutine eval_periodic_nonzero_panel_reference( &
    mesh, target, length_x, length_y, mode_layers, quadrature_order, potential, field &
    )
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: target(3), length_x, length_y
    integer(i32), intent(in) :: mode_layers, quadrature_order
    real(dp), intent(out) :: potential, field(3)
    type(panel_geometry_type) :: geometry
    type(panel_quadrature_plan_type) :: quadrature
    real(dp) :: area_xy, kx, ky, wave_number, phase, dz, decay, coefficient, sample_charge
    integer(i32) :: element, nx, ny, point, geometry_status

    if (.not. all(ieee_is_finite(target)) .or. .not. ieee_is_finite(length_x) .or. &
        .not. ieee_is_finite(length_y) .or. length_x <= 0.0_dp .or. length_y <= 0.0_dp) then
      error stop 'periodic nonzero reference requires finite target and positive periods.'
    end if
    if (mode_layers < 1_i32 .or. quadrature_order < 2_i32) then
      error stop 'periodic nonzero reference requires positive modes and quadrature order >= 2.'
    end if
    area_xy = length_x*length_y
    potential = 0.0_dp
    field = 0.0_dp
    do element = 1, mesh%nelem
      if (mesh%q_elem(element) == 0.0_dp) cycle
      call init_panel_geometry( &
        mesh%v0(:, element), mesh%v1(:, element), mesh%v2(:, element), geometry, geometry_status &
        )
      if (geometry_status /= panel_geometry_ok) then
        error stop 'periodic nonzero panel reference requires non-degenerate triangles.'
      end if
      call build_panel_duffy_quadrature(geometry, quadrature_order, quadrature)
      do point = 1, quadrature%npoint
        sample_charge = mesh%q_elem(element)*quadrature%weight(point)/geometry%area
        do nx = -mode_layers, mode_layers
          kx = 2.0_dp*pi*real(nx, dp)/length_x
          do ny = -mode_layers, mode_layers
            if (nx == 0_i32 .and. ny == 0_i32) cycle
            ky = 2.0_dp*pi*real(ny, dp)/length_y
            wave_number = sqrt(kx*kx + ky*ky)
            phase = kx*(target(1) - quadrature%position(1, point)) + &
                    ky*(target(2) - quadrature%position(2, point))
            dz = target(3) - quadrature%position(3, point)
            decay = exp(-wave_number*abs(dz))
            coefficient = sample_charge*decay/(2.0_dp*eps0*area_xy)
            potential = potential + coefficient*cos(phase)/wave_number
            field(1) = field(1) + coefficient*kx*sin(phase)/wave_number
            field(2) = field(2) + coefficient*ky*sin(phase)/wave_number
            if (dz > 0.0_dp) then
              field(3) = field(3) + coefficient*cos(phase)
            else if (dz < 0.0_dp) then
              field(3) = field(3) - coefficient*cos(phase)
            end if
          end do
        end do
      end do
    end do
  end subroutine eval_periodic_nonzero_panel_reference

end module bem_coulomb_fmm_periodic_nonzero_reference
