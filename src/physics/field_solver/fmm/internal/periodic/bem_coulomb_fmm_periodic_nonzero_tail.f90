module bem_coulomb_fmm_periodic_nonzero_tail
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi
  use bem_types, only: mesh_type
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_quadrature, only: panel_quadrature_plan_type, build_panel_duffy_quadrature
  implicit none
  private

  integer(i32), parameter, public :: periodic_nonzero_tail_ok = 0_i32
  integer(i32), parameter, public :: periodic_nonzero_tail_invalid = 1_i32

  type, public :: periodic_nonzero_tail_mode_type
    real(dp) :: kx = 0.0_dp
    real(dp) :: ky = 0.0_dp
    real(dp) :: wave_number = 0.0_dp
    real(dp) :: kappa = 0.0_dp
    real(dp) :: alpha = 0.0_dp
    real(dp) :: handoff_z = 0.0_dp
    real(dp) :: incident_cos = 0.0_dp
    real(dp) :: incident_sin = 0.0_dp
    real(dp) :: reflection = 0.0_dp
    real(dp) :: transmission = 1.0_dp
  end type periodic_nonzero_tail_mode_type

  type, public :: periodic_nonzero_tail_plan_type
    integer(i32) :: nmodes = 0_i32
    integer(i32) :: mode_layers = 0_i32
    real(dp) :: length_x = 0.0_dp
    real(dp) :: length_y = 0.0_dp
    real(dp) :: handoff_z = 0.0_dp
    real(dp) :: kappa = 0.0_dp
    real(dp) :: max_linearity = 0.0_dp
    real(dp) :: charge_abs = 0.0_dp
    real(dp) :: temperature_j = 0.0_dp
    type(periodic_nonzero_tail_mode_type), allocatable :: mode(:)
    real(dp), allocatable :: response_cos(:, :)
    real(dp), allocatable :: response_sin(:, :)
  end type periodic_nonzero_tail_plan_type

  public :: init_periodic_nonzero_tail_mode
  public :: eval_periodic_nonzero_tail_correction
  public :: periodic_nonzero_tail_linearity
  public :: build_periodic_nonzero_tail_plan
  public :: eval_periodic_nonzero_tail_plan
  public :: refresh_periodic_nonzero_tail_plan

contains

  subroutine build_periodic_nonzero_tail_plan( &
    mesh, length_x, length_y, handoff_z, kappa, mode_layers, quadrature_order, &
    charge_abs, temperature_j, plan, status &
    )
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: length_x, length_y, handoff_z, kappa, charge_abs, temperature_j
    integer(i32), intent(in) :: mode_layers, quadrature_order
    type(periodic_nonzero_tail_plan_type), intent(out) :: plan
    integer(i32), intent(out) :: status
    type(panel_geometry_type) :: geometry
    type(panel_quadrature_plan_type) :: quadrature
    real(dp) :: kx, ky, wave_number, area_xy, coefficient, phase, sample_charge
    integer(i32) :: nx, ny, index, element, point, geometry_status, mode_status

    status = periodic_nonzero_tail_invalid
    if (mesh%nelem < 1_i32 .or. length_x <= 0.0_dp .or. length_y <= 0.0_dp .or. &
        kappa < 0.0_dp .or. mode_layers < 1_i32 .or. quadrature_order < 2_i32 .or. &
        charge_abs < 0.0_dp .or. temperature_j <= 0.0_dp .or. &
        .not. all(ieee_is_finite([length_x, length_y, handoff_z, kappa, charge_abs, temperature_j]))) return
    if (handoff_z <= max(maxval(mesh%v0(3, :)), max(maxval(mesh%v1(3, :)), maxval(mesh%v2(3, :))))) return
    plan%nmodes = 2_i32*mode_layers*(mode_layers + 1_i32)
    plan%mode_layers = mode_layers
    plan%length_x = length_x
    plan%length_y = length_y
    plan%handoff_z = handoff_z
    plan%kappa = kappa
    plan%charge_abs = charge_abs
    plan%temperature_j = temperature_j
    allocate (plan%mode(plan%nmodes), plan%response_cos(plan%nmodes, mesh%nelem), &
              plan%response_sin(plan%nmodes, mesh%nelem))
    plan%response_cos = 0.0_dp
    plan%response_sin = 0.0_dp
    index = 0_i32
    do nx = 0_i32, mode_layers
      do ny = -mode_layers, mode_layers
        if (nx == 0_i32 .and. ny <= 0_i32) cycle
        index = index + 1_i32
        kx = 2.0_dp*pi*real(nx, dp)/length_x
        ky = 2.0_dp*pi*real(ny, dp)/length_y
        call init_periodic_nonzero_tail_mode( &
          kx, ky, kappa, handoff_z, 0.0_dp, 0.0_dp, plan%mode(index), mode_status &
          )
        if (mode_status /= periodic_nonzero_tail_ok) return
      end do
    end do
    area_xy = length_x*length_y
    do element = 1_i32, mesh%nelem
      call init_panel_geometry( &
        mesh%v0(:, element), mesh%v1(:, element), mesh%v2(:, element), geometry, geometry_status &
        )
      if (geometry_status /= panel_geometry_ok) return
      call build_panel_duffy_quadrature(geometry, quadrature_order, quadrature)
      do point = 1_i32, quadrature%npoint
        sample_charge = quadrature%weight(point)/geometry%area
        do index = 1_i32, plan%nmodes
          wave_number = plan%mode(index)%wave_number
          coefficient = sample_charge*exp( &
                        -wave_number*(handoff_z - quadrature%position(3, point)) &
                        )/(eps0*area_xy*wave_number)
          phase = plan%mode(index)%kx*quadrature%position(1, point) + &
                  plan%mode(index)%ky*quadrature%position(2, point)
          plan%response_cos(index, element) = plan%response_cos(index, element) + coefficient*cos(phase)
          plan%response_sin(index, element) = plan%response_sin(index, element) + coefficient*sin(phase)
        end do
      end do
    end do
    call refresh_periodic_nonzero_tail_plan(mesh%q_elem, plan, mode_status)
    if (mode_status /= periodic_nonzero_tail_ok) return
    status = periodic_nonzero_tail_ok
  end subroutine build_periodic_nonzero_tail_plan

  subroutine refresh_periodic_nonzero_tail_plan(charges, plan, status)
    real(dp), intent(in) :: charges(:)
    type(periodic_nonzero_tail_plan_type), intent(inout) :: plan
    integer(i32), intent(out) :: status
    integer(i32) :: index

    status = periodic_nonzero_tail_invalid
    if (plan%nmodes < 1_i32 .or. .not. allocated(plan%mode) .or. &
        .not. allocated(plan%response_cos) .or. .not. allocated(plan%response_sin) .or. &
        size(plan%response_cos, 2) /= size(charges) .or. .not. all(ieee_is_finite(charges))) return
    plan%max_linearity = 0.0_dp
    do index = 1_i32, plan%nmodes
      plan%mode(index)%incident_cos = dot_product(plan%response_cos(index, :), charges)
      plan%mode(index)%incident_sin = dot_product(plan%response_sin(index, :), charges)
      plan%max_linearity = max( &
                           plan%max_linearity, &
                           periodic_nonzero_tail_linearity(plan%mode(index), plan%charge_abs, plan%temperature_j) &
                           )
    end do
    status = periodic_nonzero_tail_ok
  end subroutine refresh_periodic_nonzero_tail_plan

  subroutine eval_periodic_nonzero_tail_plan(plan, position, potential, field)
    type(periodic_nonzero_tail_plan_type), intent(in) :: plan
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: potential, field(3)
    real(dp) :: mode_potential, mode_field(3)
    integer(i32) :: index

    if (plan%nmodes < 1_i32 .or. .not. allocated(plan%mode)) then
      error stop 'Periodic nonzero-tail plan is not initialized.'
    end if
    potential = 0.0_dp
    field = 0.0_dp
    do index = 1_i32, plan%nmodes
      call eval_periodic_nonzero_tail_correction(plan%mode(index), position, mode_potential, mode_field)
      potential = potential + mode_potential
      field = field + mode_field
    end do
  end subroutine eval_periodic_nonzero_tail_plan

  subroutine init_periodic_nonzero_tail_mode( &
    kx, ky, kappa, handoff_z, incident_cos, incident_sin, mode, status &
    )
    real(dp), intent(in) :: kx, ky, kappa, handoff_z, incident_cos, incident_sin
    type(periodic_nonzero_tail_mode_type), intent(out) :: mode
    integer(i32), intent(out) :: status

    mode = periodic_nonzero_tail_mode_type()
    status = periodic_nonzero_tail_invalid
    if (.not. all(ieee_is_finite([kx, ky, kappa, handoff_z, incident_cos, incident_sin])) .or. &
        kappa < 0.0_dp) return
    mode%wave_number = sqrt(kx*kx + ky*ky)
    if (mode%wave_number <= 0.0_dp) return
    mode%kx = kx
    mode%ky = ky
    mode%kappa = kappa
    mode%alpha = sqrt(mode%wave_number*mode%wave_number + kappa*kappa)
    mode%handoff_z = handoff_z
    mode%incident_cos = incident_cos
    mode%incident_sin = incident_sin
    mode%reflection = (mode%wave_number - mode%alpha)/(mode%wave_number + mode%alpha)
    mode%transmission = 2.0_dp*mode%wave_number/(mode%wave_number + mode%alpha)
    status = periodic_nonzero_tail_ok
  end subroutine init_periodic_nonzero_tail_mode

  subroutine eval_periodic_nonzero_tail_correction(mode, position, potential, field)
    type(periodic_nonzero_tail_mode_type), intent(in) :: mode
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: potential, field(3)
    real(dp) :: phase, basis, tangent_x, tangent_y, dz, factor, normal_factor

    if (mode%wave_number <= 0.0_dp .or. .not. all(ieee_is_finite(position))) then
      error stop 'Invalid periodic nonzero-tail evaluation.'
    end if
    phase = mode%kx*position(1) + mode%ky*position(2)
    basis = mode%incident_cos*cos(phase) + mode%incident_sin*sin(phase)
    tangent_x = mode%kx*(mode%incident_cos*sin(phase) - mode%incident_sin*cos(phase))
    tangent_y = mode%ky*(mode%incident_cos*sin(phase) - mode%incident_sin*cos(phase))
    dz = position(3) - mode%handoff_z
    if (dz <= 0.0_dp) then
      factor = mode%reflection*exp(mode%wave_number*dz)
      normal_factor = -mode%wave_number*factor
    else
      factor = mode%transmission*exp(-mode%alpha*dz) - exp(-mode%wave_number*dz)
      normal_factor = mode%alpha*mode%transmission*exp(-mode%alpha*dz) - &
                      mode%wave_number*exp(-mode%wave_number*dz)
    end if
    potential = factor*basis
    field = [factor*tangent_x, factor*tangent_y, normal_factor*basis]
  end subroutine eval_periodic_nonzero_tail_correction

  pure real(dp) function periodic_nonzero_tail_linearity(mode, charge_abs, temperature_j) result(eta)
    type(periodic_nonzero_tail_mode_type), intent(in) :: mode
    real(dp), intent(in) :: charge_abs, temperature_j

    if (charge_abs < 0.0_dp .or. temperature_j <= 0.0_dp) then
      eta = huge(1.0_dp)
    else
      eta = charge_abs*mode%transmission*sqrt( &
            mode%incident_cos*mode%incident_cos + mode%incident_sin*mode%incident_sin &
            )/temperature_j
    end if
  end function periodic_nonzero_tail_linearity

end module bem_coulomb_fmm_periodic_nonzero_tail
