!> 上部真空領域で周期 k/=0 panel Fourier 和をsource/target因子へ分離する。
module bem_coulomb_fmm_periodic_nonzero_upper_vacuum
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi
  use bem_types, only: mesh_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer(i32), parameter, public :: upper_vacuum_eval_ok = 0_i32
  integer(i32), parameter, public :: upper_vacuum_target_not_above_source = 1_i32

  type, public :: periodic_nonzero_upper_vacuum_plan_type
    integer(i32) :: mode_layers = 0_i32
    integer(i32) :: nmode = 0_i32
    integer(i32) :: nelem = 0_i32
    real(dp) :: length_x = 0.0_dp
    real(dp) :: length_y = 0.0_dp
    real(dp) :: source_z_max = 0.0_dp
    real(dp), allocatable :: kx(:)
    real(dp), allocatable :: ky(:)
    real(dp), allocatable :: wave_number(:)
    real(dp), allocatable :: panel_moment_real(:, :)
    real(dp), allocatable :: panel_moment_imag(:, :)
  contains
    procedure :: build => build_periodic_nonzero_upper_vacuum_plan
    procedure :: eval => eval_periodic_nonzero_upper_vacuum
  end type periodic_nonzero_upper_vacuum_plan_type

contains

  subroutine build_periodic_nonzero_upper_vacuum_plan(self, mesh, length_x, length_y, mode_layers)
    class(periodic_nonzero_upper_vacuum_plan_type), intent(out) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: length_x, length_y
    integer(i32), intent(in) :: mode_layers
    real(dp) :: source_phase, source_decay, normalized_weight
    integer(i32) :: element, mode, nx, ny, point

    if (mesh%nelem < 1_i32) error stop 'upper-vacuum Fourier plan requires a non-empty mesh.'
    if (.not. ieee_is_finite(length_x) .or. .not. ieee_is_finite(length_y) .or. &
        length_x <= 0.0_dp .or. length_y <= 0.0_dp) then
      error stop 'upper-vacuum Fourier plan requires finite positive periods.'
    end if
    if (mode_layers < 1_i32) error stop 'upper-vacuum Fourier plan requires positive mode layers.'
    if (any(mesh%panel_area <= 0.0_dp)) then
      error stop 'upper-vacuum Fourier plan requires non-degenerate panels.'
    end if

    self%mode_layers = mode_layers
    self%nmode = (2_i32*mode_layers + 1_i32)**2 - 1_i32
    self%nelem = mesh%nelem
    self%length_x = length_x
    self%length_y = length_y
    self%source_z_max = maxval(mesh%bb_max(3, :))
    allocate (self%kx(self%nmode), self%ky(self%nmode), self%wave_number(self%nmode))
    allocate (self%panel_moment_real(self%nelem, self%nmode))
    allocate (self%panel_moment_imag(self%nelem, self%nmode))
    self%panel_moment_real = 0.0_dp
    self%panel_moment_imag = 0.0_dp

    mode = 0_i32
    do nx = -mode_layers, mode_layers
      do ny = -mode_layers, mode_layers
        if (nx == 0_i32 .and. ny == 0_i32) cycle
        mode = mode + 1_i32
        self%kx(mode) = 2.0_dp*pi*real(nx, dp)/length_x
        self%ky(mode) = 2.0_dp*pi*real(ny, dp)/length_y
        self%wave_number(mode) = sqrt(self%kx(mode)**2 + self%ky(mode)**2)
        do element = 1, self%nelem
          do point = 1, 7
            normalized_weight = mesh%panel_quad_weight(point, element)/mesh%panel_area(element)
            source_phase = self%kx(mode)*mesh%panel_quad_position(1, point, element) + &
                           self%ky(mode)*mesh%panel_quad_position(2, point, element)
            source_decay = exp( &
                           self%wave_number(mode)*(mesh%panel_quad_position(3, point, element) - &
                                                   self%source_z_max) &
                           )
            self%panel_moment_real(element, mode) = self%panel_moment_real(element, mode) + &
                                                    normalized_weight*source_decay*cos(source_phase)
            self%panel_moment_imag(element, mode) = self%panel_moment_imag(element, mode) - &
                                                    normalized_weight*source_decay*sin(source_phase)
          end do
        end do
      end do
    end do
  end subroutine build_periodic_nonzero_upper_vacuum_plan

  subroutine eval_periodic_nonzero_upper_vacuum(self, charge, target, potential, field, status)
    class(periodic_nonzero_upper_vacuum_plan_type), intent(in) :: self
    real(dp), intent(in) :: charge(:)
    real(dp), intent(in) :: target(3)
    real(dp), intent(out) :: potential, field(3)
    integer(i32), intent(out) :: status
    real(dp) :: area_xy, coefficient, cosine_phase, sine_phase, target_decay
    real(dp) :: mode_real, mode_imag, source_real, source_imag, target_phase
    integer(i32) :: mode

    if (.not. all(ieee_is_finite(target))) error stop 'upper-vacuum Fourier target must be finite.'
    if (size(charge) /= self%nelem) error stop 'upper-vacuum Fourier charge size mismatch.'
    if (any(.not. ieee_is_finite(charge))) error stop 'upper-vacuum Fourier charges must be finite.'
    potential = 0.0_dp
    field = 0.0_dp
    if (target(3) <= self%source_z_max) then
      status = upper_vacuum_target_not_above_source
      return
    end if

    area_xy = self%length_x*self%length_y
    coefficient = 1.0_dp/(2.0_dp*eps0*area_xy)
    do mode = 1, self%nmode
      source_real = dot_product(charge, self%panel_moment_real(:, mode))
      source_imag = dot_product(charge, self%panel_moment_imag(:, mode))
      target_phase = self%kx(mode)*target(1) + self%ky(mode)*target(2)
      cosine_phase = cos(target_phase)
      sine_phase = sin(target_phase)
      target_decay = exp(-self%wave_number(mode)*(target(3) - self%source_z_max))
      mode_real = target_decay*( &
                  source_real*cosine_phase - source_imag*sine_phase &
                  )
      mode_imag = target_decay*( &
                  source_real*sine_phase + source_imag*cosine_phase &
                  )
      potential = potential + coefficient*mode_real/self%wave_number(mode)
      field(1) = field(1) + coefficient*self%kx(mode)*mode_imag/self%wave_number(mode)
      field(2) = field(2) + coefficient*self%ky(mode)*mode_imag/self%wave_number(mode)
      field(3) = field(3) + coefficient*mode_real
    end do
    status = upper_vacuum_eval_ok
  end subroutine eval_periodic_nonzero_upper_vacuum

end module bem_coulomb_fmm_periodic_nonzero_upper_vacuum
