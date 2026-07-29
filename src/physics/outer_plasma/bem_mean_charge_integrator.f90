!> 平面平均電荷を電流クロージャから陰的に時間発展させる数値コア。
module bem_mean_charge_integrator
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  implicit none
  private

  integer(i32), parameter, public :: mean_charge_step_ok = 0_i32
  integer(i32), parameter, public :: mean_charge_step_invalid = 1_i32
  integer(i32), parameter, public :: mean_charge_step_no_bracket = 2_i32
  integer(i32), parameter, public :: mean_charge_step_max_iterations = 3_i32
  integer(i32), parameter, public :: mean_charge_step_nonfinite = 4_i32

  !> 1 ステップのスカラー根探索設定。
  type, public :: mean_charge_step_options_type
    integer(i32) :: max_iterations = 128_i32
    integer(i32) :: max_bracket_expansions = 64_i32
    real(dp) :: absolute_potential_tolerance = 1.0e-12_dp
    real(dp) :: relative_potential_tolerance = 1.0e-12_dp
  end type mean_charge_step_options_type

  !> 平均電荷の 1 ステップ結果と数値診断。
  type, public :: mean_charge_step_result_type
    integer(i32) :: status = mean_charge_step_invalid
    integer(i32) :: iterations = 0_i32
    integer(i32) :: bracket_expansions = 0_i32
    real(dp) :: potential_v = 0.0_dp
    real(dp) :: surface_charge_density_c_m2 = 0.0_dp
    real(dp) :: current_density_a_m2 = 0.0_dp
    real(dp) :: backward_euler_residual_v = 0.0_dp
    real(dp) :: explicit_predictor_potential_v = 0.0_dp
  end type mean_charge_step_result_type

  abstract interface
    !> 平均界面電位に対する符号付き表面電流密度。
    pure function mean_current_density_model(potential_v, parameters) result(current_density_a_m2)
      import :: dp
      real(dp), intent(in) :: potential_v
      real(dp), intent(in) :: parameters(:)
      real(dp) :: current_density_a_m2
    end function mean_current_density_model
  end interface

  public :: mean_current_density_model
  public :: planar_debye_capacitance_per_area
  public :: advance_mean_charge_backward_euler

contains

  !> 線形 Debye tail の単位面積あたり容量 eps0 / lambda_D を返す。
  pure function planar_debye_capacitance_per_area(debye_length_m) result(capacitance_per_area_f_m2)
    real(dp), intent(in) :: debye_length_m
    real(dp) :: capacitance_per_area_f_m2

    if (.not. ieee_is_finite(debye_length_m) .or. debye_length_m <= 0.0_dp) then
      capacitance_per_area_f_m2 = 0.0_dp
    else
      capacitance_per_area_f_m2 = eps0/debye_length_m
    end if
  end function planar_debye_capacitance_per_area

  !> C_A d(phi)/dt = J(phi) を backward Euler で 1 ステップ進める。
  !!
  !! 平均表面電荷密度と電位は sigma = C_A phi で対応させる。電流モデルは
  !! 符号付き J(phi) だけを提供し、この routine は粒子・シース実装へ依存しない。
  pure subroutine advance_mean_charge_backward_euler( &
    surface_charge_density_old_c_m2, time_step_s, capacitance_per_area_f_m2, &
    current_model, model_parameters, result, options &
    )
    real(dp), intent(in) :: surface_charge_density_old_c_m2
    real(dp), intent(in) :: time_step_s
    real(dp), intent(in) :: capacitance_per_area_f_m2
    procedure(mean_current_density_model) :: current_model
    real(dp), intent(in) :: model_parameters(:)
    type(mean_charge_step_result_type), intent(out) :: result
    type(mean_charge_step_options_type), intent(in), optional :: options

    type(mean_charge_step_options_type) :: controls
    real(dp) :: potential_old, current_old, alpha, explicit_increment
    real(dp) :: lower, upper, residual_lower, residual_upper
    real(dp) :: midpoint, residual_midpoint, bracket_step, direction
    real(dp) :: potential_tolerance
    integer(i32) :: iteration, expansion

    result = mean_charge_step_result_type()
    controls = mean_charge_step_options_type()
    if (present(options)) controls = options

    if (.not. valid_step_inputs( &
        surface_charge_density_old_c_m2, time_step_s, capacitance_per_area_f_m2, controls &
        )) then
      result%status = mean_charge_step_invalid
      return
    end if

    potential_old = surface_charge_density_old_c_m2/capacitance_per_area_f_m2
    if (.not. ieee_is_finite(potential_old)) then
      result%status = mean_charge_step_nonfinite
      return
    end if
    current_old = current_model(potential_old, model_parameters)
    if (.not. ieee_is_finite(current_old)) then
      result%status = mean_charge_step_nonfinite
      return
    end if

    result%potential_v = potential_old
    result%surface_charge_density_c_m2 = surface_charge_density_old_c_m2
    result%current_density_a_m2 = current_old
    result%explicit_predictor_potential_v = potential_old
    if (time_step_s == 0.0_dp .or. current_old == 0.0_dp) then
      result%status = mean_charge_step_ok
      return
    end if

    alpha = time_step_s/capacitance_per_area_f_m2
    explicit_increment = alpha*current_old
    if (.not. ieee_is_finite(alpha) .or. .not. ieee_is_finite(explicit_increment) .or. &
        .not. ieee_is_finite(potential_old + explicit_increment)) then
      result%status = mean_charge_step_nonfinite
      return
    end if
    result%explicit_predictor_potential_v = potential_old + explicit_increment

    potential_tolerance = step_potential_tolerance( &
                          controls, potential_old, result%explicit_predictor_potential_v &
                          )
    if (abs(explicit_increment) <= potential_tolerance .or. &
        result%explicit_predictor_potential_v == potential_old) then
      result%backward_euler_residual_v = -explicit_increment
      result%status = mean_charge_step_ok
      return
    end if

    lower = potential_old
    residual_lower = -explicit_increment
    upper = result%explicit_predictor_potential_v
    residual_upper = backward_euler_residual( &
                     upper, potential_old, alpha, current_model, model_parameters &
                     )
    if (.not. ieee_is_finite(residual_upper)) then
      result%status = mean_charge_step_nonfinite
      return
    end if
    if (residual_upper == 0.0_dp) then
      call populate_result( &
        result, upper, capacitance_per_area_f_m2, current_model, model_parameters, &
        residual_upper, mean_charge_step_ok &
        )
      return
    end if

    if (.not. brackets_root(residual_lower, residual_upper)) then
      direction = sign(1.0_dp, explicit_increment)
      bracket_step = abs(explicit_increment)
      do expansion = 1_i32, controls%max_bracket_expansions
        if (bracket_step > 0.5_dp*huge(1.0_dp)) exit
        bracket_step = 2.0_dp*bracket_step
        upper = potential_old + direction*bracket_step
        if (.not. ieee_is_finite(upper)) exit
        residual_upper = backward_euler_residual( &
                         upper, potential_old, alpha, current_model, model_parameters &
                         )
        if (.not. ieee_is_finite(residual_upper)) then
          result%status = mean_charge_step_nonfinite
          return
        end if
        result%bracket_expansions = expansion
        if (brackets_root(residual_lower, residual_upper)) exit
      end do
      if (.not. brackets_root(residual_lower, residual_upper)) then
        result%status = mean_charge_step_no_bracket
        return
      end if
    end if

    if (lower > upper) then
      call swap_reals(lower, upper)
      call swap_reals(residual_lower, residual_upper)
    end if
    do iteration = 1_i32, controls%max_iterations
      result%iterations = iteration
      midpoint = lower + 0.5_dp*(upper - lower)
      residual_midpoint = backward_euler_residual( &
                          midpoint, potential_old, alpha, current_model, model_parameters &
                          )
      if (.not. ieee_is_finite(residual_midpoint)) then
        result%status = mean_charge_step_nonfinite
        return
      end if
      potential_tolerance = step_potential_tolerance(controls, potential_old, midpoint)
      if (abs(residual_midpoint) <= potential_tolerance .or. &
          abs(upper - lower) <= 2.0_dp*potential_tolerance) then
        call populate_result( &
          result, midpoint, capacitance_per_area_f_m2, current_model, model_parameters, &
          residual_midpoint, mean_charge_step_ok &
          )
        return
      end if
      if (brackets_root(residual_lower, residual_midpoint)) then
        upper = midpoint
      else
        lower = midpoint
        residual_lower = residual_midpoint
      end if
    end do

    call populate_result( &
      result, midpoint, capacitance_per_area_f_m2, current_model, model_parameters, &
      residual_midpoint, mean_charge_step_max_iterations &
      )
  end subroutine advance_mean_charge_backward_euler

  pure function valid_step_inputs( &
    surface_charge_density_old_c_m2, time_step_s, capacitance_per_area_f_m2, controls &
    ) result(valid)
    real(dp), intent(in) :: surface_charge_density_old_c_m2
    real(dp), intent(in) :: time_step_s
    real(dp), intent(in) :: capacitance_per_area_f_m2
    type(mean_charge_step_options_type), intent(in) :: controls
    logical :: valid

    valid = ieee_is_finite(surface_charge_density_old_c_m2) .and. &
            ieee_is_finite(time_step_s) .and. time_step_s >= 0.0_dp .and. &
            ieee_is_finite(capacitance_per_area_f_m2) .and. capacitance_per_area_f_m2 > 0.0_dp .and. &
            controls%max_iterations > 0_i32 .and. controls%max_bracket_expansions >= 0_i32 .and. &
            ieee_is_finite(controls%absolute_potential_tolerance) .and. &
            controls%absolute_potential_tolerance > 0.0_dp .and. &
            ieee_is_finite(controls%relative_potential_tolerance) .and. &
            controls%relative_potential_tolerance >= 0.0_dp
  end function valid_step_inputs

  pure function backward_euler_residual( &
    potential, potential_old, alpha, current_model, model_parameters &
    ) result(residual)
    real(dp), intent(in) :: potential, potential_old, alpha
    procedure(mean_current_density_model) :: current_model
    real(dp), intent(in) :: model_parameters(:)
    real(dp) :: residual

    residual = potential - potential_old - alpha*current_model(potential, model_parameters)
  end function backward_euler_residual

  pure function step_potential_tolerance(controls, potential_a, potential_b) result(tolerance)
    type(mean_charge_step_options_type), intent(in) :: controls
    real(dp), intent(in) :: potential_a, potential_b
    real(dp) :: tolerance

    tolerance = controls%absolute_potential_tolerance + &
                controls%relative_potential_tolerance*max(1.0_dp, abs(potential_a), abs(potential_b))
  end function step_potential_tolerance

  pure function brackets_root(residual_a, residual_b) result(brackets)
    real(dp), intent(in) :: residual_a, residual_b
    logical :: brackets

    brackets = residual_a == 0.0_dp .or. residual_b == 0.0_dp .or. &
               (residual_a < 0.0_dp .and. residual_b > 0.0_dp) .or. &
               (residual_a > 0.0_dp .and. residual_b < 0.0_dp)
  end function brackets_root

  pure subroutine populate_result( &
    result, potential, capacitance_per_area_f_m2, current_model, model_parameters, residual, status &
    )
    type(mean_charge_step_result_type), intent(inout) :: result
    real(dp), intent(in) :: potential, capacitance_per_area_f_m2, model_parameters(:), residual
    procedure(mean_current_density_model) :: current_model
    integer(i32), intent(in) :: status

    result%status = status
    result%potential_v = potential
    result%surface_charge_density_c_m2 = capacitance_per_area_f_m2*potential
    result%current_density_a_m2 = current_model(potential, model_parameters)
    result%backward_euler_residual_v = residual
  end subroutine populate_result

  pure subroutine swap_reals(a, b)
    real(dp), intent(inout) :: a, b
    real(dp) :: temporary

    temporary = a
    a = b
    b = temporary
  end subroutine swap_reals

end module bem_mean_charge_integrator
