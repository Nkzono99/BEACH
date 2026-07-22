!> Charge-driven quasi-steady Zhao outer-plasma closure.
module bem_outer_plasma_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi, qe
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid, &
                                    outer_plasma_no_physical_solution, outer_plasma_numerical_failure
  use bem_sheath_model_core, only: zhao_params_type, evaluate_zhao_density_hat, evaluate_zhao_rho_hat, &
                                   zhao_residuals_type_a, zhao_residuals_type_b, zhao_residuals_type_c, &
                                   swe_free_current_term
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer, parameter :: root_max_iterations = 60
  integer, parameter :: root_max_backtracks = 24
  integer, parameter :: rho_quadrature_panels = 256
  integer, parameter :: column_root_max_iterations = 64
  ! Automatic branch continuation samples every requested field/eta interval
  ! and halves rejected steps.  These bounds make a disconnected path fail
  ! closed instead of turning into an unbounded root search.
  integer, parameter :: column_continuation_max_attempts = 256
  integer, parameter :: column_continuation_nominal_segments = 4
  real(dp), parameter :: root_tolerance = 1.0e-9_dp
  real(dp), parameter :: column_root_relative_tolerance = 1.0e-8_dp
  real(dp), parameter :: column_eta_initial = 1.0e-6_dp
  real(dp), parameter :: column_continuation_min_relative_step = 1.0e-10_dp
  ! Maximum normalized potential or log-density change accepted as one
  ! same-branch continuation step before it is subdivided.
  real(dp), parameter :: branch_same_root_step_limit = 0.25_dp
  real(dp), parameter :: branch_degeneracy_tolerance = 1.0e-4_dp
  real(dp), parameter :: branch_continuity_tolerance = 5.0e-4_dp
  ! eta is an occupancy relative to the stationary photoelectron reference, not
  ! a probability.  A finite 16x ceiling admits order-unity transient overshoot
  ! while preventing an unbounded search from overflowing Zhao exponentials.
  real(dp), parameter :: column_eta_search_max = 16.0_dp
  real(dp), parameter :: profile_phi_end_hat = 1.0e-4_dp
  real(dp), parameter :: zero_field_tolerance_hat = 1.0e-12_dp
  integer, parameter :: atlas_max_coordinates = 4
  integer, parameter :: homotopy_max_coordinates = 4
  integer, parameter :: homotopy_max_residuals = 3
  real(dp), parameter :: atlas_default_residual_tolerance = 1.0e-12_dp
  real(dp), parameter :: homotopy_default_residual_tolerance = 1.0e-10_dp
  real(dp), parameter :: atlas_default_log_density_floor = -27.0_dp
  real(dp), parameter :: atlas_fold_tangent_tolerance = 1.0e-6_dp
  real(dp), parameter :: ab_degeneracy_probe_q = 1.0e-3_dp
  real(dp), parameter :: ab_root_condition_tolerance = 1.0e-8_dp
  real(dp), parameter :: ab_density_zero_tolerance = 1.0e-8_dp
  real(dp), parameter :: ab_limit_condition_tolerance = 1.0e-8_dp
  integer, parameter :: atlas_eval_ok = 0
  integer, parameter :: atlas_eval_physical = 1
  integer, parameter :: atlas_eval_numerical = 2

  integer(i32), parameter, public :: zhao_continuation_reason_none = 0_i32
  integer(i32), parameter, public :: zhao_continuation_reason_converged = 1_i32
  integer(i32), parameter, public :: zhao_continuation_reason_numerical_failure = 2_i32
  integer(i32), parameter, public :: zhao_continuation_reason_guard_rejected = 3_i32
  integer(i32), parameter, public :: zhao_continuation_reason_disconnected_branch = 4_i32
  integer(i32), parameter, public :: zhao_continuation_reason_nonmonotone_column = 5_i32
  integer(i32), parameter, public :: zhao_continuation_reason_physical_endpoint = 6_i32
  integer(i32), parameter, public :: zhao_continuation_reason_target_unreachable = 7_i32
  integer(i32), parameter, public :: zhao_continuation_reason_search_limit = 8_i32
  integer(i32), parameter, public :: zhao_continuation_reason_invalid_request = 9_i32
  integer(i32), parameter, public :: zhao_continuation_reason_target_bracketed = 10_i32

  type, public :: zhao_charge_root_type
    character(len=1) :: branch = ' '
    character(len=16) :: interface_side = 'lower'
    real(dp) :: phi0_v = 0.0_dp
    real(dp) :: phi_m_v = 0.0_dp
    real(dp) :: n_swe_inf_m3 = 0.0_dp
    real(dp) :: photoelectron_population_fraction = 1.0_dp
    real(dp) :: photoelectron_column_per_area = 0.0_dp
    real(dp) :: photoelectron_column_target_per_area = 0.0_dp
    real(dp) :: photoelectron_column_residual_per_area = 0.0_dp
    real(dp) :: interface_field_v_m = 0.0_dp
    real(dp) :: net_current_density_a_m2 = 0.0_dp
    real(dp) :: residual_norm = 0.0_dp
    integer(i32) :: nonlinear_iterations = 0_i32
  end type zhao_charge_root_type

  !> Full-precision telemetry for one Zhao continuation decision.
  !>
  !> This record deliberately does not reinterpret a rejected same-branch jump
  !> as a fold.  A later branch-atlas pass must establish whether the connected
  !> root curve continues before assigning that physical classification.
  type, public :: zhao_continuation_diagnostics_type
    integer(i32) :: reason_code = zhao_continuation_reason_none
    integer(i32) :: underlying_status = outer_plasma_invalid
    integer(i32) :: return_status = outer_plasma_invalid
    integer(i32) :: attempt = 0_i32
    character(len=48) :: reason = 'none'
    character(len=24) :: solver_stage = 'none'
    logical :: candidate_available = .false.
    logical :: saw_numerical_failure = .false.
    character(len=1) :: from_branch = ' '
    character(len=1) :: candidate_branch = ' '
    real(dp) :: target_field_v_m = 0.0_dp
    real(dp) :: target_eta = 0.0_dp
    real(dp) :: target_column_m2 = 0.0_dp
    real(dp) :: attempted_step = 0.0_dp
    real(dp) :: from_field_v_m = 0.0_dp
    real(dp) :: from_eta = 0.0_dp
    real(dp) :: from_phi0_v = 0.0_dp
    real(dp) :: from_phi_m_v = 0.0_dp
    real(dp) :: from_density_m3 = 0.0_dp
    real(dp) :: from_column_m2 = 0.0_dp
    real(dp) :: from_column_residual_m2 = 0.0_dp
    real(dp) :: from_root_residual = 0.0_dp
    real(dp) :: candidate_field_v_m = 0.0_dp
    real(dp) :: candidate_eta = 0.0_dp
    real(dp) :: candidate_phi0_v = 0.0_dp
    real(dp) :: candidate_phi_m_v = 0.0_dp
    real(dp) :: candidate_density_m3 = 0.0_dp
    real(dp) :: candidate_column_m2 = 0.0_dp
    real(dp) :: candidate_column_residual_m2 = 0.0_dp
    real(dp) :: candidate_root_residual = 0.0_dp
    real(dp) :: normalized_potential_jump = 0.0_dp
    real(dp) :: log_density_jump = 0.0_dp
    real(dp) :: normalized_root_jump = 0.0_dp
  end type zhao_continuation_diagnostics_type

  !> Controls for a diagnostic pseudo-arclength trace on one fixed Zhao branch.
  !>
  !> The trace is intentionally not used by the runtime closure.  It establishes
  !> whether a rejected fixed-eta step lies on a connected root curve before a
  !> later change is allowed to alter production root selection.
  type, public :: zhao_branch_atlas_options_type
    integer(i32) :: max_points = 160_i32
    integer(i32) :: max_corrector_iterations = 24_i32
    integer(i32) :: max_step_halvings = 20_i32
    integer(i32) :: eta_direction = 1_i32
    real(dp) :: initial_step = 5.0e-2_dp
    real(dp) :: minimum_step = 1.0e-7_dp
    real(dp) :: maximum_step = 1.0e-1_dp
    real(dp) :: residual_tolerance = atlas_default_residual_tolerance
    real(dp) :: seed_refinement_limit = 5.0e-1_dp
    real(dp) :: eta_min = 0.0_dp
    real(dp) :: eta_max = column_eta_search_max
    real(dp) :: log_density_floor = atlas_default_log_density_floor
  end type zhao_branch_atlas_options_type

  type, public :: zhao_branch_atlas_point_type
    type(zhao_charge_root_type) :: root
    real(dp) :: arc_length = 0.0_dp
    real(dp) :: accepted_step = 0.0_dp
    real(dp) :: tangent(atlas_max_coordinates) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: dcolumn_ds = 0.0_dp
    real(dp) :: row_rank_indicator = 0.0_dp
    real(dp) :: normalized_jump_from_previous = 0.0_dp
    integer(i32) :: corrector_iterations = 0_i32
  end type zhao_branch_atlas_point_type

  type, public :: zhao_branch_atlas_type
    type(zhao_branch_atlas_point_type), allocatable :: points(:)
    character(len=1) :: branch = ' '
    integer(i32) :: point_count = 0_i32
    integer(i32) :: termination_reason_code = zhao_continuation_reason_none
    character(len=48) :: termination_reason = 'none'
    real(dp) :: interface_field_v_m = 0.0_dp
    real(dp) :: target_column_m2 = 0.0_dp
    real(dp) :: maximum_column_m2 = 0.0_dp
    real(dp) :: log_density_floor = atlas_default_log_density_floor
    real(dp) :: seed_refinement_jump = 0.0_dp
    logical :: target_requested = .false.
    logical :: target_bracketed = .false.
    logical :: eta_fold_detected = .false.
    logical :: column_fold_detected = .false.
    logical :: seed_reanchored = .false.
  end type zhao_branch_atlas_type

  !> Controls for a diagnostic straight field-column homotopy on Zhao B/C.
  !>
  !> The homotopy coordinate prescribes both the interface field and the finite
  !> photoelectron column between two accepted time levels.  This tracer is not
  !> consulted by the runtime closure or automatic branch selection.
  type, public :: zhao_field_column_homotopy_options_type
    integer(i32) :: max_points = 512_i32
    integer(i32) :: max_corrector_iterations = 24_i32
    integer(i32) :: max_step_halvings = 20_i32
    real(dp) :: initial_step = 2.5e-2_dp
    real(dp) :: minimum_step = 1.0e-7_dp
    real(dp) :: maximum_step = 5.0e-2_dp
    real(dp) :: residual_tolerance = homotopy_default_residual_tolerance
    real(dp) :: seed_refinement_limit = 5.0e-1_dp
    real(dp) :: eta_min = 0.0_dp
    real(dp) :: eta_max = column_eta_search_max
    real(dp) :: homotopy_min = -2.5e-1_dp
    real(dp) :: homotopy_max = 1.25_dp
    real(dp) :: log_density_floor = atlas_default_log_density_floor
  end type zhao_field_column_homotopy_options_type

  type, public :: zhao_field_column_homotopy_point_type
    type(zhao_charge_root_type) :: root
    real(dp) :: homotopy_fraction = 0.0_dp
    real(dp) :: prescribed_column_m2 = 0.0_dp
    real(dp) :: normalized_column_residual = 0.0_dp
    real(dp) :: arc_length = 0.0_dp
    real(dp) :: accepted_step = 0.0_dp
    real(dp) :: tangent(homotopy_max_coordinates) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: row_rank_indicator = 0.0_dp
    real(dp) :: normalized_jump_from_previous = 0.0_dp
    integer(i32) :: corrector_iterations = 0_i32
  end type zhao_field_column_homotopy_point_type

  type, public :: zhao_field_column_homotopy_type
    type(zhao_field_column_homotopy_point_type), allocatable :: points(:)
    character(len=1) :: branch = ' '
    integer(i32) :: point_count = 0_i32
    integer(i32) :: termination_reason_code = zhao_continuation_reason_none
    character(len=48) :: termination_reason = 'none'
    real(dp) :: start_field_v_m = 0.0_dp
    real(dp) :: target_field_v_m = 0.0_dp
    real(dp) :: start_column_m2 = 0.0_dp
    real(dp) :: target_column_m2 = 0.0_dp
    real(dp) :: column_scale_m2 = 1.0_dp
    real(dp) :: seed_refinement_jump = 0.0_dp
    logical :: target_reached = .false.
    logical :: homotopy_fold_detected = .false.
    logical :: seed_reanchored = .false.
  end type zhao_field_column_homotopy_type

  !> Diagnostic chart for the singular Type-B to Type-A limit.
  !>
  !> q=sqrt(-phi_m/T_pe) regularizes the Type-A minimum coordinate.  A
  !> non-trivial Type-A tangent requires a density-zero field limit and a zero
  !> q^3 coefficient of the far-field Sagdeev residual.  These are necessary,
  !> not sufficient, connection conditions.  This record is diagnostic only
  !> and is not consulted by production branch selection.
  type, public :: zhao_ab_degeneracy_diagnostics_type
    character(len=48) :: classification = 'none'
    real(dp) :: ambient_density_ratio = 0.0_dp
    real(dp) :: photoelectron_quasineutral_term_hat = 0.0_dp
    real(dp) :: doubled_quasineutral_residual_hat = 0.0_dp
    real(dp) :: limiting_field_squared_jump_hat = 0.0_dp
    real(dp) :: b_field_squared_residual_hat = 0.0_dp
    real(dp) :: quasineutral_far_field_q3_coefficient = 0.0_dp
    real(dp) :: probe_q = ab_degeneracy_probe_q
    real(dp) :: probe_ambient_density_ratio = 0.0_dp
    real(dp) :: probe_quasineutral_residual = 0.0_dp
    real(dp) :: probe_far_field_q3_coefficient = 0.0_dp
    real(dp) :: probe_field_squared_residual_hat = 0.0_dp
    logical :: probe_available = .false.
    logical :: density_zero_limit = .false.
    logical :: b_root_field_condition_met = .false.
    logical :: limiting_field_condition_met = .false.
    logical :: far_field_tangent_condition_met = .false.
    logical :: regular_connection_conditions_met = .false.
  end type zhao_ab_degeneracy_diagnostics_type

  public :: solve_zhao_charge_root
  public :: evaluate_zhao_interface_field
  public :: build_zhao_outer_profile
  public :: solve_outer_plasma_zhao
  public :: solve_outer_plasma_zhao_column
  public :: zhao_net_current_density
  public :: write_zhao_continuation_diagnostics
  public :: trace_zhao_branch_atlas
  public :: trace_zhao_field_column_homotopy
  public :: diagnose_zhao_ab_degeneracy

contains

  subroutine diagnose_zhao_ab_degeneracy(params, b_root, diagnostics, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: b_root
    type(zhao_ab_degeneracy_diagnostics_type), intent(out) :: diagnostics
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(zhao_params_type) :: trial_params
    real(dp) :: phi0_hat, density_hat, density_for_integral_hat
    real(dp) :: photo_density_hat, ion_density_hat
    real(dp) :: a_integral, b_integral, a_field_squared, b_field_squared
    real(dp) :: probe_q, probe_phi_m_hat, probe_density_hat, probe_numerator, probe_denominator
    real(dp) :: probe_integral, probe_field_squared, target_field_hat, field_scale
    real(dp) :: x3(3), raw(3), condition_scale
    logical :: a_integral_ok, b_integral_ok, probe_integral_ok

    diagnostics = zhao_ab_degeneracy_diagnostics_type()
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. b_root%branch /= 'B' .or. &
        .not. all(ieee_is_finite([ &
                                 b_root%phi0_v, b_root%n_swe_inf_m3, &
                                 b_root%photoelectron_population_fraction, &
                                 b_root%interface_field_v_m &
                                 ])) .or. b_root%phi0_v <= 0.0_dp .or. &
        b_root%n_swe_inf_m3 < 0.0_dp .or. &
        b_root%photoelectron_population_fraction < 0.0_dp .or. &
        b_root%interface_field_v_m <= 0.0_dp) then
      diagnostics%classification = 'invalid_b_limit'
      message = 'invalid Zhao-B root for A/B degeneracy diagnostics'
      return
    end if

    trial_params = params
    trial_params%photoelectron_population_fraction = &
      b_root%photoelectron_population_fraction
    phi0_hat = b_root%phi0_v/trial_params%t_phe_ev
    density_hat = b_root%n_swe_inf_m3/trial_params%n_phe_ref_m3
    ! integrate_rho_hat rejects an exactly zero normalization even though the
    ! density formulas have a regular zero-density limit.
    density_for_integral_hat = max( &
                               density_hat, epsilon(1.0_dp)*ab_density_zero_tolerance &
                               )
    ion_density_hat = trial_params%n_swi_inf_m3/trial_params%n_phe_ref_m3
    photo_density_hat = &
      trial_params%photoelectron_population_fraction* &
      trial_params%n_phe0_m3/trial_params%n_phe_ref_m3*exp(-phi0_hat)

    diagnostics%ambient_density_ratio = density_hat
    diagnostics%photoelectron_quasineutral_term_hat = photo_density_hat
    diagnostics%doubled_quasineutral_residual_hat = &
      density_hat*(1.0_dp + erf(trial_params%u)) + &
      photo_density_hat - 2.0_dp*ion_density_hat
    ! This is the q^3 coefficient after projection onto the Type-A
    ! quasineutral curve; an off-curve state also has a q^2 term.
    diagnostics%quasineutral_far_field_q3_coefficient = &
      2.0_dp/(3.0_dp*sqrt(pi))*( &
      photo_density_hat - &
      density_hat*exp(-trial_params%u*trial_params%u)/sqrt(trial_params%tau) &
      )

    call integrate_rho_hat( &
      trial_params, 'A', 'lower', 0.0_dp, phi0_hat, phi0_hat, 0.0_dp, &
      density_for_integral_hat, a_integral, a_integral_ok &
      )
    call integrate_rho_hat( &
      trial_params, 'B', 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi0_hat, &
      density_for_integral_hat, b_integral, b_integral_ok &
      )
    if (.not. a_integral_ok .or. .not. b_integral_ok) then
      diagnostics%classification = 'limit_profile_evaluation_failed'
      status = outer_plasma_numerical_failure
      message = 'Zhao A/B limiting profile could not be evaluated'
      return
    end if
    a_field_squared = -2.0_dp*a_integral
    b_field_squared = 2.0_dp*b_integral
    diagnostics%limiting_field_squared_jump_hat = a_field_squared - b_field_squared
    field_scale = trial_params%t_phe_ev/trial_params%lambda_d_phe_ref_m
    target_field_hat = b_root%interface_field_v_m/field_scale
    diagnostics%b_field_squared_residual_hat = &
      b_field_squared - target_field_hat*target_field_hat

    probe_q = diagnostics%probe_q
    probe_phi_m_hat = -(probe_q*probe_q)
    probe_numerator = &
      2.0_dp*ion_density_hat - photo_density_hat*erfc(probe_q)
    probe_denominator = &
      1.0_dp + 2.0_dp*erf(trial_params%u) + &
      erf(probe_q/sqrt(trial_params%tau) - trial_params%u)
    if (ieee_is_finite(probe_numerator) .and. &
        ieee_is_finite(probe_denominator) .and. &
        probe_numerator > 0.0_dp .and. probe_denominator > 0.0_dp) then
      probe_density_hat = probe_numerator/probe_denominator
      diagnostics%probe_ambient_density_ratio = probe_density_hat
      x3 = [ &
           b_root%phi0_v, &
           probe_phi_m_hat*trial_params%t_phe_ev, &
           probe_density_hat*trial_params%n_phe_ref_m3 &
           ]
      call zhao_residuals_type_a(trial_params, x3, raw)
      diagnostics%probe_quasineutral_residual = raw(1)/trial_params%n_phe_ref_m3
      diagnostics%probe_far_field_q3_coefficient = raw(3)/(probe_q**3)

      call integrate_rho_hat( &
        trial_params, 'A', 'lower', probe_phi_m_hat, phi0_hat, &
        phi0_hat, probe_phi_m_hat, probe_density_hat, &
        probe_integral, probe_integral_ok &
        )
      if (probe_integral_ok) then
        probe_field_squared = -2.0_dp*probe_integral
        diagnostics%probe_field_squared_residual_hat = &
          probe_field_squared - b_field_squared
        diagnostics%probe_available = .true.
      end if
    end if

    condition_scale = max(abs(a_field_squared), abs(b_field_squared), 1.0_dp)
    diagnostics%density_zero_limit = &
      density_hat <= ab_density_zero_tolerance
    diagnostics%b_root_field_condition_met = &
      abs(diagnostics%b_field_squared_residual_hat) <= ab_root_condition_tolerance*condition_scale
    diagnostics%limiting_field_condition_met = &
      abs(diagnostics%limiting_field_squared_jump_hat) <= &
      ab_limit_condition_tolerance*condition_scale
    diagnostics%far_field_tangent_condition_met = &
      abs(diagnostics%quasineutral_far_field_q3_coefficient) <= &
      ab_limit_condition_tolerance
    diagnostics%regular_connection_conditions_met = &
      abs(diagnostics%doubled_quasineutral_residual_hat) <= ab_root_condition_tolerance .and. &
      diagnostics%density_zero_limit .and. &
      diagnostics%b_root_field_condition_met .and. &
      diagnostics%limiting_field_condition_met .and. &
      diagnostics%far_field_tangent_condition_met

    if (abs(diagnostics%doubled_quasineutral_residual_hat) > ab_root_condition_tolerance) then
      diagnostics%classification = 'b_limit_not_quasineutral'
    else if (.not. diagnostics%b_root_field_condition_met) then
      diagnostics%classification = 'b_limit_field_mismatch'
    else if (.not. diagnostics%density_zero_limit) then
      diagnostics%classification = 'nonzero_ambient_density_limit'
    else if (.not. diagnostics%limiting_field_condition_met) then
      diagnostics%classification = 'limiting_field_discontinuous'
    else if (.not. diagnostics%far_field_tangent_condition_met) then
      diagnostics%classification = 'no_regular_type_a_tangent'
    else
      diagnostics%classification = 'regular_connection_candidate'
    end if
    status = outer_plasma_ok
  end subroutine diagnose_zhao_ab_degeneracy

  subroutine write_zhao_continuation_diagnostics(unit, diagnostics, call_stage, batch_index)
    integer, intent(in) :: unit
    type(zhao_continuation_diagnostics_type), intent(in) :: diagnostics
    character(len=*), intent(in), optional :: call_stage
    integer(i32), intent(in), optional :: batch_index

    character(len=32) :: resolved_call_stage
    character(len=1) :: resolved_from_branch, resolved_candidate_branch
    integer(i32) :: resolved_batch

    resolved_call_stage = 'unknown'
    if (present(call_stage)) resolved_call_stage = trim(call_stage)
    resolved_batch = -1_i32
    if (present(batch_index)) resolved_batch = batch_index
    resolved_from_branch = diagnostics%from_branch
    if (resolved_from_branch == ' ') resolved_from_branch = '-'
    resolved_candidate_branch = diagnostics%candidate_branch
    if (resolved_candidate_branch == ' ') resolved_candidate_branch = '-'
    write (unit, '(a,a,a,a,a,i0,a,i0,a,l1)') &
      'BEACH zhao-continuation call_stage=', trim(resolved_call_stage), &
      ' solver_stage=', trim(diagnostics%solver_stage), ' batch=', resolved_batch, &
      ' reason_code=', diagnostics%reason_code, &
      ' saw_numerical_failure=', diagnostics%saw_numerical_failure
    write (unit, '(a,a,a,i0,a,i0,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3)') &
      'BEACH zhao-continuation reason=', trim(diagnostics%reason), &
      ' underlying_status=', diagnostics%underlying_status, &
      ' return_status=', diagnostics%return_status, &
      ' target_field_V_m=', diagnostics%target_field_v_m, &
      ' target_eta=', diagnostics%target_eta, &
      ' target_column_m-2=', diagnostics%target_column_m2, &
      ' attempted_step=', diagnostics%attempted_step
    write (unit, '(a,i0,a,a,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3)') &
      'BEACH zhao-continuation attempt=', diagnostics%attempt, &
      ' from_branch=', resolved_from_branch, &
      ' from_field_V_m=', diagnostics%from_field_v_m, &
      ' from_eta=', diagnostics%from_eta, &
      ' from_phi0_V=', diagnostics%from_phi0_v, &
      ' from_phi_m_V=', diagnostics%from_phi_m_v, &
      ' from_density_m-3=', diagnostics%from_density_m3, &
      ' from_column_m-2=', diagnostics%from_column_m2, &
      ' from_column_residual_m-2=', diagnostics%from_column_residual_m2
    write (unit, '(a,l1,a,a,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3)') &
      'BEACH zhao-continuation candidate_available=', diagnostics%candidate_available, &
      ' candidate_branch=', resolved_candidate_branch, &
      ' candidate_field_V_m=', diagnostics%candidate_field_v_m, &
      ' candidate_eta=', diagnostics%candidate_eta, &
      ' candidate_phi0_V=', diagnostics%candidate_phi0_v, &
      ' candidate_phi_m_V=', diagnostics%candidate_phi_m_v, &
      ' candidate_density_m-3=', diagnostics%candidate_density_m3, &
      ' candidate_column_m-2=', diagnostics%candidate_column_m2, &
      ' candidate_column_residual_m-2=', diagnostics%candidate_column_residual_m2
    write (unit, '(a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3)') &
      'BEACH zhao-continuation from_root_residual=', diagnostics%from_root_residual, &
      ' candidate_root_residual=', diagnostics%candidate_root_residual, &
      ' normalized_potential_jump=', diagnostics%normalized_potential_jump, &
      ' log_density_jump=', diagnostics%log_density_jump, &
      ' normalized_root_jump=', diagnostics%normalized_root_jump
  end subroutine write_zhao_continuation_diagnostics

  subroutine set_zhao_continuation_origin( &
    diagnostics, stage, reason_code, reason, underlying_status, from_eta, from_root &
    )
    type(zhao_continuation_diagnostics_type), intent(inout), optional :: diagnostics
    character(len=*), intent(in) :: stage, reason
    integer(i32), intent(in) :: reason_code, underlying_status
    real(dp), intent(in) :: from_eta
    type(zhao_charge_root_type), intent(in) :: from_root

    if (.not. present(diagnostics)) return
    diagnostics%reason_code = reason_code
    diagnostics%underlying_status = underlying_status
    diagnostics%return_status = underlying_status
    diagnostics%reason = trim(reason)
    diagnostics%solver_stage = trim(stage)
    diagnostics%candidate_available = .false.
    diagnostics%attempt = 0_i32
    diagnostics%attempted_step = 0.0_dp
    diagnostics%candidate_branch = ' '
    diagnostics%candidate_field_v_m = 0.0_dp
    diagnostics%candidate_eta = 0.0_dp
    diagnostics%candidate_phi0_v = 0.0_dp
    diagnostics%candidate_phi_m_v = 0.0_dp
    diagnostics%candidate_density_m3 = 0.0_dp
    diagnostics%candidate_column_m2 = 0.0_dp
    diagnostics%candidate_column_residual_m2 = 0.0_dp
    diagnostics%candidate_root_residual = 0.0_dp
    diagnostics%normalized_potential_jump = 0.0_dp
    diagnostics%log_density_jump = 0.0_dp
    diagnostics%normalized_root_jump = 0.0_dp
    diagnostics%from_branch = from_root%branch
    diagnostics%from_field_v_m = from_root%interface_field_v_m
    diagnostics%from_eta = from_eta
    diagnostics%from_phi0_v = from_root%phi0_v
    diagnostics%from_phi_m_v = from_root%phi_m_v
    diagnostics%from_density_m3 = from_root%n_swe_inf_m3
    diagnostics%from_column_m2 = from_root%photoelectron_column_per_area
    diagnostics%from_column_residual_m2 = from_root%photoelectron_column_residual_per_area
    diagnostics%from_root_residual = from_root%residual_norm
  end subroutine set_zhao_continuation_origin

  subroutine set_zhao_continuation_pair( &
    diagnostics, stage, reason_code, reason, underlying_status, params, from_eta, from_root, &
    candidate_eta, candidate_root &
    )
    type(zhao_continuation_diagnostics_type), intent(inout), optional :: diagnostics
    character(len=*), intent(in) :: stage, reason
    integer(i32), intent(in) :: reason_code, underlying_status
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: from_eta, candidate_eta
    type(zhao_charge_root_type), intent(in) :: from_root, candidate_root

    if (.not. present(diagnostics)) return
    call set_zhao_continuation_origin( &
      diagnostics, stage, reason_code, reason, underlying_status, from_eta, from_root &
      )
    diagnostics%candidate_available = .true.
    diagnostics%candidate_branch = candidate_root%branch
    diagnostics%candidate_field_v_m = candidate_root%interface_field_v_m
    diagnostics%candidate_eta = candidate_eta
    diagnostics%candidate_phi0_v = candidate_root%phi0_v
    diagnostics%candidate_phi_m_v = candidate_root%phi_m_v
    diagnostics%candidate_density_m3 = candidate_root%n_swe_inf_m3
    diagnostics%candidate_column_m2 = candidate_root%photoelectron_column_per_area
    diagnostics%candidate_column_residual_m2 = candidate_root%photoelectron_column_residual_per_area
    diagnostics%candidate_root_residual = candidate_root%residual_norm
    call zhao_root_jump_metrics( &
      params, from_root, candidate_root, diagnostics%normalized_potential_jump, &
      diagnostics%log_density_jump, diagnostics%normalized_root_jump &
      )
  end subroutine set_zhao_continuation_pair

  subroutine mark_zhao_continuation_converged(diagnostics, root)
    type(zhao_continuation_diagnostics_type), intent(inout), optional :: diagnostics
    type(zhao_charge_root_type), intent(in) :: root

    if (.not. present(diagnostics)) return
    diagnostics = zhao_continuation_diagnostics_type()
    diagnostics%reason_code = zhao_continuation_reason_converged
    diagnostics%underlying_status = outer_plasma_ok
    diagnostics%return_status = outer_plasma_ok
    diagnostics%reason = 'converged'
    diagnostics%solver_stage = 'complete'
    diagnostics%target_field_v_m = root%interface_field_v_m
    diagnostics%target_eta = root%photoelectron_population_fraction
    diagnostics%target_column_m2 = root%photoelectron_column_target_per_area
    diagnostics%candidate_available = .true.
    diagnostics%candidate_branch = root%branch
    diagnostics%candidate_field_v_m = root%interface_field_v_m
    diagnostics%candidate_eta = root%photoelectron_population_fraction
    diagnostics%candidate_phi0_v = root%phi0_v
    diagnostics%candidate_phi_m_v = root%phi_m_v
    diagnostics%candidate_density_m3 = root%n_swe_inf_m3
    diagnostics%candidate_column_m2 = root%photoelectron_column_per_area
    diagnostics%candidate_column_residual_m2 = root%photoelectron_column_residual_per_area
    diagnostics%candidate_root_residual = root%residual_norm
  end subroutine mark_zhao_continuation_converged

  subroutine zhao_root_jump_metrics(params, from_root, to_root, potential_change, density_change, jump)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: from_root, to_root
    real(dp), intent(out) :: potential_change, density_change, jump

    real(dp) :: potential_scale

    potential_change = huge(1.0_dp)
    density_change = huge(1.0_dp)
    jump = huge(1.0_dp)
    if (.not. zhao_charge_root_is_finite(from_root) .or. &
        .not. zhao_charge_root_is_finite(to_root)) return
    if (from_root%n_swe_inf_m3 <= 0.0_dp .or. to_root%n_swe_inf_m3 <= 0.0_dp) return
    potential_scale = max(params%t_phe_ev, tiny(1.0_dp))
    potential_change = max( &
                       abs(to_root%phi0_v - from_root%phi0_v), &
                       abs(to_root%phi_m_v - from_root%phi_m_v) &
                       )/potential_scale
    density_change = abs(log(to_root%n_swe_inf_m3/from_root%n_swe_inf_m3))
    jump = max(potential_change, density_change)
  end subroutine zhao_root_jump_metrics

  pure logical function zhao_charge_root_is_finite(root) result(finite)
    type(zhao_charge_root_type), intent(in) :: root

    finite = all(ieee_is_finite([ &
                                root%phi0_v, root%phi_m_v, root%n_swe_inf_m3, &
                                root%photoelectron_population_fraction, &
                                root%photoelectron_column_per_area, &
                                root%photoelectron_column_target_per_area, &
                                root%photoelectron_column_residual_per_area, &
                                root%interface_field_v_m, root%net_current_density_a_m2, root%residual_norm &
                                ]))
  end function zhao_charge_root_is_finite

  subroutine trace_zhao_branch_atlas( &
    params, interface_field_v_m, grid_points, control_length_m, seed_root, atlas, status, message, &
    target_column_m2, options &
    )
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m
    type(zhao_charge_root_type), intent(in) :: seed_root
    type(zhao_branch_atlas_type), intent(out) :: atlas
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: target_column_m2
    type(zhao_branch_atlas_options_type), intent(in), optional :: options

    type(zhao_branch_atlas_options_type) :: resolved
    type(zhao_branch_atlas_point_type), allocatable :: compact_points(:)
    type(zhao_branch_atlas_point_type) :: next_point
    type(zhao_charge_root_type) :: candidate_root
    real(dp) :: z(atlas_max_coordinates), seed_coordinates(atlas_max_coordinates)
    real(dp) :: corrected_z(atlas_max_coordinates)
    real(dp) :: predictor(atlas_max_coordinates), tangent(atlas_max_coordinates)
    real(dp) :: next_tangent(atlas_max_coordinates), anchor_tangent(atlas_max_coordinates)
    real(dp) :: residual(3), step, trial_step, row_rank, next_rank
    real(dp) :: target_tolerance, previous_target_residual, next_target_residual
    real(dp) :: column_rate_tolerance
    real(dp) :: potential_jump, density_jump, root_jump
    integer :: n, z_dimension, point_index, halving, iterations, failure_kind
    integer :: corrector_failure_kind, eta_tangent_sign, column_rate_sign, next_sign
    logical :: valid, tangent_ok, corrected, point_ok, candidate_decoded
    logical :: last_eta_lower_limit, last_eta_upper_limit

    atlas%branch = ' '
    atlas%point_count = 0_i32
    atlas%termination_reason_code = zhao_continuation_reason_none
    atlas%termination_reason = 'none'
    atlas%interface_field_v_m = 0.0_dp
    atlas%target_column_m2 = 0.0_dp
    atlas%maximum_column_m2 = 0.0_dp
    atlas%log_density_floor = atlas_default_log_density_floor
    atlas%seed_refinement_jump = 0.0_dp
    atlas%target_requested = .false.
    atlas%target_bracketed = .false.
    atlas%eta_fold_detected = .false.
    atlas%column_fold_detected = .false.
    atlas%seed_reanchored = .false.
    status = outer_plasma_invalid
    message = ''
    resolved = zhao_branch_atlas_options_type()
    if (present(options)) resolved = options
    atlas%branch = seed_root%branch
    atlas%interface_field_v_m = interface_field_v_m
    atlas%log_density_floor = resolved%log_density_floor
    if (present(target_column_m2)) then
      atlas%target_requested = .true.
      atlas%target_column_m2 = target_column_m2
    end if

    if (seed_root%branch == 'A') then
      n = 3
    else if (seed_root%branch == 'B' .or. seed_root%branch == 'C') then
      n = 2
    else
      atlas%termination_reason_code = zhao_continuation_reason_invalid_request
      atlas%termination_reason = 'invalid_seed_branch'
      message = 'Zhao branch atlas requires an A, B, or C seed root'
      return
    end if
    z_dimension = n + 1
    if (.not. valid_zhao_params(params) .or. .not. ieee_is_finite(interface_field_v_m) .or. &
        grid_points < 5_i32 .or. .not. ieee_is_finite(control_length_m) .or. &
        control_length_m <= 0.0_dp .or. resolved%max_points < 1_i32 .or. &
        resolved%max_corrector_iterations < 1_i32 .or. resolved%max_step_halvings < 0_i32 .or. &
        abs(resolved%eta_direction) /= 1_i32 .or. &
        .not. all(ieee_is_finite([ &
                                 resolved%initial_step, resolved%minimum_step, resolved%maximum_step, &
                                 resolved%residual_tolerance, resolved%seed_refinement_limit, &
                                 resolved%eta_min, resolved%eta_max, &
                                 resolved%log_density_floor &
                                 ])) .or. resolved%minimum_step <= 0.0_dp .or. &
        resolved%initial_step < resolved%minimum_step .or. &
        resolved%maximum_step < resolved%initial_step .or. &
        resolved%residual_tolerance <= 0.0_dp .or. resolved%seed_refinement_limit <= 0.0_dp .or. &
        resolved%eta_min < 0.0_dp .or. &
        resolved%eta_max <= resolved%eta_min .or. resolved%eta_max > column_eta_search_max .or. &
        resolved%log_density_floor <= -30.0_dp .or. resolved%log_density_floor >= 0.0_dp) then
      atlas%termination_reason_code = zhao_continuation_reason_invalid_request
      atlas%termination_reason = 'invalid_options'
      message = 'invalid Zhao branch-atlas request'
      return
    end if
    if (atlas%target_requested) then
      if (.not. ieee_is_finite(atlas%target_column_m2) .or. atlas%target_column_m2 < 0.0_dp) then
        atlas%termination_reason_code = zhao_continuation_reason_invalid_request
        atlas%termination_reason = 'invalid_target_column'
        message = 'invalid Zhao branch-atlas target column'
        return
      end if
    end if

    z = 0.0_dp
    call encode_unknowns( &
      params, seed_root%branch, seed_root%phi0_v, seed_root%phi_m_v, &
      seed_root%n_swe_inf_m3, z(1:3), valid &
      )
    if (.not. valid) then
      atlas%termination_reason_code = zhao_continuation_reason_invalid_request
      atlas%termination_reason = 'invalid_seed_encoding'
      message = 'Zhao branch-atlas seed root cannot be encoded'
      return
    end if
    z(z_dimension) = seed_root%photoelectron_population_fraction
    seed_coordinates = z
    if (z(z_dimension) < resolved%eta_min .or. z(z_dimension) > resolved%eta_max) then
      atlas%termination_reason_code = zhao_continuation_reason_invalid_request
      atlas%termination_reason = 'seed_eta_outside_options'
      message = 'Zhao branch-atlas seed eta lies outside its requested interval'
      return
    end if

    ! Refine the supplied runtime root at fixed eta before computing a tangent.
    anchor_tangent = 0.0_dp
    anchor_tangent(z_dimension) = 1.0_dp
    call correct_zhao_atlas_point( &
      params, seed_root%branch, interface_field_v_m, n, z, anchor_tangent, resolved, &
      corrected_z, iterations, corrected, corrector_failure_kind &
      )
    if (.not. corrected) then
      status = merge( &
               outer_plasma_no_physical_solution, outer_plasma_numerical_failure, &
               corrector_failure_kind == atlas_eval_physical &
               )
      atlas%termination_reason_code = merge( &
                                      zhao_continuation_reason_physical_endpoint, zhao_continuation_reason_numerical_failure, &
                                      corrector_failure_kind == atlas_eval_physical &
                                      )
      atlas%termination_reason = 'seed_corrector_failed'
      message = 'Zhao branch-atlas seed corrector failed'
      return
    end if
    if (maxval(abs(corrected_z(1:z_dimension) - seed_coordinates(1:z_dimension))) > &
        resolved%seed_refinement_limit) then
      status = outer_plasma_numerical_failure
      atlas%termination_reason_code = zhao_continuation_reason_guard_rejected
      atlas%termination_reason = 'seed_refinement_left_local_root'
      message = 'Zhao branch-atlas seed refinement left the local root neighborhood'
      return
    end if
    z = corrected_z
    call evaluate_zhao_atlas_system( &
      params, seed_root%branch, interface_field_v_m, n, z, residual, valid, failure_kind &
      )
    if (.not. valid .or. maxval(abs(residual(1:n))) > resolved%residual_tolerance) then
      status = outer_plasma_numerical_failure
      atlas%termination_reason_code = zhao_continuation_reason_numerical_failure
      atlas%termination_reason = 'seed_residual_not_converged'
      message = 'Zhao branch-atlas seed residual is not converged'
      return
    end if
    call zhao_atlas_tangent( &
      params, seed_root%branch, interface_field_v_m, n, z, tangent, row_rank, tangent_ok, failure_kind &
      )
    if (.not. tangent_ok) then
      status = outer_plasma_numerical_failure
      atlas%termination_reason_code = zhao_continuation_reason_numerical_failure
      atlas%termination_reason = 'seed_jacobian_rank_loss'
      message = 'Zhao branch-atlas seed Jacobian has unresolved row-rank loss'
      return
    end if
    if (tangent(z_dimension)*real(resolved%eta_direction, dp) < 0.0_dp) tangent = -tangent
    eta_tangent_sign = 0
    if (abs(tangent(z_dimension)) > atlas_fold_tangent_tolerance) then
      eta_tangent_sign = merge(1, -1, tangent(z_dimension) > 0.0_dp)
    end if
    column_rate_sign = 0

    allocate (atlas%points(resolved%max_points))
    call zhao_atlas_point_from_coordinates( &
      params, seed_root%branch, interface_field_v_m, grid_points, control_length_m, n, z, &
      atlas%target_requested, atlas%target_column_m2, atlas%points(1), status, message &
      )
    if (status /= outer_plasma_ok) then
      if (status == outer_plasma_no_physical_solution) then
        atlas%termination_reason_code = zhao_continuation_reason_physical_endpoint
      else
        atlas%termination_reason_code = zhao_continuation_reason_numerical_failure
      end if
      atlas%termination_reason = 'seed_profile_failed'
      return
    end if
    call zhao_root_jump_metrics( &
      params, seed_root, atlas%points(1)%root, potential_jump, density_jump, root_jump &
      )
    if (.not. ieee_is_finite(root_jump) .or. root_jump > resolved%seed_refinement_limit) then
      status = outer_plasma_numerical_failure
      atlas%termination_reason_code = zhao_continuation_reason_guard_rejected
      atlas%termination_reason = 'seed_refinement_root_jump'
      message = 'Zhao branch-atlas seed refinement changed the physical root too far'
      return
    end if
    atlas%seed_refinement_jump = root_jump
    atlas%seed_reanchored = root_jump > branch_same_root_step_limit
    atlas%points(1)%normalized_jump_from_previous = root_jump
    atlas%point_count = 1_i32
    atlas%points(1)%tangent = tangent
    atlas%points(1)%row_rank_indicator = row_rank
    atlas%points(1)%corrector_iterations = int(iterations, i32)
    atlas%maximum_column_m2 = atlas%points(1)%root%photoelectron_column_per_area
    target_tolerance = max( &
                       column_root_relative_tolerance*max(atlas%target_column_m2, 1.0_dp), &
                       64.0_dp*epsilon(1.0_dp)*max(atlas%target_column_m2, 1.0_dp) &
                       )
    if (atlas%target_requested) then
      previous_target_residual = &
        atlas%points(1)%root%photoelectron_column_per_area - atlas%target_column_m2
      atlas%target_bracketed = abs(previous_target_residual) <= target_tolerance
      if (atlas%target_bracketed) then
        atlas%termination_reason_code = zhao_continuation_reason_target_bracketed
        atlas%termination_reason = 'target_at_seed'
      end if
    else
      previous_target_residual = 0.0_dp
    end if

    step = resolved%initial_step
    trace_loop: do point_index = 2, int(resolved%max_points)
      if (atlas%target_bracketed) exit trace_loop
      if (zhao_atlas_coordinate_endpoint( &
          seed_root%branch, n, z, tangent, resolved, atlas%termination_reason_code, &
          atlas%termination_reason)) exit trace_loop

      trial_step = min(step, resolved%maximum_step)
      corrected = .false.
      corrector_failure_kind = atlas_eval_numerical
      do halving = 0, int(resolved%max_step_halvings)
        last_eta_lower_limit = .false.
        last_eta_upper_limit = .false.
        predictor = z + trial_step*tangent
        if (predictor(z_dimension) < resolved%eta_min) then
          last_eta_lower_limit = .true.
          corrector_failure_kind = atlas_eval_numerical
        else if (predictor(z_dimension) > resolved%eta_max) then
          last_eta_upper_limit = .true.
          corrector_failure_kind = atlas_eval_numerical
        else
          call correct_zhao_atlas_point( &
            params, seed_root%branch, interface_field_v_m, n, predictor, tangent, resolved, &
            corrected_z, iterations, corrected, corrector_failure_kind &
            )
          if (corrected) then
            if (corrected_z(z_dimension) < resolved%eta_min) then
              corrected = .false.
              last_eta_lower_limit = .true.
            else if (corrected_z(z_dimension) > resolved%eta_max) then
              corrected = .false.
              last_eta_upper_limit = .true.
            else if (sqrt(sum((corrected_z(1:z_dimension) - predictor(1:z_dimension))**2)) > &
                     2.0_dp*trial_step .or. &
                     sqrt(sum((corrected_z(1:z_dimension) - z(1:z_dimension))**2)) > &
                     2.5_dp*trial_step) then
              corrected = .false.
              corrector_failure_kind = atlas_eval_numerical
            end if
          end if
          if (corrected) then
            candidate_root = zhao_charge_root_type()
            candidate_root%branch = seed_root%branch
            candidate_root%interface_field_v_m = interface_field_v_m
            candidate_root%photoelectron_population_fraction = corrected_z(z_dimension)
            call decode_unknowns( &
              params, seed_root%branch, corrected_z(1:3), candidate_root%phi0_v, &
              candidate_root%phi_m_v, candidate_root%n_swe_inf_m3, candidate_decoded &
              )
            if (candidate_decoded) then
              call zhao_root_jump_metrics( &
                params, atlas%points(point_index - 1)%root, candidate_root, &
                potential_jump, density_jump, root_jump &
                )
              if (.not. ieee_is_finite(root_jump) .or. root_jump > branch_same_root_step_limit) then
                corrected = .false.
                corrector_failure_kind = atlas_eval_numerical
              end if
            else
              corrected = .false.
              corrector_failure_kind = atlas_eval_numerical
            end if
          end if
        end if
        if (corrected) exit
        trial_step = 0.5_dp*trial_step
        if (trial_step < resolved%minimum_step) exit
      end do
      if (.not. corrected) then
        if (last_eta_lower_limit) then
          atlas%termination_reason_code = zhao_continuation_reason_search_limit
          atlas%termination_reason = 'eta_lower_search_limit'
          status = outer_plasma_ok
          message = ''
        else if (last_eta_upper_limit) then
          atlas%termination_reason_code = zhao_continuation_reason_search_limit
          atlas%termination_reason = 'eta_upper_search_limit'
          status = outer_plasma_ok
          message = ''
        else
          atlas%termination_reason_code = zhao_continuation_reason_numerical_failure
          atlas%termination_reason = 'pseudo_arclength_corrector_failed'
          status = outer_plasma_numerical_failure
          message = 'Zhao branch-atlas pseudo-arclength corrector failed'
        end if
        exit trace_loop
      end if

      call zhao_atlas_tangent( &
        params, seed_root%branch, interface_field_v_m, n, corrected_z, next_tangent, &
        next_rank, tangent_ok, failure_kind &
        )
      if (.not. tangent_ok) then
        atlas%termination_reason_code = zhao_continuation_reason_numerical_failure
        atlas%termination_reason = 'jacobian_rank_failure'
        status = outer_plasma_numerical_failure
        message = 'Zhao branch-atlas tangent Jacobian lost numerical row rank'
        exit trace_loop
      end if
      if (dot_product(tangent(1:z_dimension), next_tangent(1:z_dimension)) < 0.0_dp) then
        next_tangent = -next_tangent
      end if
      if (abs(next_tangent(z_dimension)) > atlas_fold_tangent_tolerance) then
        next_sign = merge(1, -1, next_tangent(z_dimension) > 0.0_dp)
        if (eta_tangent_sign /= 0 .and. next_sign /= eta_tangent_sign) then
          atlas%eta_fold_detected = .true.
        end if
        eta_tangent_sign = next_sign
      end if

      call zhao_atlas_point_from_coordinates( &
        params, seed_root%branch, interface_field_v_m, grid_points, control_length_m, n, &
        corrected_z, atlas%target_requested, atlas%target_column_m2, next_point, status, message &
        )
      point_ok = status == outer_plasma_ok
      if (.not. point_ok) then
        if (status == outer_plasma_no_physical_solution) then
          atlas%termination_reason_code = zhao_continuation_reason_physical_endpoint
          atlas%termination_reason = 'profile_domain_endpoint'
          status = outer_plasma_ok
          message = ''
        else
          atlas%termination_reason_code = zhao_continuation_reason_numerical_failure
          atlas%termination_reason = 'profile_numerical_failure'
        end if
        exit trace_loop
      end if
      next_point%arc_length = atlas%points(point_index - 1)%arc_length + trial_step
      next_point%accepted_step = trial_step
      next_point%tangent = next_tangent
      next_point%row_rank_indicator = next_rank
      next_point%corrector_iterations = int(iterations, i32)
      next_point%dcolumn_ds = &
        (next_point%root%photoelectron_column_per_area - &
         atlas%points(point_index - 1)%root%photoelectron_column_per_area)/trial_step
      call zhao_root_jump_metrics( &
        params, atlas%points(point_index - 1)%root, next_point%root, &
        potential_jump, density_jump, root_jump &
        )
      next_point%normalized_jump_from_previous = root_jump
      column_rate_tolerance = column_root_relative_tolerance*max( &
                              atlas%points(point_index - 1)%root%photoelectron_column_per_area, &
                              next_point%root%photoelectron_column_per_area, 1.0_dp &
                              )/max(trial_step, resolved%minimum_step)
      if (abs(next_point%dcolumn_ds) > column_rate_tolerance) then
        next_sign = merge(1, -1, next_point%dcolumn_ds > 0.0_dp)
        if (column_rate_sign /= 0 .and. next_sign /= column_rate_sign) then
          atlas%column_fold_detected = .true.
        end if
        column_rate_sign = next_sign
      end if
      atlas%points(point_index) = next_point
      atlas%point_count = int(point_index, i32)
      atlas%maximum_column_m2 = max( &
                                atlas%maximum_column_m2, next_point%root%photoelectron_column_per_area &
                                )
      if (atlas%target_requested) then
        next_target_residual = next_point%root%photoelectron_column_per_area - atlas%target_column_m2
        atlas%target_bracketed = abs(next_target_residual) <= target_tolerance .or. &
                                 previous_target_residual*next_target_residual < 0.0_dp
        previous_target_residual = next_target_residual
        if (atlas%target_bracketed) then
          atlas%termination_reason_code = zhao_continuation_reason_target_bracketed
          atlas%termination_reason = 'target_bracketed'
        end if
      end if
      z = corrected_z
      tangent = next_tangent
      if (iterations <= 4) then
        step = min(resolved%maximum_step, 1.25_dp*trial_step)
      else if (iterations >= int(resolved%max_corrector_iterations)/2) then
        step = max(resolved%minimum_step, 0.5_dp*trial_step)
      else
        step = trial_step
      end if
    end do trace_loop

    if (atlas%termination_reason_code == zhao_continuation_reason_none) then
      atlas%termination_reason_code = zhao_continuation_reason_search_limit
      atlas%termination_reason = 'point_limit'
    end if
    if (status == outer_plasma_invalid) status = outer_plasma_ok
    if (allocated(atlas%points)) then
      allocate (compact_points(atlas%point_count))
      compact_points = atlas%points(1:atlas%point_count)
      call move_alloc(compact_points, atlas%points)
    end if
  end subroutine trace_zhao_branch_atlas

  !> Trace one fixed monotone Zhao B/C branch along a straight field-column homotopy.
  !>
  !> This diagnostic augments the encoded Zhao root with eta and lambda.  It
  !> solves the branch residuals and the finite-control-volume column residual
  !> with a pseudo-arclength corrector, then uses a fixed-lambda corrector to
  !> land exactly at lambda=1.  It does not participate in runtime root choice.
  subroutine trace_zhao_field_column_homotopy( &
    params, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
    grid_points, control_length_m, seed_root, homotopy, status, message, options &
    )
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: start_field_v_m, start_column_m2
    real(dp), intent(in) :: target_field_v_m, target_column_m2
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m
    type(zhao_charge_root_type), intent(in) :: seed_root
    type(zhao_field_column_homotopy_type), intent(out) :: homotopy
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_field_column_homotopy_options_type), intent(in), optional :: options

    type(zhao_field_column_homotopy_options_type) :: resolved
    type(zhao_field_column_homotopy_point_type), allocatable :: compact_points(:)
    type(zhao_field_column_homotopy_point_type) :: next_point
    type(zhao_charge_root_type) :: candidate_root
    real(dp) :: z(homotopy_max_coordinates), seed_coordinates(homotopy_max_coordinates)
    real(dp) :: corrected_z(homotopy_max_coordinates), target_predictor(homotopy_max_coordinates)
    real(dp) :: predictor(homotopy_max_coordinates), tangent(homotopy_max_coordinates)
    real(dp) :: next_tangent(homotopy_max_coordinates), anchor_tangent(homotopy_max_coordinates)
    real(dp) :: residual(homotopy_max_residuals), step, trial_step, accepted_distance
    real(dp) :: row_rank, next_rank, target_fraction
    real(dp) :: potential_jump, density_jump, root_jump
    integer :: n, density_coordinate, eta_coordinate, lambda_coordinate, z_dimension
    integer :: point_index, halving, iterations, target_iterations, failure_kind
    integer :: corrector_failure_kind, lambda_tangent_sign, next_sign
    logical :: valid, tangent_ok, corrected, candidate_decoded, target_landed
    logical :: last_density_limit, last_eta_lower_limit, last_eta_upper_limit
    logical :: last_homotopy_lower_limit, last_homotopy_upper_limit

    status = outer_plasma_invalid
    message = ''
    resolved = zhao_field_column_homotopy_options_type()
    if (present(options)) resolved = options
    homotopy%branch = seed_root%branch
    homotopy%start_field_v_m = start_field_v_m
    homotopy%target_field_v_m = target_field_v_m
    homotopy%start_column_m2 = start_column_m2
    homotopy%target_column_m2 = target_column_m2
    homotopy%column_scale_m2 = max( &
                               1.0_dp, start_column_m2, target_column_m2, &
                               params%n_phe_ref_m3*control_length_m &
                               )

    if (seed_root%branch == 'B' .or. seed_root%branch == 'C') then
      n = 2
    else if (seed_root%branch == 'A') then
      homotopy%termination_reason_code = zhao_continuation_reason_invalid_request
      homotopy%termination_reason = 'unsupported_type_a_homotopy'
      message = 'Zhao field-column homotopy does not yet support Type A coordinates'
      return
    else
      homotopy%termination_reason_code = zhao_continuation_reason_invalid_request
      homotopy%termination_reason = 'invalid_seed_branch'
      message = 'Zhao field-column homotopy requires a Type B or C seed root'
      return
    end if
    eta_coordinate = n + 1
    lambda_coordinate = n + 2
    z_dimension = lambda_coordinate
    density_coordinate = merge(3, 2, seed_root%branch == 'A')
    if (.not. valid_zhao_params(params) .or. &
        .not. all(ieee_is_finite([ &
                                 start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
                                 control_length_m &
                                 ])) .or. &
        start_column_m2 < 0.0_dp .or. target_column_m2 < 0.0_dp .or. &
        control_length_m <= 0.0_dp .or. grid_points < 5_i32 .or. &
        resolved%max_points < 2_i32 .or. resolved%max_corrector_iterations < 1_i32 .or. &
        resolved%max_step_halvings < 0_i32 .or. &
        .not. all(ieee_is_finite([ &
                                 resolved%initial_step, resolved%minimum_step, resolved%maximum_step, &
                                 resolved%residual_tolerance, resolved%seed_refinement_limit, &
                                 resolved%eta_min, resolved%eta_max, resolved%homotopy_min, &
                                 resolved%homotopy_max, resolved%log_density_floor &
                                 ])) .or. &
        resolved%minimum_step <= 0.0_dp .or. &
        resolved%initial_step < resolved%minimum_step .or. &
        resolved%maximum_step < resolved%initial_step .or. &
        resolved%residual_tolerance <= 0.0_dp .or. resolved%seed_refinement_limit <= 0.0_dp .or. &
        resolved%eta_min < 0.0_dp .or. resolved%eta_max <= resolved%eta_min .or. &
        resolved%eta_max > column_eta_search_max .or. resolved%homotopy_min > 0.0_dp .or. &
        resolved%homotopy_max < 1.0_dp .or. resolved%homotopy_max <= resolved%homotopy_min .or. &
        resolved%log_density_floor <= -30.0_dp .or. resolved%log_density_floor >= 0.0_dp) then
      homotopy%termination_reason_code = zhao_continuation_reason_invalid_request
      homotopy%termination_reason = 'invalid_options'
      message = 'invalid Zhao field-column homotopy request'
      return
    end if
    if (((seed_root%branch == 'A' .or. seed_root%branch == 'B') .and. &
         (start_field_v_m <= 0.0_dp .or. target_field_v_m <= 0.0_dp)) .or. &
        (seed_root%branch == 'C' .and. &
         (start_field_v_m >= 0.0_dp .or. target_field_v_m >= 0.0_dp))) then
      homotopy%termination_reason_code = zhao_continuation_reason_invalid_request
      homotopy%termination_reason = 'field_sign_crossing'
      message = 'Zhao field-column homotopy must preserve the field sign of its fixed branch'
      return
    end if
    if (params%n_phe0_m3 <= 0.0_dp) then
      homotopy%termination_reason_code = zhao_continuation_reason_invalid_request
      homotopy%termination_reason = 'degenerate_zero_photoelectron_column'
      message = 'Zhao field-column homotopy requires a nonzero photoelectron source column'
      return
    end if

    z = 0.0_dp
    call encode_unknowns( &
      params, seed_root%branch, seed_root%phi0_v, seed_root%phi_m_v, &
      seed_root%n_swe_inf_m3, z(1:3), valid &
      )
    if (.not. valid) then
      homotopy%termination_reason_code = zhao_continuation_reason_invalid_request
      homotopy%termination_reason = 'invalid_seed_encoding'
      message = 'Zhao field-column homotopy seed root cannot be encoded'
      return
    end if
    z(eta_coordinate) = seed_root%photoelectron_population_fraction
    z(lambda_coordinate) = 0.0_dp
    seed_coordinates = z
    if (z(eta_coordinate) < resolved%eta_min .or. z(eta_coordinate) > resolved%eta_max) then
      homotopy%termination_reason_code = zhao_continuation_reason_invalid_request
      homotopy%termination_reason = 'seed_eta_outside_options'
      message = 'Zhao field-column homotopy seed eta lies outside its requested interval'
      return
    end if

    ! Reconcile the recorded previous root with its exact field and queue column.
    anchor_tangent = 0.0_dp
    anchor_tangent(lambda_coordinate) = 1.0_dp
    call correct_zhao_field_column_homotopy_point( &
      params, seed_root%branch, start_field_v_m, start_column_m2, &
      target_field_v_m, target_column_m2, grid_points, control_length_m, &
      homotopy%column_scale_m2, n, z, anchor_tangent, resolved, &
      corrected_z, iterations, corrected, corrector_failure_kind &
      )
    if (.not. corrected) then
      status = merge( &
               outer_plasma_no_physical_solution, outer_plasma_numerical_failure, &
               corrector_failure_kind == atlas_eval_physical &
               )
      homotopy%termination_reason_code = merge( &
                                         zhao_continuation_reason_physical_endpoint, &
                                         zhao_continuation_reason_numerical_failure, &
                                         corrector_failure_kind == atlas_eval_physical &
                                         )
      homotopy%termination_reason = 'seed_corrector_failed'
      message = 'Zhao field-column homotopy seed corrector failed'
      return
    end if
    if (corrected_z(eta_coordinate) < resolved%eta_min .or. &
        corrected_z(eta_coordinate) > resolved%eta_max .or. &
        corrected_z(density_coordinate) < resolved%log_density_floor .or. &
        abs(corrected_z(lambda_coordinate)) > resolved%residual_tolerance) then
      status = outer_plasma_invalid
      homotopy%termination_reason_code = zhao_continuation_reason_guard_rejected
      homotopy%termination_reason = 'seed_corrector_left_options'
      message = 'Zhao field-column homotopy seed correction left its requested bounds'
      return
    end if
    if (maxval(abs(corrected_z(1:z_dimension) - seed_coordinates(1:z_dimension))) > &
        resolved%seed_refinement_limit) then
      status = outer_plasma_numerical_failure
      homotopy%termination_reason_code = zhao_continuation_reason_guard_rejected
      homotopy%termination_reason = 'seed_refinement_left_local_root'
      message = 'Zhao field-column homotopy seed refinement left the local root neighborhood'
      return
    end if
    z = corrected_z
    call evaluate_zhao_field_column_homotopy_system( &
      params, seed_root%branch, start_field_v_m, start_column_m2, &
      target_field_v_m, target_column_m2, grid_points, control_length_m, &
      homotopy%column_scale_m2, n, z, residual, valid, failure_kind &
      )
    if (.not. valid .or. maxval(abs(residual(1:n + 1))) > resolved%residual_tolerance) then
      status = outer_plasma_numerical_failure
      homotopy%termination_reason_code = zhao_continuation_reason_numerical_failure
      homotopy%termination_reason = 'seed_residual_not_converged'
      message = 'Zhao field-column homotopy seed residual is not converged'
      return
    end if
    call zhao_field_column_homotopy_tangent( &
      params, seed_root%branch, start_field_v_m, start_column_m2, &
      target_field_v_m, target_column_m2, grid_points, control_length_m, &
      homotopy%column_scale_m2, n, z, tangent, row_rank, tangent_ok, failure_kind &
      )
    if (.not. tangent_ok) then
      status = outer_plasma_numerical_failure
      homotopy%termination_reason_code = zhao_continuation_reason_numerical_failure
      homotopy%termination_reason = 'seed_jacobian_rank_loss'
      message = 'Zhao field-column homotopy seed Jacobian has unresolved row-rank loss'
      return
    end if
    if (abs(tangent(lambda_coordinate)) <= atlas_fold_tangent_tolerance) then
      status = outer_plasma_numerical_failure
      homotopy%termination_reason_code = zhao_continuation_reason_numerical_failure
      homotopy%termination_reason = 'seed_homotopy_tangent_degenerate'
      message = 'Zhao field-column homotopy has no resolved positive-lambda seed direction'
      return
    end if
    if (tangent(lambda_coordinate) < 0.0_dp) tangent = -tangent
    lambda_tangent_sign = 1

    allocate (homotopy%points(resolved%max_points))
    call zhao_field_column_homotopy_point_from_coordinates( &
      params, seed_root%branch, start_field_v_m, start_column_m2, &
      target_field_v_m, target_column_m2, grid_points, control_length_m, &
      homotopy%column_scale_m2, n, z, homotopy%points(1), status, message &
      )
    if (status /= outer_plasma_ok) then
      homotopy%termination_reason_code = merge( &
                                         zhao_continuation_reason_physical_endpoint, &
                                         zhao_continuation_reason_numerical_failure, &
                                         status == outer_plasma_no_physical_solution &
                                         )
      homotopy%termination_reason = 'seed_profile_failed'
      return
    end if
    call zhao_root_jump_metrics( &
      params, seed_root, homotopy%points(1)%root, potential_jump, density_jump, root_jump &
      )
    if (.not. ieee_is_finite(root_jump) .or. root_jump > resolved%seed_refinement_limit) then
      status = outer_plasma_numerical_failure
      homotopy%termination_reason_code = zhao_continuation_reason_guard_rejected
      homotopy%termination_reason = 'seed_refinement_root_jump'
      message = 'Zhao field-column homotopy seed refinement changed the physical root too far'
      return
    end if
    homotopy%seed_refinement_jump = root_jump
    homotopy%seed_reanchored = root_jump > branch_same_root_step_limit
    homotopy%points(1)%normalized_jump_from_previous = root_jump
    homotopy%points(1)%tangent = tangent
    homotopy%points(1)%row_rank_indicator = row_rank
    homotopy%points(1)%corrector_iterations = int(iterations, i32)
    homotopy%point_count = 1_i32

    step = resolved%initial_step
    trace_loop: do point_index = 2, int(resolved%max_points)
      if (homotopy%target_reached) exit trace_loop
      if (zhao_field_column_homotopy_coordinate_endpoint( &
          seed_root%branch, n, z, tangent, resolved, homotopy%termination_reason_code, &
          homotopy%termination_reason)) exit trace_loop

      trial_step = min(step, resolved%maximum_step)
      corrected = .false.
      target_landed = .false.
      corrector_failure_kind = atlas_eval_numerical
      do halving = 0, int(resolved%max_step_halvings)
        last_density_limit = .false.
        last_eta_lower_limit = .false.
        last_eta_upper_limit = .false.
        last_homotopy_lower_limit = .false.
        last_homotopy_upper_limit = .false.
        predictor = z + trial_step*tangent
        if (predictor(eta_coordinate) < resolved%eta_min) then
          last_eta_lower_limit = .true.
        else if (predictor(eta_coordinate) > resolved%eta_max) then
          last_eta_upper_limit = .true.
        else if (predictor(lambda_coordinate) < resolved%homotopy_min) then
          last_homotopy_lower_limit = .true.
        else if (predictor(lambda_coordinate) > resolved%homotopy_max .and. &
                 .not. (z(lambda_coordinate) < 1.0_dp .and. &
                        predictor(lambda_coordinate) >= 1.0_dp)) then
          last_homotopy_upper_limit = .true.
        else
          call correct_zhao_field_column_homotopy_point( &
            params, seed_root%branch, start_field_v_m, start_column_m2, &
            target_field_v_m, target_column_m2, grid_points, control_length_m, &
            homotopy%column_scale_m2, n, predictor, tangent, resolved, &
            corrected_z, iterations, corrected, corrector_failure_kind &
            )
          if (corrected) then
            if (corrected_z(eta_coordinate) < resolved%eta_min) then
              corrected = .false.
              last_eta_lower_limit = .true.
            else if (corrected_z(eta_coordinate) > resolved%eta_max) then
              corrected = .false.
              last_eta_upper_limit = .true.
            else if (corrected_z(lambda_coordinate) < resolved%homotopy_min) then
              corrected = .false.
              last_homotopy_lower_limit = .true.
            else if (corrected_z(lambda_coordinate) > resolved%homotopy_max .and. &
                     .not. (z(lambda_coordinate) < 1.0_dp .and. &
                            corrected_z(lambda_coordinate) >= 1.0_dp)) then
              corrected = .false.
              last_homotopy_upper_limit = .true.
            else if (sqrt(sum((corrected_z(1:z_dimension) - predictor(1:z_dimension))**2)) > &
                     2.0_dp*trial_step .or. &
                     sqrt(sum((corrected_z(1:z_dimension) - z(1:z_dimension))**2)) > &
                     2.5_dp*trial_step) then
              corrected = .false.
              corrector_failure_kind = atlas_eval_numerical
            end if
          end if
          if (corrected) then
            candidate_root = zhao_charge_root_type()
            candidate_root%branch = seed_root%branch
            candidate_root%interface_field_v_m = start_field_v_m + &
                                                 corrected_z(lambda_coordinate)* &
                                                 (target_field_v_m - start_field_v_m)
            candidate_root%photoelectron_population_fraction = corrected_z(eta_coordinate)
            call decode_unknowns( &
              params, seed_root%branch, corrected_z(1:3), candidate_root%phi0_v, &
              candidate_root%phi_m_v, candidate_root%n_swe_inf_m3, candidate_decoded &
              )
            if (candidate_decoded) then
              call zhao_root_jump_metrics( &
                params, homotopy%points(point_index - 1)%root, candidate_root, &
                potential_jump, density_jump, root_jump &
                )
              if (.not. ieee_is_finite(root_jump) .or. root_jump > branch_same_root_step_limit) then
                corrected = .false.
                corrector_failure_kind = atlas_eval_numerical
              end if
            else
              corrected = .false.
              corrector_failure_kind = atlas_eval_numerical
            end if
          end if
          if (corrected .and. z(lambda_coordinate) < 1.0_dp .and. &
              corrected_z(lambda_coordinate) >= 1.0_dp) then
            target_fraction = (1.0_dp - z(lambda_coordinate))/ &
                              (corrected_z(lambda_coordinate) - z(lambda_coordinate))
            target_predictor = z + target_fraction*(corrected_z - z)
            target_predictor(lambda_coordinate) = 1.0_dp
            anchor_tangent = 0.0_dp
            anchor_tangent(lambda_coordinate) = 1.0_dp
            call correct_zhao_field_column_homotopy_point( &
              params, seed_root%branch, start_field_v_m, start_column_m2, &
              target_field_v_m, target_column_m2, grid_points, control_length_m, &
              homotopy%column_scale_m2, n, target_predictor, anchor_tangent, resolved, &
              corrected_z, target_iterations, corrected, corrector_failure_kind &
              )
            if (corrected) then
              corrected_z(lambda_coordinate) = 1.0_dp
              if (corrected_z(eta_coordinate) < resolved%eta_min) then
                corrected = .false.
                last_eta_lower_limit = .true.
              else if (corrected_z(eta_coordinate) > resolved%eta_max) then
                corrected = .false.
                last_eta_upper_limit = .true.
              else if (corrected_z(density_coordinate) < resolved%log_density_floor) then
                corrected = .false.
                last_density_limit = .true.
              end if
            end if
            if (corrected) then
              call evaluate_zhao_field_column_homotopy_system( &
                params, seed_root%branch, start_field_v_m, start_column_m2, &
                target_field_v_m, target_column_m2, grid_points, control_length_m, &
                homotopy%column_scale_m2, n, corrected_z, residual, valid, failure_kind &
                )
              if (.not. valid .or. &
                  maxval(abs(residual(1:n + 1))) > resolved%residual_tolerance) then
                corrected = .false.
                corrector_failure_kind = failure_kind
              end if
            end if
            if (corrected) then
              candidate_root = zhao_charge_root_type()
              candidate_root%branch = seed_root%branch
              candidate_root%interface_field_v_m = target_field_v_m
              candidate_root%photoelectron_population_fraction = corrected_z(eta_coordinate)
              call decode_unknowns( &
                params, seed_root%branch, corrected_z(1:3), candidate_root%phi0_v, &
                candidate_root%phi_m_v, candidate_root%n_swe_inf_m3, candidate_decoded &
                )
              if (candidate_decoded) then
                call zhao_root_jump_metrics( &
                  params, homotopy%points(point_index - 1)%root, candidate_root, &
                  potential_jump, density_jump, root_jump &
                  )
                if (.not. ieee_is_finite(root_jump) .or. &
                    root_jump > branch_same_root_step_limit) corrected = .false.
              else
                corrected = .false.
              end if
              if (corrected) then
                iterations = iterations + target_iterations
                target_landed = .true.
              end if
            end if
          end if
        end if
        if (corrected) exit
        trial_step = 0.5_dp*trial_step
        if (trial_step < resolved%minimum_step) exit
      end do
      if (.not. corrected) then
        if (last_density_limit) then
          homotopy%termination_reason_code = zhao_continuation_reason_search_limit
          homotopy%termination_reason = 'ambient_density_floor_limit'
          status = outer_plasma_ok
          message = ''
        else if (last_eta_lower_limit) then
          homotopy%termination_reason_code = zhao_continuation_reason_search_limit
          homotopy%termination_reason = 'eta_lower_search_limit'
          status = outer_plasma_ok
          message = ''
        else if (last_eta_upper_limit) then
          homotopy%termination_reason_code = zhao_continuation_reason_search_limit
          homotopy%termination_reason = 'eta_upper_search_limit'
          status = outer_plasma_ok
          message = ''
        else if (last_homotopy_lower_limit) then
          homotopy%termination_reason_code = zhao_continuation_reason_search_limit
          homotopy%termination_reason = 'homotopy_lower_search_limit'
          status = outer_plasma_ok
          message = ''
        else if (last_homotopy_upper_limit) then
          homotopy%termination_reason_code = zhao_continuation_reason_search_limit
          homotopy%termination_reason = 'homotopy_upper_search_limit'
          status = outer_plasma_ok
          message = ''
        else
          homotopy%termination_reason_code = zhao_continuation_reason_numerical_failure
          homotopy%termination_reason = 'pseudo_arclength_corrector_failed'
          status = outer_plasma_numerical_failure
          message = 'Zhao field-column homotopy pseudo-arclength corrector failed'
        end if
        exit trace_loop
      end if

      accepted_distance = sqrt(sum((corrected_z(1:z_dimension) - z(1:z_dimension))**2))
      call zhao_field_column_homotopy_tangent( &
        params, seed_root%branch, start_field_v_m, start_column_m2, &
        target_field_v_m, target_column_m2, grid_points, control_length_m, &
        homotopy%column_scale_m2, n, corrected_z, next_tangent, &
        next_rank, tangent_ok, failure_kind &
        )
      if (.not. tangent_ok) then
        homotopy%termination_reason_code = zhao_continuation_reason_numerical_failure
        homotopy%termination_reason = 'jacobian_rank_failure'
        status = outer_plasma_numerical_failure
        message = 'Zhao field-column homotopy tangent Jacobian lost numerical row rank'
        exit trace_loop
      end if
      if (dot_product(tangent(1:z_dimension), next_tangent(1:z_dimension)) < 0.0_dp) then
        next_tangent = -next_tangent
      end if
      if (abs(next_tangent(lambda_coordinate)) > atlas_fold_tangent_tolerance) then
        next_sign = merge(1, -1, next_tangent(lambda_coordinate) > 0.0_dp)
        if (next_sign /= lambda_tangent_sign) homotopy%homotopy_fold_detected = .true.
        lambda_tangent_sign = next_sign
      end if

      call zhao_field_column_homotopy_point_from_coordinates( &
        params, seed_root%branch, start_field_v_m, start_column_m2, &
        target_field_v_m, target_column_m2, grid_points, control_length_m, &
        homotopy%column_scale_m2, n, corrected_z, next_point, status, message &
        )
      if (status /= outer_plasma_ok) then
        if (status == outer_plasma_no_physical_solution) then
          homotopy%termination_reason_code = zhao_continuation_reason_physical_endpoint
          homotopy%termination_reason = 'profile_domain_endpoint'
          status = outer_plasma_ok
          message = ''
        else
          homotopy%termination_reason_code = zhao_continuation_reason_numerical_failure
          homotopy%termination_reason = 'profile_numerical_failure'
        end if
        exit trace_loop
      end if
      next_point%arc_length = homotopy%points(point_index - 1)%arc_length + accepted_distance
      next_point%accepted_step = accepted_distance
      next_point%tangent = next_tangent
      next_point%row_rank_indicator = next_rank
      next_point%corrector_iterations = int(iterations, i32)
      call zhao_root_jump_metrics( &
        params, homotopy%points(point_index - 1)%root, next_point%root, &
        potential_jump, density_jump, root_jump &
        )
      next_point%normalized_jump_from_previous = root_jump
      homotopy%points(point_index) = next_point
      homotopy%point_count = int(point_index, i32)
      z = corrected_z
      tangent = next_tangent
      if (target_landed) then
        homotopy%target_reached = .true.
        homotopy%termination_reason_code = zhao_continuation_reason_converged
        homotopy%termination_reason = 'target_reached'
        status = outer_plasma_ok
        message = ''
        exit trace_loop
      end if
      if (iterations <= 4) then
        step = min(resolved%maximum_step, 1.25_dp*trial_step)
      else if (iterations >= int(resolved%max_corrector_iterations)/2) then
        step = max(resolved%minimum_step, 0.5_dp*trial_step)
      else
        step = trial_step
      end if
    end do trace_loop

    if (homotopy%termination_reason_code == zhao_continuation_reason_none) then
      homotopy%termination_reason_code = zhao_continuation_reason_search_limit
      homotopy%termination_reason = 'point_limit'
    end if
    if (status == outer_plasma_invalid) status = outer_plasma_ok
    if (allocated(homotopy%points)) then
      allocate (compact_points(homotopy%point_count))
      compact_points = homotopy%points(1:homotopy%point_count)
      call move_alloc(compact_points, homotopy%points)
    end if
  end subroutine trace_zhao_field_column_homotopy

  subroutine evaluate_zhao_field_column_homotopy_system( &
    params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
    grid_points, control_length_m, column_scale_m2, n, z, residual, valid, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: start_field_v_m, start_column_m2
    real(dp), intent(in) :: target_field_v_m, target_column_m2
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, column_scale_m2
    integer, intent(in) :: n
    real(dp), intent(in) :: z(homotopy_max_coordinates)
    real(dp), intent(out) :: residual(homotopy_max_residuals)
    logical, intent(out) :: valid
    integer, intent(out) :: failure_kind

    type(zhao_branch_atlas_point_type) :: fixed_point
    real(dp) :: fixed_z(atlas_max_coordinates), root_residual(3)
    real(dp) :: lambda, interface_field_v_m, prescribed_column_m2
    integer(i32) :: point_status
    character(len=256) :: point_message

    residual = 0.0_dp
    valid = .false.
    failure_kind = atlas_eval_numerical
    if (n + 2 > homotopy_max_coordinates .or. n + 1 > homotopy_max_residuals) return
    lambda = z(n + 2)
    if (.not. ieee_is_finite(lambda)) return
    interface_field_v_m = start_field_v_m + lambda*(target_field_v_m - start_field_v_m)
    prescribed_column_m2 = start_column_m2 + lambda*(target_column_m2 - start_column_m2)
    if (.not. all(ieee_is_finite([interface_field_v_m, prescribed_column_m2])) .or. &
        prescribed_column_m2 < 0.0_dp .or. column_scale_m2 <= 0.0_dp) then
      failure_kind = atlas_eval_physical
      return
    end if
    if (((branch == 'A' .or. branch == 'B') .and. interface_field_v_m <= 0.0_dp) .or. &
        (branch == 'C' .and. interface_field_v_m >= 0.0_dp)) then
      failure_kind = atlas_eval_physical
      return
    end if
    fixed_z = 0.0_dp
    fixed_z(1:n + 1) = z(1:n + 1)
    call evaluate_zhao_atlas_system( &
      params, branch, interface_field_v_m, n, fixed_z, root_residual, valid, failure_kind &
      )
    if (.not. valid) return
    call zhao_atlas_point_from_coordinates( &
      params, branch, interface_field_v_m, grid_points, control_length_m, n, fixed_z, &
      .true., prescribed_column_m2, fixed_point, point_status, point_message &
      )
    if (point_status /= outer_plasma_ok) then
      valid = .false.
      failure_kind = merge( &
                     atlas_eval_physical, atlas_eval_numerical, &
                     point_status == outer_plasma_no_physical_solution &
                     )
      return
    end if
    residual(1:n) = root_residual(1:n)
    residual(n + 1) = fixed_point%root%photoelectron_column_residual_per_area/column_scale_m2
    valid = all(ieee_is_finite(residual(1:n + 1)))
    failure_kind = merge(atlas_eval_ok, atlas_eval_numerical, valid)
  end subroutine evaluate_zhao_field_column_homotopy_system

  subroutine numerical_zhao_field_column_homotopy_jacobian( &
    params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
    grid_points, control_length_m, column_scale_m2, n, z, residual, jacobian, &
    success, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: start_field_v_m, start_column_m2
    real(dp), intent(in) :: target_field_v_m, target_column_m2
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, column_scale_m2
    integer, intent(in) :: n
    real(dp), intent(in) :: z(homotopy_max_coordinates), residual(homotopy_max_residuals)
    real(dp), intent(out) :: jacobian(homotopy_max_residuals, homotopy_max_coordinates)
    logical, intent(out) :: success
    integer, intent(out) :: failure_kind

    real(dp) :: plus_z(homotopy_max_coordinates), minus_z(homotopy_max_coordinates)
    real(dp) :: plus_residual(homotopy_max_residuals), minus_residual(homotopy_max_residuals)
    real(dp) :: h
    integer :: column, plus_kind, minus_kind, residual_count, z_dimension
    logical :: plus_valid, minus_valid

    jacobian = 0.0_dp
    success = .true.
    failure_kind = atlas_eval_ok
    residual_count = n + 1
    z_dimension = n + 2
    do column = 1, z_dimension
      h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(z(column)))
      plus_z = z
      minus_z = z
      plus_z(column) = plus_z(column) + h
      minus_z(column) = minus_z(column) - h
      call evaluate_zhao_field_column_homotopy_system( &
        params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
        grid_points, control_length_m, column_scale_m2, n, plus_z, plus_residual, &
        plus_valid, plus_kind &
        )
      call evaluate_zhao_field_column_homotopy_system( &
        params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
        grid_points, control_length_m, column_scale_m2, n, minus_z, minus_residual, &
        minus_valid, minus_kind &
        )
      if (plus_valid .and. minus_valid) then
        jacobian(1:residual_count, column) = &
          (plus_residual(1:residual_count) - minus_residual(1:residual_count))/(2.0_dp*h)
      else if (plus_valid) then
        jacobian(1:residual_count, column) = &
          (plus_residual(1:residual_count) - residual(1:residual_count))/h
      else if (minus_valid) then
        jacobian(1:residual_count, column) = &
          (residual(1:residual_count) - minus_residual(1:residual_count))/h
      else
        success = .false.
        failure_kind = merge( &
                       atlas_eval_physical, atlas_eval_numerical, &
                       plus_kind == atlas_eval_physical .and. minus_kind == atlas_eval_physical &
                       )
        return
      end if
    end do
    success = all(ieee_is_finite(jacobian(1:residual_count, 1:z_dimension)))
    if (.not. success) failure_kind = atlas_eval_numerical
  end subroutine numerical_zhao_field_column_homotopy_jacobian

  subroutine zhao_field_column_homotopy_tangent( &
    params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
    grid_points, control_length_m, column_scale_m2, n, z, tangent, row_rank_indicator, &
    success, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: start_field_v_m, start_column_m2
    real(dp), intent(in) :: target_field_v_m, target_column_m2
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, column_scale_m2
    integer, intent(in) :: n
    real(dp), intent(in) :: z(homotopy_max_coordinates)
    real(dp), intent(out) :: tangent(homotopy_max_coordinates), row_rank_indicator
    logical, intent(out) :: success
    integer, intent(out) :: failure_kind

    real(dp) :: residual(homotopy_max_residuals)
    real(dp) :: jacobian(homotopy_max_residuals, homotopy_max_coordinates)
    real(dp) :: scaled_jacobian(homotopy_max_residuals, homotopy_max_coordinates)
    real(dp) :: cofactors(homotopy_max_coordinates), minor3(3, 3)
    real(dp) :: raw_norm, row_norm
    integer :: column, source_column, minor_column, row, residual_count, z_dimension
    logical :: valid, jacobian_ok

    tangent = 0.0_dp
    row_rank_indicator = 0.0_dp
    success = .false.
    residual_count = n + 1
    z_dimension = n + 2
    call evaluate_zhao_field_column_homotopy_system( &
      params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
      grid_points, control_length_m, column_scale_m2, n, z, residual, valid, failure_kind &
      )
    if (.not. valid) return
    call numerical_zhao_field_column_homotopy_jacobian( &
      params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
      grid_points, control_length_m, column_scale_m2, n, z, residual, jacobian, &
      jacobian_ok, failure_kind &
      )
    if (.not. jacobian_ok) return
    scaled_jacobian = jacobian
    do row = 1, residual_count
      row_norm = sqrt(sum(jacobian(row, 1:z_dimension)**2))
      if (.not. ieee_is_finite(row_norm) .or. row_norm <= tiny(1.0_dp)) return
      scaled_jacobian(row, 1:z_dimension) = jacobian(row, 1:z_dimension)/row_norm
    end do
    cofactors = 0.0_dp
    do column = 1, z_dimension
      minor_column = 0
      minor3 = 0.0_dp
      do source_column = 1, z_dimension
        if (source_column == column) cycle
        minor_column = minor_column + 1
        minor3(:, minor_column) = scaled_jacobian(1:3, source_column)
      end do
      cofactors(column) = (-1.0_dp)**(column + 1)*determinant3(minor3)
    end do
    raw_norm = sqrt(sum(cofactors(1:z_dimension)**2))
    if (.not. ieee_is_finite(raw_norm) .or. raw_norm <= tiny(1.0_dp)) return
    row_rank_indicator = raw_norm
    if (row_rank_indicator <= 1.0e-12_dp) return
    tangent(1:z_dimension) = cofactors(1:z_dimension)/raw_norm
    failure_kind = atlas_eval_ok
    success = .true.
  end subroutine zhao_field_column_homotopy_tangent

  subroutine correct_zhao_field_column_homotopy_point( &
    params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
    grid_points, control_length_m, column_scale_m2, n, predictor, tangent, options, &
    corrected_z, iterations, success, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: start_field_v_m, start_column_m2
    real(dp), intent(in) :: target_field_v_m, target_column_m2
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, column_scale_m2
    integer, intent(in) :: n
    real(dp), intent(in) :: predictor(homotopy_max_coordinates)
    real(dp), intent(in) :: tangent(homotopy_max_coordinates)
    type(zhao_field_column_homotopy_options_type), intent(in) :: options
    real(dp), intent(out) :: corrected_z(homotopy_max_coordinates)
    integer, intent(out) :: iterations
    logical, intent(out) :: success
    integer, intent(out) :: failure_kind

    real(dp) :: z(homotopy_max_coordinates), residual(homotopy_max_residuals)
    real(dp) :: jacobian(homotopy_max_residuals, homotopy_max_coordinates)
    real(dp) :: system(homotopy_max_coordinates, homotopy_max_coordinates)
    real(dp) :: rhs(homotopy_max_coordinates), delta(homotopy_max_coordinates)
    real(dp) :: trial_z(homotopy_max_coordinates), trial_residual(homotopy_max_residuals)
    real(dp) :: norm, trial_norm, arc_residual, trial_arc_residual, backtrack_step
    integer :: iteration, backtrack, residual_count, z_dimension, trial_kind
    logical :: valid, jacobian_ok, linear_ok, trial_valid

    residual_count = n + 1
    z_dimension = n + 2
    z = predictor
    corrected_z = predictor
    success = .false.
    failure_kind = atlas_eval_numerical
    iterations = 0
    call evaluate_zhao_field_column_homotopy_system( &
      params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
      grid_points, control_length_m, column_scale_m2, n, z, residual, valid, failure_kind &
      )
    if (.not. valid) return
    arc_residual = dot_product(tangent(1:z_dimension), z(1:z_dimension) - predictor(1:z_dimension))
    norm = max(maxval(abs(residual(1:residual_count))), abs(arc_residual))
    do iteration = 0, int(options%max_corrector_iterations)
      if (norm <= options%residual_tolerance) then
        success = .true.
        exit
      end if
      if (iteration == int(options%max_corrector_iterations)) exit
      call numerical_zhao_field_column_homotopy_jacobian( &
        params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
        grid_points, control_length_m, column_scale_m2, n, z, residual, jacobian, &
        jacobian_ok, failure_kind &
        )
      if (.not. jacobian_ok) exit
      system = 0.0_dp
      rhs = 0.0_dp
      system(1:residual_count, 1:z_dimension) = jacobian(1:residual_count, 1:z_dimension)
      system(z_dimension, 1:z_dimension) = tangent(1:z_dimension)
      rhs(1:residual_count) = -residual(1:residual_count)
      rhs(z_dimension) = -arc_residual
      call solve_zhao_homotopy_system(system, rhs, z_dimension, delta, linear_ok)
      if (.not. linear_ok) exit
      backtrack_step = 1.0_dp
      do backtrack = 1, root_max_backtracks
        trial_z = z + backtrack_step*delta
        call evaluate_zhao_field_column_homotopy_system( &
          params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
          grid_points, control_length_m, column_scale_m2, n, trial_z, trial_residual, &
          trial_valid, trial_kind &
          )
        if (trial_valid) then
          trial_arc_residual = dot_product( &
                               tangent(1:z_dimension), &
                               trial_z(1:z_dimension) - predictor(1:z_dimension) &
                               )
          trial_norm = max(maxval(abs(trial_residual(1:residual_count))), abs(trial_arc_residual))
          if (trial_norm < norm) then
            z = trial_z
            residual = trial_residual
            arc_residual = trial_arc_residual
            norm = trial_norm
            failure_kind = atlas_eval_ok
            exit
          end if
        end if
        backtrack_step = 0.5_dp*backtrack_step
      end do
      if (backtrack > root_max_backtracks) exit
    end do
    iterations = iteration
    corrected_z = z
    if (.not. success) failure_kind = atlas_eval_numerical
  end subroutine correct_zhao_field_column_homotopy_point

  subroutine solve_zhao_homotopy_system(a_in, b_in, n, x, success)
    real(dp), intent(in) :: a_in(homotopy_max_coordinates, homotopy_max_coordinates)
    real(dp), intent(in) :: b_in(homotopy_max_coordinates)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(homotopy_max_coordinates)
    logical, intent(out) :: success

    real(dp) :: a(homotopy_max_coordinates, homotopy_max_coordinates)
    real(dp) :: b(homotopy_max_coordinates), row_scale, factor, pivot_value, tmp
    integer :: i, j, k, pivot

    a = a_in
    b = b_in
    x = 0.0_dp
    success = .false.
    do i = 1, n
      row_scale = maxval(abs(a(i, 1:n)))
      if (.not. ieee_is_finite(row_scale) .or. row_scale <= tiny(1.0_dp)) return
      a(i, 1:n) = a(i, 1:n)/row_scale
      b(i) = b(i)/row_scale
    end do
    do k = 1, n
      pivot = k
      do i = k + 1, n
        if (abs(a(i, k)) > abs(a(pivot, k))) pivot = i
      end do
      if (.not. ieee_is_finite(a(pivot, k)) .or. abs(a(pivot, k)) <= 1.0e-13_dp) return
      if (pivot /= k) then
        do j = k, n
          tmp = a(k, j)
          a(k, j) = a(pivot, j)
          a(pivot, j) = tmp
        end do
        tmp = b(k)
        b(k) = b(pivot)
        b(pivot) = tmp
      end if
      pivot_value = a(k, k)
      do i = k + 1, n
        factor = a(i, k)/pivot_value
        a(i, k:n) = a(i, k:n) - factor*a(k, k:n)
        b(i) = b(i) - factor*b(k)
      end do
    end do
    do i = n, 1, -1
      x(i) = b(i)
      do j = i + 1, n
        x(i) = x(i) - a(i, j)*x(j)
      end do
      x(i) = x(i)/a(i, i)
    end do
    success = all(ieee_is_finite(x(1:n)))
  end subroutine solve_zhao_homotopy_system

  subroutine zhao_field_column_homotopy_point_from_coordinates( &
    params, branch, start_field_v_m, start_column_m2, target_field_v_m, target_column_m2, &
    grid_points, control_length_m, column_scale_m2, n, z, point, status, message &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: start_field_v_m, start_column_m2
    real(dp), intent(in) :: target_field_v_m, target_column_m2
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, column_scale_m2
    integer, intent(in) :: n
    real(dp), intent(in) :: z(homotopy_max_coordinates)
    type(zhao_field_column_homotopy_point_type), intent(out) :: point
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(zhao_branch_atlas_point_type) :: fixed_point
    real(dp) :: fixed_z(atlas_max_coordinates), interface_field_v_m, prescribed_column_m2

    point%homotopy_fraction = z(n + 2)
    interface_field_v_m = start_field_v_m + &
                          point%homotopy_fraction*(target_field_v_m - start_field_v_m)
    prescribed_column_m2 = start_column_m2 + &
                           point%homotopy_fraction*(target_column_m2 - start_column_m2)
    point%prescribed_column_m2 = prescribed_column_m2
    fixed_z = 0.0_dp
    fixed_z(1:n + 1) = z(1:n + 1)
    call zhao_atlas_point_from_coordinates( &
      params, branch, interface_field_v_m, grid_points, control_length_m, n, fixed_z, &
      .true., prescribed_column_m2, fixed_point, status, message &
      )
    if (status /= outer_plasma_ok) return
    point%root = fixed_point%root
    point%normalized_column_residual = &
      point%root%photoelectron_column_residual_per_area/column_scale_m2
  end subroutine zhao_field_column_homotopy_point_from_coordinates

  logical function zhao_field_column_homotopy_coordinate_endpoint( &
    branch, n, z, tangent, options, reason_code, reason &
    ) result(at_endpoint)
    character(len=1), intent(in) :: branch
    integer, intent(in) :: n
    real(dp), intent(in) :: z(homotopy_max_coordinates), tangent(homotopy_max_coordinates)
    type(zhao_field_column_homotopy_options_type), intent(in) :: options
    integer(i32), intent(out) :: reason_code
    character(len=*), intent(out) :: reason

    integer :: density_coordinate, eta_coordinate, lambda_coordinate

    at_endpoint = .false.
    reason_code = zhao_continuation_reason_none
    reason = 'none'
    density_coordinate = merge(3, 2, branch == 'A')
    eta_coordinate = n + 1
    lambda_coordinate = n + 2
    if (z(density_coordinate) <= options%log_density_floor .and. &
        tangent(density_coordinate) < 0.0_dp) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'ambient_density_floor_limit'
    else if (z(eta_coordinate) <= options%eta_min + options%minimum_step .and. &
             tangent(eta_coordinate) < 0.0_dp) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'eta_lower_search_limit'
    else if (z(eta_coordinate) >= options%eta_max - options%minimum_step .and. &
             tangent(eta_coordinate) > 0.0_dp) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'eta_upper_search_limit'
    else if (z(lambda_coordinate) <= options%homotopy_min + options%minimum_step .and. &
             tangent(lambda_coordinate) < 0.0_dp) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'homotopy_lower_search_limit'
    else if (z(lambda_coordinate) >= options%homotopy_max - options%minimum_step .and. &
             tangent(lambda_coordinate) > 0.0_dp .and. &
             .not. (z(lambda_coordinate) < 1.0_dp .and. options%homotopy_max >= 1.0_dp)) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'homotopy_upper_search_limit'
    end if
  end function zhao_field_column_homotopy_coordinate_endpoint

  subroutine evaluate_zhao_atlas_system( &
    params, branch, interface_field_v_m, n, z, residual, valid, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: interface_field_v_m
    integer, intent(in) :: n
    real(dp), intent(in) :: z(atlas_max_coordinates)
    real(dp), intent(out) :: residual(3)
    logical, intent(out) :: valid
    integer, intent(out) :: failure_kind

    type(zhao_params_type) :: trial_params
    real(dp) :: phi0_v, phi_m_v, density_m3, field_scale, target_field_hat
    logical :: decoded

    residual = 0.0_dp
    valid = .false.
    failure_kind = atlas_eval_numerical
    if (n + 1 > atlas_max_coordinates) return
    if (.not. ieee_is_finite(z(n + 1)) .or. z(n + 1) < 0.0_dp .or. &
        z(n + 1) > column_eta_search_max) then
      failure_kind = atlas_eval_numerical
      return
    end if
    call decode_unknowns( &
      params, branch, z(1:3), phi0_v, phi_m_v, density_m3, decoded &
      )
    if (.not. decoded) then
      failure_kind = atlas_eval_numerical
      return
    end if
    trial_params = params
    trial_params%photoelectron_population_fraction = z(n + 1)
    field_scale = trial_params%t_phe_ev/trial_params%lambda_d_phe_ref_m
    target_field_hat = interface_field_v_m/field_scale
    call eval_charge_residual( &
      trial_params, branch, target_field_hat, z(1:3), residual, valid &
      )
    if (valid) then
      failure_kind = atlas_eval_ok
    else if (.not. ion_accessible(trial_params, max(phi0_v/trial_params%t_phe_ev, 0.0_dp))) then
      failure_kind = atlas_eval_physical
    end if
  end subroutine evaluate_zhao_atlas_system

  subroutine numerical_zhao_atlas_jacobian( &
    params, branch, interface_field_v_m, n, z, residual, jacobian, success, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: interface_field_v_m
    integer, intent(in) :: n
    real(dp), intent(in) :: z(atlas_max_coordinates), residual(3)
    real(dp), intent(out) :: jacobian(3, atlas_max_coordinates)
    logical, intent(out) :: success
    integer, intent(out) :: failure_kind

    real(dp) :: plus_z(atlas_max_coordinates), minus_z(atlas_max_coordinates)
    real(dp) :: plus_residual(3), minus_residual(3), h
    integer :: column, plus_kind, minus_kind, z_dimension
    logical :: plus_valid, minus_valid

    jacobian = 0.0_dp
    success = .true.
    failure_kind = atlas_eval_ok
    z_dimension = n + 1
    do column = 1, z_dimension
      h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(z(column)))
      plus_z = z
      minus_z = z
      plus_z(column) = plus_z(column) + h
      minus_z(column) = minus_z(column) - h
      call evaluate_zhao_atlas_system( &
        params, branch, interface_field_v_m, n, plus_z, plus_residual, plus_valid, plus_kind &
        )
      call evaluate_zhao_atlas_system( &
        params, branch, interface_field_v_m, n, minus_z, minus_residual, minus_valid, minus_kind &
        )
      if (plus_valid .and. minus_valid) then
        jacobian(1:n, column) = (plus_residual(1:n) - minus_residual(1:n))/(2.0_dp*h)
      else if (plus_valid) then
        jacobian(1:n, column) = (plus_residual(1:n) - residual(1:n))/h
      else if (minus_valid) then
        jacobian(1:n, column) = (residual(1:n) - minus_residual(1:n))/h
      else
        success = .false.
        failure_kind = merge( &
                       atlas_eval_physical, atlas_eval_numerical, &
                       plus_kind == atlas_eval_physical .and. minus_kind == atlas_eval_physical &
                       )
        return
      end if
    end do
    success = all(ieee_is_finite(jacobian(1:n, 1:z_dimension)))
    if (.not. success) failure_kind = atlas_eval_numerical
  end subroutine numerical_zhao_atlas_jacobian

  subroutine zhao_atlas_tangent( &
    params, branch, interface_field_v_m, n, z, tangent, row_rank_indicator, success, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: interface_field_v_m
    integer, intent(in) :: n
    real(dp), intent(in) :: z(atlas_max_coordinates)
    real(dp), intent(out) :: tangent(atlas_max_coordinates), row_rank_indicator
    logical, intent(out) :: success
    integer, intent(out) :: failure_kind

    real(dp) :: residual(3), jacobian(3, atlas_max_coordinates)
    real(dp) :: scaled_jacobian(3, atlas_max_coordinates)
    real(dp) :: cofactors(atlas_max_coordinates), minor(3, 3), raw_norm, row_norm
    integer :: column, source_column, minor_column, row, z_dimension
    logical :: valid, jacobian_ok

    tangent = 0.0_dp
    row_rank_indicator = 0.0_dp
    success = .false.
    z_dimension = n + 1
    call evaluate_zhao_atlas_system( &
      params, branch, interface_field_v_m, n, z, residual, valid, failure_kind &
      )
    if (.not. valid) then
      failure_kind = atlas_eval_numerical
      return
    end if
    call numerical_zhao_atlas_jacobian( &
      params, branch, interface_field_v_m, n, z, residual, jacobian, jacobian_ok, failure_kind &
      )
    if (.not. jacobian_ok) return
    scaled_jacobian = jacobian
    do row = 1, n
      row_norm = sqrt(sum(jacobian(row, 1:z_dimension)**2))
      if (.not. ieee_is_finite(row_norm) .or. row_norm <= tiny(1.0_dp)) return
      scaled_jacobian(row, 1:z_dimension) = jacobian(row, 1:z_dimension)/row_norm
    end do
    cofactors = 0.0_dp
    if (n == 2) then
      cofactors(1) = scaled_jacobian(1, 2)*scaled_jacobian(2, 3) - &
                     scaled_jacobian(1, 3)*scaled_jacobian(2, 2)
      cofactors(2) = scaled_jacobian(1, 3)*scaled_jacobian(2, 1) - &
                     scaled_jacobian(1, 1)*scaled_jacobian(2, 3)
      cofactors(3) = scaled_jacobian(1, 1)*scaled_jacobian(2, 2) - &
                     scaled_jacobian(1, 2)*scaled_jacobian(2, 1)
    else
      do column = 1, 4
        minor = 0.0_dp
        minor_column = 0
        do source_column = 1, 4
          if (source_column == column) cycle
          minor_column = minor_column + 1
          minor(1:3, minor_column) = scaled_jacobian(1:3, source_column)
        end do
        cofactors(column) = (-1.0_dp)**(column + 1)*determinant3(minor)
      end do
    end if
    raw_norm = sqrt(sum(cofactors(1:z_dimension)**2))
    if (.not. ieee_is_finite(raw_norm) .or. raw_norm <= tiny(1.0_dp)) return
    row_rank_indicator = raw_norm
    if (.not. ieee_is_finite(row_rank_indicator) .or. row_rank_indicator <= 1.0e-12_dp) return
    tangent(1:z_dimension) = cofactors(1:z_dimension)/raw_norm
    failure_kind = atlas_eval_ok
    success = .true.
  end subroutine zhao_atlas_tangent

  subroutine correct_zhao_atlas_point( &
    params, branch, interface_field_v_m, n, predictor, tangent, options, &
    corrected_z, iterations, success, failure_kind &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: interface_field_v_m
    integer, intent(in) :: n
    real(dp), intent(in) :: predictor(atlas_max_coordinates), tangent(atlas_max_coordinates)
    type(zhao_branch_atlas_options_type), intent(in) :: options
    real(dp), intent(out) :: corrected_z(atlas_max_coordinates)
    integer, intent(out) :: iterations
    logical, intent(out) :: success
    integer, intent(out) :: failure_kind

    real(dp) :: z(atlas_max_coordinates), residual(3), jacobian(3, atlas_max_coordinates)
    real(dp) :: system(atlas_max_coordinates, atlas_max_coordinates)
    real(dp) :: rhs(atlas_max_coordinates), delta(atlas_max_coordinates)
    real(dp) :: trial_z(atlas_max_coordinates), trial_residual(3)
    real(dp) :: norm, trial_norm, arc_residual, trial_arc_residual, backtrack_step
    integer :: iteration, backtrack, z_dimension, trial_kind
    logical :: valid, jacobian_ok, linear_ok, trial_valid

    z_dimension = n + 1
    z = predictor
    corrected_z = predictor
    success = .false.
    failure_kind = atlas_eval_numerical
    iterations = 0
    call evaluate_zhao_atlas_system( &
      params, branch, interface_field_v_m, n, z, residual, valid, failure_kind &
      )
    if (.not. valid) return
    arc_residual = dot_product(tangent(1:z_dimension), z(1:z_dimension) - predictor(1:z_dimension))
    norm = max(maxval(abs(residual(1:n))), abs(arc_residual))
    do iteration = 0, int(options%max_corrector_iterations)
      if (norm <= options%residual_tolerance) then
        success = .true.
        exit
      end if
      if (iteration == int(options%max_corrector_iterations)) exit
      call numerical_zhao_atlas_jacobian( &
        params, branch, interface_field_v_m, n, z, residual, jacobian, &
        jacobian_ok, failure_kind &
        )
      if (.not. jacobian_ok) exit
      system = 0.0_dp
      rhs = 0.0_dp
      system(1:n, 1:z_dimension) = jacobian(1:n, 1:z_dimension)
      system(z_dimension, 1:z_dimension) = tangent(1:z_dimension)
      rhs(1:n) = -residual(1:n)
      rhs(z_dimension) = -arc_residual
      call solve_zhao_atlas_system(system, rhs, z_dimension, delta, linear_ok)
      if (.not. linear_ok) exit
      backtrack_step = 1.0_dp
      do backtrack = 1, root_max_backtracks
        trial_z = z + backtrack_step*delta
        call evaluate_zhao_atlas_system( &
          params, branch, interface_field_v_m, n, trial_z, trial_residual, &
          trial_valid, trial_kind &
          )
        if (trial_valid) then
          trial_arc_residual = dot_product( &
                               tangent(1:z_dimension), trial_z(1:z_dimension) - predictor(1:z_dimension) &
                               )
          trial_norm = max(maxval(abs(trial_residual(1:n))), abs(trial_arc_residual))
          if (trial_norm < norm) then
            z = trial_z
            residual = trial_residual
            arc_residual = trial_arc_residual
            norm = trial_norm
            failure_kind = atlas_eval_ok
            exit
          end if
        end if
        backtrack_step = 0.5_dp*backtrack_step
      end do
      if (backtrack > root_max_backtracks) exit
    end do
    iterations = iteration
    corrected_z = z
    if (.not. success) failure_kind = atlas_eval_numerical
  end subroutine correct_zhao_atlas_point

  subroutine solve_zhao_atlas_system(a_in, b_in, n, x, success)
    real(dp), intent(in) :: a_in(atlas_max_coordinates, atlas_max_coordinates)
    real(dp), intent(in) :: b_in(atlas_max_coordinates)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(atlas_max_coordinates)
    logical, intent(out) :: success

    real(dp) :: a(atlas_max_coordinates, atlas_max_coordinates)
    real(dp) :: b(atlas_max_coordinates), row_scale, factor, pivot_value, tmp
    integer :: i, j, k, pivot

    a = a_in
    b = b_in
    x = 0.0_dp
    success = .false.
    do i = 1, n
      row_scale = maxval(abs(a(i, 1:n)))
      if (.not. ieee_is_finite(row_scale) .or. row_scale <= tiny(1.0_dp)) return
      a(i, 1:n) = a(i, 1:n)/row_scale
      b(i) = b(i)/row_scale
    end do
    do k = 1, n
      pivot = k
      do i = k + 1, n
        if (abs(a(i, k)) > abs(a(pivot, k))) pivot = i
      end do
      if (.not. ieee_is_finite(a(pivot, k)) .or. abs(a(pivot, k)) <= 1.0e-13_dp) return
      if (pivot /= k) then
        do j = k, n
          tmp = a(k, j)
          a(k, j) = a(pivot, j)
          a(pivot, j) = tmp
        end do
        tmp = b(k)
        b(k) = b(pivot)
        b(pivot) = tmp
      end if
      pivot_value = a(k, k)
      do i = k + 1, n
        factor = a(i, k)/pivot_value
        a(i, k:n) = a(i, k:n) - factor*a(k, k:n)
        b(i) = b(i) - factor*b(k)
      end do
    end do
    do i = n, 1, -1
      x(i) = b(i)
      do j = i + 1, n
        x(i) = x(i) - a(i, j)*x(j)
      end do
      x(i) = x(i)/a(i, i)
    end do
    success = all(ieee_is_finite(x(1:n)))
  end subroutine solve_zhao_atlas_system

  subroutine zhao_atlas_point_from_coordinates( &
    params, branch, interface_field_v_m, grid_points, control_length_m, n, z, &
    target_requested, target_column_m2, point, status, message &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: interface_field_v_m
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m
    integer, intent(in) :: n
    real(dp), intent(in) :: z(atlas_max_coordinates)
    logical, intent(in) :: target_requested
    real(dp), intent(in) :: target_column_m2
    type(zhao_branch_atlas_point_type), intent(out) :: point
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(zhao_params_type) :: trial_params
    type(outer_plasma_state_type) :: state
    real(dp) :: residual(3)
    logical :: decoded, valid
    integer :: failure_kind

    point%root = zhao_charge_root_type()
    point%arc_length = 0.0_dp
    point%accepted_step = 0.0_dp
    point%tangent = 0.0_dp
    point%dcolumn_ds = 0.0_dp
    point%row_rank_indicator = 0.0_dp
    point%normalized_jump_from_previous = 0.0_dp
    point%corrector_iterations = 0_i32
    trial_params = params
    trial_params%photoelectron_population_fraction = z(n + 1)
    point%root%branch = branch
    point%root%interface_field_v_m = interface_field_v_m
    point%root%photoelectron_population_fraction = z(n + 1)
    call decode_unknowns( &
      trial_params, branch, z(1:3), point%root%phi0_v, point%root%phi_m_v, &
      point%root%n_swe_inf_m3, decoded &
      )
    if (.not. decoded) then
      status = outer_plasma_no_physical_solution
      message = 'Zhao branch-atlas point left the physical root domain'
      return
    end if
    if (branch == 'A') then
      point%root%interface_side = 'lower'
    else
      point%root%interface_side = 'monotonic'
    end if
    call evaluate_zhao_atlas_system( &
      params, branch, interface_field_v_m, n, z, residual, valid, failure_kind &
      )
    if (.not. valid) then
      status = outer_plasma_numerical_failure
      message = 'Zhao branch-atlas point residual could not be evaluated'
      return
    end if
    point%root%residual_norm = maxval(abs(residual(1:n)))
    point%root%net_current_density_a_m2 = zhao_net_current_density(trial_params, point%root)
    call build_zhao_outer_profile( &
      trial_params, point%root, grid_points, state, status, message, &
      control_length_m=control_length_m &
      )
    if (status /= outer_plasma_ok) return
    point%root%photoelectron_column_per_area = state%photoelectron_column_per_area
    if (target_requested) then
      point%root%photoelectron_column_target_per_area = target_column_m2
      point%root%photoelectron_column_residual_per_area = &
        point%root%photoelectron_column_per_area - target_column_m2
    end if
  end subroutine zhao_atlas_point_from_coordinates

  logical function zhao_atlas_coordinate_endpoint( &
    branch, n, z, tangent, options, reason_code, reason &
    ) result(at_endpoint)
    character(len=1), intent(in) :: branch
    integer, intent(in) :: n
    real(dp), intent(in) :: z(atlas_max_coordinates), tangent(atlas_max_coordinates)
    type(zhao_branch_atlas_options_type), intent(in) :: options
    integer(i32), intent(out) :: reason_code
    character(len=*), intent(out) :: reason

    integer :: density_coordinate, z_dimension

    at_endpoint = .false.
    reason_code = zhao_continuation_reason_none
    reason = 'none'
    z_dimension = n + 1
    density_coordinate = merge(3, 2, branch == 'A')
    if (z(density_coordinate) <= options%log_density_floor .and. tangent(density_coordinate) < 0.0_dp) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'ambient_density_floor_limit'
    else if (z(z_dimension) <= options%eta_min + options%minimum_step .and. &
             tangent(z_dimension) < 0.0_dp) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'eta_lower_search_limit'
    else if (z(z_dimension) >= options%eta_max - options%minimum_step .and. &
             tangent(z_dimension) > 0.0_dp) then
      at_endpoint = .true.
      reason_code = zhao_continuation_reason_search_limit
      reason = 'eta_upper_search_limit'
    end if
  end function zhao_atlas_coordinate_endpoint

  pure real(dp) function determinant3(matrix) result(value)
    real(dp), intent(in) :: matrix(3, 3)

    value = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
            matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
  end function determinant3

  subroutine evaluate_zhao_interface_field(params, root, interface_field_v_m, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    real(dp), intent(out) :: interface_field_v_m
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: phi0_hat, phi_m_hat, density_hat, integral, field_squared, sign_field
    logical :: integral_ok

    interface_field_v_m = 0.0_dp
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. root%n_swe_inf_m3 <= 0.0_dp) then
      message = 'invalid Zhao root for interface-field evaluation'
      return
    end if
    if (root%branch == '0') then
      status = outer_plasma_ok
      return
    end if
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    density_hat = root%n_swe_inf_m3/params%n_phe_ref_m3
    select case (root%branch)
    case ('A')
      call integrate_rho_hat( &
        params, 'A', 'lower', phi_m_hat, phi0_hat, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      field_squared = -2.0_dp*integral
      sign_field = 1.0_dp
    case ('B', 'C')
      call integrate_rho_hat( &
        params, root%branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      field_squared = 2.0_dp*integral
      if (root%branch == 'B') then
        sign_field = 1.0_dp
      else
        sign_field = -1.0_dp
      end if
    case default
      message = 'unknown Zhao branch for interface-field evaluation'
      return
    end select
    if (.not. integral_ok .or. field_squared < -1.0e-10_dp) then
      status = outer_plasma_no_physical_solution
      message = 'Zhao root has no real interface field'
      return
    end if
    interface_field_v_m = sign_field*sqrt(max(0.0_dp, field_squared))* &
                          params%t_phe_ev/params%lambda_d_phe_ref_m
    if (.not. ieee_is_finite(interface_field_v_m)) then
      status = outer_plasma_numerical_failure
      message = 'non-finite Zhao interface field'
      return
    end if
    status = outer_plasma_ok
  end subroutine evaluate_zhao_interface_field

  subroutine solve_outer_plasma_zhao( &
    model, params, interface_field_v_m, grid_points, state, root, status, message, initial_root, &
    control_length_m, allow_transient_bootstrap &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    integer(i32), intent(in) :: grid_points
    type(outer_plasma_state_type), intent(out) :: state
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root
    real(dp), intent(in), optional :: control_length_m
    logical, intent(in), optional :: allow_transient_bootstrap

    logical :: bootstrap_allowed

    bootstrap_allowed = .false.
    if (present(allow_transient_bootstrap)) bootstrap_allowed = allow_transient_bootstrap
    if (present(initial_root)) then
      call solve_zhao_charge_root( &
        model, params, interface_field_v_m, root, status, message, initial_root=initial_root, &
        allow_transient_bootstrap=bootstrap_allowed &
        )
    else
      call solve_zhao_charge_root( &
        model, params, interface_field_v_m, root, status, message, &
        allow_transient_bootstrap=bootstrap_allowed &
        )
    end if
    if (status /= outer_plasma_ok) then
      state = outer_plasma_state_type()
      state%model = 'kinetic_1d'
      state%applicability_status = status
      return
    end if
    if (present(control_length_m)) then
      call build_zhao_outer_profile( &
        params, root, grid_points, state, status, message, control_length_m=control_length_m &
        )
    else
      call build_zhao_outer_profile(params, root, grid_points, state, status, message)
    end if
    if (status == outer_plasma_ok) then
      root%photoelectron_column_per_area = state%photoelectron_column_per_area
      root%photoelectron_column_target_per_area = state%photoelectron_column_target_per_area
      root%photoelectron_column_residual_per_area = state%photoelectron_column_residual_per_area
    end if
  end subroutine solve_outer_plasma_zhao

  subroutine solve_outer_plasma_zhao_column( &
    model, params, interface_field_v_m, grid_points, control_length_m, target_column_m2, &
    state, root, status, message, initial_root, diagnostics &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, target_column_m2
    type(outer_plasma_state_type), intent(out) :: state
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root
    type(zhao_continuation_diagnostics_type), intent(out), optional :: diagnostics

    type(outer_plasma_state_type) :: lower_state, trial_state
    type(zhao_charge_root_type) :: lower_root, upper_root, trial_root, current_root
    real(dp) :: lower_eta, upper_eta, trial_eta, current_eta
    real(dp) :: lower_column, upper_column, trial_column, current_column
    real(dp) :: trial_residual, current_residual, column_tolerance
    integer :: iteration
    integer(i32) :: first_trial_status
    logical :: have_bracket
    character(len=256) :: first_trial_message
    type(zhao_continuation_diagnostics_type) :: first_trial_diagnostics

    state = outer_plasma_state_type()
    root = zhao_charge_root_type()
    if (present(diagnostics)) then
      diagnostics = zhao_continuation_diagnostics_type()
      diagnostics%target_field_v_m = interface_field_v_m
      diagnostics%target_column_m2 = target_column_m2
    end if
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. .not. ieee_is_finite(interface_field_v_m) .or. &
        .not. ieee_is_finite(control_length_m) .or. control_length_m <= 0.0_dp .or. &
        .not. ieee_is_finite(target_column_m2) .or. target_column_m2 < 0.0_dp .or. &
        grid_points < 5_i32) then
      message = 'invalid Zhao photoelectron-column closure request'
      if (present(diagnostics)) then
        diagnostics%reason_code = zhao_continuation_reason_invalid_request
        diagnostics%reason = 'invalid_request'
        diagnostics%solver_stage = 'request'
        diagnostics%underlying_status = status
        diagnostics%return_status = status
      end if
      return
    end if
    if (.not. automatic_zhao_model(model)) then
      message = 'Zhao photoelectron-column closure requires zhao_branch=auto'
      if (present(diagnostics)) then
        diagnostics%reason_code = zhao_continuation_reason_invalid_request
        diagnostics%reason = 'fixed_branch_request'
        diagnostics%solver_stage = 'request'
        diagnostics%underlying_status = status
        diagnostics%return_status = status
      end if
      return
    end if
    if (present(initial_root)) then
      if (initial_root%branch == '0') then
        message = 'Zhao photoelectron-column closure cannot continue from transient branch 0'
        call set_zhao_continuation_origin( &
          diagnostics, 'request', zhao_continuation_reason_invalid_request, 'transient_anchor', &
          status, initial_root%photoelectron_population_fraction, initial_root &
          )
        return
      end if
    end if

    column_tolerance = max( &
                       column_root_relative_tolerance*max(target_column_m2, 1.0_dp), &
                       64.0_dp*epsilon(1.0_dp)*max( &
                       target_column_m2, params%n_phe_ref_m3*control_length_m, 1.0_dp &
                       ) &
                       )
    have_bracket = .false.
    if (target_column_m2 == 0.0_dp .or. .not. present(initial_root)) then
      call evaluate_column_candidate( &
        model, params, interface_field_v_m, grid_points, control_length_m, 0.0_dp, &
        target_column_m2, lower_state, lower_root, status, message &
        )
      if (status /= outer_plasma_ok) then
        if (status /= outer_plasma_numerical_failure) status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'photoelectron-column closure has no eta=0 Zhao anchor: '//trim(message)
        if (present(diagnostics)) then
          diagnostics%reason_code = merge( &
                                    zhao_continuation_reason_numerical_failure, zhao_continuation_reason_physical_endpoint, &
                                    status == outer_plasma_numerical_failure &
                                    )
          if (status == outer_plasma_numerical_failure) then
            diagnostics%reason = 'eta_zero_numerical_failure'
          else
            diagnostics%reason = 'eta_zero_physical_endpoint'
          end if
          diagnostics%solver_stage = 'eta_anchor'
          diagnostics%underlying_status = status
          diagnostics%return_status = status
        end if
        return
      end if
      if (target_column_m2 == 0.0_dp) then
        state = lower_state
        root = lower_root
        call mark_zhao_continuation_converged(diagnostics, root)
        return
      end if
      current_eta = 0.0_dp
      current_column = lower_root%photoelectron_column_per_area
      current_residual = lower_root%photoelectron_column_residual_per_area
      current_root = lower_root
      if (current_residual >= 0.0_dp) then
        status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'positive target column is not bracketed above the eta=0 Zhao state'
        call set_zhao_continuation_origin( &
          diagnostics, 'eta_anchor', zhao_continuation_reason_target_unreachable, &
          'target_below_eta_zero', status, current_eta, current_root &
          )
        return
      end if
    else
      current_eta = initial_root%photoelectron_population_fraction
      if (.not. ieee_is_finite(current_eta) .or. current_eta < 0.0_dp .or. &
          current_eta > column_eta_search_max) then
        status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'previous Zhao eta lies outside the connected column-search interval'
        call set_zhao_continuation_origin( &
          diagnostics, 'eta_anchor', zhao_continuation_reason_search_limit, &
          'previous_eta_outside_search', status, current_eta, initial_root &
          )
        return
      end if
      call continue_column_candidate_in_field( &
        model, params, initial_root%interface_field_v_m, current_eta, initial_root, &
        interface_field_v_m, grid_points, control_length_m, target_column_m2, &
        column_tolerance, trial_state, trial_root, status, message, diagnostics &
        )
      if (status /= outer_plasma_ok) then
        if (status /= outer_plasma_numerical_failure) status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'photoelectron-column continuation could not update the previous Zhao root: '// &
                  trim(message)
        return
      end if
      current_root = trial_root
      current_column = trial_root%photoelectron_column_per_area
      current_residual = trial_root%photoelectron_column_residual_per_area
      if (abs(current_residual) <= column_tolerance) then
        state = trial_state
        root = trial_root
        call mark_zhao_continuation_converged(diagnostics, root)
        return
      end if
    end if

    if (current_residual < 0.0_dp) then
      lower_eta = current_eta
      lower_column = current_column
      lower_root = current_root
      do iteration = 1, column_root_max_iterations
        if (current_eta >= column_eta_search_max) exit
        trial_eta = min( &
                    column_eta_search_max, max(column_eta_initial, 2.0_dp*current_eta) &
                    )
        call continue_column_candidate_in_eta( &
          model, params, interface_field_v_m, current_eta, current_root, trial_eta, &
          grid_points, control_length_m, target_column_m2, column_tolerance, &
          trial_state, trial_root, status, message, diagnostics &
          )
        if (status /= outer_plasma_ok) then
          if (status /= outer_plasma_numerical_failure) status = outer_plasma_no_physical_solution
          state%applicability_status = status
          message = 'photoelectron-column Zhao path ended while increasing eta: '//trim(message)
          return
        end if
        trial_eta = trial_root%photoelectron_population_fraction
        trial_column = trial_root%photoelectron_column_per_area
        trial_residual = trial_root%photoelectron_column_residual_per_area
        if (abs(trial_residual) <= column_tolerance) then
          state = trial_state
          root = trial_root
          call mark_zhao_continuation_converged(diagnostics, root)
          return
        else if (trial_residual > 0.0_dp) then
          upper_eta = trial_eta
          upper_column = trial_column
          upper_root = trial_root
          have_bracket = .true.
          exit
        end if
        current_eta = trial_eta
        current_column = trial_column
        current_residual = trial_residual
        current_root = trial_root
        lower_eta = current_eta
        lower_column = current_column
        lower_root = current_root
      end do
    else
      upper_eta = current_eta
      upper_column = current_column
      upper_root = current_root
      do iteration = 1, column_root_max_iterations
        if (current_eta == 0.0_dp) exit
        if (current_eta <= column_eta_initial) then
          trial_eta = 0.0_dp
        else
          trial_eta = 0.5_dp*current_eta
        end if
        call continue_column_candidate_in_eta( &
          model, params, interface_field_v_m, current_eta, current_root, trial_eta, &
          grid_points, control_length_m, target_column_m2, column_tolerance, &
          trial_state, trial_root, status, message, diagnostics &
          )
        if (status /= outer_plasma_ok) then
          if (status /= outer_plasma_numerical_failure) status = outer_plasma_no_physical_solution
          state%applicability_status = status
          message = 'photoelectron-column Zhao path ended while decreasing eta: '//trim(message)
          return
        end if
        trial_eta = trial_root%photoelectron_population_fraction
        trial_column = trial_root%photoelectron_column_per_area
        trial_residual = trial_root%photoelectron_column_residual_per_area
        if (abs(trial_residual) <= column_tolerance) then
          state = trial_state
          root = trial_root
          call mark_zhao_continuation_converged(diagnostics, root)
          return
        else if (trial_residual < 0.0_dp) then
          lower_eta = trial_eta
          lower_column = trial_column
          lower_root = trial_root
          have_bracket = .true.
          exit
        end if
        current_eta = trial_eta
        current_column = trial_column
        current_residual = trial_residual
        current_root = trial_root
        upper_eta = current_eta
        upper_column = current_column
        upper_root = current_root
      end do
    end if

    if (.not. have_bracket) then
      status = outer_plasma_no_physical_solution
      state%applicability_status = status
      if (current_eta >= column_eta_search_max - 64.0_dp*epsilon(1.0_dp)*column_eta_search_max) then
        write (message, '(a,es12.4,a,es12.4)') &
          'photoelectron-column target ', target_column_m2, &
          ' m^-2 was not bracketed before eta search limit ', column_eta_search_max
        call set_zhao_continuation_origin( &
          diagnostics, 'eta_search', zhao_continuation_reason_search_limit, &
          'eta_upper_search_limit', status, current_eta, current_root &
          )
      else if (current_eta <= 64.0_dp*epsilon(1.0_dp)) then
        message = 'photoelectron-column target lies below the connected Zhao path at eta=0'
        call set_zhao_continuation_origin( &
          diagnostics, 'eta_search', zhao_continuation_reason_target_unreachable, &
          'target_below_eta_zero', status, current_eta, current_root &
          )
      else
        message = 'photoelectron-column eta search reached its iteration limit without a bracket'
        call set_zhao_continuation_origin( &
          diagnostics, 'eta_search', zhao_continuation_reason_search_limit, &
          'eta_iteration_search_limit', status, current_eta, current_root &
          )
      end if
      return
    end if

    do iteration = 1, column_root_max_iterations
      trial_eta = 0.5_dp*(lower_eta + upper_eta)
      call continue_column_candidate_in_eta( &
        model, params, interface_field_v_m, lower_eta, lower_root, trial_eta, &
        grid_points, control_length_m, target_column_m2, column_tolerance, &
        trial_state, trial_root, status, message, diagnostics &
        )
      first_trial_status = status
      first_trial_message = message
      if (present(diagnostics)) first_trial_diagnostics = diagnostics
      if (status /= outer_plasma_ok) then
        call continue_column_candidate_in_eta( &
          model, params, interface_field_v_m, upper_eta, upper_root, trial_eta, &
          grid_points, control_length_m, target_column_m2, column_tolerance, &
          trial_state, trial_root, status, message, diagnostics &
          )
      end if
      if (status /= outer_plasma_ok) then
        if (status == outer_plasma_numerical_failure .or. &
            first_trial_status == outer_plasma_numerical_failure) then
          status = outer_plasma_numerical_failure
          if (first_trial_status == outer_plasma_numerical_failure) then
            message = 'photoelectron-column Zhao bracket has an unresolved numerical failure: '// &
                      trim(first_trial_message)
            if (present(diagnostics)) diagnostics = first_trial_diagnostics
          else
            message = 'photoelectron-column Zhao bracket has an unresolved numerical failure: '// &
                      trim(message)
          end if
        else
          status = outer_plasma_no_physical_solution
          message = 'photoelectron-column Zhao path lost its physical root inside the eta bracket'
        end if
        state%applicability_status = status
        return
      end if
      trial_eta = trial_root%photoelectron_population_fraction
      trial_column = trial_root%photoelectron_column_per_area
      if (trial_column < lower_column - column_tolerance .or. &
          trial_column > upper_column + column_tolerance) then
        status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'photoelectron-column Zhao path is non-monotone inside the eta bracket'
        call set_zhao_continuation_pair( &
          diagnostics, 'eta_bracket', zhao_continuation_reason_nonmonotone_column, &
          'nonmonotone_column', status, params, lower_eta, lower_root, trial_eta, trial_root &
          )
        return
      end if
      trial_residual = trial_root%photoelectron_column_residual_per_area
      if (abs(trial_residual) <= column_tolerance) then
        state = trial_state
        root = trial_root
        call mark_zhao_continuation_converged(diagnostics, root)
        return
      end if
      if (upper_eta - lower_eta <= 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, trial_eta)) then
        status = outer_plasma_numerical_failure
        state = trial_state
        root = trial_root
        state%applicability_status = status
        write (message, '(a,es12.4,a,es12.4)') &
          'photoelectron-column eta bracket collapsed before convergence; target=', &
          target_column_m2, ', residual=', trial_residual
        call set_zhao_continuation_pair( &
          diagnostics, 'eta_bracket', zhao_continuation_reason_numerical_failure, &
          'collapsed_eta_bracket', status, params, lower_eta, lower_root, trial_eta, trial_root &
          )
        return
      end if
      if (trial_residual > 0.0_dp) then
        upper_eta = trial_eta
        upper_column = trial_column
        upper_root = trial_root
      else
        lower_eta = trial_eta
        lower_column = trial_column
        lower_state = trial_state
        lower_root = trial_root
      end if
    end do

    status = outer_plasma_numerical_failure
    state = trial_state
    root = trial_root
    state%applicability_status = status
    write (message, '(a,es12.4,a,es12.4)') &
      'photoelectron-column eta solve reached its iteration limit; target=', &
      target_column_m2, ', residual=', trial_residual
    call set_zhao_continuation_pair( &
      diagnostics, 'eta_bracket', zhao_continuation_reason_search_limit, &
      'eta_root_iteration_limit', status, params, lower_eta, lower_root, trial_eta, trial_root &
      )
  end subroutine solve_outer_plasma_zhao_column

  subroutine evaluate_column_candidate( &
    model, params, interface_field_v_m, grid_points, control_length_m, eta, target_column_m2, &
    state, root, status, message, initial_root &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, eta, target_column_m2
    type(outer_plasma_state_type), intent(out) :: state
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root

    type(zhao_params_type) :: trial_params

    trial_params = params
    trial_params%photoelectron_population_fraction = eta
    if (present(initial_root)) then
      call solve_outer_plasma_zhao( &
        model, trial_params, interface_field_v_m, grid_points, state, root, status, message, &
        initial_root=initial_root, control_length_m=control_length_m &
        )
    else
      call solve_outer_plasma_zhao( &
        model, trial_params, interface_field_v_m, grid_points, state, root, status, message, &
        control_length_m=control_length_m &
        )
    end if
    if (status /= outer_plasma_ok) return
    root%photoelectron_population_fraction = eta
    root%photoelectron_column_target_per_area = target_column_m2
    root%photoelectron_column_residual_per_area = &
      root%photoelectron_column_per_area - target_column_m2
    state%photoelectron_population_fraction = eta
    state%photoelectron_column_target_per_area = target_column_m2
    state%photoelectron_column_residual_per_area = &
      state%photoelectron_column_per_area - target_column_m2
  end subroutine evaluate_column_candidate

  subroutine continue_column_candidate_in_field( &
    model, params, start_field_v_m, eta, start_root, target_field_v_m, grid_points, &
    control_length_m, target_column_m2, column_tolerance, state, root, status, message, diagnostics &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: start_field_v_m, eta, target_field_v_m
    type(zhao_charge_root_type), intent(in) :: start_root
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, target_column_m2, column_tolerance
    type(outer_plasma_state_type), intent(out) :: state
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_continuation_diagnostics_type), intent(inout), optional :: diagnostics

    type(outer_plasma_state_type) :: waypoint_state
    type(zhao_charge_root_type) :: waypoint_root
    logical :: crosses_zero

    state = outer_plasma_state_type()
    root = zhao_charge_root_type()
    status = outer_plasma_no_physical_solution
    message = ''
    if (.not. ieee_is_finite(start_field_v_m) .or. &
        .not. ieee_is_finite(target_field_v_m)) then
      status = outer_plasma_invalid
      message = 'photoelectron-column field continuation has a non-finite endpoint'
      call set_zhao_continuation_origin( &
        diagnostics, 'field_path', zhao_continuation_reason_invalid_request, &
        'nonfinite_field_endpoint', status, eta, start_root &
        )
      if (present(diagnostics)) then
        diagnostics%target_field_v_m = target_field_v_m
        diagnostics%target_eta = eta
        diagnostics%target_column_m2 = target_column_m2
      end if
      return
    end if
    if (.not. automatic_zhao_model(model)) then
      call evaluate_column_candidate( &
        model, params, target_field_v_m, grid_points, control_length_m, eta, &
        target_column_m2, state, root, status, message, start_root &
        )
      if (status == outer_plasma_ok) then
        call validate_column_continuation_step( &
          model, params, eta, start_root, eta, root, column_tolerance, status, message, diagnostics &
          )
      end if
      return
    end if

    crosses_zero = (start_field_v_m < 0.0_dp .and. target_field_v_m > 0.0_dp) .or. &
                   (start_field_v_m > 0.0_dp .and. target_field_v_m < 0.0_dp)
    if (.not. crosses_zero) then
      call continue_column_candidate_field_segment( &
        model, params, start_field_v_m, eta, start_root, target_field_v_m, grid_points, &
        control_length_m, target_column_m2, column_tolerance, state, root, status, message, diagnostics &
        )
      return
    end if

    call continue_column_candidate_field_segment( &
      model, params, start_field_v_m, eta, start_root, 0.0_dp, grid_points, &
      control_length_m, target_column_m2, column_tolerance, waypoint_state, waypoint_root, &
      status, message, diagnostics &
      )
    if (status /= outer_plasma_ok) then
      if (present(diagnostics)) diagnostics%target_field_v_m = target_field_v_m
      message = 'photoelectron-column field continuation could not reach E=0: '//trim(message)
      return
    end if
    call continue_column_candidate_field_segment( &
      model, params, 0.0_dp, eta, waypoint_root, target_field_v_m, grid_points, &
      control_length_m, target_column_m2, column_tolerance, state, root, status, message, diagnostics &
      )
  end subroutine continue_column_candidate_in_field

  subroutine continue_column_candidate_field_segment( &
    model, params, start_field_v_m, eta, start_root, target_field_v_m, grid_points, &
    control_length_m, target_column_m2, column_tolerance, state, root, status, message, diagnostics &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: start_field_v_m, eta, target_field_v_m
    type(zhao_charge_root_type), intent(in) :: start_root
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, target_column_m2, column_tolerance
    type(outer_plasma_state_type), intent(out) :: state
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_continuation_diagnostics_type), intent(inout), optional :: diagnostics

    type(outer_plasma_state_type) :: candidate_state
    type(zhao_charge_root_type) :: candidate_root, current_root
    real(dp) :: current_field, trial_field, remaining, direction
    real(dp) :: total_distance, maximum_step, minimum_step, step, scale
    integer :: attempt
    logical :: saw_numerical_failure, candidate_solved
    character(len=256) :: last_message

    state = outer_plasma_state_type()
    root = zhao_charge_root_type()
    status = outer_plasma_no_physical_solution
    message = ''
    last_message = ''
    saw_numerical_failure = .false.
    if (present(diagnostics)) then
      diagnostics%target_field_v_m = target_field_v_m
      diagnostics%target_eta = eta
      diagnostics%target_column_m2 = target_column_m2
    end if
    current_field = start_field_v_m
    current_root = start_root
    total_distance = abs(target_field_v_m - start_field_v_m)
    scale = max( &
            abs(start_field_v_m), abs(target_field_v_m), &
            params%t_phe_ev/params%lambda_d_phe_ref_m, tiny(1.0_dp) &
            )
    minimum_step = max( &
                   column_continuation_min_relative_step*scale, &
                   64.0_dp*epsilon(1.0_dp)*scale &
                   )
    if (total_distance == 0.0_dp) then
      call evaluate_column_candidate( &
        model, params, target_field_v_m, grid_points, control_length_m, eta, &
        target_column_m2, state, root, status, message, current_root &
        )
      if (status == outer_plasma_ok) then
        call validate_column_continuation_step( &
          model, params, eta, current_root, eta, root, column_tolerance, status, message, diagnostics &
          )
      end if
      if (status /= outer_plasma_ok) then
        if (present(diagnostics)) then
          if (diagnostics%reason_code == zhao_continuation_reason_none) then
            call set_zhao_continuation_origin( &
              diagnostics, 'field_continuation', &
              merge(zhao_continuation_reason_numerical_failure, &
                    zhao_continuation_reason_physical_endpoint, &
                    status == outer_plasma_numerical_failure), &
              'unchanged_field_update_failed', status, eta, current_root &
              )
          end if
          diagnostics%attempted_step = 0.0_dp
        end if
        if (status /= outer_plasma_numerical_failure) status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'photoelectron-column root update failed at unchanged field: '//trim(message)
      end if
      return
    end if

    maximum_step = total_distance/real(column_continuation_nominal_segments, dp)
    step = maximum_step
    do attempt = 1, column_continuation_max_attempts
      remaining = target_field_v_m - current_field
      direction = merge(1.0_dp, -1.0_dp, remaining > 0.0_dp)
      if (step >= abs(remaining)) then
        trial_field = target_field_v_m
      else
        trial_field = current_field + direction*step
      end if
      if (trial_field == current_field) exit

      call evaluate_column_candidate( &
        model, params, trial_field, grid_points, control_length_m, eta, &
        target_column_m2, candidate_state, candidate_root, status, message, current_root &
        )
      candidate_solved = status == outer_plasma_ok
      if (status == outer_plasma_ok) then
        call validate_column_continuation_step( &
          model, params, eta, current_root, eta, candidate_root, &
          column_tolerance, status, message, diagnostics &
          )
      end if
      if (status == outer_plasma_ok) then
        current_field = trial_field
        current_root = candidate_root
        if (current_field == target_field_v_m) then
          state = candidate_state
          root = candidate_root
          return
        end if
        step = min(maximum_step, 2.0_dp*step)
      else
        if (present(diagnostics)) then
          if (.not. candidate_solved) then
            call set_zhao_continuation_origin( &
              diagnostics, 'field_continuation', &
              merge(zhao_continuation_reason_numerical_failure, &
                    zhao_continuation_reason_physical_endpoint, &
                    status == outer_plasma_numerical_failure), &
              'field_candidate_failed', status, eta, current_root &
              )
          end if
          diagnostics%attempt = int(attempt, i32)
          diagnostics%attempted_step = abs(trial_field - current_field)
          diagnostics%target_field_v_m = target_field_v_m
          diagnostics%target_eta = eta
          diagnostics%target_column_m2 = target_column_m2
        end if
        saw_numerical_failure = saw_numerical_failure .or. status == outer_plasma_numerical_failure
        last_message = message
        if (step <= minimum_step) exit
        step = 0.5_dp*step
      end if
    end do

    if (saw_numerical_failure) then
      status = outer_plasma_numerical_failure
    else
      status = outer_plasma_no_physical_solution
    end if
    state%applicability_status = status
    message = 'bounded Zhao field continuation exhausted step-halving: '//trim(last_message)
    if (present(diagnostics)) then
      diagnostics%saw_numerical_failure = saw_numerical_failure
      diagnostics%return_status = status
    end if
  end subroutine continue_column_candidate_field_segment

  subroutine continue_column_candidate_in_eta( &
    model, params, interface_field_v_m, start_eta, start_root, target_eta, grid_points, &
    control_length_m, target_column_m2, column_tolerance, state, root, status, message, diagnostics &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m, start_eta, target_eta
    type(zhao_charge_root_type), intent(in) :: start_root
    integer(i32), intent(in) :: grid_points
    real(dp), intent(in) :: control_length_m, target_column_m2, column_tolerance
    type(outer_plasma_state_type), intent(out) :: state
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_continuation_diagnostics_type), intent(inout), optional :: diagnostics

    type(outer_plasma_state_type) :: candidate_state
    type(zhao_charge_root_type) :: candidate_root, current_root
    real(dp) :: current_eta, trial_eta, remaining, direction
    real(dp) :: total_distance, maximum_step, minimum_step, step, scale
    integer :: attempt
    logical :: saw_numerical_failure, candidate_solved
    character(len=256) :: last_message

    state = outer_plasma_state_type()
    root = zhao_charge_root_type()
    status = outer_plasma_no_physical_solution
    message = ''
    last_message = ''
    saw_numerical_failure = .false.
    if (present(diagnostics)) then
      diagnostics%target_field_v_m = interface_field_v_m
      diagnostics%target_eta = target_eta
      diagnostics%target_column_m2 = target_column_m2
    end if
    if (.not. ieee_is_finite(start_eta) .or. .not. ieee_is_finite(target_eta) .or. &
        start_eta < 0.0_dp .or. target_eta < 0.0_dp) then
      message = 'photoelectron-column eta continuation has an invalid endpoint'
      call set_zhao_continuation_origin( &
        diagnostics, 'eta_continuation', zhao_continuation_reason_invalid_request, &
        'invalid_eta_endpoint', status, start_eta, start_root &
        )
      return
    end if
    if (.not. automatic_zhao_model(model)) then
      call evaluate_column_candidate( &
        model, params, interface_field_v_m, grid_points, control_length_m, target_eta, &
        target_column_m2, state, root, status, message, start_root &
        )
      if (status == outer_plasma_ok) then
        call validate_column_continuation_step( &
          model, params, start_eta, start_root, target_eta, root, &
          column_tolerance, status, message, diagnostics &
          )
      end if
      if (status /= outer_plasma_ok .and. present(diagnostics)) then
        if (diagnostics%reason_code == zhao_continuation_reason_none) then
          call set_zhao_continuation_origin( &
            diagnostics, 'eta_continuation', &
            merge(zhao_continuation_reason_numerical_failure, &
                  zhao_continuation_reason_physical_endpoint, &
                  status == outer_plasma_numerical_failure), &
            'fixed_branch_candidate_failed', status, start_eta, start_root &
            )
        end if
        diagnostics%attempted_step = abs(target_eta - start_eta)
      end if
      return
    end if

    current_eta = start_eta
    current_root = start_root
    total_distance = abs(target_eta - start_eta)
    scale = max(abs(start_eta), abs(target_eta), 1.0_dp)
    minimum_step = max( &
                   column_continuation_min_relative_step*scale, &
                   64.0_dp*epsilon(1.0_dp)*scale &
                   )
    if (total_distance == 0.0_dp) then
      call evaluate_column_candidate( &
        model, params, interface_field_v_m, grid_points, control_length_m, target_eta, &
        target_column_m2, state, root, status, message, current_root &
        )
      if (status == outer_plasma_ok) then
        call validate_column_continuation_step( &
          model, params, current_eta, current_root, target_eta, root, &
          column_tolerance, status, message, diagnostics &
          )
      end if
      if (status /= outer_plasma_ok) then
        if (present(diagnostics)) then
          if (diagnostics%reason_code == zhao_continuation_reason_none) then
            call set_zhao_continuation_origin( &
              diagnostics, 'eta_continuation', &
              merge(zhao_continuation_reason_numerical_failure, &
                    zhao_continuation_reason_physical_endpoint, &
                    status == outer_plasma_numerical_failure), &
              'unchanged_eta_update_failed', status, current_eta, current_root &
              )
          end if
          diagnostics%attempted_step = 0.0_dp
        end if
        if (status /= outer_plasma_numerical_failure) status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'photoelectron-column root update failed at unchanged eta: '//trim(message)
      end if
      return
    end if

    maximum_step = total_distance/real(column_continuation_nominal_segments, dp)
    step = maximum_step
    do attempt = 1, column_continuation_max_attempts
      remaining = target_eta - current_eta
      direction = merge(1.0_dp, -1.0_dp, remaining > 0.0_dp)
      if (step >= abs(remaining)) then
        trial_eta = target_eta
      else
        trial_eta = current_eta + direction*step
      end if
      if (trial_eta == current_eta) exit

      call evaluate_column_candidate( &
        model, params, interface_field_v_m, grid_points, control_length_m, trial_eta, &
        target_column_m2, candidate_state, candidate_root, status, message, current_root &
        )
      candidate_solved = status == outer_plasma_ok
      if (status == outer_plasma_ok) then
        call validate_column_continuation_step( &
          model, params, current_eta, current_root, trial_eta, candidate_root, &
          column_tolerance, status, message, diagnostics &
          )
      end if
      if (status == outer_plasma_ok) then
        if (column_target_is_bracketed(current_root, candidate_root, column_tolerance)) then
          state = candidate_state
          root = candidate_root
          return
        end if
        current_eta = trial_eta
        current_root = candidate_root
        if (current_eta == target_eta) then
          state = candidate_state
          root = candidate_root
          return
        end if
        step = min(maximum_step, 2.0_dp*step)
      else
        if (present(diagnostics)) then
          if (.not. candidate_solved) then
            call set_zhao_continuation_origin( &
              diagnostics, 'eta_continuation', &
              merge(zhao_continuation_reason_numerical_failure, &
                    zhao_continuation_reason_physical_endpoint, &
                    status == outer_plasma_numerical_failure), &
              'eta_candidate_failed', status, current_eta, current_root &
              )
          end if
          diagnostics%attempt = int(attempt, i32)
          diagnostics%attempted_step = abs(trial_eta - current_eta)
          diagnostics%target_field_v_m = interface_field_v_m
          diagnostics%target_eta = target_eta
          diagnostics%target_column_m2 = target_column_m2
        end if
        saw_numerical_failure = saw_numerical_failure .or. status == outer_plasma_numerical_failure
        last_message = message
        if (step <= minimum_step) exit
        step = 0.5_dp*step
      end if
    end do

    if (saw_numerical_failure) then
      status = outer_plasma_numerical_failure
    else
      status = outer_plasma_no_physical_solution
    end if
    state%applicability_status = status
    message = 'bounded Zhao eta continuation exhausted step-halving: '//trim(last_message)
    if (present(diagnostics)) then
      diagnostics%saw_numerical_failure = saw_numerical_failure
      diagnostics%return_status = status
    end if
  end subroutine continue_column_candidate_in_eta

  logical function column_target_is_bracketed(from_root, to_root, column_tolerance) result(bracketed)
    type(zhao_charge_root_type), intent(in) :: from_root, to_root
    real(dp), intent(in) :: column_tolerance

    real(dp) :: from_residual, to_residual

    from_residual = from_root%photoelectron_column_residual_per_area
    to_residual = to_root%photoelectron_column_residual_per_area
    bracketed = abs(to_residual) <= column_tolerance .or. &
                (from_residual < 0.0_dp .and. to_residual > 0.0_dp) .or. &
                (from_residual > 0.0_dp .and. to_residual < 0.0_dp)
  end function column_target_is_bracketed

  subroutine validate_column_continuation_step( &
    model, params, from_eta, from_root, to_eta, to_root, column_tolerance, status, message, diagnostics &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: from_eta, to_eta, column_tolerance
    type(zhao_charge_root_type), intent(in) :: from_root, to_root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_continuation_diagnostics_type), intent(inout), optional :: diagnostics

    real(dp) :: from_column, to_column, monotonic_tolerance

    status = outer_plasma_ok
    message = ''
    from_column = from_root%photoelectron_column_per_area
    to_column = to_root%photoelectron_column_per_area
    if (.not. ieee_is_finite(from_eta) .or. .not. ieee_is_finite(to_eta) .or. &
        .not. zhao_charge_root_is_finite(from_root) .or. &
        .not. zhao_charge_root_is_finite(to_root)) then
      status = outer_plasma_numerical_failure
      message = 'photoelectron-column Zhao continuation contains a non-finite state'
      call set_zhao_continuation_pair( &
        diagnostics, 'step_validation', zhao_continuation_reason_numerical_failure, &
        'nonfinite_state', status, params, from_eta, from_root, to_eta, to_root &
        )
      return
    end if

    if (from_root%branch == '0' .or. to_root%branch == '0') then
      status = outer_plasma_no_physical_solution
      message = 'photoelectron-column continuation cannot use transient branch 0'
      call set_zhao_continuation_pair( &
        diagnostics, 'step_validation', zhao_continuation_reason_physical_endpoint, &
        'transient_branch', status, params, from_eta, from_root, to_eta, to_root &
        )
      return
    end if
    if (automatic_zhao_model(model)) then
      if (from_root%branch /= to_root%branch) then
        if (.not. zhao_branch_transition_is_continuous(params, from_root, to_root)) then
          status = outer_plasma_no_physical_solution
          write (message, '(a,a,a,a)') &
            'automatic Zhao continuation rejected a disconnected ', from_root%branch, &
            '->', to_root%branch
          call set_zhao_continuation_pair( &
            diagnostics, 'step_validation', zhao_continuation_reason_disconnected_branch, &
            'disconnected_branch', status, params, from_eta, from_root, to_eta, to_root &
            )
          return
        end if
      else if (.not. zhao_same_branch_step_is_bounded(params, from_root, to_root)) then
        status = outer_plasma_no_physical_solution
        message = 'automatic Zhao continuation step moved too far on one root branch'
        call set_zhao_continuation_pair( &
          diagnostics, 'step_validation', zhao_continuation_reason_guard_rejected, &
          'same_branch_jump_guard', status, params, from_eta, from_root, to_eta, to_root &
          )
        return
      end if
    end if

    monotonic_tolerance = max( &
                          column_tolerance, &
                          64.0_dp*epsilon(1.0_dp)*max(abs(from_column), abs(to_column), 1.0_dp) &
                          )
    if (to_eta > from_eta .and. to_column < from_column - monotonic_tolerance) then
      status = outer_plasma_no_physical_solution
      message = 'photoelectron-column Zhao path decreases while eta increases'
      call set_zhao_continuation_pair( &
        diagnostics, 'step_validation', zhao_continuation_reason_nonmonotone_column, &
        'column_decreases_with_eta', status, params, from_eta, from_root, to_eta, to_root &
        )
    else if (to_eta < from_eta .and. to_column > from_column + monotonic_tolerance) then
      status = outer_plasma_no_physical_solution
      message = 'photoelectron-column Zhao path increases while eta decreases'
      call set_zhao_continuation_pair( &
        diagnostics, 'step_validation', zhao_continuation_reason_nonmonotone_column, &
        'column_increases_as_eta_decreases', status, params, from_eta, from_root, to_eta, to_root &
        )
    end if
  end subroutine validate_column_continuation_step

  logical function automatic_zhao_model(model) result(is_automatic)
    character(len=*), intent(in) :: model

    character(len=16) :: normalized_model

    normalized_model = trim(lower_ascii(model))
    is_automatic = normalized_model == 'auto' .or. normalized_model == 'zhao_auto'
  end function automatic_zhao_model

  logical function zhao_same_branch_step_is_bounded(params, from_root, to_root) result(bounded)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: from_root, to_root

    real(dp) :: potential_change, density_change, jump

    bounded = .false.
    if (from_root%branch /= to_root%branch) return
    if (index('ABC', from_root%branch) == 0) return
    if (from_root%n_swe_inf_m3 <= 0.0_dp .or. to_root%n_swe_inf_m3 <= 0.0_dp) return

    call zhao_root_jump_metrics(params, from_root, to_root, potential_change, density_change, jump)
    bounded = ieee_is_finite(potential_change) .and. ieee_is_finite(density_change) .and. &
              jump <= branch_same_root_step_limit
  end function zhao_same_branch_step_is_bounded

  logical function zhao_branch_transition_is_continuous(params, from_root, to_root) result(continuous)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: from_root, to_root

    character(len=2) :: transition
    real(dp) :: potential_scale, density_scale
    logical :: density_is_continuous

    continuous = .false.
    if (index('ABC', from_root%branch) == 0 .or. index('ABC', to_root%branch) == 0) return
    continuous = from_root%branch == to_root%branch
    if (continuous) return

    potential_scale = max(params%t_phe_ev, tiny(1.0_dp))
    density_scale = max( &
                    abs(from_root%n_swe_inf_m3), abs(to_root%n_swe_inf_m3), &
                    params%n_phe_ref_m3, 1.0_dp &
                    )
    density_is_continuous = &
      abs(from_root%n_swe_inf_m3 - to_root%n_swe_inf_m3) <= &
      branch_continuity_tolerance*density_scale
    if (.not. density_is_continuous) return

    transition = from_root%branch//to_root%branch
    select case (transition)
    case ('AB')
      continuous = abs(from_root%phi_m_v) <= branch_degeneracy_tolerance*potential_scale .and. &
                   abs(from_root%phi0_v - to_root%phi0_v) <= &
                   branch_continuity_tolerance*potential_scale
    case ('BA')
      continuous = abs(to_root%phi_m_v) <= branch_degeneracy_tolerance*potential_scale .and. &
                   abs(from_root%phi0_v - to_root%phi0_v) <= &
                   branch_continuity_tolerance*potential_scale
    case ('BC', 'CB')
      continuous = max(abs(from_root%phi0_v), abs(to_root%phi0_v)) <= &
        branch_degeneracy_tolerance*potential_scale .and. &
        abs(from_root%phi0_v - to_root%phi0_v) <= &
        branch_continuity_tolerance*potential_scale
    case ('AC')
      continuous = max( &
                   abs(from_root%phi0_v), abs(from_root%phi_m_v), abs(to_root%phi0_v) &
                   ) <= branch_degeneracy_tolerance*potential_scale
    case ('CA')
      continuous = max( &
                   abs(from_root%phi0_v), abs(to_root%phi0_v), abs(to_root%phi_m_v) &
                   ) <= branch_degeneracy_tolerance*potential_scale
    case default
      continuous = .false.
    end select
  end function zhao_branch_transition_is_continuous

  subroutine solve_zhao_charge_root( &
    model, params, interface_field_v_m, root, status, message, initial_root, allow_transient_bootstrap &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root
    logical, intent(in), optional :: allow_transient_bootstrap

    character(len=16) :: normalized_model
    character(len=1) :: order(3), candidate
    type(zhao_charge_root_type) :: trial_root
    real(dp) :: field_scale, target_field_hat, degenerate_density_hat
    integer :: candidate_count, candidate_index, warm_index
    character(len=1) :: first_candidate
    logical :: bootstrap_allowed, saw_numerical_failure
    character(len=256) :: numerical_message

    root = zhao_charge_root_type()
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. .not. ieee_is_finite(interface_field_v_m)) then
      message = 'invalid Zhao charge-driven parameters'
      return
    end if
    root%photoelectron_population_fraction = params%photoelectron_population_fraction
    bootstrap_allowed = .false.
    if (present(allow_transient_bootstrap)) bootstrap_allowed = allow_transient_bootstrap
    field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
    target_field_hat = interface_field_v_m/field_scale
    normalized_model = trim(lower_ascii(model))
    call branch_order(normalized_model, params, target_field_hat, order, candidate_count, status, message)
    if (status /= outer_plasma_ok) return

    if (abs(target_field_hat) <= zero_field_tolerance_hat) then
      if (normalized_model /= 'auto' .and. normalized_model /= 'zhao_auto' .and. &
          normalized_model /= 'b' .and. normalized_model /= 'zhao_b') then
        status = outer_plasma_no_physical_solution
        message = 'requested Zhao branch does not contain the zero-field degenerate state'
        return
      end if
      degenerate_density_hat = ( &
                               2.0_dp*params%n_swi_inf_m3/params%n_phe_ref_m3 - &
                               params%photoelectron_population_fraction*params%n_phe0_m3/params%n_phe_ref_m3 &
                               )/(1.0_dp + erf(params%u))
      root%phi0_v = 0.0_dp
      root%phi_m_v = 0.0_dp
      root%interface_field_v_m = 0.0_dp
      root%residual_norm = 0.0_dp
      root%nonlinear_iterations = 0_i32
      if (degenerate_density_hat > 0.0_dp .and. ieee_is_finite(degenerate_density_hat)) then
        root%branch = 'B'
        root%interface_side = 'monotonic'
        root%n_swe_inf_m3 = degenerate_density_hat*params%n_phe_ref_m3
        message = 'zero-field degenerate Zhao-B state'
      else if ((normalized_model == 'auto' .or. normalized_model == 'zhao_auto') .and. bootstrap_allowed) then
        root%branch = '0'
        root%interface_side = 'bootstrap'
        root%n_swe_inf_m3 = params%n_swi_inf_m3
        message = 'zero-field transient bootstrap; photoelectron outer population is not yet formed'
      else
        status = outer_plasma_no_physical_solution
        message = 'zero-field Zhao state has no positive ambient-electron density and bootstrap is disabled'
        return
      end if
      root%net_current_density_a_m2 = zhao_net_current_density(params, root)
      status = outer_plasma_ok
      return
    end if

    if (present(initial_root)) then
      if ((normalized_model == 'auto' .or. normalized_model == 'zhao_auto') .and. &
          index('ABC', initial_root%branch) > 0) then
        do warm_index = 1, candidate_count
          if (order(warm_index) /= initial_root%branch) cycle
          first_candidate = order(1)
          order(1) = order(warm_index)
          order(warm_index) = first_candidate
          exit
        end do
      end if
    end if
    saw_numerical_failure = .false.
    numerical_message = ''
    do candidate_index = 1, candidate_count
      candidate = order(candidate_index)
      if (present(initial_root)) then
        call solve_one_branch( &
          params, candidate, interface_field_v_m, trial_root, status, message, initial_root=initial_root &
          )
      else
        call solve_one_branch(params, candidate, interface_field_v_m, trial_root, status, message)
      end if
      if (status == outer_plasma_ok) then
        root = trial_root
        return
      end if
      if (status == outer_plasma_numerical_failure) then
        saw_numerical_failure = .true.
        numerical_message = message
      end if
    end do
    if (saw_numerical_failure) then
      status = outer_plasma_numerical_failure
      message = 'charge-driven Zhao branch search has an unresolved numerical failure: '//trim(numerical_message)
    else
      status = outer_plasma_no_physical_solution
      message = 'no charge-driven Zhao branch satisfies the prescribed interface field'
    end if
  end subroutine solve_zhao_charge_root

  subroutine branch_order(model, params, target_field_hat, order, count, status, message)
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: target_field_hat
    character(len=1), intent(out) :: order(3)
    integer, intent(out) :: count
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    order = ' '
    count = 0
    status = outer_plasma_ok
    message = ''
    select case (trim(model))
    case ('zhao_a', 'a')
      order(1) = 'A'
      count = 1
    case ('zhao_b', 'b')
      order(1) = 'B'
      count = 1
    case ('zhao_c', 'c')
      order(1) = 'C'
      count = 1
    case ('zhao_auto', 'auto')
      if (target_field_hat > 0.0_dp) then
        order = ['A', 'B', 'C']
      else
        order = ['C', 'A', 'B']
      end if
      if (params%alpha_rad*180.0_dp/pi < 20.0_dp .and. target_field_hat > 0.0_dp) then
        order = ['B', 'A', 'C']
      end if
      count = 3
    case default
      status = outer_plasma_invalid
      message = 'unknown Zhao charge-driven model'
    end select
  end subroutine branch_order

  subroutine solve_one_branch( &
    params, branch, interface_field_v_m, root, status, message, initial_root &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: interface_field_v_m
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root

    real(dp) :: target_field_hat, field_scale
    real(dp) :: guesses(3, 8), y(3), best_y(3), norm, best_norm
    integer :: guess_count, guess_index, iterations, best_iterations
    logical :: guess_valid, success

    root = zhao_charge_root_type()
    root%branch = branch
    root%interface_field_v_m = interface_field_v_m
    root%photoelectron_population_fraction = params%photoelectron_population_fraction
    status = outer_plasma_no_physical_solution
    message = ''
    field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
    target_field_hat = interface_field_v_m/field_scale
    if ((branch == 'A' .or. branch == 'B') .and. target_field_hat <= 0.0_dp) then
      message = 'Zhao A/B lower interface requires a positive normal field'
      return
    end if
    if (branch == 'C' .and. target_field_hat >= 0.0_dp) then
      message = 'Zhao C interface requires a negative normal field'
      return
    end if

    call make_branch_guesses(params, branch, guesses, guess_count)
    if (present(initial_root)) then
      if (initial_root%branch == branch) then
        call encode_unknowns( &
          params, branch, initial_root%phi0_v, initial_root%phi_m_v, &
          initial_root%n_swe_inf_m3, y, guess_valid &
          )
        if (guess_valid) then
          guesses(:, 2:guess_count + 1) = guesses(:, 1:guess_count)
          guesses(:, 1) = y
          guess_count = guess_count + 1
        end if
      end if
    end if

    best_norm = huge(1.0_dp)
    best_y = guesses(:, 1)
    best_iterations = 0
    do guess_index = 1, guess_count
      call newton_branch( &
        params, branch, target_field_hat, guesses(:, guess_index), y, norm, iterations, success &
        )
      if (norm < best_norm) then
        best_norm = norm
        best_y = y
        best_iterations = iterations
      end if
      if (success) exit
    end do
    if (.not. success) then
      status = outer_plasma_numerical_failure
      message = 'charge-driven Zhao root Newton solve did not converge'
      return
    end if

    call decode_unknowns( &
      params, branch, y, root%phi0_v, root%phi_m_v, root%n_swe_inf_m3, guess_valid &
      )
    if (.not. guess_valid) then
      status = outer_plasma_numerical_failure
      message = 'charge-driven Zhao root decoded to an invalid state'
      return
    end if
    if (branch == 'A') then
      root%interface_side = 'lower'
    else
      root%interface_side = 'monotonic'
    end if
    root%residual_norm = norm
    root%nonlinear_iterations = int(iterations, i32)
    root%net_current_density_a_m2 = zhao_net_current_density(params, root)
    status = outer_plasma_ok
    message = ''
  end subroutine solve_one_branch

  subroutine make_branch_guesses(params, branch, guesses, count)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(out) :: guesses(3, 8)
    integer, intent(out) :: count

    real(dp) :: phi0, phi_m, density
    logical :: valid

    guesses = 0.0_dp
    select case (branch)
    case ('A')
      count = 5
      call encode_unknowns(params, branch, 3.6_dp, -0.5_dp, 8.2e6_dp, guesses(:, 1), valid)
      call encode_unknowns(params, branch, 2.8_dp, -0.3_dp, 8.0e6_dp, guesses(:, 2), valid)
      call encode_unknowns(params, branch, 4.5_dp, -0.8_dp, 8.4e6_dp, guesses(:, 3), valid)
      phi0 = 1.3_dp*params%t_phe_ev
      phi_m = -0.3_dp*params%t_phe_ev
      density = max(0.9_dp*params%n_swi_inf_m3, tiny(1.0_dp))
      call encode_unknowns(params, branch, phi0, phi_m, density, guesses(:, 4), valid)
      call encode_unknowns( &
        params, branch, 0.8_dp*params%t_phe_ev, -0.1_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 5), valid &
        )
    case ('B')
      count = 7
      call encode_unknowns(params, branch, 1.3_dp, 1.3_dp, 7.0e6_dp, guesses(:, 1), valid)
      call encode_unknowns(params, branch, 0.8_dp, 0.8_dp, 6.5e6_dp, guesses(:, 2), valid)
      call encode_unknowns(params, branch, 2.0_dp, 2.0_dp, 7.8e6_dp, guesses(:, 3), valid)
      call encode_unknowns( &
        params, branch, 0.6_dp*params%t_phe_ev, 0.6_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 4), valid &
        )
      call encode_unknowns( &
        params, branch, 0.2_dp*params%t_phe_ev, 0.2_dp*params%t_phe_ev, &
        0.8_dp*params%n_swi_inf_m3, guesses(:, 5), valid &
        )
      density = max( &
                (2.0_dp*params%n_swi_inf_m3 - &
                 params%photoelectron_population_fraction*params%n_phe0_m3)/(1.0_dp + erf(params%u)), &
                0.5_dp*params%n_swi_inf_m3 &
                )
      call encode_unknowns( &
        params, branch, 0.02_dp*params%t_phe_ev, 0.02_dp*params%t_phe_ev, &
        density, guesses(:, 6), valid &
        )
      call encode_unknowns( &
        params, branch, 0.002_dp*params%t_phe_ev, 0.002_dp*params%t_phe_ev, &
        density, guesses(:, 7), valid &
        )
    case ('C')
      count = 7
      call encode_unknowns(params, branch, -0.5_dp, -0.5_dp, 6.0e6_dp, guesses(:, 1), valid)
      call encode_unknowns(params, branch, -2.0_dp, -2.0_dp, 7.0e6_dp, guesses(:, 2), valid)
      call encode_unknowns(params, branch, -5.0_dp, -5.0_dp, 8.0e6_dp, guesses(:, 3), valid)
      call encode_unknowns(params, branch, -10.0_dp, -10.0_dp, 8.2e6_dp, guesses(:, 4), valid)
      call encode_unknowns(params, branch, -15.0_dp, -15.0_dp, 8.5e6_dp, guesses(:, 5), valid)
      call encode_unknowns( &
        params, branch, -1.0_dp*params%t_phe_ev, -1.0_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 6), valid &
        )
      call encode_unknowns( &
        params, branch, -3.0_dp*params%t_phe_ev, -3.0_dp*params%t_phe_ev, &
        0.9_dp*params%n_swi_inf_m3, guesses(:, 7), valid &
        )
    case default
      count = 0
    end select
  end subroutine make_branch_guesses

  subroutine encode_unknowns(params, branch, phi0_v, phi_m_v, density_m3, y, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, density_m3
    real(dp), intent(out) :: y(3)
    logical, intent(out) :: valid

    y = 0.0_dp
    valid = density_m3 > 0.0_dp .and. params%n_phe_ref_m3 > 0.0_dp .and. params%t_phe_ev > 0.0_dp
    if (.not. valid) return
    select case (branch)
    case ('A')
      valid = phi0_v > 0.0_dp .and. phi_m_v < 0.0_dp
      if (.not. valid) return
      y(1) = log(phi0_v/params%t_phe_ev)
      y(2) = log(-phi_m_v/params%t_phe_ev)
      y(3) = log(density_m3/params%n_phe_ref_m3)
    case ('B')
      valid = phi0_v > 0.0_dp
      if (.not. valid) return
      y(1) = log(phi0_v/params%t_phe_ev)
      y(2) = log(density_m3/params%n_phe_ref_m3)
    case ('C')
      valid = phi0_v < 0.0_dp
      if (.not. valid) return
      y(1) = log(-phi0_v/params%t_phe_ev)
      y(2) = log(density_m3/params%n_phe_ref_m3)
    case default
      valid = .false.
    end select
    valid = valid .and. all(ieee_is_finite(y))
  end subroutine encode_unknowns

  subroutine decode_unknowns(params, branch, y, phi0_v, phi_m_v, density_m3, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: y(3)
    real(dp), intent(out) :: phi0_v, phi_m_v, density_m3
    logical, intent(out) :: valid

    phi0_v = 0.0_dp
    phi_m_v = 0.0_dp
    density_m3 = 0.0_dp
    valid = all(ieee_is_finite(y))
    if (.not. valid) return
    if (y(1) < -50.0_dp .or. y(1) > log(200.0_dp)) then
      valid = .false.
      return
    end if
    select case (branch)
    case ('A')
      if (y(2) < -50.0_dp .or. y(2) > log(200.0_dp) .or. &
          y(3) < -30.0_dp .or. y(3) > log(1.0e6_dp)) then
        valid = .false.
        return
      end if
      phi0_v = params%t_phe_ev*exp(y(1))
      phi_m_v = -params%t_phe_ev*exp(y(2))
      density_m3 = params%n_phe_ref_m3*exp(y(3))
    case ('B')
      if (y(2) < -30.0_dp .or. y(2) > log(1.0e6_dp)) then
        valid = .false.
        return
      end if
      phi0_v = params%t_phe_ev*exp(y(1))
      phi_m_v = phi0_v
      density_m3 = params%n_phe_ref_m3*exp(y(2))
    case ('C')
      if (y(2) < -30.0_dp .or. y(2) > log(1.0e6_dp)) then
        valid = .false.
        return
      end if
      phi0_v = -params%t_phe_ev*exp(y(1))
      phi_m_v = phi0_v
      density_m3 = params%n_phe_ref_m3*exp(y(2))
    case default
      valid = .false.
    end select
    valid = valid .and. all(ieee_is_finite([phi0_v, phi_m_v, density_m3]))
  end subroutine decode_unknowns

  subroutine newton_branch(params, branch, target_field_hat, y0, y_out, final_norm, iterations, success)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y0(3)
    real(dp), intent(out) :: y_out(3), final_norm
    integer, intent(out) :: iterations
    logical, intent(out) :: success

    real(dp) :: y(3), f(3), jac(3, 3), delta(3), trial(3), trial_f(3)
    real(dp) :: norm, trial_norm, step
    integer :: n, iteration, backtrack
    logical :: valid, jacobian_ok, linear_ok, trial_valid

    if (branch == 'A') then
      n = 3
    else
      n = 2
    end if
    y = y0
    call eval_charge_residual(params, branch, target_field_hat, y, f, valid)
    if (.not. valid) then
      y_out = y
      final_norm = huge(1.0_dp)
      iterations = 0
      success = .false.
      return
    end if
    norm = maxval(abs(f(1:n)))
    success = .false.
    do iteration = 0, root_max_iterations
      if (norm <= root_tolerance) then
        success = .true.
        exit
      end if
      if (iteration == root_max_iterations) exit
      call numerical_jacobian(params, branch, target_field_hat, y, f, n, jac, jacobian_ok)
      if (.not. jacobian_ok) exit
      call solve_small_system(jac, -f, n, delta, linear_ok)
      if (.not. linear_ok) exit
      step = 1.0_dp
      do backtrack = 1, root_max_backtracks
        trial = y + step*delta
        call eval_charge_residual(params, branch, target_field_hat, trial, trial_f, trial_valid)
        if (trial_valid) then
          trial_norm = maxval(abs(trial_f(1:n)))
          if (trial_norm < norm) then
            y = trial
            f = trial_f
            norm = trial_norm
            exit
          end if
        end if
        step = 0.5_dp*step
      end do
      if (backtrack > root_max_backtracks) exit
    end do
    y_out = y
    final_norm = norm
    iterations = iteration
  end subroutine newton_branch

  subroutine numerical_jacobian(params, branch, target_field_hat, y, f0, n, jac, success)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y(3), f0(3)
    integer, intent(in) :: n
    real(dp), intent(out) :: jac(3, 3)
    logical, intent(out) :: success

    real(dp) :: yp(3), ym(3), fp(3), fm(3), h
    integer :: column
    logical :: plus_valid, minus_valid

    jac = 0.0_dp
    success = .true.
    do column = 1, n
      h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(y(column)))
      yp = y
      ym = y
      yp(column) = yp(column) + h
      ym(column) = ym(column) - h
      call eval_charge_residual(params, branch, target_field_hat, yp, fp, plus_valid)
      call eval_charge_residual(params, branch, target_field_hat, ym, fm, minus_valid)
      if (plus_valid .and. minus_valid) then
        jac(1:n, column) = (fp(1:n) - fm(1:n))/(2.0_dp*h)
      else if (plus_valid) then
        jac(1:n, column) = (fp(1:n) - f0(1:n))/h
      else if (minus_valid) then
        jac(1:n, column) = (f0(1:n) - fm(1:n))/h
      else
        success = .false.
        return
      end if
    end do
    success = all(ieee_is_finite(jac(1:n, 1:n)))
  end subroutine numerical_jacobian

  subroutine solve_small_system(a_in, b_in, n, x, success)
    real(dp), intent(in) :: a_in(3, 3), b_in(3)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(3)
    logical, intent(out) :: success

    real(dp) :: a(3, 3), b(3), factor, pivot_value, tmp
    integer :: i, j, k, pivot

    a = a_in
    b = b_in
    x = 0.0_dp
    success = .false.
    do k = 1, n
      pivot = k
      do i = k + 1, n
        if (abs(a(i, k)) > abs(a(pivot, k))) pivot = i
      end do
      if (.not. ieee_is_finite(a(pivot, k)) .or. abs(a(pivot, k)) <= 1.0e-14_dp) return
      if (pivot /= k) then
        do j = k, n
          tmp = a(k, j)
          a(k, j) = a(pivot, j)
          a(pivot, j) = tmp
        end do
        tmp = b(k)
        b(k) = b(pivot)
        b(pivot) = tmp
      end if
      pivot_value = a(k, k)
      do i = k + 1, n
        factor = a(i, k)/pivot_value
        a(i, k:n) = a(i, k:n) - factor*a(k, k:n)
        b(i) = b(i) - factor*b(k)
      end do
    end do
    do i = n, 1, -1
      x(i) = b(i)
      do j = i + 1, n
        x(i) = x(i) - a(i, j)*x(j)
      end do
      x(i) = x(i)/a(i, i)
    end do
    success = all(ieee_is_finite(x(1:n)))
  end subroutine solve_small_system

  subroutine eval_charge_residual(params, branch, target_field_hat, y, residual, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y(3)
    real(dp), intent(out) :: residual(3)
    logical, intent(out) :: valid

    real(dp) :: phi0_v, phi_m_v, density_m3, phi0_hat, phi_m_hat, density_hat
    real(dp) :: raw(3), integral, field_squared, field_scale
    real(dp) :: x3(3), x2(2)
    logical :: integral_ok

    residual = 0.0_dp
    call decode_unknowns(params, branch, y, phi0_v, phi_m_v, density_m3, valid)
    if (.not. valid) return
    phi0_hat = phi0_v/params%t_phe_ev
    phi_m_hat = phi_m_v/params%t_phe_ev
    density_hat = density_m3/params%n_phe_ref_m3
    if (.not. ion_accessible(params, max(phi0_hat, 0.0_dp))) then
      valid = .false.
      return
    end if

    select case (branch)
    case ('A')
      x3 = [phi0_v, phi_m_v, density_m3]
      call zhao_residuals_type_a(params, x3, raw)
      call integrate_rho_hat( &
        params, branch, 'lower', phi_m_hat, phi0_hat, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) then
        valid = .false.
        return
      end if
      field_squared = -2.0_dp*integral
      if (field_squared < -1.0e-10_dp) then
        valid = .false.
        return
      end if
      field_scale = max(1.0_dp, target_field_hat*target_field_hat)
      residual(1) = raw(1)/params%n_phe_ref_m3
      residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_scale
      residual(3) = raw(3)
    case ('B')
      x2 = [phi0_v, density_m3]
      call zhao_residuals_type_b(params, x2, raw(1:2))
      call integrate_rho_hat( &
        params, branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) then
        valid = .false.
        return
      end if
      field_squared = 2.0_dp*integral
      if (field_squared < -1.0e-10_dp) then
        valid = .false.
        return
      end if
      field_scale = max(1.0_dp, target_field_hat*target_field_hat)
      residual(1) = raw(1)/params%n_phe_ref_m3
      residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_scale
    case ('C')
      x2 = [phi0_v, density_m3]
      call zhao_residuals_type_c(params, x2, raw(1:2))
      call integrate_rho_hat( &
        params, branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) then
        valid = .false.
        return
      end if
      field_squared = 2.0_dp*integral
      if (field_squared < -1.0e-10_dp) then
        valid = .false.
        return
      end if
      field_scale = max(1.0_dp, target_field_hat*target_field_hat)
      residual(1) = raw(1)/params%n_phe_ref_m3
      residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_scale
    case default
      valid = .false.
      return
    end select
    valid = all(ieee_is_finite(residual))
  end subroutine eval_charge_residual

  subroutine integrate_rho_hat( &
    params, branch, side, lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat, integral, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat
    real(dp), intent(out) :: integral
    logical, intent(out) :: success

    real(dp) :: t, phi_hat, jacobian, rho_hat, summand, weight, h
    integer :: point

    integral = 0.0_dp
    success = .false.
    if (rho_quadrature_panels < 2 .or. mod(rho_quadrature_panels, 2) /= 0) return
    if (.not. all(ieee_is_finite([ &
                                 lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat &
                                 ])) .or. density_hat <= 0.0_dp) return
    h = 1.0_dp/real(rho_quadrature_panels, dp)
    do point = 0, rho_quadrature_panels
      t = real(point, dp)*h
      phi_hat = lower_phi_hat + (upper_phi_hat - lower_phi_hat)*sin(0.5_dp*pi*t)**2
      jacobian = (upper_phi_hat - lower_phi_hat)*0.5_dp*pi*sin(pi*t)
      if (.not. ion_accessible(params, phi_hat)) return
      call evaluate_zhao_rho_hat( &
        params, branch, side, phi_hat, phi0_hat, phi_m_hat, density_hat, rho_hat &
        )
      if (.not. ieee_is_finite(rho_hat)) return
      summand = rho_hat*jacobian
      if (point == 0 .or. point == rho_quadrature_panels) then
        weight = 1.0_dp
      else if (mod(point, 2) == 0) then
        weight = 2.0_dp
      else
        weight = 4.0_dp
      end if
      integral = integral + weight*summand
    end do
    integral = integral*h/3.0_dp
    success = ieee_is_finite(integral)
  end subroutine integrate_rho_hat

  subroutine build_zhao_outer_profile( &
    params, root, grid_points, state, status, message, control_length_m &
    )
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    integer(i32), intent(in) :: grid_points
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: control_length_m

    real(dp), allocatable :: phi_hat(:), field_hat(:), rho_hat(:)
    real(dp) :: phi0_hat, phi_m_hat, density_hat, electron_current, ion_current, photo_current
    integer :: point
    logical :: success

    state = outer_plasma_state_type()
    state%model = 'kinetic_1d'
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. grid_points < 5_i32 .or. &
        .not. ieee_is_finite(root%interface_field_v_m) .or. root%n_swe_inf_m3 <= 0.0_dp) then
      state%applicability_status = status
      message = 'invalid Zhao profile request'
      return
    end if
    if (present(control_length_m)) then
      if (.not. ieee_is_finite(control_length_m) .or. control_length_m <= 0.0_dp) then
        state%applicability_status = status
        message = 'Zhao profile control length must be finite and positive'
        return
      end if
    end if

    allocate (phi_hat(grid_points), field_hat(grid_points), rho_hat(grid_points))
    allocate (state%z(grid_points), state%potential(grid_points), &
              state%field(grid_points), state%charge_density(grid_points))
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    density_hat = root%n_swe_inf_m3/params%n_phe_ref_m3
    if (root%branch == '0' .or. &
        (abs(phi0_hat) <= zero_field_tolerance_hat .and. &
         abs(root%interface_field_v_m) <= &
         zero_field_tolerance_hat*params%t_phe_ev/params%lambda_d_phe_ref_m)) then
      phi_hat = 0.0_dp
      field_hat = 0.0_dp
      rho_hat = 0.0_dp
      do point = 1, grid_points
        state%z(point) = params%lambda_d_phe_ref_m*real(point - 1, dp)/real(grid_points - 1, dp)
      end do
    else
      select case (root%branch)
      case ('A')
        call build_type_a_profile( &
          params, phi0_hat, phi_m_hat, density_hat, phi_hat, field_hat, rho_hat, state%z, success &
          )
      case ('B', 'C')
        call build_monotonic_profile( &
          params, root%branch, phi0_hat, density_hat, phi_hat, field_hat, rho_hat, state%z, success &
          )
      case default
        success = .false.
      end select
      if (.not. success) then
        status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'failed to reconstruct the charge-driven Zhao profile'
        return
      end if
    end if
    if (present(control_length_m)) then
      call resample_profile_to_control_length( &
        control_length_m, state%z, phi_hat, field_hat, rho_hat, success &
        )
      if (.not. success) then
        status = outer_plasma_numerical_failure
        state%applicability_status = status
        message = 'failed to resample the Zhao profile onto its finite control volume'
        return
      end if
      if (root%branch == 'A') then
        point = minloc(phi_hat, dim=1)
        if (point > 1 .and. point < grid_points) phi_hat(point) = phi_m_hat
      end if
      call evaluate_resampled_zhao_state(params, root, phi_hat, field_hat, rho_hat, success)
      if (.not. success) then
        status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'finite Zhao control volume does not contain a valid resolved profile'
        return
      end if
    end if

    state%profile_n = grid_points
    state%kinetic_closure = 'zhao_charge_driven'
    state%zhao_branch = root%branch
    state%interface_z = 0.0_dp
    state%interface_potential = root%phi0_v
    state%zhao_phi0 = root%phi0_v
    state%zhao_phi_minimum = root%phi_m_v
    state%zhao_electron_density_infinity = root%n_swe_inf_m3
    state%photoelectron_population_fraction = params%photoelectron_population_fraction
    state%infinity_potential = 0.0_dp
    state%debye_length = params%lambda_d_phe_ref_m
    state%interface_field = root%interface_field_v_m
    state%nonlinear_iterations = root%nonlinear_iterations
    state%nonlinear_residual = root%residual_norm
    state%potential = params%t_phe_ev*phi_hat
    state%field = (params%t_phe_ev/params%lambda_d_phe_ref_m)*field_hat
    state%field(1) = root%interface_field_v_m
    state%charge_density = qe*params%n_phe_ref_m3*rho_hat
    state%integrated_charge_per_area = trapz(state%z, state%charge_density)
    call evaluate_profile_photoelectron_column(params, root, state, success)
    if (.not. success) then
      status = outer_plasma_numerical_failure
      state%applicability_status = status
      message = 'failed to integrate the resolved Zhao photoelectron column'
      return
    end if
    call zhao_current_components(params, root, electron_current, ion_current, photo_current)
    state%electron_current_density = electron_current
    state%ion_current_density = ion_current
    state%photoelectron_current_density = photo_current
    state%total_current_density = electron_current + ion_current + photo_current
    state%ready = .true.
    state%applicability_status = outer_plasma_ok
    status = outer_plasma_ok
  end subroutine build_zhao_outer_profile

  subroutine resample_profile_to_control_length( &
    control_length_m, z, phi_hat, field_hat, rho_hat, success &
    )
    real(dp), intent(in) :: control_length_m
    real(dp), intent(inout) :: z(:), phi_hat(:), field_hat(:), rho_hat(:)
    logical, intent(out) :: success

    real(dp), allocatable :: source_z(:), source_phi(:), source_field(:), source_rho(:)
    real(dp) :: query_z, weight
    integer :: point, source_point, n

    success = .false.
    n = size(z)
    if (n < 2 .or. size(phi_hat) /= n .or. size(field_hat) /= n .or. size(rho_hat) /= n .or. &
        .not. ieee_is_finite(control_length_m) .or. control_length_m <= 0.0_dp) return
    if (.not. all(ieee_is_finite(z)) .or. .not. all(z(2:) > z(:n - 1))) return

    allocate (source_z(n), source_phi(n), source_field(n), source_rho(n))
    source_z = z
    source_phi = phi_hat
    source_field = field_hat
    source_rho = rho_hat
    source_point = 1
    do point = 1, n
      query_z = control_length_m*real(point - 1, dp)/real(n - 1, dp)
      z(point) = query_z
      if (query_z > source_z(n)) then
        ! The native reconstruction ends at |phi_hat|=1e-4.  Beyond it the
        ! finite control volume uses the analytic infinity asymptote.
        phi_hat(point) = 0.0_dp
        field_hat(point) = 0.0_dp
        rho_hat(point) = 0.0_dp
        cycle
      end if
      do while (source_point < n - 1 .and. source_z(source_point + 1) < query_z)
        source_point = source_point + 1
      end do
      weight = (query_z - source_z(source_point))/ &
               (source_z(source_point + 1) - source_z(source_point))
      phi_hat(point) = (1.0_dp - weight)*source_phi(source_point) + &
                       weight*source_phi(source_point + 1)
      field_hat(point) = (1.0_dp - weight)*source_field(source_point) + &
                         weight*source_field(source_point + 1)
      rho_hat(point) = (1.0_dp - weight)*source_rho(source_point) + &
                       weight*source_rho(source_point + 1)
    end do
    z(n) = control_length_m
    success = all(ieee_is_finite(z)) .and. all(ieee_is_finite(phi_hat)) .and. &
              all(ieee_is_finite(field_hat)) .and. all(ieee_is_finite(rho_hat))
  end subroutine resample_profile_to_control_length

  subroutine evaluate_resampled_zhao_state(params, root, phi_hat, field_hat, rho_hat, success)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    real(dp), intent(in) :: phi_hat(:)
    real(dp), intent(out) :: field_hat(:), rho_hat(:)
    logical, intent(out) :: success

    real(dp) :: phi0_hat, phi_m_hat, density_hat, integral, field_squared
    integer :: point, minimum_point, n
    character(len=16) :: side
    logical :: integral_ok

    success = .false.
    n = size(phi_hat)
    if (n < 3 .or. size(field_hat) /= n .or. size(rho_hat) /= n .or. &
        .not. all(ieee_is_finite(phi_hat))) return
    field_hat = 0.0_dp
    rho_hat = 0.0_dp
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    density_hat = root%n_swe_inf_m3/params%n_phe_ref_m3
    if (root%branch == '0' .or. &
        (abs(phi0_hat) <= zero_field_tolerance_hat .and. &
         abs(root%interface_field_v_m) <= &
         zero_field_tolerance_hat*params%t_phe_ev/params%lambda_d_phe_ref_m)) then
      success = maxval(abs(phi_hat)) <= zero_field_tolerance_hat
      return
    end if

    minimum_point = minloc(phi_hat, dim=1)
    if (root%branch == 'A' .and. (minimum_point <= 1 .or. minimum_point >= n)) return
    do point = 1, n
      if (abs(phi_hat(point)) <= 64.0_dp*epsilon(1.0_dp)) cycle
      select case (root%branch)
      case ('A')
        if (point <= minimum_point) then
          side = 'lower'
        else
          side = 'upper'
        end if
        call integrate_rho_hat( &
          params, 'A', side, phi_m_hat, phi_hat(point), phi0_hat, phi_m_hat, density_hat, &
          integral, integral_ok &
          )
        field_squared = -2.0_dp*integral
        if (.not. integral_ok .or. field_squared < -1.0e-9_dp) return
        if (point <= minimum_point) then
          field_hat(point) = sqrt(max(0.0_dp, field_squared))
        else
          field_hat(point) = -sqrt(max(0.0_dp, field_squared))
        end if
      case ('B', 'C')
        side = 'monotonic'
        call integrate_rho_hat( &
          params, root%branch, side, phi_hat(point), 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
          integral, integral_ok &
          )
        field_squared = 2.0_dp*integral
        if (.not. integral_ok .or. field_squared < -1.0e-9_dp) return
        if (root%branch == 'B') then
          field_hat(point) = sqrt(max(0.0_dp, field_squared))
        else
          field_hat(point) = -sqrt(max(0.0_dp, field_squared))
        end if
      case default
        return
      end select
      if (.not. ion_accessible(params, phi_hat(point))) return
      call evaluate_zhao_rho_hat( &
        params, root%branch, side, phi_hat(point), phi0_hat, phi_m_hat, density_hat, rho_hat(point) &
        )
      if (.not. ieee_is_finite(rho_hat(point))) return
    end do
    success = all(ieee_is_finite(field_hat)) .and. all(ieee_is_finite(rho_hat))
  end subroutine evaluate_resampled_zhao_state

  subroutine evaluate_profile_photoelectron_column(params, root, state, success)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    type(outer_plasma_state_type), intent(inout) :: state
    logical, intent(out) :: success

    real(dp), allocatable :: photoelectron_density(:)
    real(dp) :: phi_hat, phi0_hat, phi_m_hat, density_hat
    real(dp) :: n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat
    integer :: point, minimum_point
    character(len=16) :: side

    success = .false.
    state%photoelectron_column_per_area = 0.0_dp
    if (root%branch == '0') then
      success = .true.
      return
    end if
    if (.not. allocated(state%z) .or. .not. allocated(state%potential) .or. &
        size(state%z) /= size(state%potential) .or. size(state%z) < 2) return

    allocate (photoelectron_density(size(state%z)))
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    density_hat = root%n_swe_inf_m3/params%n_phe_ref_m3
    minimum_point = minloc(state%potential, dim=1)
    do point = 1, size(state%z)
      if (root%branch == 'A') then
        if (point <= minimum_point) then
          side = 'lower'
        else
          side = 'upper'
        end if
      else
        side = 'monotonic'
      end if
      phi_hat = state%potential(point)/params%t_phe_ev
      if (.not. ion_accessible(params, phi_hat)) return
      call evaluate_zhao_density_hat( &
        params, root%branch, side, phi_hat, phi0_hat, phi_m_hat, density_hat, &
        n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat &
        )
      photoelectron_density(point) = params%n_phe_ref_m3*(n_phe_f_hat + n_phe_c_hat)
      if (.not. ieee_is_finite(photoelectron_density(point)) .or. &
          photoelectron_density(point) < 0.0_dp) return
    end do
    state%photoelectron_column_per_area = trapz(state%z, photoelectron_density)
    success = ieee_is_finite(state%photoelectron_column_per_area) .and. &
              state%photoelectron_column_per_area >= 0.0_dp
  end subroutine evaluate_profile_photoelectron_column

  subroutine build_monotonic_profile( &
    params, branch, phi0_hat, density_hat, phi_hat, field_hat, rho_hat, z, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_hat, density_hat
    real(dp), intent(out) :: phi_hat(:), field_hat(:), rho_hat(:), z(:)
    logical, intent(out) :: success

    real(dp) :: x, phi_end_hat, integral, field_squared, phi_m_hat
    integer :: point, n
    logical :: integral_ok

    success = .false.
    n = size(phi_hat)
    if (size(field_hat) /= n .or. size(rho_hat) /= n .or. size(z) /= n) return
    phi_m_hat = phi0_hat
    if (branch == 'B') then
      phi_end_hat = profile_phi_end_hat
      if (phi_end_hat >= phi0_hat) phi_end_hat = 0.5_dp*phi0_hat
    else if (branch == 'C') then
      phi_end_hat = -profile_phi_end_hat
      if (phi_end_hat <= phi0_hat) phi_end_hat = 0.5_dp*phi0_hat
    else
      return
    end if
    do point = 1, n
      x = real(point - 1, dp)/real(n - 1, dp)
      phi_hat(point) = phi0_hat + (phi_end_hat - phi0_hat)*sin(0.5_dp*pi*x)**2
      call integrate_rho_hat( &
        params, branch, 'monotonic', phi_hat(point), 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) return
      field_squared = 2.0_dp*integral
      if (field_squared < -1.0e-9_dp) return
      if (branch == 'B') then
        field_hat(point) = sqrt(max(0.0_dp, field_squared))
      else
        field_hat(point) = -sqrt(max(0.0_dp, field_squared))
      end if
      if (.not. ion_accessible(params, phi_hat(point))) return
      call evaluate_zhao_rho_hat( &
        params, branch, 'monotonic', phi_hat(point), phi0_hat, phi_m_hat, density_hat, rho_hat(point) &
        )
    end do
    call integrate_profile_z(params%lambda_d_phe_ref_m, phi_hat, field_hat, z, success)
  end subroutine build_monotonic_profile

  subroutine build_type_a_profile( &
    params, phi0_hat, phi_m_hat, density_hat, phi_hat, field_hat, rho_hat, z, success &
    )
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: phi0_hat, phi_m_hat, density_hat
    real(dp), intent(out) :: phi_hat(:), field_hat(:), rho_hat(:), z(:)
    logical, intent(out) :: success

    real(dp) :: x, phi_end_hat, integral, field_squared
    integer :: point, upper_point, n, n_lower, n_upper
    logical :: integral_ok

    success = .false.
    n = size(phi_hat)
    if (size(field_hat) /= n .or. size(rho_hat) /= n .or. size(z) /= n) return
    n_lower = max(3, (n + 1)/3)
    n_upper = n - n_lower + 1
    if (n_upper < 3) return
    phi_end_hat = -profile_phi_end_hat
    if (phi_end_hat <= phi_m_hat) phi_end_hat = 0.5_dp*phi_m_hat

    do point = 1, n_lower
      x = real(point - 1, dp)/real(n_lower - 1, dp)
      phi_hat(point) = phi_m_hat + (phi0_hat - phi_m_hat)*cos(0.5_dp*pi*x)**2
      call integrate_rho_hat( &
        params, 'A', 'lower', phi_m_hat, phi_hat(point), phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) return
      field_squared = -2.0_dp*integral
      if (field_squared < -1.0e-9_dp) return
      field_hat(point) = sqrt(max(0.0_dp, field_squared))
      call evaluate_zhao_rho_hat( &
        params, 'A', 'lower', phi_hat(point), phi0_hat, phi_m_hat, density_hat, rho_hat(point) &
        )
    end do
    do upper_point = 2, n_upper
      point = n_lower + upper_point - 1
      x = real(upper_point - 1, dp)/real(n_upper - 1, dp)
      phi_hat(point) = phi_m_hat + (phi_end_hat - phi_m_hat)*sin(0.5_dp*pi*x)**2
      call integrate_rho_hat( &
        params, 'A', 'upper', phi_m_hat, phi_hat(point), phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) return
      field_squared = -2.0_dp*integral
      if (field_squared < -1.0e-9_dp) return
      field_hat(point) = -sqrt(max(0.0_dp, field_squared))
      call evaluate_zhao_rho_hat( &
        params, 'A', 'upper', phi_hat(point), phi0_hat, phi_m_hat, density_hat, rho_hat(point) &
        )
    end do
    call integrate_profile_z(params%lambda_d_phe_ref_m, phi_hat, field_hat, z, success)
  end subroutine build_type_a_profile

  subroutine integrate_profile_z(lambda_d, phi_hat, field_hat, z, success)
    real(dp), intent(in) :: lambda_d, phi_hat(:), field_hat(:)
    real(dp), intent(out) :: z(:)
    logical, intent(out) :: success

    real(dp) :: field_average, delta_phi
    integer :: point

    success = .false.
    if (lambda_d <= 0.0_dp .or. size(phi_hat) /= size(field_hat) .or. size(z) /= size(phi_hat)) return
    z(1) = 0.0_dp
    do point = 2, size(phi_hat)
      delta_phi = abs(phi_hat(point) - phi_hat(point - 1))
      field_average = 0.5_dp*(abs(field_hat(point)) + abs(field_hat(point - 1)))
      if (field_average <= tiny(1.0_dp) .and. delta_phi > tiny(1.0_dp)) return
      z(point) = z(point - 1) + lambda_d*delta_phi/max(field_average, tiny(1.0_dp))
    end do
    success = all(ieee_is_finite(z)) .and. all(z(2:) > z(:size(z) - 1))
  end subroutine integrate_profile_z

  real(dp) function zhao_net_current_density(params, root) result(current_density)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root

    real(dp) :: electron_current, ion_current, photo_current

    call zhao_current_components(params, root, electron_current, ion_current, photo_current)
    current_density = electron_current + ion_current + photo_current
  end function zhao_net_current_density

  subroutine zhao_current_components(params, root, electron_current, ion_current, photo_current)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    real(dp), intent(out) :: electron_current, ion_current, photo_current

    real(dp) :: a_swe, electron_term, ion_term, photo_term, scale
    real(dp) :: phi0_hat, phi_m_hat
    character(len=1) :: branch

    electron_current = 0.0_dp
    ion_current = 0.0_dp
    photo_current = 0.0_dp
    if (.not. valid_zhao_params(params) .or. root%n_swe_inf_m3 <= 0.0_dp) return
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    branch = root%branch
    if (branch == '0') branch = 'B'
    ion_term = params%n_swi_inf_m3*sqrt( &
               2.0_dp*pi*params%t_swe_ev/params%t_phe_ev*params%m_e_kg/params%m_i_kg &
               )*params%mach
    select case (branch)
    case ('A')
      a_swe = sqrt(max(0.0_dp, -phi_m_hat/params%tau)) - params%u
      electron_term = swe_free_current_term(params, root%n_swe_inf_m3, a_swe)
      photo_term = params%n_phe0_m3*exp(phi_m_hat - phi0_hat)
    case ('B')
      electron_term = swe_free_current_term(params, root%n_swe_inf_m3, -params%u)
      photo_term = params%n_phe0_m3*exp(-phi0_hat)
    case ('C')
      a_swe = sqrt(max(0.0_dp, -phi0_hat/params%tau)) - params%u
      electron_term = swe_free_current_term(params, root%n_swe_inf_m3, a_swe)
      photo_term = params%n_phe0_m3
    case default
      return
    end select
    scale = qe*params%v_phe_th_mps/(2.0_dp*sqrt(pi))
    electron_current = -scale*electron_term
    ion_current = scale*ion_term
    ! This diagnostic is the tracked full surface source.  Partial outer-cloud
    ! occupancy changes density closure, not the emitted current ledger.
    photo_current = scale*photo_term
  end subroutine zhao_current_components

  pure logical function ion_accessible(params, phi_hat) result(accessible)
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: phi_hat

    accessible = params%tau > 0.0_dp .and. params%mach > 0.0_dp .and. &
                 1.0_dp - 2.0_dp*phi_hat/(params%tau*params%mach*params%mach) > 0.0_dp
  end function ion_accessible

  pure logical function valid_zhao_params(params) result(valid)
    type(zhao_params_type), intent(in) :: params

    valid = all(ieee_is_finite([ &
                               params%n_swi_inf_m3, params%n_phe_ref_m3, params%n_phe0_m3, &
                               params%photoelectron_population_fraction, &
                               params%t_swe_ev, params%t_phe_ev, params%m_i_kg, params%m_e_kg, &
                               params%v_phe_th_mps, params%mach, params%u, params%tau, &
                               params%lambda_d_phe_ref_m &
                               ])) .and. &
            params%n_swi_inf_m3 > 0.0_dp .and. params%n_phe_ref_m3 > 0.0_dp .and. &
            params%photoelectron_population_fraction >= 0.0_dp .and. &
            params%t_swe_ev > 0.0_dp .and. params%t_phe_ev > 0.0_dp .and. &
            params%m_i_kg > 0.0_dp .and. params%m_e_kg > 0.0_dp .and. &
            params%v_phe_th_mps > 0.0_dp .and. params%mach > 0.0_dp .and. &
            params%tau > 0.0_dp .and. params%lambda_d_phe_ref_m > 0.0_dp
  end function valid_zhao_params

  pure real(dp) function trapz(x, y) result(integral)
    real(dp), intent(in) :: x(:), y(:)
    integer :: point

    integral = 0.0_dp
    do point = 1, size(x) - 1
      integral = integral + 0.5_dp*(y(point) + y(point + 1))*(x(point + 1) - x(point))
    end do
  end function trapz

end module bem_outer_plasma_zhao
