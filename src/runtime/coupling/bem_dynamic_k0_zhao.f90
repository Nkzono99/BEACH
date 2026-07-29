!> 実測界面エネルギー分布とZhao profileを結ぶ強光電子用mean-charge更新。
module bem_dynamic_k0_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi
  use bem_dynamic_k0_mean, only: &
    dynamic_k0_ok, dynamic_k0_invalid, dynamic_k0_numerical_failure, dynamic_k0_step_type
  use bem_outer_plasma_kinetic, only: &
    kinetic_outer_plasma_options_type, continue_outer_plasma_zhao_connected_root, &
    materialize_outer_plasma_zhao_root, zhao_charge_root_type, &
    zhao_charge_root_barrier_energy, zhao_connected_path_certificate_type
  use bem_outer_plasma_types, only: &
    outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid, &
    outer_plasma_no_physical_solution, outer_plasma_numerical_failure
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: dynamic_zhao_nonmonotone_barrier = 3_i32
  integer(i32), parameter, public :: dynamic_zhao_no_physical_root = 4_i32
  real(dp), parameter :: frozen_interface_shift_limit = 0.25_dp
  real(dp), parameter :: frozen_source_log_shift_limit = 0.25_dp
  real(dp), parameter :: frozen_field_shift_limit = 0.25_dp
  real(dp), parameter :: frozen_barrier_shift_limit = 0.25_dp
  integer(i32), parameter :: marginal_bisection_limit = 96_i32

  !> charge-weighted empirical interface normal-energy distribution.
  type, public :: measured_interface_energy_distribution_type
    logical :: ready = .false.
    integer(i32) :: group_count = 0_i32
    real(dp) :: total_charge_c = 0.0_dp
    real(dp), allocatable :: energy_j(:)
    real(dp), allocatable :: group_charge_c(:)
    real(dp), allocatable :: cumulative_charge_c(:)
  end type measured_interface_energy_distribution_type

  public :: build_measured_interface_energy_distribution
  public :: advance_dynamic_k0_zhao
  public :: zhao_profile_barrier_energy
  public :: measured_sample_escape_fraction

contains

  !> rank 0へ集めた各rayの(E_n, |dq|)を降順groupへ正規化する。
  subroutine build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
    real(dp), intent(in) :: sample_energy_j(:), sample_charge_c(:)
    type(measured_interface_energy_distribution_type), intent(out) :: distribution
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp), allocatable :: sorted_energy(:), sorted_charge(:)
    real(dp), allocatable :: grouped_energy(:), grouped_charge(:), cumulative(:)
    integer :: i, group_count

    distribution = measured_interface_energy_distribution_type()
    status = dynamic_k0_invalid
    message = ''
    if (size(sample_energy_j) /= size(sample_charge_c)) then
      message = 'measured Zhao energy and charge sample sizes differ'
      return
    end if
    if (size(sample_energy_j) < 1) then
      message = 'measured Zhao energy distribution has no interface samples'
      return
    end if
    if (any(.not. ieee_is_finite(sample_energy_j)) .or. &
        any(.not. ieee_is_finite(sample_charge_c)) .or. &
        any(sample_energy_j < 0.0_dp) .or. any(sample_charge_c <= 0.0_dp)) then
      message = 'measured Zhao energy distribution contains an invalid sample'
      return
    end if

    allocate (sorted_energy(size(sample_energy_j)), sorted_charge(size(sample_charge_c)))
    sorted_energy = sample_energy_j
    sorted_charge = sample_charge_c
    call stable_sort_energy_descending(sorted_energy, sorted_charge)

    allocate (grouped_energy(size(sorted_energy)), grouped_charge(size(sorted_charge)))
    group_count = 0
    do i = 1, size(sorted_energy)
      if (group_count == 0) then
        group_count = 1
        grouped_energy(group_count) = sorted_energy(i)
        grouped_charge(group_count) = sorted_charge(i)
      else if (sorted_energy(i) /= grouped_energy(group_count)) then
        group_count = group_count + 1
        grouped_energy(group_count) = sorted_energy(i)
        grouped_charge(group_count) = sorted_charge(i)
      else
        grouped_charge(group_count) = grouped_charge(group_count) + sorted_charge(i)
      end if
    end do
    allocate (distribution%energy_j(group_count), distribution%group_charge_c(group_count))
    distribution%energy_j = grouped_energy(:group_count)
    distribution%group_charge_c = grouped_charge(:group_count)
    allocate (cumulative(0:group_count))
    cumulative(0) = 0.0_dp
    do i = 1, group_count
      cumulative(i) = cumulative(i - 1) + distribution%group_charge_c(i)
    end do
    if (.not. all(ieee_is_finite(cumulative)) .or. cumulative(group_count) <= 0.0_dp) then
      distribution = measured_interface_energy_distribution_type()
      message = 'measured Zhao cumulative interface charge is invalid'
      return
    end if
    call move_alloc(cumulative, distribution%cumulative_charge_c)
    distribution%group_count = int(group_count, i32)
    distribution%total_charge_c = distribution%cumulative_charge_c(group_count)
    distribution%ready = .true.
    status = dynamic_k0_ok
  end subroutine build_measured_interface_energy_distribution

  !> Zhao branch上でQ=Q_base+Q_escape[B(Q)]を経験分布の一般化根として解く。
  subroutine advance_dynamic_k0_zhao( &
    kinetic_options, initial_outer, lower_boundary_model, area_xy_m2, &
    surface_charge_before_c, surface_charge_base_c, time_step_s, distribution, &
    step, accepted_outer, effective_source_scale, message &
    )
    type(kinetic_outer_plasma_options_type), intent(in) :: kinetic_options
    type(outer_plasma_state_type), intent(in) :: initial_outer
    character(len=*), intent(in) :: lower_boundary_model
    real(dp), intent(in) :: area_xy_m2, surface_charge_before_c, surface_charge_base_c, time_step_s
    type(measured_interface_energy_distribution_type), intent(in) :: distribution
    type(dynamic_k0_step_type), intent(out) :: step
    type(outer_plasma_state_type), intent(out) :: accepted_outer
    real(dp), intent(out) :: effective_source_scale
    character(len=*), intent(out) :: message

    type(kinetic_outer_plasma_options_type) :: trial_options
    type(outer_plasma_state_type) :: accepted_state
    type(zhao_charge_root_type) :: source_anchor_root, trial_root, lower_root, upper_root
    type(zhao_charge_root_type) :: accepted_root
    type(zhao_charge_root_type) :: marginal_lower_root, marginal_upper_root
    real(dp) :: boundary_factor, reference_current_density, outward_current_density
    real(dp) :: trial_charge, trial_barrier, lower_charge, upper_charge
    real(dp) :: lower_barrier, upper_barrier, target_energy, midpoint_charge, midpoint_barrier
    real(dp) :: escaped_charge, empirical_escaped_charge, marginal_fraction, raw_marginal_fraction
    real(dp) :: accepted_charge, accepted_barrier, accepted_marginal_fraction
    real(dp) :: charge_tolerance, barrier_tolerance, endpoint_barrier_tolerance
    real(dp) :: anchor_charge, photoelectron_voltage_scale
    real(dp) :: reference_charge, reference_escaped_charge, initial_barrier
    real(dp) :: source_log_shift, normalized_field_shift, normalized_barrier_shift
    real(dp) :: initial_ambient_electron_barrier, trial_ambient_electron_barrier
    real(dp) :: normalized_ambient_potential_shift, normalized_ambient_electron_barrier_shift
    real(dp) :: normalized_ambient_ion_potential_shift, ambient_density_log_shift
    real(dp) :: ion_energy_scale, ion_drift_energy, ion_scaled_speed, sqrt_half_ion_mass
    real(dp) :: resolved_barrier, resolved_barrier_tolerance
    integer(i32) :: candidate_index, iteration, trial_status
    integer(i32) :: group_index, lower_index, upper_index, midpoint_index
    character(len=256) :: trial_message
    logical :: lower_predicate, upper_predicate, midpoint_predicate

    step = dynamic_k0_step_type()
    accepted_outer = outer_plasma_state_type()
    accepted_state = outer_plasma_state_type()
    effective_source_scale = 0.0_dp
    message = ''
    if (trim(lower_ascii(kinetic_options%kinetic_closure)) /= 'zhao_charge_driven') then
      message = 'dynamic Zhao mean requires kinetic_closure=zhao_charge_driven'
      return
    end if
    if (.not. initial_outer%ready .or. index('ABC', initial_outer%zhao_branch) == 0 .or. &
        .not. ieee_is_finite(initial_outer%photoelectron_source_scale) .or. &
        .not. ieee_is_finite(initial_outer%debye_length) .or. &
        .not. ieee_is_finite(initial_outer%zhao_phi_minimum) .or. &
        .not. ieee_is_finite(initial_outer%zhao_electron_density_infinity) .or. &
        initial_outer%photoelectron_source_scale <= 0.0_dp .or. initial_outer%debye_length <= 0.0_dp) then
      message = 'dynamic Zhao mean requires a committed A/B/C branch anchor'
      return
    end if
    if (initial_outer%zhao_electron_density_infinity <= 0.0_dp) then
      message = 'dynamic Zhao mean requires a positive committed ambient electron density'
      return
    end if
    select case (trim(lower_ascii(lower_boundary_model)))
    case ('e_bottom_zero')
      boundary_factor = 1.0_dp
    case ('symmetric_vacuum')
      boundary_factor = 2.0_dp
    case default
      message = 'dynamic Zhao mean received an unsupported lower boundary model'
      return
    end select
    if (.not. distribution%ready .or. distribution%group_count < 1_i32) then
      message = 'dynamic Zhao mean received an incomplete measured energy distribution'
      return
    end if
    if (.not. allocated(distribution%energy_j) .or. &
        .not. allocated(distribution%group_charge_c) .or. &
        .not. allocated(distribution%cumulative_charge_c)) then
      message = 'dynamic Zhao mean received an unallocated measured energy distribution'
      return
    end if
    if (size(distribution%energy_j) /= distribution%group_count .or. &
        size(distribution%group_charge_c) /= distribution%group_count) then
      message = 'dynamic Zhao mean received inconsistent measured energy group sizes'
      return
    end if
    if (lbound(distribution%cumulative_charge_c, 1) /= 0 .or. &
        ubound(distribution%cumulative_charge_c, 1) /= distribution%group_count) then
      message = 'dynamic Zhao mean received inconsistent cumulative energy groups'
      return
    end if
    if (any(.not. ieee_is_finite(distribution%energy_j)) .or. &
        any(.not. ieee_is_finite(distribution%group_charge_c)) .or. &
        any(.not. ieee_is_finite(distribution%cumulative_charge_c)) .or. &
        any(distribution%energy_j < 0.0_dp) .or. &
        any(distribution%group_charge_c <= 0.0_dp)) then
      message = 'dynamic Zhao mean received invalid measured energy groups'
      return
    end if
    if (distribution%group_count > 1_i32) then
      if (any(distribution%energy_j(1:distribution%group_count - 1_i32) <= &
              distribution%energy_j(2:distribution%group_count))) then
        message = 'dynamic Zhao measured energy groups are not strictly descending'
        return
      end if
    end if
    if (distribution%cumulative_charge_c(0) /= 0.0_dp .or. &
        any(distribution%cumulative_charge_c(1:distribution%group_count) <= &
            distribution%cumulative_charge_c(0:distribution%group_count - 1_i32)) .or. &
        distribution%cumulative_charge_c(distribution%group_count) /= distribution%total_charge_c) then
      message = 'dynamic Zhao cumulative measured charge is inconsistent'
      return
    end if
    if (.not. all(ieee_is_finite([ &
                                 area_xy_m2, surface_charge_before_c, surface_charge_base_c, time_step_s, &
                                 distribution%total_charge_c, kinetic_options%electron_charge, &
                                 kinetic_options%electron_temperature_j, kinetic_options%photoelectron_charge, &
                                 kinetic_options%photoelectron_mass, kinetic_options%photoelectron_temperature_j, &
                                 kinetic_options%photoelectron_reference_density, kinetic_options%zhao_alpha_deg, &
                                 kinetic_options%ion_charge, kinetic_options%ion_mass, &
                                 kinetic_options%ion_temperature_j, kinetic_options%ion_drift_infinity &
                                 ])) .or. area_xy_m2 <= 0.0_dp .or. time_step_s <= 0.0_dp .or. &
        distribution%total_charge_c <= 0.0_dp .or. kinetic_options%electron_charge >= 0.0_dp .or. &
        kinetic_options%electron_temperature_j <= 0.0_dp .or. kinetic_options%photoelectron_charge >= 0.0_dp .or. &
        kinetic_options%photoelectron_mass <= 0.0_dp .or. kinetic_options%photoelectron_temperature_j <= 0.0_dp .or. &
        kinetic_options%photoelectron_reference_density <= 0.0_dp .or. kinetic_options%ion_charge <= 0.0_dp .or. &
        kinetic_options%ion_mass <= 0.0_dp .or. kinetic_options%ion_temperature_j < 0.0_dp .or. &
        kinetic_options%ion_drift_infinity <= 0.0_dp) then
      message = 'dynamic Zhao mean received invalid geometry, source, ambient, or photoelectron parameters'
      return
    end if

    ! Validate the frozen-ion response scale before any Zhao continuation.
    ! Evaluate 0.5*m*u**2 without forming u**2 first: an overflowing
    ! intermediate must fail closed instead of turning dU_i/E_i into zero.
    sqrt_half_ion_mass = sqrt(0.5_dp*kinetic_options%ion_mass)
    if (.not. ieee_is_finite(sqrt_half_ion_mass) .or. sqrt_half_ion_mass <= 0.0_dp .or. &
        abs(kinetic_options%ion_drift_infinity) > &
        sqrt(huge(1.0_dp))/sqrt_half_ion_mass) then
      step%status = dynamic_zhao_no_physical_root
      message = 'dynamic Zhao ambient ion energy scale is not representable'
      return
    end if
    ion_scaled_speed = sqrt_half_ion_mass*abs(kinetic_options%ion_drift_infinity)
    if (.not. ieee_is_finite(ion_scaled_speed) .or. &
        ion_scaled_speed > sqrt(huge(1.0_dp))) then
      step%status = dynamic_zhao_no_physical_root
      message = 'dynamic Zhao ambient ion energy scale is not representable'
      return
    end if
    ion_drift_energy = ion_scaled_speed*ion_scaled_speed
    ion_energy_scale = max(kinetic_options%ion_temperature_j, ion_drift_energy)
    if (.not. ieee_is_finite(ion_energy_scale) .or. ion_energy_scale <= 0.0_dp) then
      step%status = dynamic_zhao_no_physical_root
      message = 'dynamic Zhao ambient ion energy scale is invalid'
      return
    end if

    reference_current_density = abs(kinetic_options%photoelectron_charge)* &
                                kinetic_options%photoelectron_reference_density* &
                                sin(kinetic_options%zhao_alpha_deg*pi/180.0_dp)* &
                                sqrt(2.0_dp*kinetic_options%photoelectron_temperature_j/ &
                                     kinetic_options%photoelectron_mass)/(2.0_dp*sqrt(pi))
    outward_current_density = distribution%total_charge_c/(area_xy_m2*time_step_s)
    if (.not. ieee_is_finite(reference_current_density) .or. reference_current_density <= 0.0_dp .or. &
        .not. ieee_is_finite(outward_current_density) .or. outward_current_density <= 0.0_dp) then
      message = 'dynamic Zhao measured source normalization is invalid'
      return
    end if
    effective_source_scale = outward_current_density/reference_current_density
    if (.not. ieee_is_finite(effective_source_scale) .or. effective_source_scale <= 0.0_dp) then
      message = 'dynamic Zhao effective source scale is invalid'
      return
    end if
    if (abs(log(effective_source_scale) - log(initial_outer%photoelectron_source_scale)) <= &
        64.0_dp*epsilon(1.0_dp)) then
      effective_source_scale = initial_outer%photoelectron_source_scale
    end if

    trial_options = kinetic_options
    trial_options%photoelectron_source_scale = effective_source_scale
    trial_options%photoelectron_population_fraction = 1.0_dp
    trial_options%photoelectron_column_closure_enabled = .false.
    trial_options%zhao_branch = lower_ascii(initial_outer%zhao_branch)
    charge_tolerance = max( &
                       4096.0_dp*epsilon(1.0_dp)*max( &
                       abs(surface_charge_base_c), abs(surface_charge_before_c), &
                       distribution%total_charge_c, tiny(1.0_dp)), &
                       1.0e-12_dp*distribution%total_charge_c &
                       )
    barrier_tolerance = max( &
                        4096.0_dp*epsilon(1.0_dp)*max( &
                        maxval(distribution%energy_j), kinetic_options%photoelectron_temperature_j, tiny(1.0_dp)), &
                        1.0e-10_dp*max(maxval(distribution%energy_j), tiny(1.0_dp)), &
                        1.0e-8_dp*kinetic_options%photoelectron_temperature_j &
                        )
    endpoint_barrier_tolerance = max( &
                                 8192.0_dp*epsilon(1.0_dp)*max( &
                                 maxval(distribution%energy_j), &
                                 kinetic_options%photoelectron_temperature_j, tiny(1.0_dp) &
                                 ), &
                                 1.0e-11_dp*kinetic_options%photoelectron_temperature_j &
                                 )
    if (any(abs( &
            distribution%cumulative_charge_c(1:distribution%group_count) - &
            distribution%cumulative_charge_c(0:distribution%group_count - 1_i32) - &
            distribution%group_charge_c &
            ) > charge_tolerance)) then
      message = 'dynamic Zhao cumulative measured charge increments are inconsistent'
      return
    end if
    anchor_charge = boundary_factor*eps0*area_xy_m2*initial_outer%interface_field
    if (.not. ieee_is_finite(anchor_charge) .or. &
        abs(surface_charge_before_c - anchor_charge) > charge_tolerance) then
      message = 'dynamic Zhao committed charge and branch-anchor field are inconsistent'
      return
    end if

    ! Build one common anchor on the final measured-source slice. The old
    ! profile's empirical classification supplies only a nearby reference Q;
    ! every actual candidate below is then traced at fixed final source scale.
    initial_barrier = zhao_profile_barrier_energy(initial_outer, kinetic_options%photoelectron_charge)
    if (.not. ieee_is_finite(initial_barrier)) then
      message = 'dynamic Zhao committed profile has no finite reference barrier'
      return
    end if
    reference_escaped_charge = 0.0_dp
    do group_index = 1_i32, distribution%group_count
      if (distribution%energy_j(group_index) > initial_barrier) then
        reference_escaped_charge = reference_escaped_charge + distribution%group_charge_c(group_index)
      end if
    end do
    reference_charge = surface_charge_base_c + reference_escaped_charge
    call continue_zhao_root_from_state( &
      trial_options, initial_outer, &
      reference_charge/(boundary_factor*eps0*area_xy_m2), &
      .false., source_anchor_root, trial_status, trial_message &
      )
    step%iterations = step%iterations + 1_i32
    if (trial_status /= dynamic_k0_ok) then
      step%status = trial_status
      write (message, '(a,es12.4,a,es12.4,2a)') &
        'dynamic Zhao final-source anchor failed at Q_C=', reference_charge, &
        ' E_I_V_m=', reference_charge/(boundary_factor*eps0*area_xy_m2), ': ', trim(trial_message)
      return
    end if

    ! The two endpoint paths from one final-source root cover the complete
    ! empirical charge interval. Their adaptive connected-path guards provide
    ! the same finite-precision monotonicity certificate used by the former
    ! all-breakpoint sweep. The order predicate can then be located in O(log M).
    candidate_index = -1_i32
    accepted_charge = 0.0_dp
    accepted_barrier = 0.0_dp
    accepted_marginal_fraction = -1.0_dp

    lower_charge = surface_charge_base_c
    upper_charge = surface_charge_base_c + distribution%total_charge_c
    if (.not. all(ieee_is_finite([lower_charge, upper_charge])) .or. &
        upper_charge <= lower_charge) then
      step%status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao empirical charge interval is not representable'
      return
    end if
    call evaluate_zhao_trial( &
      trial_options, source_anchor_root, boundary_factor, area_xy_m2, lower_charge, &
      lower_root, lower_barrier, trial_status, trial_message &
      )
    step%iterations = step%iterations + 1_i32
    if (trial_status /= dynamic_k0_ok) then
      step%status = trial_status
      message = 'dynamic Zhao lower endpoint certificate failed: '//trim(trial_message)
      return
    end if
    call evaluate_zhao_trial( &
      trial_options, source_anchor_root, boundary_factor, area_xy_m2, upper_charge, &
      upper_root, upper_barrier, trial_status, trial_message &
      )
    step%iterations = step%iterations + 1_i32
    if (trial_status /= dynamic_k0_ok) then
      step%status = trial_status
      message = 'dynamic Zhao upper endpoint certificate failed: '//trim(trial_message)
      return
    end if
    if (upper_barrier < lower_barrier - barrier_tolerance) then
      step%status = dynamic_zhao_nonmonotone_barrier
      message = 'dynamic Zhao endpoint barriers violate charge ordering'
      return
    end if

    lower_index = 0_i32
    upper_index = distribution%group_count
    lower_predicate = lower_barrier >= distribution%energy_j(1)
    upper_predicate = upper_barrier >= distribution%energy_j(distribution%group_count)

    if (lower_predicate) then
      candidate_index = 0_i32
      accepted_charge = lower_charge
      accepted_barrier = lower_barrier
      accepted_root = lower_root
    else if (.not. upper_predicate) then
      candidate_index = distribution%group_count
      accepted_charge = upper_charge
      accepted_barrier = upper_barrier
      accepted_root = upper_root
      ! At Q_M the all-escape solution and a unit-weight marginal solution
      ! carry exactly the same charge. Preserve the latter only when the
      ! analytic root barrier and measured endpoint agree to roundoff. Do not
      ! apply this tolerance to interior order predicates: doing so can hide a
      ! genuinely escaping, very small macro-charge group.
      if (abs(upper_barrier - distribution%energy_j(distribution%group_count)) <= &
          endpoint_barrier_tolerance) then
        accepted_marginal_fraction = 1.0_dp
      end if
    else
      ! P_k=[B(Q_k)>=E_{k+1}] for k<M and
      ! P_M=[B(Q_M)>=E_M] is monotone on the certified charge interval.
      ! Retain a false/true integer bracket and find its first true index.
      do while (upper_index - lower_index > 1_i32)
        midpoint_index = lower_index + (upper_index - lower_index)/2_i32
        midpoint_charge = surface_charge_base_c + &
                          distribution%cumulative_charge_c(midpoint_index)
        call evaluate_zhao_trial( &
          trial_options, source_anchor_root, boundary_factor, area_xy_m2, &
          midpoint_charge, trial_root, midpoint_barrier, trial_status, trial_message &
          )
        step%iterations = step%iterations + 1_i32
        if (trial_status /= dynamic_k0_ok) then
          step%status = trial_status
          write (message, '(a,i0,2a)') &
            'dynamic Zhao order-statistic path failed at group ', midpoint_index, &
            ': ', trim(trial_message)
          return
        end if
        if (midpoint_barrier < lower_barrier - barrier_tolerance .or. &
            midpoint_barrier > upper_barrier + barrier_tolerance) then
          step%status = dynamic_zhao_nonmonotone_barrier
          write (message, '(a,i0)') &
            'dynamic Zhao order-statistic barrier ordering failed at group ', midpoint_index
          return
        end if
        midpoint_predicate = midpoint_barrier >= &
          distribution%energy_j(midpoint_index + 1_i32)
        if (midpoint_predicate) then
          upper_index = midpoint_index
          upper_charge = midpoint_charge
          upper_barrier = midpoint_barrier
          upper_root = trial_root
        else
          lower_index = midpoint_index
          lower_charge = midpoint_charge
          lower_barrier = midpoint_barrier
          lower_root = trial_root
        end if
      end do

      group_index = upper_index
      target_energy = distribution%energy_j(group_index)
      if (group_index < distribution%group_count .and. &
          upper_barrier < target_energy) then
        candidate_index = group_index
        accepted_charge = upper_charge
        accepted_barrier = upper_barrier
        accepted_root = upper_root
      else
        ! Equality belongs to the return side. The false/true bracket gives
        ! B(Q_{k-1})<E_k<=B(Q_k); theta=0 remains the preceding pure
        ! solution, while a theta=1 endpoint is retained as marginal.
        if (lower_barrier >= target_energy .or. upper_barrier < target_energy) then
          step%status = dynamic_k0_numerical_failure
          message = 'dynamic Zhao marginal order-statistic bracket is inconsistent'
          return
        end if
        marginal_lower_root = lower_root
        marginal_upper_root = upper_root
        do iteration = 1_i32, marginal_bisection_limit
          if (abs(upper_charge - lower_charge) <= charge_tolerance) exit
          midpoint_charge = lower_charge + 0.5_dp*(upper_charge - lower_charge)
          call evaluate_zhao_trial( &
            trial_options, marginal_lower_root, boundary_factor, area_xy_m2, &
            midpoint_charge, trial_root, midpoint_barrier, trial_status, trial_message &
            )
          step%iterations = step%iterations + 1_i32
          if (trial_status /= dynamic_k0_ok) then
            step%status = trial_status
            message = 'dynamic Zhao marginal order-statistic path failed: '//trim(trial_message)
            return
          end if
          if (midpoint_barrier < lower_barrier - barrier_tolerance .or. &
              midpoint_barrier > upper_barrier + barrier_tolerance) then
            step%status = dynamic_zhao_nonmonotone_barrier
            message = 'dynamic Zhao marginal barrier ordering failed'
            return
          end if
          if (midpoint_barrier >= target_energy) then
            upper_charge = midpoint_charge
            upper_barrier = midpoint_barrier
            marginal_upper_root = trial_root
          else
            lower_charge = midpoint_charge
            lower_barrier = midpoint_barrier
            marginal_lower_root = trial_root
          end if
        end do
        if (abs(upper_charge - lower_charge) > charge_tolerance) then
          step%status = dynamic_zhao_no_physical_root
          message = 'dynamic Zhao marginal root did not converge in charge'
          return
        end if
        ! Keep the upper return-side endpoint even when the lower endpoint has
        ! a slightly smaller finite-precision energy residual.
        trial_charge = upper_charge
        trial_barrier = upper_barrier
        trial_root = marginal_upper_root
        if (abs(trial_barrier - target_energy) > barrier_tolerance) then
          step%status = dynamic_zhao_no_physical_root
          write (message, '(a,3(a,es12.4))') &
            'dynamic Zhao marginal root did not meet its energy residual:', &
            ' barrier_J=', trial_barrier, ' target_J=', target_energy, &
            ' tolerance_J=', barrier_tolerance
          return
        end if
        raw_marginal_fraction = ( &
                                trial_charge - surface_charge_base_c - &
                                distribution%cumulative_charge_c(group_index - 1_i32) &
                                )/distribution%group_charge_c(group_index)
        if (raw_marginal_fraction < &
            -charge_tolerance/distribution%group_charge_c(group_index) .or. &
            raw_marginal_fraction > &
            1.0_dp + charge_tolerance/distribution%group_charge_c(group_index)) then
          step%status = dynamic_zhao_no_physical_root
          message = 'dynamic Zhao marginal charge fraction left its empirical group'
          return
        end if
        candidate_index = group_index
        accepted_charge = trial_charge
        accepted_barrier = trial_barrier
        accepted_root = trial_root
        accepted_marginal_fraction = min(1.0_dp, max(0.0_dp, raw_marginal_fraction))
      end if
    end if

    if (candidate_index < 0_i32 .or. &
        candidate_index > distribution%group_count) then
      step%status = dynamic_zhao_no_physical_root
      message = 'dynamic Zhao order-statistic search found no generalized charge root'
      return
    end if

    trial_charge = accepted_charge
    trial_barrier = accepted_barrier
    marginal_fraction = accepted_marginal_fraction
    if (marginal_fraction >= 0.0_dp) then
      escaped_charge = distribution%cumulative_charge_c(candidate_index - 1_i32) + &
                       marginal_fraction*distribution%group_charge_c(candidate_index)
      step%marginal_photoelectron_energy_j = distribution%energy_j(candidate_index)
    else
      escaped_charge = distribution%cumulative_charge_c(candidate_index)
    end if
    if (abs(surface_charge_base_c + escaped_charge - trial_charge) > charge_tolerance) then
      step%status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao selected CDF charge does not reproduce its root charge'
      return
    end if

    photoelectron_voltage_scale = kinetic_options%photoelectron_temperature_j/ &
                                  abs(kinetic_options%photoelectron_charge)
    source_log_shift = abs( &
                       log(effective_source_scale) - &
                       log(initial_outer%photoelectron_source_scale) &
                       )
    normalized_field_shift = abs( &
                             accepted_root%interface_field_v_m - initial_outer%interface_field &
                             )*initial_outer%debye_length/photoelectron_voltage_scale
    normalized_barrier_shift = abs(trial_barrier - initial_barrier)/ &
                               kinetic_options%photoelectron_temperature_j
    if (.not. ieee_is_finite(accepted_root%n_swe_inf_m3) .or. &
        accepted_root%n_swe_inf_m3 <= 0.0_dp) then
      step%status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao accepted root has invalid ambient electron density'
      return
    end if
    initial_ambient_electron_barrier = max( &
                                       0.0_dp, &
                                       kinetic_options%electron_charge*initial_outer%interface_potential, &
                                       kinetic_options%electron_charge*initial_outer%zhao_phi_minimum &
                                       )
    trial_ambient_electron_barrier = max( &
                                     0.0_dp, &
                                     kinetic_options%electron_charge*accepted_root%phi0_v, &
                                     kinetic_options%electron_charge*accepted_root%phi_m_v &
                                     )
    normalized_ambient_potential_shift = abs( &
                                         kinetic_options%electron_charge* &
                                         (accepted_root%phi0_v - initial_outer%interface_potential) &
                                         )/kinetic_options%electron_temperature_j
    normalized_ambient_electron_barrier_shift = abs( &
                                                trial_ambient_electron_barrier - &
                                                initial_ambient_electron_barrier &
                                                )/kinetic_options%electron_temperature_j
    normalized_ambient_ion_potential_shift = abs( &
                                             kinetic_options%ion_charge* &
                                             (accepted_root%phi0_v - initial_outer%interface_potential) &
                                             )/ion_energy_scale
    ambient_density_log_shift = abs( &
                                log(accepted_root%n_swe_inf_m3) - &
                                log(initial_outer%zhao_electron_density_infinity) &
                                )
    ! A common interface-potential shift changes the infinity-to-interface
    ! mapping of the ambient reservoir, whose response scale is Te/|qe|.
    ! Its profile cutoff, density amplitude, and drifting-ion energy are
    ! independent frozen-cohort coordinates. Emitted photoelectrons instead
    ! respond to the profile-relative barrier and field deformation, which
    ! retain their Tpe and lambda_D,pe scales.
    if (.not. all(ieee_is_finite([ &
                                 normalized_ambient_potential_shift, &
                                 normalized_ambient_electron_barrier_shift, &
                                 normalized_ambient_ion_potential_shift, ambient_density_log_shift &
                                 ])) .or. &
        normalized_ambient_potential_shift > frozen_interface_shift_limit .or. &
        normalized_ambient_electron_barrier_shift > frozen_interface_shift_limit .or. &
        normalized_ambient_ion_potential_shift > frozen_interface_shift_limit .or. &
        ambient_density_log_shift > frozen_interface_shift_limit) then
      step%status = dynamic_zhao_no_physical_root
      write (message, '(a,7(a,es11.3))') &
        'dynamic Zhao frozen ambient cohort trust failure:', &
        ' dphi/Te=', &
        normalized_ambient_potential_shift, &
        ' dBe/Te=', normalized_ambient_electron_barrier_shift, &
        ' dphi/Ei=', normalized_ambient_ion_potential_shift, &
        ' dlogne=', ambient_density_log_shift, &
        ' lim=', frozen_interface_shift_limit, &
        ' phi0=', initial_outer%interface_potential, &
        ' phi1=', accepted_root%phi0_v
      return
    end if
    if (.not. all(ieee_is_finite([ &
                                 source_log_shift, normalized_field_shift, normalized_barrier_shift &
                                 ])) .or. &
        source_log_shift > frozen_source_log_shift_limit .or. &
        normalized_field_shift > frozen_field_shift_limit .or. &
        normalized_barrier_shift > frozen_barrier_shift_limit) then
      step%status = dynamic_zhao_no_physical_root
      write (message, '(a,3(a,es12.4))') &
        'dynamic Zhao frozen interface cohort left its trust region:', &
        ' delta_log_source=', source_log_shift, &
        ' lambda_delta_E_over_Tpe=', normalized_field_shift, &
        ' delta_barrier_over_Tpe=', normalized_barrier_shift
      return
    end if
    if (abs(boundary_factor*eps0*area_xy_m2*accepted_root%interface_field_v_m - trial_charge) > &
        charge_tolerance) then
      step%status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao accepted root does not reproduce its prescribed charge'
      return
    end if

    trial_options%interface_field = accepted_root%interface_field_v_m
    call materialize_outer_plasma_zhao_root( &
      trial_options, accepted_root, accepted_state, trial_status, trial_message &
      )
    if (trial_status /= outer_plasma_ok) then
      select case (trial_status)
      case (outer_plasma_invalid)
        step%status = dynamic_k0_invalid
      case (outer_plasma_no_physical_solution)
        step%status = dynamic_zhao_no_physical_root
      case default
        step%status = dynamic_k0_numerical_failure
      end select
      message = 'dynamic Zhao accepted profile materialization failed: '//trim(trial_message)
      return
    end if
    resolved_barrier = zhao_profile_barrier_energy( &
                       accepted_state, kinetic_options%photoelectron_charge &
                       )
    resolved_barrier_tolerance = max( &
                                 8192.0_dp*epsilon(1.0_dp)*max( &
                                 abs(resolved_barrier), abs(trial_barrier), &
                                 kinetic_options%photoelectron_temperature_j, tiny(1.0_dp) &
                                 ), &
                                 1.0e-11_dp*kinetic_options%photoelectron_temperature_j &
                                 )
    if (.not. accepted_state%ready .or. accepted_state%zhao_branch /= accepted_root%branch .or. &
        .not. ieee_is_finite(resolved_barrier) .or. &
        abs(resolved_barrier - trial_barrier) > resolved_barrier_tolerance .or. &
        abs(accepted_state%interface_potential - accepted_root%phi0_v) > &
        256.0_dp*epsilon(1.0_dp)*max( &
        abs(accepted_state%interface_potential), abs(accepted_root%phi0_v), 1.0_dp &
        ) .or. &
        abs(boundary_factor*eps0*area_xy_m2*accepted_state%interface_field - trial_charge) > &
        charge_tolerance) then
      step%status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao accepted analytic root and resolved profile disagree'
      return
    end if
    accepted_state%photoelectron_current_density = escaped_charge/(area_xy_m2*time_step_s)
    accepted_state%total_current_density = &
      accepted_state%electron_current_density + accepted_state%ion_current_density + &
      accepted_state%photoelectron_current_density

    step%status = dynamic_k0_ok
    step%surface_charge_before_c = surface_charge_before_c
    step%surface_charge_after_c = surface_charge_base_c + escaped_charge
    step%interface_potential_before_v = initial_outer%interface_potential
    step%interface_potential_after_v = accepted_state%interface_potential
    step%interface_field_after_v_m = accepted_state%interface_field
    step%photoelectron_escape_fraction = escaped_charge/distribution%total_charge_c
    step%photoelectron_return_fraction = 1.0_dp - step%photoelectron_escape_fraction
    step%photoelectron_source_charge_c = distribution%total_charge_c
    step%photoelectron_barrier_energy_j = trial_barrier
    step%photoelectron_energy_max_j = distribution%energy_j(1)
    step%marginal_photoelectron_escape_fraction = marginal_fraction
    step%zhao_effective_source_scale = effective_source_scale
    step%zhao_branch = accepted_state%zhao_branch
    empirical_escaped_charge = 0.0_dp
    do group_index = 1_i32, distribution%group_count
      empirical_escaped_charge = empirical_escaped_charge + &
                                 distribution%group_charge_c(group_index)* &
                                 measured_sample_escape_fraction(distribution%energy_j(group_index), step)
    end do
    step%backward_euler_residual_charge_c = &
      step%surface_charge_after_c - surface_charge_base_c - empirical_escaped_charge
    if (.not. all(ieee_is_finite([ &
                                 step%surface_charge_after_c, step%interface_potential_after_v, &
                                 step%interface_field_after_v_m, step%photoelectron_escape_fraction, &
                                 step%photoelectron_return_fraction, step%photoelectron_barrier_energy_j, &
                                 step%photoelectron_source_charge_c, step%photoelectron_energy_max_j, &
                                 step%zhao_effective_source_scale, &
                                 step%backward_euler_residual_charge_c &
                                 ])) .or. index('ABC', step%zhao_branch) == 0 .or. &
        abs(step%backward_euler_residual_charge_c) > charge_tolerance) then
      step%status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao accepted result failed its charge or branch invariant'
      return
    end if
    accepted_outer = accepted_state
  end subroutine advance_dynamic_k0_zhao

  !> profile全体とinfinity endpointから必要な法線運動エネルギーを返す。
  pure function zhao_profile_barrier_energy(state, charge) result(barrier_j)
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: charge
    real(dp) :: barrier_j

    barrier_j = ieee_value(0.0_dp, ieee_quiet_nan)
    if (.not. state%ready .or. state%profile_n < 2_i32) return
    if (.not. allocated(state%potential)) return
    if (size(state%potential) /= state%profile_n) return
    if (.not. ieee_is_finite(charge) .or. charge >= 0.0_dp .or. &
        .not. ieee_is_finite(state%interface_potential) .or. &
        .not. ieee_is_finite(state%infinity_potential) .or. &
        any(.not. ieee_is_finite(state%potential))) return
    barrier_j = max( &
                0.0_dp, maxval(charge*(state%potential - state%interface_potential)), &
                charge*(state%infinity_potential - state%interface_potential) &
                )
  end function zhao_profile_barrier_energy

  !> accepted generalized rootに対する個々の実測rayのescape weightを返す。
  pure function measured_sample_escape_fraction(energy_j, step) result(fraction)
    real(dp), intent(in) :: energy_j
    type(dynamic_k0_step_type), intent(in) :: step
    real(dp) :: fraction

    fraction = 0.0_dp
    if (.not. ieee_is_finite(energy_j) .or. energy_j < 0.0_dp .or. &
        step%status /= dynamic_k0_ok) return
    if (step%marginal_photoelectron_escape_fraction >= 0.0_dp .and. &
        energy_j == step%marginal_photoelectron_energy_j) then
      fraction = step%marginal_photoelectron_escape_fraction
    else if (energy_j > step%photoelectron_barrier_energy_j) then
      fraction = 1.0_dp
    end if
  end function measured_sample_escape_fraction

  subroutine evaluate_zhao_trial( &
    options, initial_root, boundary_factor, area_xy_m2, surface_charge_c, &
    root, barrier_j, status, message &
    )
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(zhao_charge_root_type), intent(in) :: initial_root
    real(dp), intent(in) :: boundary_factor, area_xy_m2, surface_charge_c
    type(zhao_charge_root_type), intent(out) :: root
    real(dp), intent(out) :: barrier_j
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: target_field

    root = zhao_charge_root_type()
    barrier_j = ieee_value(0.0_dp, ieee_quiet_nan)
    status = dynamic_k0_numerical_failure
    message = ''
    target_field = surface_charge_c/(boundary_factor*eps0*area_xy_m2)
    if (.not. ieee_is_finite(target_field)) then
      message = 'dynamic Zhao candidate interface field is non-finite'
      return
    end if
    call continue_zhao_root_from_root( &
      options, initial_root, target_field, .true., root, status, message &
      )
    if (status /= dynamic_k0_ok) return
    call validate_zhao_trial_root( &
      initial_root%branch, options%photoelectron_charge, root, barrier_j, status, message &
      )
  end subroutine evaluate_zhao_trial

  subroutine validate_zhao_trial_root( &
    initial_branch, charge, root, barrier_j, status, message &
    )
    character(len=1), intent(in) :: initial_branch
    real(dp), intent(in) :: charge
    type(zhao_charge_root_type), intent(in) :: root
    real(dp), intent(out) :: barrier_j
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    barrier_j = ieee_value(0.0_dp, ieee_quiet_nan)
    status = dynamic_k0_numerical_failure
    message = ''
    if (root%branch /= initial_branch .or. index('ABC', root%branch) == 0) then
      status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao candidate left its committed branch'
      return
    end if
    barrier_j = zhao_charge_root_barrier_energy(root, charge)
    if (.not. ieee_is_finite(barrier_j) .or. barrier_j < 0.0_dp) then
      status = dynamic_k0_numerical_failure
      message = 'dynamic Zhao candidate produced an invalid analytic root barrier'
      return
    end if
    status = dynamic_k0_ok
  end subroutine validate_zhao_trial_root

  !> source scale と prescribed field を同じroot sheet上で追う。
  subroutine continue_zhao_root_from_state( &
    options, initial_state, target_field, require_monotone_barrier, &
    root, status, message &
    )
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_state_type), intent(in) :: initial_state
    real(dp), intent(in) :: target_field
    logical, intent(in) :: require_monotone_barrier
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(kinetic_outer_plasma_options_type) :: trial_options
    type(zhao_connected_path_certificate_type) :: certificate
    integer(i32) :: outer_status
    character(len=256) :: outer_message

    root = zhao_charge_root_type()
    status = dynamic_k0_invalid
    message = ''
    if (.not. initial_state%ready .or. index('ABC', initial_state%zhao_branch) == 0 .or. &
        .not. all(ieee_is_finite([ &
                                 target_field, options%photoelectron_source_scale, &
                                 initial_state%interface_field, initial_state%photoelectron_source_scale &
                                 ])) .or. options%photoelectron_source_scale <= 0.0_dp .or. &
        initial_state%photoelectron_source_scale <= 0.0_dp) then
      message = 'dynamic Zhao continuation received an invalid committed state or target'
      return
    end if

    trial_options = options
    trial_options%zhao_branch = lower_ascii(initial_state%zhao_branch)
    trial_options%interface_field = target_field
    call continue_outer_plasma_zhao_connected_root( &
      trial_options, initial_state, require_monotone_barrier, root, certificate, &
      outer_status, outer_message &
      )
    call map_connected_zhao_status(certificate, outer_status, outer_message, status, message)
  end subroutine continue_zhao_root_from_state

  !> fixed source slice上のcertified rootから次のfield rootを追う。
  subroutine continue_zhao_root_from_root( &
    options, initial_root, target_field, require_monotone_barrier, &
    root, status, message &
    )
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(zhao_charge_root_type), intent(in) :: initial_root
    real(dp), intent(in) :: target_field
    logical, intent(in) :: require_monotone_barrier
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(kinetic_outer_plasma_options_type) :: trial_options
    type(zhao_connected_path_certificate_type) :: certificate
    integer(i32) :: outer_status
    character(len=256) :: outer_message

    root = zhao_charge_root_type()
    status = dynamic_k0_invalid
    message = ''
    if (index('ABC', initial_root%branch) == 0 .or. &
        .not. all(ieee_is_finite([ &
                                 target_field, options%photoelectron_source_scale, &
                                 initial_root%interface_field_v_m &
                                 ])) .or. options%photoelectron_source_scale <= 0.0_dp) then
      message = 'dynamic Zhao root continuation received an invalid anchor or target'
      return
    end if

    trial_options = options
    trial_options%zhao_branch = lower_ascii(initial_root%branch)
    trial_options%interface_field = target_field
    call continue_outer_plasma_zhao_connected_root( &
      trial_options, initial_root, require_monotone_barrier, root, certificate, &
      outer_status, outer_message &
      )
    call map_connected_zhao_status(certificate, outer_status, outer_message, status, message)
  end subroutine continue_zhao_root_from_root

  subroutine map_connected_zhao_status(certificate, outer_status, outer_message, status, message)
    type(zhao_connected_path_certificate_type), intent(in) :: certificate
    integer(i32), intent(in) :: outer_status
    character(len=*), intent(in) :: outer_message
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = dynamic_k0_invalid
    message = ''
    select case (outer_status)
    case (outer_plasma_ok)
      status = dynamic_k0_ok
    case (outer_plasma_no_physical_solution)
      if (certificate%nonmonotone_barrier_detected) then
        status = dynamic_zhao_nonmonotone_barrier
      else
        status = dynamic_zhao_no_physical_root
      end if
      message = trim(certificate%reason)//': '//trim(outer_message)
    case (outer_plasma_invalid)
      status = dynamic_k0_invalid
      message = 'invalid Zhao connected-path request: '//trim(outer_message)
    case default
      status = dynamic_k0_numerical_failure
      message = 'numerical Zhao connected-path failure: '//trim(outer_message)
    end select
  end subroutine map_connected_zhao_status

  subroutine stable_sort_energy_descending(energy, charge)
    real(dp), intent(inout) :: energy(:), charge(:)
    real(dp), allocatable :: work_energy(:), work_charge(:)
    integer :: width, left, middle, right, i, j, k, n

    n = size(energy)
    if (size(charge) /= n .or. n <= 1) return
    allocate (work_energy(n), work_charge(n))
    width = 1
    do while (width < n)
      left = 1
      do while (left <= n)
        middle = min(left + width, n + 1)
        right = min(left + 2*width - 1, n)
        i = left
        j = middle
        k = left
        do while (i < middle .and. j <= right)
          if (energy(i) >= energy(j)) then
            work_energy(k) = energy(i)
            work_charge(k) = charge(i)
            i = i + 1
          else
            work_energy(k) = energy(j)
            work_charge(k) = charge(j)
            j = j + 1
          end if
          k = k + 1
        end do
        do while (i < middle)
          work_energy(k) = energy(i)
          work_charge(k) = charge(i)
          i = i + 1
          k = k + 1
        end do
        do while (j <= right)
          work_energy(k) = energy(j)
          work_charge(k) = charge(j)
          j = j + 1
          k = k + 1
        end do
        left = left + 2*width
      end do
      energy = work_energy
      charge = work_charge
      if (width > n/2) exit
      width = 2*width
    end do
  end subroutine stable_sort_energy_descending

end module bem_dynamic_k0_zhao
