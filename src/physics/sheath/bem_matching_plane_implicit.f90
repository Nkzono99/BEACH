!> Backward-Euler displacement solve for a matching-plane response, including continuation and MPI broadcast.
module bem_matching_plane_implicit
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_matching_plane_response_provider, only: matching_plane_response_provider_type, &
                                                  matching_plane_provider_ok, matching_plane_provider_invalid_argument, &
                                                  matching_plane_provider_no_physical_solution, &
                                                  matching_plane_provider_numerical_failure
  use bem_matching_plane_zhao, only: matching_plane_zhao_root_seed_type
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  implicit none
  private
  public :: solve_matching_implicit_zero_mode
contains

  subroutine solve_matching_implicit_zero_mode( &
    provider, mpi, displacement_before, displacement_seed, duration, displacement_bounded, &
    displacement_min, displacement_max, &
    displacement_scale, search_direction, feedback_reference, electron_charge, ion_charge, photoelectron_active, &
    photoelectron_charge, root_before, root_after, displacement_after, response_after &
    )
    type(matching_plane_response_provider_type), intent(inout) :: provider
    type(mpi_context), intent(in) :: mpi
    logical, intent(in) :: displacement_bounded
    real(dp), intent(in) :: displacement_before, displacement_seed, duration
    real(dp), intent(in) :: displacement_min, displacement_max, displacement_scale
    integer(i32), intent(in) :: search_direction
    real(dp), intent(in) :: feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    logical, intent(in) :: photoelectron_active
    type(matching_plane_zhao_root_seed_type), intent(in) :: root_before
    type(matching_plane_zhao_root_seed_type), intent(out) :: root_after
    real(dp), intent(out) :: displacement_after, response_after(6)

    real(dp) :: result_packet(7)
    integer(i32) :: status, status_packet(1)
    character(len=512) :: message

    displacement_after = 0.0_dp
    response_after = 0.0_dp
    root_after = matching_plane_zhao_root_seed_type()
    result_packet = 0.0_dp
    status = matching_plane_provider_ok
    message = ''
    if (mpi_is_root(mpi)) then
      call solve_matching_implicit_zero_mode_local( &
        provider, displacement_before, displacement_seed, duration, displacement_bounded, &
        displacement_min, displacement_max, &
        displacement_scale, search_direction, feedback_reference, electron_charge, ion_charge, photoelectron_active, &
        photoelectron_charge, root_before, root_after, displacement_after, response_after, status, message &
        )
      if (status == matching_plane_provider_ok .and. len_trim(message) > 0) then
        write (error_unit, '(a)') trim(message)
        flush (error_unit)
      end if
      if (status == matching_plane_provider_ok) result_packet = [displacement_after, response_after]
    end if
    status_packet = [status]
    call mpi_bcast_i32_array(mpi, status_packet, 0_i32)
    status = status_packet(1)
    if (status /= matching_plane_provider_ok) then
      if (mpi_is_root(mpi)) then
        write (error_unit, '(a)') trim(message)
        flush (error_unit)
      end if
      error stop 128
    end if
    call mpi_bcast_real_dp_array(mpi, result_packet, 0_i32)
    displacement_after = result_packet(1)
    response_after = result_packet(2:7)
  end subroutine solve_matching_implicit_zero_mode

  subroutine solve_matching_implicit_zero_mode_local( &
    provider, displacement_before, displacement_seed, duration, displacement_bounded, &
    displacement_min, displacement_max, &
    displacement_scale, search_direction, feedback_reference, electron_charge, ion_charge, photoelectron_active, &
    photoelectron_charge, root_before, root_after, displacement_after, response_after, status, message &
    )
    type(matching_plane_response_provider_type), intent(inout) :: provider
    logical, intent(in) :: displacement_bounded
    real(dp), intent(in) :: displacement_before, displacement_seed, duration
    real(dp), intent(in) :: displacement_min, displacement_max, displacement_scale
    integer(i32), intent(in) :: search_direction
    real(dp), intent(in) :: feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    logical, intent(in) :: photoelectron_active
    type(matching_plane_zhao_root_seed_type), intent(in) :: root_before
    type(matching_plane_zhao_root_seed_type), intent(out) :: root_after
    real(dp), intent(out) :: displacement_after, response_after(6)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: lower, upper, candidate, step, invalid_candidate, denominator, guard
    real(dp) :: lower_residual, upper_residual, candidate_residual
    real(dp) :: displacement_tolerance, residual_tolerance
    real(dp) :: lower_response(6), upper_response(6), candidate_response(6)
    type(matching_plane_zhao_root_seed_type) :: lower_root, upper_root, candidate_root, evaluation_root
    real(dp) :: current_density
    integer(i32), parameter :: online_expansion_count = 64_i32
    integer(i32), parameter :: online_initial_scan_count = 256_i32
    real(dp), parameter :: online_initial_scan_spacing = 1.0_dp/32.0_dp
    integer(i32) :: iteration, boundary_iteration, rejected_status
    character(len=512) :: evaluation_message
    logical :: bracketed, lower_candidate_brackets, have_valid_point
    logical :: saw_numerical_candidate
    logical :: boundary_failure_numerical

    displacement_after = 0.0_dp
    response_after = 0.0_dp
    root_after = matching_plane_zhao_root_seed_type()
    lower_root = matching_plane_zhao_root_seed_type()
    upper_root = matching_plane_zhao_root_seed_type()
    candidate_root = matching_plane_zhao_root_seed_type()
    evaluation_root = matching_plane_zhao_root_seed_type()
    status = matching_plane_provider_ok
    message = ''
    saw_numerical_candidate = .false.
    displacement_tolerance = 128.0_dp*epsilon(1.0_dp)*max( &
                             displacement_scale, abs(displacement_before), tiny(1.0_dp) &
                             )
    residual_tolerance = displacement_tolerance
    if (.not. displacement_bounded) then
      residual_tolerance = max( &
                           residual_tolerance, &
                           sqrt(epsilon(1.0_dp))*max(displacement_scale, abs(displacement_before)) &
                           )
    end if

    if (displacement_bounded) then
      lower = displacement_min
      upper = displacement_max
      call evaluate_matching_implicit_residual_local( &
        provider, lower, displacement_before, duration, feedback_reference, &
        electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
        root_before, lower_root, lower_residual, lower_response, current_density, status, evaluation_message &
        )
      if (status /= matching_plane_provider_ok) then
        message = 'implicit matching-plane lower endpoint failed: '//trim(evaluation_message)
        return
      end if
      evaluation_root = root_before
      if (root_before%valid .and. lower_root%valid) evaluation_root = lower_root
      call evaluate_matching_implicit_residual_local( &
        provider, upper, displacement_before, duration, feedback_reference, &
        electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
        evaluation_root, upper_root, upper_residual, upper_response, current_density, status, evaluation_message &
        )
      if (status /= matching_plane_provider_ok) then
        message = 'implicit matching-plane upper endpoint failed: '//trim(evaluation_message)
        return
      end if
      bracketed = residuals_bracket_zero(lower_residual, upper_residual)
      if (.not. bracketed) then
        status = matching_plane_provider_no_physical_solution
        message = 'implicit matching-plane zero-mode root is not bracketed by the response table.'
        return
      end if
    else
      lower = displacement_seed
      call evaluate_matching_implicit_residual_local( &
        provider, lower, displacement_before, duration, feedback_reference, &
        electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
        root_before, lower_root, lower_residual, lower_response, current_density, status, evaluation_message &
        )
      bracketed = .false.
      have_valid_point = status == matching_plane_provider_ok
      if (have_valid_point) then
        if (abs(lower_residual) <= residual_tolerance) then
          displacement_after = lower
          response_after = lower_response
          root_after = lower_root
          return
        end if
        ! At the seed, the residual gives a local endpoint correction scale.
        ! Probe that local scale first and only then expand geometrically.
        step = min(displacement_scale, max(abs(lower_residual), displacement_tolerance))
      else
        if ((status /= matching_plane_provider_no_physical_solution .and. &
             status /= matching_plane_provider_numerical_failure) .or. &
            search_direction == 0_i32) then
          message = 'implicit matching-plane Zhao starting point failed: '//trim(evaluation_message)
          return
        end if
        ! Explicit A/B/C can have a narrow certified interval separated from
        ! zero by points where the finite-start Zhao solve is inconclusive.
        ! Scan the natural displacement scale without bracketing across such a
        ! gap; geometric powers can skip the whole Type-A interval.
        saw_numerical_candidate = status == matching_plane_provider_numerical_failure
        status = matching_plane_provider_ok
        step = online_initial_scan_spacing*displacement_scale
        have_valid_point = .false.
        do iteration = 1_i32, online_initial_scan_count
          candidate = real(search_direction*iteration, dp)*step
          evaluation_root = root_before
          if (root_before%valid .and. lower_root%valid) evaluation_root = lower_root
          call evaluate_matching_implicit_residual_local( &
            provider, candidate, displacement_before, duration, feedback_reference, &
            electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
            evaluation_root, candidate_root, candidate_residual, candidate_response, current_density, &
            status, evaluation_message &
            )
          if ((status == matching_plane_provider_no_physical_solution .or. &
               status == matching_plane_provider_numerical_failure) .and. have_valid_point) then
            invalid_candidate = candidate
            rejected_status = status
            call recover_matching_continuation_substep_local( &
              provider, lower, invalid_candidate, rejected_status, displacement_before, duration, feedback_reference, &
              electron_charge, ion_charge, photoelectron_active, photoelectron_charge, evaluation_root, &
              displacement_tolerance, candidate, candidate_root, candidate_residual, candidate_response, &
              current_density, status, evaluation_message &
              )
            if (status /= matching_plane_provider_ok) then
              message = 'implicit matching-plane Zhao initial-scan subdivision failed: '//trim(evaluation_message)
              return
            end if
          end if
          if (status == matching_plane_provider_ok) then
            if (abs(candidate_residual) <= residual_tolerance) then
              displacement_after = candidate
              response_after = candidate_response
              root_after = candidate_root
              return
            end if
            if (have_valid_point .and. residuals_bracket_zero(lower_residual, candidate_residual)) then
              if (candidate < lower) then
                upper = lower
                upper_residual = lower_residual
                upper_response = lower_response
                upper_root = lower_root
                lower = candidate
                lower_residual = candidate_residual
                lower_response = candidate_response
                lower_root = candidate_root
              else
                upper = candidate
                upper_residual = candidate_residual
                upper_response = candidate_response
                upper_root = candidate_root
              end if
              bracketed = .true.
              exit
            end if
            lower = candidate
            lower_residual = candidate_residual
            lower_response = candidate_response
            lower_root = candidate_root
            have_valid_point = .true.
          else if (status == matching_plane_provider_no_physical_solution .or. &
                   status == matching_plane_provider_numerical_failure) then
            saw_numerical_candidate = saw_numerical_candidate .or. &
                                      status == matching_plane_provider_numerical_failure
            have_valid_point = .false.
            status = matching_plane_provider_ok
          else
            message = 'implicit matching-plane Zhao initial signed scan failed: '//trim(evaluation_message)
            return
          end if
        end do
        if (.not. bracketed) then
          if (saw_numerical_candidate) then
            status = matching_plane_provider_numerical_failure
          else
            status = matching_plane_provider_no_physical_solution
          end if
          message = 'implicit matching-plane Zhao root was not bracketed by the signed natural-scale scan.'
          return
        end if
      end if

      do iteration = 1_i32, online_expansion_count
        if (bracketed) exit
        if (have_valid_point) then
          if (lower_residual < 0.0_dp) then
            candidate = lower + step
          else
            candidate = lower - step
          end if
        else
          candidate = real(search_direction, dp)*step
        end if
        if (.not. ieee_is_finite(candidate)) then
          status = matching_plane_provider_numerical_failure
          message = 'implicit matching-plane Zhao bracket expansion overflowed.'
          return
        end if
        evaluation_root = root_before
        if (root_before%valid .and. lower_root%valid) evaluation_root = lower_root
        call evaluate_matching_implicit_residual_local( &
          provider, candidate, displacement_before, duration, feedback_reference, &
          electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
          evaluation_root, candidate_root, candidate_residual, candidate_response, current_density, &
          status, evaluation_message &
          )
        if (status == matching_plane_provider_ok) then
          if (abs(candidate_residual) <= residual_tolerance) then
            displacement_after = candidate
            response_after = candidate_response
            root_after = candidate_root
            return
          end if
          if (.not. have_valid_point) then
            lower = candidate
            lower_residual = candidate_residual
            lower_response = candidate_response
            lower_root = candidate_root
            have_valid_point = .true.
          else if (residuals_bracket_zero(lower_residual, candidate_residual)) then
            if (candidate < lower) then
              upper = lower
              upper_residual = lower_residual
              upper_response = lower_response
              upper_root = lower_root
              lower = candidate
              lower_residual = candidate_residual
              lower_response = candidate_response
              lower_root = candidate_root
            else
              upper = candidate
              upper_residual = candidate_residual
              upper_response = candidate_response
              upper_root = candidate_root
            end if
            bracketed = .true.
            exit
          else
            lower = candidate
            lower_residual = candidate_residual
            lower_response = candidate_response
            lower_root = candidate_root
          end if
        else if ((status == matching_plane_provider_no_physical_solution .or. &
                  status == matching_plane_provider_numerical_failure) .and. have_valid_point) then
          ! Do not skip a root merely because the geometric probe crossed the
          ! branch boundary.  Approach the invalid endpoint from the last valid
          ! point and look for a sign change without extrapolating the response.
          invalid_candidate = candidate
          boundary_failure_numerical = status == matching_plane_provider_numerical_failure
          status = matching_plane_provider_ok
          do boundary_iteration = 1_i32, online_expansion_count
            candidate = 0.5_dp*lower + 0.5_dp*invalid_candidate
            if (abs(candidate - lower) <= displacement_tolerance) exit
            evaluation_root = root_before
            if (root_before%valid .and. lower_root%valid) evaluation_root = lower_root
            call evaluate_matching_implicit_residual_local( &
              provider, candidate, displacement_before, duration, feedback_reference, &
              electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
              evaluation_root, candidate_root, candidate_residual, candidate_response, current_density, &
              status, evaluation_message &
              )
            if (status == matching_plane_provider_ok) then
              if (abs(candidate_residual) <= residual_tolerance) then
                displacement_after = candidate
                response_after = candidate_response
                root_after = candidate_root
                return
              end if
              if (residuals_bracket_zero(lower_residual, candidate_residual)) then
                if (candidate < lower) then
                  upper = lower
                  upper_residual = lower_residual
                  upper_response = lower_response
                  upper_root = lower_root
                  lower = candidate
                  lower_residual = candidate_residual
                  lower_response = candidate_response
                  lower_root = candidate_root
                else
                  upper = candidate
                  upper_residual = candidate_residual
                  upper_response = candidate_response
                  upper_root = candidate_root
                end if
                bracketed = .true.
                exit
              end if
              lower = candidate
              lower_residual = candidate_residual
              lower_response = candidate_response
              lower_root = candidate_root
            else if (status == matching_plane_provider_no_physical_solution .or. &
                     status == matching_plane_provider_numerical_failure) then
              boundary_failure_numerical = boundary_failure_numerical .or. &
                                           status == matching_plane_provider_numerical_failure
              invalid_candidate = candidate
              status = matching_plane_provider_ok
            else
              message = 'implicit matching-plane Zhao branch-boundary search failed: '//trim(evaluation_message)
              return
            end if
          end do
          if (bracketed) exit
          if (boundary_failure_numerical) then
            status = matching_plane_provider_numerical_failure
          else
            status = matching_plane_provider_no_physical_solution
          end if
          message = 'implicit matching-plane Zhao branch ended before the backward-Euler root.'
          return
        else if (status /= matching_plane_provider_no_physical_solution) then
          message = 'implicit matching-plane Zhao bracket expansion failed: '//trim(evaluation_message)
          return
        else
          status = matching_plane_provider_ok
        end if

        if (step > huge(step)/2.0_dp) then
          status = matching_plane_provider_numerical_failure
          message = 'implicit matching-plane Zhao bracket expansion exceeded the numeric range.'
          return
        end if
        step = 2.0_dp*step
      end do
      if (.not. bracketed) then
        status = matching_plane_provider_no_physical_solution
        message = 'implicit matching-plane Zhao root was not bracketed after automatic expansion.'
        return
      end if
    end if

    if (abs(lower_residual) <= residual_tolerance) then
      displacement_after = lower
      response_after = lower_response
      root_after = lower_root
      return
    end if
    if (abs(upper_residual) <= residual_tolerance) then
      displacement_after = upper
      response_after = upper_response
      root_after = upper_root
      return
    end if

    do iteration = 1_i32, online_expansion_count
      candidate = 0.5_dp*lower + 0.5_dp*upper
      if (.not. displacement_bounded) then
        denominator = upper_residual - lower_residual
        if (ieee_is_finite(denominator) .and. denominator /= 0.0_dp) then
          candidate = (lower*upper_residual - upper*lower_residual)/denominator
        end if
        guard = 0.25_dp*(upper - lower)
        if (.not. ieee_is_finite(candidate) .or. candidate <= lower + guard .or. candidate >= upper - guard) then
          candidate = 0.5_dp*lower + 0.5_dp*upper
        end if
      end if
      if (.not. root_before%valid) then
        evaluation_root = root_before
      else if (abs(candidate - lower) <= abs(upper - candidate)) then
        evaluation_root = lower_root
      else
        evaluation_root = upper_root
      end if
      call evaluate_matching_implicit_residual_local( &
        provider, candidate, displacement_before, duration, feedback_reference, &
        electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
        evaluation_root, candidate_root, candidate_residual, candidate_response, current_density, &
        status, evaluation_message &
        )
      if (status == matching_plane_provider_no_physical_solution .or. &
          status == matching_plane_provider_numerical_failure) then
        invalid_candidate = candidate
        rejected_status = status
        if (abs(candidate - lower) <= abs(upper - candidate)) then
          call recover_matching_continuation_substep_local( &
            provider, lower, invalid_candidate, rejected_status, displacement_before, duration, feedback_reference, &
            electron_charge, ion_charge, photoelectron_active, photoelectron_charge, evaluation_root, &
            displacement_tolerance, candidate, candidate_root, candidate_residual, candidate_response, &
            current_density, status, evaluation_message &
            )
        else
          call recover_matching_continuation_substep_local( &
            provider, upper, invalid_candidate, rejected_status, displacement_before, duration, feedback_reference, &
            electron_charge, ion_charge, photoelectron_active, photoelectron_charge, evaluation_root, &
            displacement_tolerance, candidate, candidate_root, candidate_residual, candidate_response, &
            current_density, status, evaluation_message &
            )
        end if
      end if
      if (status /= matching_plane_provider_ok) then
        message = 'implicit matching-plane bracket refinement failed: '//trim(evaluation_message)
        return
      end if
      if (abs(candidate_residual) <= residual_tolerance) then
        displacement_after = candidate
        response_after = candidate_response
        root_after = candidate_root
        return
      end if
      lower_candidate_brackets = residuals_bracket_zero(lower_residual, candidate_residual)
      if (lower_candidate_brackets) then
        upper = candidate
        upper_residual = candidate_residual
        upper_response = candidate_response
        upper_root = candidate_root
      else
        lower = candidate
        lower_residual = candidate_residual
        lower_response = candidate_response
        lower_root = candidate_root
      end if
      if (upper - lower <= displacement_tolerance) exit
    end do
    if (abs(lower_residual) <= abs(upper_residual)) then
      displacement_after = lower
      response_after = lower_response
      root_after = lower_root
    else
      displacement_after = upper
      response_after = upper_response
      root_after = upper_root
    end if
    if (min(abs(lower_residual), abs(upper_residual)) > 8.0_dp*residual_tolerance) then
      write (message, '(a,es12.4,a,es12.4,a,es12.4)') &
        'WARNING: implicit matching-plane bracket refinement accepted the finite best endpoint: residual=', &
        min(abs(lower_residual), abs(upper_residual)), ', tolerance=', 8.0_dp*residual_tolerance, &
        ', bracket_width=', upper - lower
    end if
  end subroutine solve_matching_implicit_zero_mode_local

  subroutine recover_matching_continuation_substep_local( &
    provider, anchor_displacement, rejected_displacement, rejected_status, &
    displacement_before, duration, feedback_reference, &
    electron_charge, ion_charge, photoelectron_active, photoelectron_charge, anchor_root, displacement_tolerance, &
    recovered_displacement, recovered_root, residual, response, current_density, status, message &
    )
    type(matching_plane_response_provider_type), intent(inout) :: provider
    real(dp), intent(in) :: anchor_displacement, rejected_displacement, displacement_before, duration
    integer(i32), intent(in) :: rejected_status
    real(dp), intent(in) :: feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    logical, intent(in) :: photoelectron_active
    type(matching_plane_zhao_root_seed_type), intent(in) :: anchor_root
    real(dp), intent(in) :: displacement_tolerance
    real(dp), intent(out) :: recovered_displacement, residual, response(6), current_density
    type(matching_plane_zhao_root_seed_type), intent(out) :: recovered_root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer(i32), parameter :: subdivision_count = 128_i32
    real(dp) :: invalid_displacement
    integer(i32) :: iteration
    logical :: saw_numerical_failure

    recovered_displacement = anchor_displacement
    recovered_root = matching_plane_zhao_root_seed_type()
    residual = 0.0_dp
    response = 0.0_dp
    current_density = 0.0_dp
    status = rejected_status
    message = ''
    invalid_displacement = rejected_displacement
    saw_numerical_failure = rejected_status == matching_plane_provider_numerical_failure

    do iteration = 1_i32, subdivision_count
      recovered_displacement = 0.5_dp*anchor_displacement + 0.5_dp*invalid_displacement
      if (abs(recovered_displacement - anchor_displacement) <= displacement_tolerance) exit
      call evaluate_matching_implicit_residual_local( &
        provider, recovered_displacement, displacement_before, duration, feedback_reference, &
        electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
        anchor_root, recovered_root, residual, response, current_density, status, message &
        )
      if (status == matching_plane_provider_ok) return
      if (status /= matching_plane_provider_no_physical_solution .and. &
          status /= matching_plane_provider_numerical_failure) return
      saw_numerical_failure = saw_numerical_failure .or. status == matching_plane_provider_numerical_failure
      invalid_displacement = recovered_displacement
    end do

    if (saw_numerical_failure) then
      status = matching_plane_provider_numerical_failure
    else
      status = matching_plane_provider_no_physical_solution
    end if
    message = 'no valid Type-A Zhao response was found within the continuation subdivision tolerance.'
  end subroutine recover_matching_continuation_substep_local

  subroutine evaluate_matching_implicit_residual_local( &
    provider, displacement, displacement_before, duration, feedback_reference, &
    electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
    root_before, root_after, residual, response, current_density, status, message &
    )
    type(matching_plane_response_provider_type), intent(inout) :: provider
    real(dp), intent(in) :: displacement, displacement_before, duration, feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    logical, intent(in) :: photoelectron_active
    type(matching_plane_zhao_root_seed_type), intent(in) :: root_before
    type(matching_plane_zhao_root_seed_type), intent(out) :: root_after
    real(dp), intent(out) :: residual, response(6), current_density
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: input(5), escape_fraction, escape_flux, barrier_energy_ev

    input = [displacement, feedback_reference]
    call provider%evaluate_local( &
      input, response, status, message, continuation_seed=root_before, continuation_candidate=root_after &
      )
    residual = 0.0_dp
    current_density = 0.0_dp
    if (status /= matching_plane_provider_ok) return
    escape_fraction = 0.0_dp
    escape_flux = 0.0_dp
    current_density = electron_charge*response(2) + ion_charge*response(3)
    if (photoelectron_active) then
      barrier_energy_ev = response(1) - response(6)
      if (.not. ieee_is_finite(barrier_energy_ev) .or. barrier_energy_ev < 0.0_dp) then
        status = matching_plane_provider_invalid_argument
        message = 'implicit matching-plane PE barrier or reference energy is invalid.'
        return
      end if
      if (feedback_reference(1) > 0.0_dp) then
        if (.not. ieee_is_finite(feedback_reference(2)) .or. feedback_reference(2) <= 0.0_dp) then
          status = matching_plane_provider_invalid_argument
          message = 'positive implicit matching-plane PE flux requires positive mean energy.'
          return
        end if
        escape_fraction = exp(-barrier_energy_ev/feedback_reference(2))
        escape_flux = feedback_reference(1)*escape_fraction
      end if
      current_density = current_density - photoelectron_charge*escape_flux
    end if
    residual = displacement - displacement_before - duration*current_density
    if (.not. all(ieee_is_finite([escape_fraction, escape_flux, current_density, residual]))) then
      status = matching_plane_provider_numerical_failure
      message = 'implicit matching-plane zero-mode residual is not finite.'
    end if
  end subroutine evaluate_matching_implicit_residual_local

  pure logical function residuals_bracket_zero(first, second) result(bracketed)
    real(dp), intent(in) :: first, second

    bracketed = (first <= 0.0_dp .and. second >= 0.0_dp) .or. &
                (first >= 0.0_dp .and. second <= 0.0_dp)
  end function residuals_bracket_zero

end module bem_matching_plane_implicit
