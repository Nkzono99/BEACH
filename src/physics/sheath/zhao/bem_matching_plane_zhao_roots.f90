!> Zhao の分岐候補の列挙、物理解の選択、accepted endpoint からの継続。
!! 数値解法と物理式は numerics / physics に委譲し、選択順と縮退判定をここに集める。
submodule(bem_matching_plane_zhao) bem_matching_plane_zhao_roots
  implicit none

  real(dp), parameter :: zero_field_tolerance_hat = 1.0e-12_dp
  real(dp), parameter :: root_cluster_tolerance = 1.0e-6_dp
  real(dp), parameter :: energy_tie_tolerance = 1.0e-6_dp
  real(dp), parameter :: continuation_root_jump_limit = 0.25_dp
  real(dp), parameter :: continuation_distance_tie_tolerance = 1.0e-6_dp

contains

  module procedure solve_matching_root

  type(zhao_matching_root_type) :: roots(3*matching_root_seed_count)
  real(dp) :: field_scale, target_field_hat
  integer :: root_count
  logical :: saw_nonphysical, saw_numerical

  root = zhao_matching_root_type()
  field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
  target_field_hat = interface_field_v_m/field_scale
  if (.not. all(ieee_is_finite([field_scale, target_field_hat])) .or. field_scale <= 0.0_dp) then
    status = matching_plane_zhao_numerical_failure
    message = 'matching-plane Zhao field normalization is invalid.'
    return
  end if
  call collect_matching_roots(model, params, target_field_hat, roots, root_count, &
                              saw_nonphysical, saw_numerical, status, message)
  if (status /= matching_plane_zhao_ok) return
  if (saw_numerical) then
    status = matching_plane_zhao_numerical_failure
    message = 'matching-plane Zhao candidate profile could not be certified numerically.'
  else if (root_count == 1) then
    root = roots(1)
    if (trim(root_selection) == 'minimum_energy') then
      call evaluate_root_potential_energy(params, root, status, message)
    end if
  else if (root_count > 1) then
    if (trim(root_selection) == 'minimum_energy') then
      call select_minimum_energy_root(params, roots, root_count, root, status, message)
    else
      status = matching_plane_zhao_ambiguous_solution
      message = 'matching-plane Zhao search found multiple physical roots.'
    end if
  else if (saw_nonphysical) then
    status = matching_plane_zhao_no_physical_solution
    message = 'all located Zhao candidates have nonphysical connecting profiles.'
  else
    status = matching_plane_zhao_numerical_failure
    message = 'matching-plane Zhao branch search did not converge.'
  end if
  end procedure solve_matching_root

  subroutine collect_matching_roots(model, params, target_field_hat, roots, root_count, &
                                    saw_nonphysical, saw_numerical, status, message)
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: target_field_hat
    type(zhao_matching_root_type), intent(out) :: roots(3*matching_root_seed_count)
    integer, intent(out) :: root_count
    logical, intent(out) :: saw_nonphysical, saw_numerical
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=1) :: order(3)
    type(zhao_matching_root_type) :: branch_roots(matching_root_seed_count)
    integer :: branch_count, branch_index, count, index, previous
    logical :: nonphysical, numerical, duplicate, any_compatible

    roots = zhao_matching_root_type()
    root_count = 0
    saw_nonphysical = .false.
    saw_numerical = .false.
    any_compatible = .false.
    call matching_branch_order(model, target_field_hat, order, branch_count, status, message)
    if (status /= matching_plane_zhao_ok) return
    do branch_index = 1, branch_count
      call collect_matching_branch_roots(params, order(branch_index), target_field_hat, branch_roots, &
                                         count, nonphysical, numerical, status, message)
      if (status == matching_plane_zhao_no_physical_solution) cycle
      if (status /= matching_plane_zhao_ok) return
      any_compatible = .true.
      saw_nonphysical = saw_nonphysical .or. nonphysical
      saw_numerical = saw_numerical .or. numerical
      do index = 1, count
        duplicate = .false.
        do previous = 1, root_count
          if (matching_roots_equivalent(params, branch_roots(index), roots(previous))) duplicate = .true.
          ! At zero field a vanishing C endpoint is the analytic flat B state.
          ! The field-squared Newton tolerance gives O(sqrt(tol)) voltage uncertainty.
          if (abs(target_field_hat) <= zero_field_tolerance_hat .and. &
              (branch_roots(index)%phi0_v == 0.0_dp .or. roots(previous)%phi0_v == 0.0_dp)) then
            if (max(abs(branch_roots(index)%phi0_v), abs(roots(previous)%phi0_v)) &
                < 1.0e-4_dp*params%t_phe_ev) then
              duplicate = .true.
              if (branch_roots(index)%phi0_v == 0.0_dp) roots(previous) = branch_roots(index)
            end if
          end if
        end do
        if (duplicate) cycle
        root_count = root_count + 1
        roots(root_count) = branch_roots(index)
      end do
    end do
    status = matching_plane_zhao_ok
    message = ''
    if (.not. any_compatible) then
      status = matching_plane_zhao_no_physical_solution
      message = 'requested Zhao branches are incompatible with the prescribed field or upstream state.'
    end if
  end subroutine collect_matching_roots

  subroutine matching_branch_order(model, target_field_hat, order, count, status, message)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: target_field_hat
    character(len=1), intent(out) :: order(3)
    integer, intent(out) :: count
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    order = ' '
    count = 0
    status = matching_plane_zhao_ok
    message = ''
    select case (trim(model))
    case ('a')
      order(1) = 'A'
      count = 1
    case ('b')
      order(1) = 'B'
      count = 1
    case ('c')
      order(1) = 'C'
      count = 1
    case ('auto')
      if (abs(target_field_hat) <= zero_field_tolerance_hat) then
        order = ['B', 'C', 'A']
      else if (target_field_hat > 0.0_dp) then
        order = ['A', 'B', 'C']
      else
        order = ['C', 'A', 'B']
      end if
      count = 3
    case default
      status = matching_plane_zhao_invalid_argument
      message = 'unknown matching-plane Zhao branch.'
    end select
  end subroutine matching_branch_order

  subroutine collect_matching_branch_roots( &
    params, branch, target_field_hat, unique_roots, unique_count, &
    saw_nonphysical_profile, saw_numerical_profile_failure, status, message &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat
    type(zhao_matching_root_type), intent(out) :: unique_roots(matching_root_seed_count)
    integer, intent(out) :: unique_count
    logical, intent(out) :: saw_nonphysical_profile, saw_numerical_profile_failure
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: guesses(3, matching_root_seed_count), y(3), norm, ion, free, reflected, photo, captured
    type(zhao_matching_root_type) :: candidate_root, candidate_roots(matching_root_seed_count)
    integer :: guess_count, guess_index, iterations, root_index
    integer(i32) :: profile_status
    logical :: success, compatible, duplicate_root, candidate_valid(matching_root_seed_count)
    logical :: guess_nonphysical_profile(matching_root_seed_count), guess_numerical_profile_failure(matching_root_seed_count)
    character(len=512) :: profile_message

    unique_roots = zhao_matching_root_type()
    unique_count = 0
    saw_nonphysical_profile = .false.
    saw_numerical_profile_failure = .false.
    status = matching_plane_zhao_no_physical_solution
    message = ''
    compatible = (branch == 'C' .and. target_field_hat <= 0.0_dp) .or. &
                 (branch == 'B' .and. target_field_hat >= 0.0_dp) .or. &
                 (branch == 'A' .and. target_field_hat > 0.0_dp)
    if ((branch == 'A' .or. branch == 'C') .and. params%u > 0.0_dp) then
      message = 'positive electron drift admits no A/C connection to a strict field-free infinity.'
      return
    end if
    if (.not. compatible) then
      message = 'Zhao branch and matching-plane field signs are incompatible.'
      return
    end if

    if (branch == 'B' .and. abs(target_field_hat) <= zero_field_tolerance_hat) then
      call evaluate_zhao_density_hat(params, 'B', 'monotonic', 0.0_dp, 0.0_dp, 0.0_dp, &
                                     1.0_dp, ion, free, reflected, photo, captured)
      if (ion <= photo + captured .or. free + reflected <= 0.0_dp) then
        message = 'zero-field Zhao-B state has no positive ambient electron density.'
        return
      end if
      unique_count = 1
      unique_roots(1)%branch = 'B'
      unique_roots(1)%ambient_electron_density_m3 = &
        params%n_phe_ref_m3*(ion - photo - captured)/(free + reflected)
      unique_roots(1)%residual_norm = 0.0_dp
      unique_roots(1)%minimum_field_squared_hat = 0.0_dp
      unique_roots(1)%potential_energy_j_m2 = 0.0_dp
      status = matching_plane_zhao_ok
      return
    end if
    call make_matching_branch_guesses(params, branch, target_field_hat, guesses, guess_count)
    candidate_roots = zhao_matching_root_type()
    candidate_valid = .false.
    guess_nonphysical_profile = .false.
    guess_numerical_profile_failure = .false.
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp shared(params,branch,target_field_hat,guesses,guess_count,candidate_roots,candidate_valid) &
    !$omp shared(guess_nonphysical_profile,guess_numerical_profile_failure) &
    !$omp private(guess_index,y,norm,iterations,success,candidate_root,profile_status,profile_message)
    do guess_index = 1, guess_count
      call newton_matching_branch( &
        params, branch, target_field_hat, guesses(:, guess_index), y, norm, iterations, success &
        )
      if (.not. success) cycle
      candidate_root = zhao_matching_root_type()
      candidate_root%branch = branch
      call decode_matching_unknowns( &
        params, branch, y, candidate_root%phi0_v, candidate_root%phi_m_v, &
        candidate_root%ambient_electron_density_m3, success &
        )
      if (.not. success) cycle
      candidate_root%residual_norm = norm
      candidate_root%nonlinear_iterations = int(iterations, i32)
      call validate_matching_root_profile( &
        params, candidate_root, target_field_hat, profile_status, profile_message &
        )
      if (profile_status == matching_plane_zhao_no_physical_solution) then
        guess_nonphysical_profile(guess_index) = .true.
        cycle
      else if (profile_status /= matching_plane_zhao_ok) then
        guess_numerical_profile_failure(guess_index) = .true.
        cycle
      end if

      candidate_roots(guess_index) = candidate_root
      candidate_valid(guess_index) = .true.
    end do
    !$omp end parallel do

    ! Candidate completion order must not affect root clustering or tie breaking.
    do guess_index = 1, guess_count
      saw_nonphysical_profile = saw_nonphysical_profile .or. guess_nonphysical_profile(guess_index)
      saw_numerical_profile_failure = &
        saw_numerical_profile_failure .or. guess_numerical_profile_failure(guess_index)
      if (.not. candidate_valid(guess_index)) cycle
      candidate_root = candidate_roots(guess_index)
      duplicate_root = .false.
      do root_index = 1, unique_count
        if (.not. matching_roots_equivalent(params, candidate_root, unique_roots(root_index))) cycle
        duplicate_root = .true.
        if (candidate_root%residual_norm < unique_roots(root_index)%residual_norm) then
          unique_roots(root_index) = candidate_root
        end if
        exit
      end do
      if (.not. duplicate_root) then
        unique_count = unique_count + 1
        unique_roots(unique_count) = candidate_root
      end if
    end do
    status = matching_plane_zhao_ok
  end subroutine collect_matching_branch_roots

  module procedure solve_matching_continuation

  type(zhao_matching_root_type) :: candidate_root, roots(3*matching_root_seed_count)
  real(dp) :: field_scale, target_field_hat, seed_y(3), candidate_y(3), norm
  real(dp) :: candidate_jump, nearest_jump, second_nearest_jump
  integer :: iterations, root_count, root_index, nearest_index
  integer(i32) :: profile_status
  logical :: seed_valid, success, candidate_valid, branch_allowed
  logical :: saw_nonphysical_profile, saw_numerical_profile_failure
  character(len=512) :: profile_message

  root = zhao_matching_root_type()
  fallback_used = .false.
  root_jump = huge(1.0_dp)
  status = matching_plane_zhao_invalid_argument
  message = ''
  if (.not. all(ieee_is_finite([seed%phi0_v, seed%phi_m_v, seed%ambient_electron_density_m3])) .or. &
      seed%ambient_electron_density_m3 <= 0.0_dp .or. index('ABC', seed%branch) == 0) then
    message = 'matching-plane Zhao continuation seed is invalid.'
    return
  end if
  field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
  target_field_hat = interface_field_v_m/field_scale
  if (.not. all(ieee_is_finite([field_scale, target_field_hat])) .or. field_scale <= 0.0_dp) then
    status = matching_plane_zhao_numerical_failure
    message = 'matching-plane Zhao field normalization is invalid.'
    return
  end if
  branch_allowed = trim(model) == 'auto' .or. trim(model) == lower_ascii(seed%branch)
  branch_allowed = branch_allowed .and. &
                   ((seed%branch == 'A' .and. target_field_hat > 0.0_dp .and. params%u <= 0.0_dp) .or. &
                    (seed%branch == 'B' .and. target_field_hat > 0.0_dp) .or. &
                    (seed%branch == 'C' .and. target_field_hat <= 0.0_dp .and. params%u <= 0.0_dp))
  call encode_matching_unknowns(params, seed%branch, seed%phi0_v, seed%phi_m_v, &
                                seed%ambient_electron_density_m3, seed_y, seed_valid)
  if (seed_valid .and. branch_allowed) then
    call newton_matching_branch(params, seed%branch, target_field_hat, seed_y, &
                                candidate_y, norm, iterations, success)
    if (success) then
      candidate_root = zhao_matching_root_type()
      candidate_root%branch = seed%branch
      call decode_matching_unknowns(params, seed%branch, candidate_y, candidate_root%phi0_v, &
                                    candidate_root%phi_m_v, candidate_root%ambient_electron_density_m3, candidate_valid)
      if (candidate_valid) then
        candidate_root%residual_norm = norm
        candidate_root%nonlinear_iterations = int(iterations, i32)
        call validate_matching_root_profile(params, candidate_root, target_field_hat, profile_status, profile_message)
        if (profile_status == matching_plane_zhao_ok) then
          root_jump = matching_seed_distance(params, seed, candidate_root)
          if (root_jump <= continuation_root_jump_limit) then
            root = candidate_root
            status = matching_plane_zhao_ok
            return
          end if
        end if
      end if
    end if
  end if

  fallback_used = .true.
  call collect_matching_roots(model, params, target_field_hat, roots, root_count, &
                              saw_nonphysical_profile, saw_numerical_profile_failure, status, message)
  if (status /= matching_plane_zhao_ok) return
  if (saw_numerical_profile_failure) then
    status = matching_plane_zhao_numerical_failure
    message = 'matching-plane Zhao continuation could not certify a candidate profile numerically.'
    return
  end if
  nearest_jump = huge(1.0_dp)
  second_nearest_jump = huge(1.0_dp)
  nearest_index = 0
  do root_index = 1, root_count
    candidate_jump = matching_seed_distance(params, seed, roots(root_index))
    if (candidate_jump < nearest_jump) then
      second_nearest_jump = nearest_jump
      nearest_jump = candidate_jump
      nearest_index = root_index
    else if (candidate_jump < second_nearest_jump) then
      second_nearest_jump = candidate_jump
    end if
  end do
  root_jump = nearest_jump
  if (nearest_index > 0 .and. &
      abs(second_nearest_jump - nearest_jump) <= &
      continuation_distance_tie_tolerance*max(1.0_dp, nearest_jump)) then
    status = matching_plane_zhao_ambiguous_solution
    message = 'matching-plane Zhao continuation found indistinguishable nearest roots.'
  else if (nearest_index > 0) then
    root = roots(nearest_index)
    status = matching_plane_zhao_ok
    message = ''
  else if (saw_nonphysical_profile) then
    status = matching_plane_zhao_no_physical_solution
    message = 'all located Zhao continuation candidates have nonphysical connecting profiles.'
  else
    status = matching_plane_zhao_numerical_failure
    message = 'matching-plane Zhao continuation multistart did not converge.'
  end if
  end procedure solve_matching_continuation

  pure real(dp) function matching_seed_distance(params, seed, root) result(distance)
    type(zhao_params_type), intent(in) :: params
    type(matching_plane_zhao_root_seed_type), intent(in) :: seed
    type(zhao_matching_root_type), intent(in) :: root
    distance = max(abs(root%phi0_v - seed%phi0_v)/params%t_phe_ev, &
                   abs(min(root%phi_m_v, 0.0_dp) - min(seed%phi_m_v, 0.0_dp))/params%t_phe_ev, &
                   abs(log(root%ambient_electron_density_m3/seed%ambient_electron_density_m3)))
  end function matching_seed_distance

  module procedure select_minimum_energy_root

  type(zhao_matching_root_type) :: candidates(size(roots))
  real(dp) :: energy_scale
  integer :: candidate_index, best_index

  root = zhao_matching_root_type()
  candidates = roots
  status = matching_plane_zhao_numerical_failure
  message = ''
  if (root_count < 1 .or. root_count > size(roots)) then
    message = 'matching-plane Zhao minimum-energy selection received an invalid candidate count.'
    return
  end if
  do candidate_index = 1, root_count
    call evaluate_root_potential_energy(params, candidates(candidate_index), status, message)
    if (status /= matching_plane_zhao_ok) return
  end do
  best_index = 1
  do candidate_index = 2, root_count
    if (candidates(candidate_index)%potential_energy_j_m2 < &
        candidates(best_index)%potential_energy_j_m2) best_index = candidate_index
  end do
  do candidate_index = 1, root_count
    if (candidate_index == best_index) cycle
    energy_scale = max( &
                   abs(candidates(best_index)%potential_energy_j_m2), &
                   abs(candidates(candidate_index)%potential_energy_j_m2), tiny(1.0_dp) &
                   )
    if (abs( &
        candidates(candidate_index)%potential_energy_j_m2 - &
        candidates(best_index)%potential_energy_j_m2 &
        ) <= energy_tie_tolerance*energy_scale) then
      status = matching_plane_zhao_ambiguous_solution
      message = 'matching-plane Zhao minimum-energy candidates are numerically tied.'
      return
    end if
  end do
  root = candidates(best_index)
  status = matching_plane_zhao_ok
  message = ''
  end procedure select_minimum_energy_root

  pure logical function matching_roots_equivalent(params, first, second) result(equivalent)
    type(zhao_params_type), intent(in) :: params
    type(zhao_matching_root_type), intent(in) :: first, second

    real(dp) :: first_phi0_hat, second_phi0_hat, first_phi_m_hat, second_phi_m_hat
    real(dp) :: log_density_ratio

    equivalent = .false.
    if (first%branch /= second%branch) return
    if (min(first%ambient_electron_density_m3, second%ambient_electron_density_m3) <= 0.0_dp) return
    first_phi0_hat = first%phi0_v/params%t_phe_ev
    second_phi0_hat = second%phi0_v/params%t_phe_ev
    first_phi_m_hat = first%phi_m_v/params%t_phe_ev
    second_phi_m_hat = second%phi_m_v/params%t_phe_ev
    log_density_ratio = log(first%ambient_electron_density_m3/second%ambient_electron_density_m3)
    if (.not. all(ieee_is_finite([ &
                                 first_phi0_hat, second_phi0_hat, first_phi_m_hat, second_phi_m_hat, &
                                 log_density_ratio &
                                 ]))) return
    equivalent = abs(first_phi0_hat - second_phi0_hat) <= &
      root_cluster_tolerance*max(1.0_dp, abs(first_phi0_hat), abs(second_phi0_hat)) .and. &
      abs(first_phi_m_hat - second_phi_m_hat) <= &
      root_cluster_tolerance*max(1.0_dp, abs(first_phi_m_hat), abs(second_phi_m_hat)) .and. &
      abs(log_density_ratio) <= root_cluster_tolerance
  end function matching_roots_equivalent

end submodule bem_matching_plane_zhao_roots
