!> 電荷差分の確定反映、表面電流補正、バッチ電荷台帳を実装する。
submodule(bem_simulator) bem_simulator_charge
  use, intrinsic :: iso_fortran_env, only: error_unit
  use bem_charge_ledger, only: checked_accumulate_charge, finite_charge_sum
  implicit none
  real(dp), parameter :: neutral_return_max_unresolved_fraction = 0.05_dp
contains

  module procedure commit_batch_charge
  real(dp) :: norm_dq, norm_q
  if (.not. all(ieee_is_finite(mesh%q_elem))) then
    error stop 'committed surface charge is not finite before batch update.'
  end if
  workspace%q_before = mesh%q_elem
  if (workspace%charge_candidate_ready) then
    mesh%q_elem = workspace%candidate_charge
  else
    workspace%dq = sum(workspace%dq_thread, dim=2) + sum(workspace%photo_emission_dq, dim=2)
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
    call validate_finite_charge_addition( &
      mesh%q_elem, workspace%dq, 'batch surface-charge update' &
      )
    mesh%q_elem = mesh%q_elem + workspace%dq
  end if
  if (.not. all(ieee_is_finite(mesh%q_elem))) then
    error stop 'batch surface-charge update produced a non-finite charge.'
  end if
  call apply_surface_model_charge_relaxation(mesh, external_e, field_bc_mode=field_bc_mode)
  if (.not. all(ieee_is_finite(mesh%q_elem))) then
    error stop 'surface-model charge relaxation produced a non-finite charge.'
  end if
  call validate_finite_charge_addition( &
    mesh%q_elem, -workspace%q_before, 'batch surface-charge difference' &
    )
  workspace%dq = mesh%q_elem - workspace%q_before
  if (.not. all(ieee_is_finite(workspace%dq))) then
    error stop 'batch surface-charge difference is not finite.'
  end if
  norm_dq = stable_l2_norm(workspace%dq)
  norm_q = stable_l2_norm(mesh%q_elem)
  rel = finite_nonnegative_ratio(norm_dq, max(norm_q, q_floor))
  workspace%charge_candidate_ready = .false.
  end procedure commit_batch_charge

  module procedure prepare_adaptive_charge_candidate
  workspace%dq = sum(workspace%dq_thread, dim=2) + sum(workspace%photo_emission_dq, dim=2)
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
  call validate_finite_charge_addition( &
    mesh%q_elem, workspace%dq, 'adaptive candidate surface-charge update' &
    )
  workspace%candidate_charge = mesh%q_elem + workspace%dq
  if (.not. all(ieee_is_finite(workspace%candidate_charge))) then
    error stop 'adaptive candidate charge is not finite.'
  end if
  workspace%charge_candidate_ready = .true.
  end procedure prepare_adaptive_charge_candidate

  !> Finiteな二項の加算が表現可能範囲を越えないことを、演算前に検証する。
  subroutine validate_finite_charge_addition(base, increment, operation)
    real(dp), intent(in) :: base(:), increment(:)
    character(len=*), intent(in) :: operation
    integer :: i

    if (.not. all(ieee_is_finite(base)) .or. .not. all(ieee_is_finite(increment))) then
      error stop trim(operation)//' contains a non-finite operand.'
    end if
    do i = 1, size(base)
      if (increment(i) > 0.0_dp .and. base(i) > huge(base(i)) - increment(i)) then
        error stop trim(operation)//' overflowed.'
      end if
      if (increment(i) < 0.0_dp .and. base(i) < -huge(base(i)) - increment(i)) then
        error stop trim(operation)//' overflowed.'
      end if
    end do
  end subroutine validate_finite_charge_addition

  !> 二乗の中間overflowを避けて有限vectorのL2 normを返す。
  pure real(dp) function stable_l2_norm(values) result(norm)
    real(dp), intent(in) :: values(:)
    real(dp) :: scale, unit_norm

    if (size(values) == 0) then
      norm = 0.0_dp
      return
    end if
    scale = maxval(abs(values))
    if (scale == 0.0_dp) then
      norm = 0.0_dp
      return
    end if
    unit_norm = sqrt(sum((values/scale)*(values/scale)))
    if (scale > huge(norm)/unit_norm) then
      norm = huge(norm)
    else
      norm = scale*unit_norm
      if (.not. ieee_is_finite(norm)) norm = huge(norm)
    end if
  end function stable_l2_norm

  !> 非負の有限比をoverflow時は最大有限値へ飽和させる。
  pure real(dp) function finite_nonnegative_ratio(numerator, denominator) result(ratio)
    real(dp), intent(in) :: numerator, denominator

    if (numerator == 0.0_dp) then
      ratio = 0.0_dp
    else if (denominator < 1.0_dp .and. numerator > huge(ratio)*denominator) then
      ratio = huge(ratio)
    else
      ratio = numerator/denominator
      if (.not. ieee_is_finite(ratio)) ratio = huge(ratio)
    end if
  end function finite_nonnegative_ratio

  module procedure apply_neutral_return_surface_closure
  integer(i32) :: i, species_idx, elem_idx, n, terminal_count
  integer(i64) :: escaped_count, soft_count, invalid_count
  real(dp) :: macro_charge, emitted_charge, absorbed_charge, unresolved_charge
  real(dp) :: escaped_charge, soft_charge, invalid_charge
  real(dp) :: charge_scale, charge_tolerance
  real(dp) :: weight_scale, correction_charge, unresolved_fraction
  logical :: has_neutral_return

  n = app%n_particle_species
  has_neutral_return = .false.
  do species_idx = 1_i32, n
    if (.not. app%particle_species(species_idx)%enabled) cycle
    if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == 'neutral_return') then
      has_neutral_return = .true.
    end if
  end do
  if (.not. has_neutral_return) return
  if (.not. app%sim%use_box) error stop 'neutral_return requires a finite box.'

  workspace%neutral_return_charge_values = 0.0_dp
  workspace%neutral_return_terminal_counts = 0_i64
  do i = 1_i32, pcls_batch%n
    species_idx = pcls_batch%species_id(i)
    if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    if (.not. ieee_is_finite(macro_charge) .or. macro_charge >= 0.0_dp .or. i > fresh_particle_count) then
      workspace%neutral_return_charge_values(5*n + species_idx) = &
        workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
      workspace%neutral_return_terminal_counts(2*n + species_idx) = &
        workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
      cycle
    end if
    workspace%neutral_return_charge_values(species_idx) = &
      workspace%neutral_return_charge_values(species_idx) + macro_charge
    terminal_count = merge(1_i32, 0_i32, workspace%absorbed_flag(i)) + &
                     merge(1_i32, 0_i32, workspace%escaped_boundary_flag(i)) + &
                     merge(1_i32, 0_i32, workspace%soft_discarded_boundary_flag(i)) + &
                     merge(1_i32, 0_i32, pcls_batch%alive(i))
    if (terminal_count /= 1_i32) then
      workspace%neutral_return_charge_values(5*n + species_idx) = &
        workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
      workspace%neutral_return_terminal_counts(2*n + species_idx) = &
        workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
    else if (workspace%absorbed_flag(i)) then
      elem_idx = workspace%absorbed_element(i)
      if (elem_idx >= 1_i32 .and. elem_idx <= size(workspace%dq_thread, 1)) then
        workspace%neutral_return_charge_values(n + species_idx) = &
          workspace%neutral_return_charge_values(n + species_idx) + macro_charge
      else
        workspace%neutral_return_charge_values(5*n + species_idx) = &
          workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
      end if
    else if (workspace%escaped_boundary_flag(i)) then
      workspace%neutral_return_charge_values(3*n + species_idx) = &
        workspace%neutral_return_charge_values(3*n + species_idx) + macro_charge
      workspace%neutral_return_terminal_counts(species_idx) = &
        workspace%neutral_return_terminal_counts(species_idx) + 1_i64
    else if (workspace%soft_discarded_boundary_flag(i)) then
      workspace%neutral_return_charge_values(4*n + species_idx) = &
        workspace%neutral_return_charge_values(4*n + species_idx) + macro_charge
      workspace%neutral_return_terminal_counts(n + species_idx) = &
        workspace%neutral_return_terminal_counts(n + species_idx) + 1_i64
    else
      workspace%neutral_return_charge_values(2*n + species_idx) = &
        workspace%neutral_return_charge_values(2*n + species_idx) + macro_charge
    end if
  end do

  call mpi_allreduce_sum_real_dp_array(mpi, workspace%neutral_return_charge_values)
  call mpi_allreduce_sum_i64_array(mpi, workspace%neutral_return_terminal_counts)
  workspace%neutral_return_emitted_charge = workspace%neutral_return_charge_values(1:n)
  workspace%neutral_return_absorbed_charge = workspace%neutral_return_charge_values(n + 1:2*n)
  workspace%neutral_return_unresolved_charge = workspace%neutral_return_charge_values(2*n + 1:3*n)

  do species_idx = 1_i32, n
    if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
    emitted_charge = workspace%neutral_return_emitted_charge(species_idx)
    absorbed_charge = workspace%neutral_return_absorbed_charge(species_idx)
    unresolved_charge = workspace%neutral_return_unresolved_charge(species_idx)
    escaped_charge = workspace%neutral_return_charge_values(3*n + species_idx)
    soft_charge = workspace%neutral_return_charge_values(4*n + species_idx)
    invalid_charge = workspace%neutral_return_charge_values(5*n + species_idx)
    escaped_count = workspace%neutral_return_terminal_counts(species_idx)
    soft_count = workspace%neutral_return_terminal_counts(n + species_idx)
    invalid_count = workspace%neutral_return_terminal_counts(2*n + species_idx)
    if (escaped_count > 0_i64 .or. soft_count > 0_i64 .or. invalid_count > 0_i64) then
      write (error_unit, '(a,i0,3(a,i0),3(a,es13.5))') &
        'neutral_return has unsupported terminal outcome for species ', species_idx, &
        ': escaped=', escaped_count, ' soft=', soft_count, ' invalid=', invalid_count, &
        ' escaped_C=', escaped_charge, ' soft_C=', soft_charge, ' invalid_C=', invalid_charge
      error stop 'neutral_return terminal outcome is unsupported.'
    end if
    charge_scale = max(abs(emitted_charge), abs(absorbed_charge), abs(unresolved_charge), tiny(1.0_dp))
    charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*charge_scale
    if (abs(emitted_charge) <= charge_tolerance) cycle
    if (emitted_charge >= 0.0_dp .or. absorbed_charge >= -charge_tolerance) then
      error stop 'neutral_return charge signs are invalid.'
    end if
    weight_scale = emitted_charge/absorbed_charge
    correction_charge = emitted_charge - absorbed_charge
    unresolved_fraction = unresolved_charge/emitted_charge
    if (.not. all(ieee_is_finite([weight_scale, correction_charge, unresolved_fraction])) .or. &
        unresolved_fraction > neutral_return_max_unresolved_fraction + sqrt(epsilon(1.0_dp))) then
      error stop 'neutral_return unresolved fraction exceeds the applicability limit.'
    end if
    workspace%neutral_return_weight_scale(species_idx) = max(1.0_dp, weight_scale)
    workspace%neutral_return_correction(species_idx) = correction_charge
    workspace%neutral_return_unresolved_fraction(species_idx) = max(0.0_dp, unresolved_fraction)
  end do

  do i = 1_i32, fresh_particle_count
    species_idx = pcls_batch%species_id(i)
    if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
    if (.not. workspace%absorbed_flag(i)) cycle
    elem_idx = workspace%absorbed_element(i)
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                       (workspace%neutral_return_weight_scale(species_idx) - 1.0_dp)*macro_charge
  end do
  end procedure apply_neutral_return_surface_closure

  module procedure apply_fixed_surface_current_closure
  integer(i32) :: i, species_idx, elem_idx, n
  real(dp) :: macro_charge, raw_charge, target_current, target_charge, weight_scale, correction
  real(dp) :: charge_tolerance
  logical :: has_fixed_current

  n = app%n_particle_species
  has_fixed_current = .false.
  do species_idx = 1_i32, n
    if (.not. app%particle_species(species_idx)%enabled) cycle
    if (fixed_current_species_active(app, current_model, species_idx)) then
      has_fixed_current = .true.
    end if
  end do
  if (.not. has_fixed_current) return

  workspace%fixed_current_charge_values = 0.0_dp
  do i = 1_i32, fresh_particle_count
    species_idx = pcls_batch%species_id(i)
    if (.not. fixed_current_species_active(app, current_model, species_idx)) cycle
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    if (workspace%absorbed_flag(i)) then
      workspace%fixed_current_charge_values(species_idx) = &
        workspace%fixed_current_charge_values(species_idx) + macro_charge
    else if (workspace%escaped_boundary_flag(i)) then
      workspace%fixed_current_charge_values(2*n + species_idx) = &
        workspace%fixed_current_charge_values(2*n + species_idx) + macro_charge
    end if
  end do
  do species_idx = 1_i32, n
    if (.not. fixed_current_species_active(app, current_model, species_idx)) cycle
    workspace%fixed_current_charge_values(n + species_idx) = sum(workspace%photo_emission_dq(:, species_idx))
  end do
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%fixed_current_charge_values)

  do species_idx = 1_i32, n
    if (.not. fixed_current_species_active(app, current_model, species_idx)) cycle
    if (app%particle_species(species_idx)%has_target_absorbed_current_a .or. &
        current_model%has_absorbed_target(species_idx)) then
      raw_charge = workspace%fixed_current_charge_values(species_idx)
      if (current_model%has_absorbed_target(species_idx)) then
        target_current = current_model%absorbed_current_a(species_idx)
      else
        target_current = app%particle_species(species_idx)%target_absorbed_current_a
      end if
      if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
        error stop 'fixed_current absorbed raw/current/duration value is not finite.'
      end if
      target_charge = checked_fixed_target_charge( &
                      target_current, app%sim%batch_duration, 'fixed_current absorbed' &
                      )
      charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*max(abs(raw_charge), abs(target_charge), tiny(1.0_dp))
      if (abs(raw_charge) <= charge_tolerance) then
        if (abs(target_charge) > charge_tolerance) then
          error stop 'fixed_current cannot map a nonzero absorbed target onto an empty raw channel.'
        end if
        weight_scale = 1.0_dp
        correction = 0.0_dp
      else
        weight_scale = target_charge/raw_charge
        correction = target_charge - raw_charge
      end if
      if (.not. ieee_is_finite(weight_scale) .or. weight_scale < 0.0_dp) then
        error stop 'fixed_current absorbed scale is invalid.'
      end if
      workspace%fixed_absorbed_target_charge(species_idx) = target_charge
      workspace%fixed_absorbed_weight_scale(species_idx) = weight_scale
      workspace%fixed_current_correction(species_idx) = &
        workspace%fixed_current_correction(species_idx) + correction
      do i = 1_i32, fresh_particle_count
        if (pcls_batch%species_id(i) /= species_idx .or. .not. workspace%absorbed_flag(i)) cycle
        elem_idx = workspace%absorbed_element(i)
        macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
        workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                           (weight_scale - 1.0_dp)*macro_charge
      end do
    end if

    if (app%particle_species(species_idx)%has_target_emission_current_a .or. &
        current_model%has_emission_target(species_idx)) then
      raw_charge = workspace%fixed_current_charge_values(n + species_idx)
      if (current_model%has_emission_target(species_idx)) then
        target_current = current_model%emission_current_a(species_idx)
      else
        target_current = app%particle_species(species_idx)%target_emission_current_a
      end if
      if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
        error stop 'fixed_current emission raw/current/duration value is not finite.'
      end if
      target_charge = checked_fixed_target_charge( &
                      target_current, app%sim%batch_duration, 'fixed_current emission' &
                      )
      charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*max(abs(raw_charge), abs(target_charge), tiny(1.0_dp))
      if (abs(raw_charge) <= charge_tolerance) then
        if (abs(target_charge) > charge_tolerance) then
          error stop 'fixed_current cannot map a nonzero emission target onto an empty raw channel.'
        end if
        weight_scale = 1.0_dp
        correction = 0.0_dp
      else
        weight_scale = target_charge/raw_charge
        correction = target_charge - raw_charge
      end if
      if (.not. ieee_is_finite(weight_scale) .or. weight_scale < 0.0_dp) then
        error stop 'fixed_current emission scale is invalid.'
      end if
      workspace%fixed_emission_target_charge(species_idx) = target_charge
      workspace%fixed_emission_weight_scale(species_idx) = weight_scale
      workspace%fixed_current_correction(species_idx) = &
        workspace%fixed_current_correction(species_idx) + correction
      workspace%photo_emission_dq(:, species_idx) = &
        weight_scale*workspace%photo_emission_dq(:, species_idx)
    end if

    if (current_model%has_escape_target(species_idx)) then
      raw_charge = workspace%fixed_current_charge_values(2*n + species_idx)
      target_current = current_model%escaped_particle_current_a(species_idx)
      if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
        error stop 'fixed_current escape raw/current/duration value is not finite.'
      end if
      target_charge = checked_fixed_target_charge( &
                      target_current, app%sim%batch_duration, 'fixed_current escape' &
                      )
      if (target_charge /= 0.0_dp .and. &
          sign(1.0_dp, target_charge) /= sign(1.0_dp, app%particle_species(species_idx)%q_particle)) then
        error stop 'fixed_current escape target sign must match the escaped particle charge.'
      end if
      workspace%fixed_escape_target_charge(species_idx) = target_charge
      workspace%fixed_escape_correction(species_idx) = target_charge - raw_charge
      if (.not. ieee_is_finite(workspace%fixed_escape_correction(species_idx))) then
        error stop 'fixed_current escape correction is not finite.'
      end if
    end if
  end do
  end procedure apply_fixed_surface_current_closure

  !> finiteな電流とdurationの積を、丸め境界を含めて有限かつ非zero underflowなしに変換する。
  real(dp) function checked_fixed_target_charge(target_current, duration, context) result(target_charge)
    real(dp), intent(in) :: target_current, duration
    character(len=*), intent(in) :: context

    if (duration > 1.0_dp .and. abs(target_current) > huge(target_charge)/duration) then
      error stop trim(context)//' target charge overflowed for this batch duration.'
    end if
    target_charge = target_current*duration
    if (.not. ieee_is_finite(target_charge)) then
      error stop trim(context)//' target charge overflowed for this batch duration.'
    end if
    if (target_current /= 0.0_dp .and. target_charge == 0.0_dp) then
      error stop trim(context)//' target charge underflowed for this batch duration.'
    end if
  end function checked_fixed_target_charge

  logical function fixed_current_species_active(app, current_model, species_idx) result(active)
    type(app_config), intent(in) :: app
    type(surface_closure_contract_type), intent(in) :: current_model
    integer(i32), intent(in) :: species_idx

    active = trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == 'fixed_current' .or. &
             current_model%has_absorbed_target(species_idx) .or. &
             current_model%has_emission_target(species_idx) .or. &
             current_model%has_escape_target(species_idx)
  end function fixed_current_species_active

  module procedure record_batch_initial_charge
  integer(i32) :: i, species_idx
  real(dp) :: macro_charge
  do i = 1_i32, fresh_particle_count
    species_idx = pcls_batch%species_id(i)
    macro_charge = checked_macro_charge(pcls_batch%q(i), pcls_batch%w(i), 'initial batch charge ledger')
    if (trim(lower_ascii(app%particle_species(species_idx)%source_mode)) == 'photo_raycast') then
      call checked_accumulate_charge( &
        ledger%emitted_from_surface(species_idx), macro_charge, 'local emitted charge ledger' &
        )
      ledger%emitted_count(species_idx) = ledger%emitted_count(species_idx) + 1_i64
    else
      call checked_accumulate_charge( &
        ledger%injected_from_remote(species_idx), macro_charge, 'local injected charge ledger' &
        )
      ledger%injected_count(species_idx) = ledger%injected_count(species_idx) + 1_i64
    end if
  end do
  end procedure record_batch_initial_charge

  module procedure record_batch_outcome_charge
  integer(i32) :: i, species_idx
  real(dp) :: macro_charge
  do i = 1_i32, pcls_batch%n
    species_idx = pcls_batch%species_id(i)
    macro_charge = checked_macro_charge(pcls_batch%q(i), pcls_batch%w(i), 'outcome batch charge ledger')
    if (absorbed_flag(i)) then
      call checked_accumulate_charge( &
        ledger%absorbed_on_surface(species_idx), macro_charge, 'local absorbed charge ledger' &
        )
      ledger%absorbed_count(species_idx) = ledger%absorbed_count(species_idx) + 1_i64
    else if (escaped_boundary_flag(i)) then
      call checked_accumulate_charge( &
        ledger%escaped_to_infinity(species_idx), macro_charge, 'local escaped charge ledger' &
        )
      ledger%escaped_count(species_idx) = ledger%escaped_count(species_idx) + 1_i64
    else if (soft_discarded_boundary_flag(i) .or. pcls_batch%alive(i)) then
      call checked_accumulate_charge( &
        ledger%discarded_unresolved(species_idx), macro_charge, 'local discarded charge ledger' &
        )
      ledger%discarded_unresolved_count(species_idx) = ledger%discarded_unresolved_count(species_idx) + 1_i64
    end if
  end do
  end procedure record_batch_outcome_charge

  module procedure reduce_charge_ledger_fluxes
  integer(i32) :: n
  n = ledger%nspecies
  workspace%ledger_charge_values = [ &
                                   ledger%injected_from_remote, ledger%emitted_from_surface, ledger%absorbed_on_surface, &
                                   ledger%escaped_to_infinity, ledger%discarded_unresolved &
                                   ]
  if (.not. all(ieee_is_finite(workspace%ledger_charge_values))) then
    error stop 'local batch charge ledger contains non-finite fluxes before MPI reduction.'
  end if
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%ledger_charge_values)
  if (.not. all(ieee_is_finite(workspace%ledger_charge_values))) then
    error stop 'global batch charge ledger overflowed during MPI reduction.'
  end if
  ledger%injected_from_remote = workspace%ledger_charge_values(1:n)
  ledger%emitted_from_surface = workspace%ledger_charge_values(n + 1:2*n)
  ledger%absorbed_on_surface = workspace%ledger_charge_values(2*n + 1:3*n)
  ledger%escaped_to_infinity = workspace%ledger_charge_values(3*n + 1:4*n)
  ledger%discarded_unresolved = workspace%ledger_charge_values(4*n + 1:5*n)
  workspace%ledger_count_values = [ &
                                  ledger%injected_count, ledger%emitted_count, ledger%absorbed_count, ledger%escaped_count, &
                                  ledger%discarded_unresolved_count &
                                  ]
  if (any(workspace%ledger_count_values < 0_i64)) then
    error stop 'local batch charge ledger contains invalid counts before MPI reduction.'
  end if
  call mpi_allreduce_sum_i64_array(mpi, workspace%ledger_count_values)
  if (any(workspace%ledger_count_values < 0_i64)) then
    error stop 'global batch charge ledger count overflowed during MPI reduction.'
  end if
  ledger%injected_count = workspace%ledger_count_values(1:n)
  ledger%emitted_count = workspace%ledger_count_values(n + 1:2*n)
  ledger%absorbed_count = workspace%ledger_count_values(2*n + 1:3*n)
  ledger%escaped_count = workspace%ledger_count_values(3*n + 1:4*n)
  ledger%discarded_unresolved_count = workspace%ledger_count_values(4*n + 1:5*n)
  end procedure reduce_charge_ledger_fluxes

  !> finiteな粒子電荷とweightから、overflow/zero-underflowを拒否してmacro chargeを返す。
  real(dp) function checked_macro_charge(particle_charge, particle_weight, context) result(macro_charge)
    real(dp), intent(in) :: particle_charge, particle_weight
    character(len=*), intent(in) :: context

    if (.not. ieee_is_finite(particle_charge) .or. .not. ieee_is_finite(particle_weight)) then
      error stop trim(context)//' has a non-finite particle charge or weight.'
    end if
    if (abs(particle_weight) > 1.0_dp .and. abs(particle_charge) > huge(macro_charge)/abs(particle_weight)) then
      error stop trim(context)//' macro charge overflowed.'
    end if
    macro_charge = particle_charge*particle_weight
    if (.not. ieee_is_finite(macro_charge)) error stop trim(context)//' macro charge is not finite.'
    if (particle_charge /= 0.0_dp .and. particle_weight /= 0.0_dp .and. macro_charge == 0.0_dp) then
      error stop trim(context)//' macro charge underflowed to zero.'
    end if
  end function checked_macro_charge

end submodule bem_simulator_charge
