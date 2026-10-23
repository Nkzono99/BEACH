!> 電荷差分の確定反映、表面電流補正、バッチ電荷台帳を実装する。
submodule(bem_simulator) bem_simulator_charge
  use, intrinsic :: iso_fortran_env, only: error_unit
  use bem_charge_ledger, only: checked_accumulate_charge, finite_charge_sum
  use bem_runtime_modes, only: surface_charge_closure_id, surface_charge_closure_fixed_current, &
                               surface_charge_closure_neutral_return, closure_mode_none, closure_mode_target, &
                               closure_mode_tracked_source
  use bem_surface_closure_contract, only: species_closure_mode
  use bem_simulator_workspace, only: fixed_charge_absorbed, fixed_charge_emission, fixed_charge_escaped, &
                                     fixed_charge_barrier_return, neutral_charge_emitted, neutral_charge_absorbed, &
                                     neutral_charge_unresolved, neutral_charge_escaped, neutral_charge_soft, &
                                     neutral_charge_invalid, neutral_count_escaped, neutral_count_soft, &
                                     neutral_count_invalid
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
  ! adaptive candidate は prepare 側で既に等電位化済み。
  if (conductor%active() .and. .not. workspace%charge_candidate_ready) then
    call conductor%relax(mesh, sim, snapshot)
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
  rel = finite_nonnegative_ratio(norm_dq, max(norm_q, sim%q_floor))
  workspace%charge_candidate_ready = .false.
  end procedure commit_batch_charge

  module procedure prepare_adaptive_charge_candidate
  real(dp), allocatable :: committed_charge(:)
  workspace%dq = sum(workspace%dq_thread, dim=2) + sum(workspace%photo_emission_dq, dim=2)
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
  call validate_finite_charge_addition( &
    mesh%q_elem, workspace%dq, 'adaptive candidate surface-charge update' &
    )
  workspace%candidate_charge = mesh%q_elem + workspace%dq
  if (.not. all(ieee_is_finite(workspace%candidate_charge))) then
    error stop 'adaptive candidate charge is not finite.'
  end if
  if (conductor%active()) then
    ! 電位ステップの判定は再配分後の候補電荷で行う。snapshot は確定電荷の状態へ戻す。
    committed_charge = mesh%q_elem
    mesh%q_elem = workspace%candidate_charge
    call conductor%relax(mesh, sim, snapshot)
    workspace%candidate_charge = mesh%q_elem
    mesh%q_elem = committed_charge
    call snapshot%refresh(mesh)
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
  logical :: neutral_return_species(app%n_particle_species)

  n = app%n_particle_species
  ! 種ごとの closure をバッチの先頭で一度だけ解決し、粒子ループでは文字列を比べない。
  do species_idx = 1_i32, n
    neutral_return_species(species_idx) = &
      surface_charge_closure_id(app%particle_species(species_idx)%surface_charge_closure) == &
      surface_charge_closure_neutral_return
  end do
  has_neutral_return = .false.
  do species_idx = 1_i32, n
    if (.not. app%particle_species(species_idx)%enabled) cycle
    if (neutral_return_species(species_idx)) then
      has_neutral_return = .true.
    end if
  end do
  if (.not. has_neutral_return) return
  if (.not. app%sim%use_box) error stop 'neutral_return requires a finite box.'

  workspace%neutral_return_charge_values = 0.0_dp
  workspace%neutral_return_terminal_counts = 0_i64
  do i = 1_i32, pcls_batch%n
    species_idx = pcls_batch%species_id(i)
    if (.not. neutral_return_species(species_idx)) cycle
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    if (.not. ieee_is_finite(macro_charge) .or. macro_charge >= 0.0_dp .or. i > fresh_particle_count) then
      call count_invalid_terminal()
      cycle
    end if
    workspace%neutral_return_charge_values(species_idx, neutral_charge_emitted) = &
      workspace%neutral_return_charge_values(species_idx, neutral_charge_emitted) + macro_charge
    terminal_count = merge(1_i32, 0_i32, workspace%absorbed_flag(i)) + &
                     merge(1_i32, 0_i32, workspace%escaped_boundary_flag(i)) + &
                     merge(1_i32, 0_i32, workspace%soft_discarded_boundary_flag(i)) + &
                     merge(1_i32, 0_i32, pcls_batch%alive(i))
    if (terminal_count /= 1_i32) then
      call count_invalid_terminal()
    else if (workspace%absorbed_flag(i)) then
      elem_idx = workspace%absorbed_element(i)
      if (elem_idx >= 1_i32 .and. elem_idx <= size(workspace%dq_thread, 1)) then
        workspace%neutral_return_charge_values(species_idx, neutral_charge_absorbed) = &
          workspace%neutral_return_charge_values(species_idx, neutral_charge_absorbed) + macro_charge
      else
        call count_invalid_terminal()
      end if
    else if (workspace%escaped_boundary_flag(i)) then
      workspace%neutral_return_charge_values(species_idx, neutral_charge_escaped) = &
        workspace%neutral_return_charge_values(species_idx, neutral_charge_escaped) + macro_charge
      workspace%neutral_return_terminal_counts(species_idx, neutral_count_escaped) = &
        workspace%neutral_return_terminal_counts(species_idx, neutral_count_escaped) + 1_i64
    else if (workspace%soft_discarded_boundary_flag(i)) then
      workspace%neutral_return_charge_values(species_idx, neutral_charge_soft) = &
        workspace%neutral_return_charge_values(species_idx, neutral_charge_soft) + macro_charge
      workspace%neutral_return_terminal_counts(species_idx, neutral_count_soft) = &
        workspace%neutral_return_terminal_counts(species_idx, neutral_count_soft) + 1_i64
    else
      workspace%neutral_return_charge_values(species_idx, neutral_charge_unresolved) = &
        workspace%neutral_return_charge_values(species_idx, neutral_charge_unresolved) + macro_charge
    end if
  end do

  call allreduce_sum_real_matrix(mpi, workspace%neutral_return_charge_values)
  call allreduce_sum_i64_matrix(mpi, workspace%neutral_return_terminal_counts)
  workspace%neutral_return_emitted_charge = workspace%neutral_return_charge_values(:, neutral_charge_emitted)
  workspace%neutral_return_absorbed_charge = workspace%neutral_return_charge_values(:, neutral_charge_absorbed)
  workspace%neutral_return_unresolved_charge = workspace%neutral_return_charge_values(:, neutral_charge_unresolved)

  do species_idx = 1_i32, n
    if (.not. neutral_return_species(species_idx)) cycle
    emitted_charge = workspace%neutral_return_emitted_charge(species_idx)
    absorbed_charge = workspace%neutral_return_absorbed_charge(species_idx)
    unresolved_charge = workspace%neutral_return_unresolved_charge(species_idx)
    escaped_charge = workspace%neutral_return_charge_values(species_idx, neutral_charge_escaped)
    soft_charge = workspace%neutral_return_charge_values(species_idx, neutral_charge_soft)
    invalid_charge = workspace%neutral_return_charge_values(species_idx, neutral_charge_invalid)
    escaped_count = workspace%neutral_return_terminal_counts(species_idx, neutral_count_escaped)
    soft_count = workspace%neutral_return_terminal_counts(species_idx, neutral_count_soft)
    invalid_count = workspace%neutral_return_terminal_counts(species_idx, neutral_count_invalid)
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
    if (.not. neutral_return_species(species_idx)) cycle
    if (.not. workspace%absorbed_flag(i)) cycle
    elem_idx = workspace%absorbed_element(i)
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                       (workspace%neutral_return_weight_scale(species_idx) - 1.0_dp)*macro_charge
  end do
contains

  !> 終端が一意でない、または電荷・吸収要素が不正な粒子を invalid として数える。
  subroutine count_invalid_terminal()
    workspace%neutral_return_charge_values(species_idx, neutral_charge_invalid) = &
      workspace%neutral_return_charge_values(species_idx, neutral_charge_invalid) + macro_charge
    workspace%neutral_return_terminal_counts(species_idx, neutral_count_invalid) = &
      workspace%neutral_return_terminal_counts(species_idx, neutral_count_invalid) + 1_i64
  end subroutine count_invalid_terminal
  end procedure apply_neutral_return_surface_closure

  module procedure apply_fixed_surface_current_closure
  integer(i32) :: species_idx, n, mode
  real(dp) :: partner_charge, reflected_escape, reflected_return, folded_return_charge, charge_sign
  logical :: tracked_source
  logical :: fixed_current_species(app%n_particle_species)

  n = app%n_particle_species
  do species_idx = 1_i32, n
    fixed_current_species(species_idx) = &
      surface_charge_closure_id(app%particle_species(species_idx)%surface_charge_closure) == &
      surface_charge_closure_fixed_current
  end do
  workspace%fixed_escape_closed_charge = 0.0_dp
  workspace%fixed_escape_unclosed_charge = 0.0_dp
  call tally_raw_charge()

  do species_idx = 1_i32, n
    mode = species_closure_mode(current_model, species_idx, fixed_current_species(species_idx))
    if (mode == closure_mode_none) cycle
    tracked_source = mode == closure_mode_tracked_source
    charge_sign = sign(1.0_dp, app%particle_species(species_idx)%q_particle)
    partner_charge = 0.0_dp
    reflected_escape = 0.0_dp
    reflected_return = 0.0_dp
    if (mode /= closure_mode_target) call attribute_tracked_escape()
    folded_return_charge = 0.0_dp
    if (app%particle_species(species_idx)%has_target_absorbed_current_a .or. &
        current_model%has_absorbed_target(species_idx)) call close_absorbed_channel()
    if (app%particle_species(species_idx)%has_target_emission_current_a .or. &
        current_model%has_emission_target(species_idx)) call close_emission_channel()
    if (current_model%has_escape_target(species_idx)) call close_escape_channel()
  end do
contains

  !> 種ごとの raw 電荷（吸収・放出反作用・脱出・外部障壁で戻った電荷）を数え、rank 合計する。
  subroutine tally_raw_charge()
    integer(i32) :: i, s
    real(dp) :: macro_charge

    workspace%fixed_current_charge_values = 0.0_dp
    do i = 1_i32, fresh_particle_count
      s = pcls_batch%species_id(i)
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (workspace%absorbed_flag(i)) then
        workspace%fixed_current_charge_values(s, fixed_charge_absorbed) = &
          workspace%fixed_current_charge_values(s, fixed_charge_absorbed) + macro_charge
      else if (workspace%escaped_boundary_flag(i)) then
        workspace%fixed_current_charge_values(s, fixed_charge_escaped) = &
          workspace%fixed_current_charge_values(s, fixed_charge_escaped) + macro_charge
      end if
    end do
    do s = 1_i32, n
      workspace%fixed_current_charge_values(s, fixed_charge_emission) = sum(workspace%photo_emission_dq(:, s))
      ! 外部障壁で戻った粒子の電荷。H を外向きに通過したが上流へ抜けられず、セルへ戻された分。
      workspace%fixed_current_charge_values(s, fixed_charge_barrier_return) = &
        app%particle_species(s)%q_particle*sum(workspace%matching_plane_moments_thread(3, s, :))
    end do
    call allreduce_sum_real_matrix(mpi, workspace%fixed_current_charge_values)
  end subroutine tally_raw_charge

  !> 外部根の脱出束が partner species（追跡する SEE）や反射して H を抜けた上流電子と合算の値なら、それらが追跡で
  !! 脱出した電荷（帯電符号の大きさ）を根の脱出束まで差し引き、この species の脱出目標と帰還目標をその分だけ縮める。
  !! 超過は収支残差として残す。
  subroutine attribute_tracked_escape()
    integer(i32) :: partner_idx, reflected_idx
    real(dp) :: partner_escape_charge, total_escape_charge

    partner_escape_charge = 0.0_dp
    partner_idx = 0_i32
    reflected_idx = 0_i32
    if (allocated(current_model%escape_partner_species)) partner_idx = current_model%escape_partner_species(species_idx)
    if (allocated(current_model%reflected_electron_species)) then
      reflected_idx = current_model%reflected_electron_species(species_idx)
    end if
    if (partner_idx > 0_i32) then
      partner_escape_charge = max(0.0_dp, -workspace%fixed_current_charge_values(partner_idx, fixed_charge_escaped))
    end if
    if (reflected_idx > 0_i32) then
      reflected_escape = max(0.0_dp, -workspace%fixed_current_charge_values(reflected_idx, fixed_charge_escaped))
      reflected_return = max(0.0_dp, -workspace%fixed_current_charge_values(reflected_idx, fixed_charge_barrier_return))
    end if
    if ((partner_idx > 0_i32 .or. reflected_idx > 0_i32) .and. current_model%has_escape_target(species_idx)) then
      total_escape_charge = max(0.0_dp, -checked_fixed_target_charge( &
                                current_model%escaped_particle_current_a(species_idx), app%sim%batch_duration, &
                                'fixed_current escape'))
      partner_charge = min(partner_escape_charge + reflected_escape, total_escape_charge)
      workspace%fixed_escape_closed_charge(species_idx) = partner_charge
      workspace%fixed_escape_unclosed_charge(species_idx) = partner_escape_charge + reflected_escape - partner_charge
    end if
  end subroutine attribute_tracked_escape

  !> 吸収（帰還）channel の重みを目標へ合わせる。追跡した源の種（tracked_source）は raw channel が空でも止めず、
  !! 帰還目標を放出の反作用電荷へ畳む（folded_return_charge）。
  subroutine close_absorbed_channel()
    integer(i32) :: i, elem_idx
    real(dp) :: raw_charge, target_current, target_charge, weight_scale, correction, charge_tolerance, macro_charge

    raw_charge = workspace%fixed_current_charge_values(species_idx, fixed_charge_absorbed)
    if (current_model%has_absorbed_target(species_idx)) then
      target_current = current_model%absorbed_current_a(species_idx)
    else
      target_current = app%particle_species(species_idx)%target_absorbed_current_a
    end if
    if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
      error stop 'fixed_current absorbed raw/current/duration value is not finite.'
    end if
    target_charge = checked_fixed_target_charge(target_current, app%sim%batch_duration, 'fixed_current absorbed')
    if (tracked_source) then
      ! 追跡した源の種では、外部障壁で戻った反射電子は電子の raw 吸収に含まれるので、帰還目標から外す。
      target_charge = target_charge + reflected_return
    else
      target_charge = target_charge - partner_charge
    end if
    if (allocated(current_model%absorbed_target_credits_tracked_escape)) then
      if (current_model%absorbed_target_credits_tracked_escape(species_idx)) then
        ! 外部障壁を越えて脱出した反射電子は壁に届かないので、根の流入から差し引く（吸収目標は帯電符号で小さくなる）。
        target_charge = target_charge - workspace%fixed_current_charge_values(species_idx, fixed_charge_escaped)
      end if
    end if
    ! 差し引きで目標の符号が粒子電荷と逆になれば、吸収なしまで縮める。
    if (target_charge*charge_sign < 0.0_dp) target_charge = 0.0_dp
    charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*max(abs(raw_charge), abs(target_charge), tiny(1.0_dp))
    if (abs(raw_charge) <= charge_tolerance) then
      if (abs(target_charge) > charge_tolerance) then
        if (.not. tracked_source) then
          error stop 'fixed_current cannot map a nonzero absorbed target onto an empty raw channel.'
        end if
        ! 生の帰還がないバッチでは帰還目標を放出の反作用電荷へ畳む（放出 − 帰還 = 脱出は変わらない）。
        folded_return_charge = target_charge
        target_charge = 0.0_dp
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
    workspace%fixed_channel_active(1, species_idx) = .true.
    workspace%fixed_current_correction(species_idx) = &
      workspace%fixed_current_correction(species_idx) + correction
    do i = 1_i32, fresh_particle_count
      if (pcls_batch%species_id(i) /= species_idx .or. .not. workspace%absorbed_flag(i)) cycle
      elem_idx = workspace%absorbed_element(i)
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                         (weight_scale - 1.0_dp)*macro_charge
    end do
  end subroutine close_absorbed_channel

  !> 放出 channel の反作用電荷を目標へ合わせる。追跡した源の種は、放出がないバッチでは目標を載せず未閉として記録する。
  subroutine close_emission_channel()
    real(dp) :: raw_charge, target_current, target_charge, weight_scale, correction, charge_tolerance

    raw_charge = workspace%fixed_current_charge_values(species_idx, fixed_charge_emission)
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
                    ) + folded_return_charge
    if (tracked_source) then
      ! 根の外側源には反射電子の外向き束も入っているので、追跡した反射電子の脱出と帰還を放出目標から外す。
      ! 放出の反作用電荷は粒子電荷と逆符号。差し引きで符号が入れ替われば、放出なしまで縮める。
      target_charge = target_charge - (reflected_escape + reflected_return)
      if (target_charge*charge_sign > 0.0_dp) then
        workspace%fixed_escape_unclosed_charge(species_idx) = &
          workspace%fixed_escape_unclosed_charge(species_idx) + abs(target_charge)
        target_charge = 0.0_dp
      end if
    end if
    charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*max(abs(raw_charge), abs(target_charge), tiny(1.0_dp))
    if (abs(raw_charge) <= charge_tolerance) then
      if (abs(target_charge) > charge_tolerance) then
        if (.not. tracked_source) then
          error stop 'fixed_current cannot map a nonzero emission target onto an empty raw channel.'
        end if
        ! 放出がないバッチには目標を載せられない。根の脱出分が収支に残るので未閉として記録する。
        workspace%fixed_escape_unclosed_charge(species_idx) = abs(target_charge)
        target_charge = 0.0_dp
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
    workspace%fixed_channel_active(2, species_idx) = .true.
    workspace%fixed_current_correction(species_idx) = &
      workspace%fixed_current_correction(species_idx) + correction
    workspace%photo_emission_dq(:, species_idx) = &
      weight_scale*workspace%photo_emission_dq(:, species_idx)
  end subroutine close_emission_channel

  !> 脱出 channel の目標（partner・反射電子に充てた分を除いた根の脱出）と raw との差を補正として記録する。
  subroutine close_escape_channel()
    real(dp) :: raw_charge, target_current, target_charge

    raw_charge = workspace%fixed_current_charge_values(species_idx, fixed_charge_escaped)
    target_current = current_model%escaped_particle_current_a(species_idx)
    if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
      error stop 'fixed_current escape raw/current/duration value is not finite.'
    end if
    target_charge = checked_fixed_target_charge( &
                    target_current, app%sim%batch_duration, 'fixed_current escape' &
                    ) + partner_charge
    if (target_charge /= 0.0_dp .and. &
        sign(1.0_dp, target_charge) /= sign(1.0_dp, app%particle_species(species_idx)%q_particle)) then
      error stop 'fixed_current escape target sign must match the escaped particle charge.'
    end if
    workspace%fixed_escape_target_charge(species_idx) = target_charge
    workspace%fixed_escape_correction(species_idx) = target_charge - raw_charge
    workspace%fixed_channel_active(3, species_idx) = .true.
    if (tracked_source) then
      workspace%fixed_escape_closed_charge(species_idx) = &
        max(0.0_dp, abs(target_charge) - workspace%fixed_escape_unclosed_charge(species_idx))
    end if
    if (.not. ieee_is_finite(workspace%fixed_escape_correction(species_idx))) then
      error stop 'fixed_current escape correction is not finite.'
    end if
  end subroutine close_escape_channel
  end procedure apply_fixed_surface_current_closure

  !> species × channel の 2 次元配列を rank 合計する。列優先の並びは旧 1 次元配列と同じなので MPI の buffer も同じ。
  subroutine allreduce_sum_real_matrix(mpi, values)
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(inout) :: values(:, :)
    real(dp) :: flat(size(values))

    flat = reshape(values, [size(values)])
    call mpi_allreduce_sum_real_dp_array(mpi, flat)
    values = reshape(flat, shape(values))
  end subroutine allreduce_sum_real_matrix

  !> 64bit 整数の species × kind 配列を rank 合計する。
  subroutine allreduce_sum_i64_matrix(mpi, values)
    type(mpi_context), intent(in) :: mpi
    integer(i64), intent(inout) :: values(:, :)
    integer(i64) :: flat(size(values))

    flat = reshape(values, [size(values)])
    call mpi_allreduce_sum_i64_array(mpi, flat)
    values = reshape(flat, shape(values))
  end subroutine allreduce_sum_i64_matrix

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

  module procedure record_batch_initial_charge
  integer(i32) :: i, species_idx
  real(dp) :: macro_charge
  do i = 1_i32, fresh_particle_count
    species_idx = pcls_batch%species_id(i)
    macro_charge = checked_macro_charge(pcls_batch%q(i), pcls_batch%w(i), 'initial batch charge ledger')
    if (pcls_batch%source_element(i) > 0_i32) then
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
