!> チェックポイントの統計・電荷・台帳の読み込みと形式検証。
submodule(bem_restart) bem_restart_records
  use bem_kinds, only: dp, i64
  use bem_matching_plane_summary, only: matching_plane_summary_key_count, read_matching_plane_summary_key, &
                                        require_matching_plane_summary_keys, read_matching_summary_real
  implicit none
contains

  !> Read one complete numeric CSV record before list-directed conversion. Null fields, slash termination,
  !! repetition syntax and cross-record reads must never reuse a previous checkpoint value.
  module procedure read_checkpoint_csv_record
  integer :: length, first, last, column, comma

  read (unit, '(a)', advance='no', size=length, iostat=ios) line
  if (ios == iostat_end .and. length == 0) return
  if (ios /= iostat_eor) then
    ios = 1
    return
  end if
  ios = 1
  first = 1
  do column = 1, columns
    comma = index(line(first:length), ',')
    if (column < columns) then
      if (comma == 0) return
      last = first + comma - 2
    else
      if (comma /= 0) return
      last = length
    end if
    if (.not. is_decimal_real_token(trim(adjustl(line(first:last))))) return
    first = last + 2
  end do
  ios = 0
  end procedure read_checkpoint_csv_record

  !> `summary.txt` を読み込み、必須キーの存在と要素数整合を検証する。
  !! 欠落キーやメッシュ要素数不一致は再開不能として停止する。
  !! @param[in] path `summary.txt` のファイルパス。
  !! @param[in] expected_nelem 現在メッシュの要素数（整合性検証に使用）。
  !! @param[out] stats 復元した統計値。
  module procedure load_summary_file

  integer :: u, ios, pos
  integer(i32) :: mesh_nelem, saved_world_size, summary_schema_version
  ! The outer photoelectron spectrum is one long line of bin fluxes.
  character(len=8192) :: line
  character(len=64) :: key
  character(len=8192) :: value
  logical :: found_mesh, found_processed, found_absorbed, found_escaped
  logical :: found_batches, found_rel, found_world_size
  logical :: found_matching(matching_plane_summary_key_count)

  stats = sim_stats()
  mesh_nelem = -1_i32
  saved_world_size = 1_i32
  summary_schema_version = -1_i32
  found_mesh = .false.
  found_processed = .false.
  found_absorbed = .false.
  found_escaped = .false.
  found_batches = .false.
  found_rel = .false.
  found_world_size = .false.
  found_matching = .false.

  open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'Failed to open summary.txt for resume.'

  do
    read (u, '(A)', iostat=ios) line
    if (ios /= 0) exit
    line = trim(adjustl(line))
    if (len_trim(line) == 0) cycle
    pos = index(line, '=')
    if (pos <= 0) cycle

    key = trim(adjustl(line(:pos - 1)))
    value = trim(adjustl(line(pos + 1:)))

    select case (trim(key))
    case ('checkpoint_schema_version')
      read (value, *) summary_schema_version
    case ('mesh_nelem')
      read (value, *) mesh_nelem
      found_mesh = .true.
    case ('mpi_world_size')
      read (value, *) saved_world_size
      found_world_size = .true.
    case ('processed_particles')
      read (value, *) stats%processed_particles
      found_processed = .true.
    case ('absorbed')
      read (value, *) stats%absorbed
      found_absorbed = .true.
    case ('escaped')
      read (value, *) stats%escaped
      found_escaped = .true.
    case ('batches')
      read (value, *) stats%batches
      found_batches = .true.
    case ('escaped_boundary')
      read (value, *) stats%escaped_boundary
    case ('survived_max_step')
      read (value, *) stats%survived_max_step
    case ('multiple_box_events_retry_attempted')
      read (value, *) stats%multiple_box_events_retry_attempted
    case ('multiple_box_events_retry_resolved')
      read (value, *) stats%multiple_box_events_retry_resolved
    case ('adaptive_substep_halved_steps')
      read (value, *) stats%adaptive_substep_halved_steps
    case ('adaptive_substeps')
      read (value, *) stats%adaptive_substeps
    case ('adaptive_substep_depth_limited_steps')
      read (value, *) stats%adaptive_substep_depth_limited_steps
    case ('adaptive_grown_steps')
      read (value, *) stats%adaptive_grown_steps
    case ('particle_pool_claims')
      read (value, *) stats%particle_pool_claims
    case ('particle_pool_claim_seconds')
      read (value, *) stats%particle_pool_claim_seconds
    case ('multiple_box_events_soft_discarded')
      read (value, *) stats%multiple_box_events_soft_discarded
    case ('multiple_box_events_soft_discarded_abs_charge_C')
      read (value, *) stats%multiple_box_events_soft_discarded_abs_charge
    case ('last_rel_change')
      read (value, *) stats%last_rel_change
      found_rel = .true.
    case ('simulated_time_s')
      read (value, *) stats%simulated_time
    case ('adaptive_nonzero_mode_rejected_trials')
      read (value, *) stats%adaptive_nonzero_mode%rejected_trials
    case ('adaptive_nonzero_mode_last_batch_duration_s')
      read (value, *) stats%adaptive_nonzero_mode%last_batch_duration
    case ('adaptive_nonzero_mode_last_potential_step_V')
      read (value, *) stats%adaptive_nonzero_mode%last_potential_step
    case ('adaptive_nonzero_mode_omp_threads')
      read (value, *) stats%adaptive_nonzero_mode%omp_threads
    case ('secondary_emission_hits')
      read (value, *) stats%secondary_emission_hits
    case ('secondary_emission_count')
      read (value, *) stats%secondary_emission_count
    case ('secondary_generation_capped_count')
      read (value, *) stats%secondary_generation_capped_count
    case ('secondary_emitted_charge_C')
      call read_matching_summary_real(value, key, stats%secondary_emitted_charge_c)
    case ('secondary_yield_sum')
      call read_matching_summary_real(value, key, stats%secondary_yield_sum)
    case default
      call read_matching_plane_summary_key(key, value, stats%matching_plane, found_matching)
    end select
  end do
  close (u)

  if (.not. (found_mesh .and. found_processed .and. found_absorbed .and. &
             found_escaped .and. found_batches .and. found_rel)) then
    error stop 'Resume checkpoint summary is missing required keys.'
  end if
  if (mesh_nelem /= expected_nelem) then
    error stop 'Resume checkpoint mesh element count does not match current mesh.'
  end if
  if (summary_schema_version >= 9_i32) call require_matching_plane_summary_keys(found_matching)
  call validate_summary_stats(stats)
  if (present(expected_world_size)) then
    if (.not. found_world_size .and. expected_world_size > 1_i32) then
      error stop 'Resume checkpoint summary is missing mpi_world_size.'
    end if
    if (max(1_i32, expected_world_size) /= saved_world_size) then
      error stop 'Resume checkpoint mpi_world_size does not match current MPI world size.'
    end if
  end if
  end procedure load_summary_file

  !> `charges.csv` を読み込み、各要素の電荷をメッシュへ復元する。
  !! 行重複や要素数不足を検出した場合は停止する。
  !! @param[in] path `charges.csv` のファイルパス。
  !! @param[inout] mesh 要素電荷 `q_elem` を復元値で上書きするメッシュ。
  module procedure load_charge_file

  integer :: u, ios
  integer(i32) :: elem_idx, n_loaded
  real(dp) :: charge
  character(len=512) :: header, record
  logical, allocatable :: seen(:)

  if (.not. allocated(mesh%q_elem)) error stop 'Mesh charges are not allocated.'

  allocate (seen(mesh%nelem))
  seen = .false.
  mesh%q_elem = 0.0d0
  n_loaded = 0_i32

  open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'Failed to open charges.csv for resume.'

  read (u, '(A)', iostat=ios) header
  if (ios /= 0) error stop 'Failed to read charges.csv header.'

  do
    call read_checkpoint_csv_record(u, record, 2, ios)
    if (ios < 0) exit
    if (ios == 0) read (record, *, iostat=ios) elem_idx, charge
    if (ios /= 0) error stop 'Failed to parse charges.csv during resume.'
    if (elem_idx < 1_i32 .or. elem_idx > mesh%nelem) then
      error stop 'Resume checkpoint charges.csv has an invalid element index.'
    end if
    if (seen(elem_idx)) then
      error stop 'Resume checkpoint charges.csv contains duplicate element rows.'
    end if
    if (.not. ieee_is_finite(charge)) then
      error stop 'Resume checkpoint charges.csv contains non-finite charge values.'
    end if
    seen(elem_idx) = .true.
    mesh%q_elem(elem_idx) = charge
    n_loaded = n_loaded + 1_i32
  end do
  close (u)

  if (n_loaded /= mesh%nelem) then
    error stop 'Resume checkpoint charges.csv does not match the current mesh.'
  end if
  end procedure load_charge_file

  !> summary stock と species 別 CSV flux/count から累積 charge ledger を復元する。
  module procedure load_charge_ledger_checkpoint
  integer :: u, ios, pos, columns
  integer(i32) :: nspecies, batch_count, row_batch, species_idx, loaded
  integer(i64) :: count_values(5)
  real(dp) :: charge_values(18), fixed_charge_values(13), legacy_charge_values(8), stock_values(6)
  character(len=512) :: line
  character(len=2048) :: header
  character(len=2048) :: record
  character(len=96) :: key
  character(len=256) :: value
  logical :: found_nspecies, found_batch, found_stocks(6), has_fixed_current_columns, has_current_budget_columns
  logical, allocatable :: seen(:)

  nspecies = 0_i32
  batch_count = 0_i32
  found_nspecies = .false.
  found_batch = .false.
  found_stocks = .false.
  stock_values = 0.0_dp
  open (newunit=u, file=trim(summary_path), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'Failed to open summary.txt for charge ledger resume.'
  do
    read (u, '(A)', iostat=ios) line
    if (ios /= 0) exit
    pos = index(line, '=')
    if (pos <= 0) cycle
    key = trim(adjustl(line(:pos - 1)))
    value = trim(adjustl(line(pos + 1:)))
    select case (trim(key))
    case ('charge_ledger_nspecies')
      read (value, *, iostat=ios) nspecies
      found_nspecies = ios == 0
    case ('charge_ledger_batch_count')
      read (value, *, iostat=ios) batch_count
      found_batch = ios == 0
    case ('charge_ledger_surface_charge_before_C')
      read (value, *, iostat=ios) stock_values(1)
      found_stocks(1) = ios == 0
    case ('charge_ledger_surface_charge_after_C')
      read (value, *, iostat=ios) stock_values(2)
      found_stocks(2) = ios == 0
    case ('charge_ledger_local_flight_charge_before_C')
      read (value, *, iostat=ios) stock_values(3)
      found_stocks(3) = ios == 0
    case ('charge_ledger_local_flight_charge_after_C')
      read (value, *, iostat=ios) stock_values(4)
      found_stocks(4) = ios == 0
    case ('charge_ledger_unresolved_stock_before_C')
      read (value, *, iostat=ios) stock_values(5)
      found_stocks(5) = ios == 0
    case ('charge_ledger_unresolved_stock_after_C')
      read (value, *, iostat=ios) stock_values(6)
      found_stocks(6) = ios == 0
    end select
  end do
  close (u)
  if (.not. found_nspecies .or. .not. found_batch .or. .not. all(found_stocks)) then
    error stop 'Resume checkpoint charge ledger summary is incomplete.'
  end if
  if (nspecies < 1_i32 .or. batch_count < 0_i32) then
    error stop 'Resume checkpoint charge ledger dimensions are invalid.'
  end if
  if (any(.not. ieee_is_finite(stock_values))) then
    error stop 'Resume checkpoint charge ledger stocks must be finite.'
  end if

  call ledger%init(nspecies)
  ledger%batch_count = batch_count
  ledger%surface_charge_before = stock_values(1)
  ledger%surface_charge_after = stock_values(2)
  ledger%local_flight_charge_before = stock_values(3)
  ledger%local_flight_charge_after = stock_values(4)
  ledger%unresolved_stock_before = stock_values(5)
  ledger%unresolved_stock_after = stock_values(6)
  allocate (seen(nspecies))
  seen = .false.
  loaded = 0_i32
  open (newunit=u, file=trim(ledger_path), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'Failed to open charge_ledger.csv for resume.'
  read (u, '(A)', iostat=ios) header
  if (ios /= 0) error stop 'Failed to read charge_ledger.csv header.'
  has_fixed_current_columns = index(header, 'fixed_current_correction_C') > 0
  has_current_budget_columns = index(header, 'fixed_escape_correction_C') > 0
  columns = 15
  if (has_fixed_current_columns) columns = 20
  if (has_current_budget_columns) columns = 25
  do
    call read_checkpoint_csv_record(u, record, columns, ios)
    if (ios < 0) exit
    if (ios > 0) error stop 'Failed to parse charge_ledger.csv during resume.'
    charge_values = [ &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                    0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp &
                    ]
    if (has_current_budget_columns) then
      read (record, *, iostat=ios) row_batch, species_idx, charge_values, count_values
    else if (has_fixed_current_columns) then
      read (record, *, iostat=ios) row_batch, species_idx, fixed_charge_values, count_values
      if (ios == 0) then
        charge_values(1:13) = fixed_charge_values
        charge_values(14) = fixed_charge_values(9)
        charge_values(15) = fixed_charge_values(11)
      end if
    else
      read (record, *, iostat=ios) row_batch, species_idx, legacy_charge_values, count_values
      if (ios == 0) charge_values(1:8) = legacy_charge_values
    end if
    if (ios /= 0) error stop 'Failed to parse charge_ledger.csv during resume.'
    if (row_batch /= batch_count) error stop 'Resume charge ledger batch count mismatch.'
    if (species_idx < 1_i32 .or. species_idx > nspecies) error stop 'Resume charge ledger species index is invalid.'
    if (seen(species_idx)) error stop 'Resume charge ledger contains duplicate species rows.'
    if (any(.not. ieee_is_finite(charge_values))) error stop 'Resume charge ledger charges must be finite.'
    if (charge_values(14) /= charge_values(9) .or. charge_values(15) /= charge_values(11) .or. &
        charge_values(17) /= charge_values(16)) then
      error stop 'Resume charge ledger applied-charge aliases must match their target charges.'
    end if
    if (any(count_values < 0_i64)) error stop 'Resume charge ledger counts must be nonnegative.'
    seen(species_idx) = .true.
    ledger%injected_from_remote(species_idx) = charge_values(1)
    ledger%emitted_from_surface(species_idx) = charge_values(2)
    ledger%absorbed_on_surface(species_idx) = charge_values(3)
    ledger%escaped_to_infinity(species_idx) = charge_values(4)
    ledger%discarded_unresolved(species_idx) = charge_values(5)
    ledger%neutral_return_correction(species_idx) = charge_values(6)
    ledger%neutral_return_weight_scale(species_idx) = charge_values(7)
    ledger%neutral_return_unresolved_fraction(species_idx) = charge_values(8)
    ! Applied-charge CSV columns are backward-compatible aliases.  Only the
    ! target charge is retained in memory because the closure applies it exactly.
    ledger%fixed_absorbed_target_charge(species_idx) = charge_values(9)
    ledger%fixed_absorbed_weight_scale(species_idx) = charge_values(10)
    ledger%fixed_emission_target_charge(species_idx) = charge_values(11)
    ledger%fixed_emission_weight_scale(species_idx) = charge_values(12)
    ledger%fixed_current_correction(species_idx) = charge_values(13)
    ledger%fixed_escape_target_charge(species_idx) = charge_values(16)
    ledger%fixed_escape_correction(species_idx) = charge_values(18)
    ledger%injected_count(species_idx) = count_values(1)
    ledger%emitted_count(species_idx) = count_values(2)
    ledger%absorbed_count(species_idx) = count_values(3)
    ledger%escaped_count(species_idx) = count_values(4)
    ledger%discarded_unresolved_count(species_idx) = count_values(5)
    loaded = loaded + 1_i32
  end do
  close (u)
  if (loaded /= nspecies) error stop 'Resume charge ledger species rows are incomplete.'
  end procedure load_charge_ledger_checkpoint

  !> summary.txt から復元した統計値が壊れていないことを検証する。
  subroutine validate_summary_stats(stats)
    type(sim_stats), intent(in) :: stats
    real(dp) :: matching_budget_scale

    if (stats%processed_particles < 0_i64) error stop 'Resume checkpoint processed_particles must be >= 0.'
    if (stats%absorbed < 0_i64) error stop 'Resume checkpoint absorbed must be >= 0.'
    if (stats%escaped < 0_i64) error stop 'Resume checkpoint escaped must be >= 0.'
    if (stats%escaped_boundary < 0_i64) error stop 'Resume checkpoint escaped_boundary must be >= 0.'
    if (stats%survived_max_step < 0_i64) error stop 'Resume checkpoint survived_max_step must be >= 0.'
    if (stats%multiple_box_events_retry_attempted < 0_i64 .or. &
        stats%multiple_box_events_retry_resolved < 0_i64 .or. &
        stats%multiple_box_events_retry_resolved > stats%multiple_box_events_retry_attempted) then
      error stop 'Resume checkpoint multiple_box_events retry counters are inconsistent.'
    end if
    if (stats%multiple_box_events_soft_discarded < 0_i64) then
      error stop 'Resume checkpoint multiple_box_events_soft_discarded must be >= 0.'
    end if
    if (stats%adaptive_grown_steps < 0_i64) error stop 'Resume checkpoint adaptive_grown_steps must be >= 0.'
    if (stats%adaptive_substep_halved_steps < 0_i64 .or. stats%adaptive_substeps < 0_i64 .or. &
        stats%adaptive_substep_depth_limited_steps < 0_i64 .or. &
        stats%adaptive_substep_halved_steps > stats%adaptive_substeps .or. &
        stats%adaptive_substep_depth_limited_steps > stats%adaptive_substep_halved_steps) then
      error stop 'Resume checkpoint adaptive_substep counters are inconsistent.'
    end if
    if (.not. ieee_is_finite(stats%multiple_box_events_soft_discarded_abs_charge) .or. &
        stats%multiple_box_events_soft_discarded_abs_charge < 0.0_dp) then
      error stop 'Resume checkpoint multiple_box_events_soft_discarded_abs_charge_C must be finite and >= 0.'
    end if
    if (stats%batches < 0_i32) error stop 'Resume checkpoint batches must be >= 0.'
    if (.not. ieee_is_finite(stats%last_rel_change) .or. stats%last_rel_change < 0.0d0) then
      error stop 'Resume checkpoint last_rel_change must be finite and >= 0.'
    end if
    if (.not. ieee_is_finite(stats%simulated_time) .or. stats%simulated_time < 0.0_dp) then
      error stop 'Resume checkpoint simulated_time_s must be finite and >= 0.'
    end if
    if (stats%adaptive_nonzero_mode%rejected_trials < 0_i64) then
      error stop 'Resume checkpoint adaptive rejected-trial count must be >= 0.'
    end if
    if (.not. ieee_is_finite(stats%adaptive_nonzero_mode%last_batch_duration) .or. &
        stats%adaptive_nonzero_mode%last_batch_duration < 0.0_dp) then
      error stop 'Resume checkpoint adaptive last batch duration must be finite and >= 0.'
    end if
    if (.not. ieee_is_finite(stats%adaptive_nonzero_mode%last_potential_step) .or. &
        stats%adaptive_nonzero_mode%last_potential_step < 0.0_dp) then
      error stop 'Resume checkpoint adaptive last potential step must be finite and >= 0.'
    end if
    if (stats%adaptive_nonzero_mode%omp_threads < 0_i32) then
      error stop 'Resume checkpoint adaptive OpenMP thread count must be >= 0.'
    end if
    if (stats%secondary_emission_hits < 0_i64 .or. stats%secondary_emission_count < 0_i64 .or. &
        stats%secondary_generation_capped_count < 0_i64 .or. &
        .not. ieee_is_finite(stats%secondary_yield_sum) .or. stats%secondary_yield_sum < 0.0_dp .or. &
        .not. ieee_is_finite(stats%secondary_emitted_charge_c) .or. stats%secondary_emitted_charge_c > 0.0_dp) then
      error stop 'Resume checkpoint secondary-emission diagnostics are invalid.'
    end if
    if (stats%matching_plane%state_valid) then
      if (.not. all(ieee_is_finite([ &
                                   stats%matching_plane%displacement_c_m2, stats%matching_plane%phi_v, &
                                   stats%matching_plane%response, stats%matching_plane%feedback, &
                                   stats%matching_plane%photoelectron_return_flux_m2_s, &
                                   stats%matching_plane%photoelectron_escape_flux_m2_s, stats%matching_plane%residual, &
                                   stats%matching_plane%secondary_share, stats%matching_plane%secondary_budget_residual_a_m2 &
                                   ]))) then
        error stop 'Resume checkpoint matching-plane state must be finite.'
      end if
      if (any(stats%matching_plane%feedback < 0.0_dp) .or. &
          any(stats%matching_plane%response(2:3) < 0.0_dp) .or. &
          stats%matching_plane%photoelectron_return_flux_m2_s < 0.0_dp .or. &
          stats%matching_plane%photoelectron_escape_flux_m2_s < 0.0_dp .or. &
          stats%matching_plane%iterations <= 0_i32 .or. stats%matching_plane%residual < 0.0_dp .or. &
          stats%matching_plane%secondary_share < 0.0_dp .or. stats%matching_plane%secondary_share > 1.0_dp) then
        error stop 'Resume checkpoint matching-plane state is outside its physical range.'
      end if
      matching_budget_scale = max( &
                              1.0_dp, stats%matching_plane%feedback(1), &
                              stats%matching_plane%photoelectron_return_flux_m2_s, &
                              stats%matching_plane%photoelectron_escape_flux_m2_s &
                              )
      if (abs( &
          stats%matching_plane%feedback(1) - &
          stats%matching_plane%photoelectron_return_flux_m2_s - &
          stats%matching_plane%photoelectron_escape_flux_m2_s &
          ) > sqrt(epsilon(1.0_dp))*matching_budget_scale) then
        error stop 'Resume checkpoint matching-plane photoelectron budget is inconsistent.'
      end if
      if (.not. all(ieee_is_finite(stats%matching_plane%photoelectron_spectrum_m2_s)) .or. &
          any(stats%matching_plane%photoelectron_spectrum_m2_s < 0.0_dp)) then
        error stop 'Resume checkpoint matching-plane photoelectron spectrum must be finite and >= 0.'
      end if
      if (any(stats%matching_plane%photoelectron_spectrum_m2_s > 0.0_dp) .and. &
          abs(sum(stats%matching_plane%photoelectron_spectrum_m2_s) - stats%matching_plane%feedback(1)) > &
          sqrt(epsilon(1.0_dp))*matching_budget_scale) then
        error stop 'Resume checkpoint matching-plane photoelectron spectrum does not match its outward flux.'
      end if
    end if
  end subroutine validate_summary_stats

end submodule bem_restart_records
