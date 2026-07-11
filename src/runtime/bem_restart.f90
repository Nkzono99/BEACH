!> チェックポイントファイルの保存/復元を扱う補助モジュール。
module bem_restart
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: sim_stats, mesh_type, injection_state
  use bem_app_config_types, only: app_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_restart_state_type
  use bem_model_fingerprint, only: model_fingerprint, mesh_fingerprint, species_fingerprint
  use bem_physics_config_types, only: validate_phase0_physics_config, physics_config_ok
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  implicit none

  private
  integer(i32), parameter, public :: restart_contract_ok = 0_i32
  integer(i32), parameter, public :: restart_contract_mismatch = 1_i32
  integer(i32), parameter, public :: restart_contract_unsupported_schema = 2_i32
  integer(i32), parameter, public :: restart_contract_malformed = 3_i32
  public :: load_restart_checkpoint
  public :: validate_restart_contract
  public :: write_rng_state_file
  public :: write_macro_residuals_file
  public :: restart_rng_state_path
  public :: restart_macro_residual_path

contains

  !> 既存チェックポイントディレクトリから統計・要素電荷・乱数状態を復元する。
  !! @param[in] out_dir チェックポイントを探索するディレクトリ。
  !! @param[inout] mesh 現在のメッシュ。`q_elem` を復元値で上書きする。
  !! @param[out] stats 復元された統計値。
  !! @param[out] has_restart 復元可能なチェックポイントが存在したか。
  !! @param[inout] state 種別ごとのマクロ粒子残差（指定時のみ復元）。
  subroutine load_restart_checkpoint( &
    out_dir, mesh, stats, has_restart, state, mpi_rank, mpi_size, mpi, require_checkpoint, app, charge_ledger, &
    electrostatic_state &
    )
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(inout) :: mesh
    type(sim_stats), intent(out) :: stats
    logical, intent(out) :: has_restart
    type(injection_state), intent(inout), optional :: state
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    logical, intent(in), optional :: require_checkpoint
    type(app_config), intent(in), optional :: app
    type(charge_ledger_type), intent(inout), optional :: charge_ledger
    type(electrostatic_restart_state_type), intent(out), optional :: electrostatic_state

    character(len=1024) :: summary_path, charges_path, rng_path, residual_path, ledger_path
    character(len=256) :: contract_message
    logical :: has_summary, has_charges, has_rng, has_residual, has_legacy_residual, has_ledger
    logical :: must_have_checkpoint
    integer(i32) :: local_rank, world_size, contract_status

    stats = sim_stats()
    if (present(electrostatic_state)) electrostatic_state = electrostatic_restart_state_type()
    has_restart = .false.
    must_have_checkpoint = .false.
    if (present(require_checkpoint)) must_have_checkpoint = require_checkpoint
    call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'load_restart_checkpoint')

    summary_path = trim(out_dir)//'/summary.txt'
    charges_path = trim(out_dir)//'/charges.csv'
    rng_path = restart_rng_state_path(trim(out_dir), mpi_rank=local_rank, mpi_size=world_size)
    residual_path = restart_macro_residual_path(trim(out_dir), mpi_rank=local_rank, mpi_size=world_size)
    ledger_path = trim(out_dir)//'/charge_ledger.csv'

    inquire (file=trim(summary_path), exist=has_summary)
    inquire (file=trim(charges_path), exist=has_charges)
    inquire (file=trim(rng_path), exist=has_rng)
    inquire (file=trim(residual_path), exist=has_residual)
    inquire (file=trim(ledger_path), exist=has_ledger)
    call detect_legacy_ranked_residuals(trim(out_dir), local_rank, world_size, has_legacy_residual, mpi)

    if (has_legacy_residual) then
      error stop 'Resume checkpoint contains legacy rank-local macro residual files; use one global macro_residuals.csv.'
    end if

    if (.not. has_summary .and. .not. has_charges .and. .not. has_rng) then
      if (must_have_checkpoint) error stop 'Resume requested but checkpoint files are missing in checkpoint directory.'
      return
    end if

    if (.not. (has_summary .and. has_charges .and. has_rng)) then
      error stop 'Resume requested but checkpoint files are incomplete in checkpoint directory.'
    end if

    if (present(app)) then
      call validate_restart_contract(trim(summary_path), mesh, app, contract_status, contract_message)
      if (contract_status /= restart_contract_ok) then
        error stop 'Resume checkpoint contract mismatch: '//trim(contract_message)
      end if
    end if

    call load_summary_file(trim(summary_path), mesh%nelem, stats, expected_world_size=world_size)
    if (present(electrostatic_state)) then
      call load_electrostatic_state(trim(summary_path), electrostatic_state)
      if (present(app)) then
        if (trim(app%periodic2%zero_mode_policy) == 'exclude_k0' .and. .not. electrostatic_state%outer_ready) then
          error stop 'Resume checkpoint is missing the required split-periodic outer state.'
        end if
      end if
    end if
    call load_charge_file(trim(charges_path), mesh)
    if (present(charge_ledger)) then
      if (has_ledger) then
        call load_charge_ledger_checkpoint(trim(summary_path), trim(ledger_path), charge_ledger)
        if (charge_ledger%batch_count /= stats%batches) then
          error stop 'Resume charge ledger batch count does not match summary statistics.'
        end if
        if (present(app)) then
          if (charge_ledger%nspecies /= app%n_particle_species) then
            error stop 'Resume charge ledger species count does not match current config.'
          end if
        end if
      else if (present(app) .and. app%n_particle_species > 0_i32) then
        call charge_ledger%init(app%n_particle_species)
        charge_ledger%batch_count = stats%batches
        charge_ledger%surface_charge_before = sum(mesh%q_elem)
        charge_ledger%surface_charge_after = sum(mesh%q_elem)
      end if
    end if
    call restore_rng_state(trim(rng_path))
    if (present(state)) then
      if (allocated(state%macro_residual)) state%macro_residual = 0.0d0
      if (has_residual .and. allocated(state%macro_residual)) then
        if (.not. present(mpi) .or. local_rank == 0_i32) call load_macro_residual_file(trim(residual_path), state)
        if (present(mpi)) call mpi_bcast_real_dp_array(mpi, state%macro_residual, 0_i32)
      end if
    end if
    has_restart = .true.
  end subroutine load_restart_checkpoint

  subroutine load_electrostatic_state(path, state)
    character(len=*), intent(in) :: path
    type(electrostatic_restart_state_type), intent(out) :: state
    integer :: u, ios, pos
    character(len=512) :: line
    character(len=64) :: key
    character(len=256) :: value
    logical :: found_potential, found_batch

    state = electrostatic_restart_state_type()
    found_potential = .false.
    found_batch = .false.
    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary.txt for electrostatic restart state.'
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      pos = index(line, '=')
      if (pos <= 0) cycle
      key = trim(adjustl(line(:pos - 1)))
      value = trim(adjustl(line(pos + 1:)))
      select case (trim(key))
      case ('electrostatic_split_periodic_active')
        read (value, *, iostat=ios) state%outer_ready
      case ('interface_potential_V')
        read (value, *, iostat=ios) state%outer_interface_potential
        found_potential = ios == 0
      case ('last_outer_update_batch')
        read (value, *, iostat=ios) state%last_outer_update_batch
        found_batch = ios == 0
      end select
      if (ios /= 0) error stop 'Malformed electrostatic restart state in summary.txt.'
    end do
    close (u)
    if (state%outer_ready .and. .not. (found_potential .and. found_batch)) then
      error stop 'Incomplete electrostatic restart state in summary.txt.'
    end if
  end subroutine load_electrostatic_state

  !> schema v2 fingerprint を現在の ordered model/mesh/species contract と照合する。
  subroutine validate_restart_contract(path, mesh, app, status, message)
    character(len=*), intent(in) :: path
    type(mesh_type), intent(in) :: mesh
    type(app_config), intent(in) :: app
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer :: u, ios, pos
    integer(i32) :: schema_version
    character(len=512) :: line
    character(len=64) :: key
    character(len=256) :: value
    character(len=16) :: saved_model, saved_mesh, saved_species
    logical :: found_schema, found_model, found_mesh, found_species
    integer(i32) :: physics_status
    character(len=256) :: physics_message

    status = restart_contract_ok
    message = ''
    schema_version = -1_i32
    saved_model = ''
    saved_mesh = ''
    saved_species = ''
    found_schema = .false.
    found_model = .false.
    found_mesh = .false.
    found_species = .false.

    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      status = restart_contract_malformed
      message = 'cannot open summary.txt'
      return
    end if
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      pos = index(line, '=')
      if (pos <= 0) cycle
      key = trim(adjustl(line(:pos - 1)))
      value = trim(adjustl(line(pos + 1:)))
      select case (trim(key))
      case ('checkpoint_schema_version')
        read (value, *, iostat=ios) schema_version
        if (ios /= 0) then
          close (u)
          status = restart_contract_malformed
          message = 'invalid checkpoint_schema_version'
          return
        end if
        found_schema = .true.
      case ('model_fingerprint')
        saved_model = trim(value)
        found_model = .true.
      case ('mesh_fingerprint')
        saved_mesh = trim(value)
        found_mesh = .true.
      case ('species_fingerprint')
        saved_species = trim(value)
        found_species = .true.
      end select
    end do
    close (u)

    ! Legacy checkpoints predate fingerprints and are accepted only for implemented Phase 0 point-source modes.
    if (.not. found_schema) then
      call validate_phase0_physics_config( &
        app%field, app%periodic2, app%panel, app%outer_plasma, app%coupling, physics_status, physics_message &
        )
      if (physics_status /= physics_config_ok) then
        status = restart_contract_mismatch
        message = 'legacy checkpoint is incompatible with this physics model'
      end if
      return
    end if
    if (schema_version /= 2_i32) then
      status = restart_contract_unsupported_schema
      message = 'unsupported checkpoint schema version'
      return
    end if
    if (.not. (found_model .and. found_mesh .and. found_species)) then
      status = restart_contract_malformed
      message = 'schema v2 summary is missing fingerprints'
      return
    end if
    if (saved_model /= model_fingerprint(app)) then
      status = restart_contract_mismatch
      message = 'model fingerprint differs'
    else if (saved_mesh /= mesh_fingerprint(mesh)) then
      status = restart_contract_mismatch
      message = 'mesh fingerprint differs'
    else if (saved_species /= species_fingerprint(app)) then
      status = restart_contract_mismatch
      message = 'species fingerprint differs'
    end if
  end subroutine validate_restart_contract

  !> 現在の Fortran 乱数状態をファイルへ保存する。
  !! @param[in] out_dir 出力ディレクトリ。
  subroutine write_rng_state_file(out_dir, mpi_rank, mpi_size, mpi)
    character(len=*), intent(in) :: out_dir
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi

    character(len=1024) :: path
    integer :: n, u, ios, i
    integer, allocatable :: seed(:)
    integer(i32) :: local_rank, world_size

    call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'write_rng_state_file')
    call random_seed(size=n)
    allocate (seed(n))
    call random_seed(get=seed)

    path = restart_rng_state_path(trim(out_dir), mpi_rank=local_rank, mpi_size=world_size)
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open rng_state.txt.'

    write (u, '(i0)') n
    do i = 1, n
      write (u, '(i0)') seed(i)
    end do
    close (u)
  end subroutine write_rng_state_file

  !> マクロ粒子残差を `macro_residuals.csv` として保存する。
    !! @param[in] out_dir 出力ディレクトリ。
    !! @param[in] state 種別ごとのマクロ粒子残差を保持した注入状態。
  subroutine write_macro_residuals_file(out_dir, state, mpi_rank, mpi_size, mpi)
    character(len=*), intent(in) :: out_dir
    type(injection_state), intent(in) :: state
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi

    character(len=1024) :: path
    integer :: u, ios, i
    integer(i32) :: local_rank, world_size

    if (.not. allocated(state%macro_residual)) return

    call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'write_macro_residuals_file')
    if (local_rank /= 0_i32) return

    path = restart_macro_residual_path(trim(out_dir), mpi_rank=local_rank, mpi_size=world_size)
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open macro_residuals.csv.'

    write (u, '(a)') 'species_idx,residual'
    do i = 1, size(state%macro_residual)
      write (u, '(i0,a,es24.16)') i, ',', state%macro_residual(i)
    end do
    close (u)
  end subroutine write_macro_residuals_file

  !> `summary.txt` を読み込み、必須キーの存在と要素数整合を検証する。
  !! 欠落キーやメッシュ要素数不一致は再開不能として停止する。
  !! @param[in] path `summary.txt` のファイルパス。
  !! @param[in] expected_nelem 現在メッシュの要素数（整合性検証に使用）。
  !! @param[out] stats 復元した統計値。
  subroutine load_summary_file(path, expected_nelem, stats, expected_world_size)
    character(len=*), intent(in) :: path
    integer(i32), intent(in) :: expected_nelem
    type(sim_stats), intent(out) :: stats
    integer(i32), intent(in), optional :: expected_world_size

    integer :: u, ios, pos
    integer(i32) :: mesh_nelem, saved_world_size
    character(len=512) :: line
    character(len=64) :: key
    character(len=256) :: value
    logical :: found_mesh, found_processed, found_absorbed, found_escaped
    logical :: found_batches, found_rel, found_world_size

    stats = sim_stats()
    mesh_nelem = -1_i32
    saved_world_size = 1_i32
    found_mesh = .false.
    found_processed = .false.
    found_absorbed = .false.
    found_escaped = .false.
    found_batches = .false.
    found_rel = .false.
    found_world_size = .false.

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
      case ('last_rel_change')
        read (value, *) stats%last_rel_change
        found_rel = .true.
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
    call validate_summary_stats(stats)
    if (present(expected_world_size)) then
      if (.not. found_world_size .and. expected_world_size > 1_i32) then
        error stop 'Resume checkpoint summary is missing mpi_world_size.'
      end if
      if (max(1_i32, expected_world_size) /= saved_world_size) then
        error stop 'Resume checkpoint mpi_world_size does not match current MPI world size.'
      end if
    end if
  end subroutine load_summary_file

  !> `charges.csv` を読み込み、各要素の電荷をメッシュへ復元する。
  !! 行重複や要素数不足を検出した場合は停止する。
  !! @param[in] path `charges.csv` のファイルパス。
  !! @param[inout] mesh 要素電荷 `q_elem` を復元値で上書きするメッシュ。
  subroutine load_charge_file(path, mesh)
    character(len=*), intent(in) :: path
    type(mesh_type), intent(inout) :: mesh

    integer :: u, ios
    integer(i32) :: elem_idx, n_loaded
    real(dp) :: charge
    character(len=512) :: header
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
      read (u, *, iostat=ios) elem_idx, charge
      if (ios < 0) exit
      if (ios > 0) error stop 'Failed to parse charges.csv during resume.'
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
  end subroutine load_charge_file

  !> summary stock と species 別 CSV flux/count から累積 charge ledger を復元する。
  subroutine load_charge_ledger_checkpoint(summary_path, ledger_path, ledger)
    character(len=*), intent(in) :: summary_path, ledger_path
    type(charge_ledger_type), intent(inout) :: ledger
    integer :: u, ios, pos
    integer(i32) :: nspecies, batch_count, row_batch, species_idx, loaded
    integer(i64) :: count_values(5)
    real(dp) :: charge_values(7), stock_values(8)
    character(len=512) :: line, header
    character(len=96) :: key
    character(len=256) :: value
    logical :: found_nspecies, found_batch, found_stocks(8)
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
      case ('charge_ledger_outer_flight_charge_before_C')
        read (value, *, iostat=ios) stock_values(5)
        found_stocks(5) = ios == 0
      case ('charge_ledger_outer_flight_charge_after_C')
        read (value, *, iostat=ios) stock_values(6)
        found_stocks(6) = ios == 0
      case ('charge_ledger_unresolved_stock_before_C')
        read (value, *, iostat=ios) stock_values(7)
        found_stocks(7) = ios == 0
      case ('charge_ledger_unresolved_stock_after_C')
        read (value, *, iostat=ios) stock_values(8)
        found_stocks(8) = ios == 0
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
    ledger%outer_flight_charge_before = stock_values(5)
    ledger%outer_flight_charge_after = stock_values(6)
    ledger%unresolved_stock_before = stock_values(7)
    ledger%unresolved_stock_after = stock_values(8)
    allocate (seen(nspecies))
    seen = .false.
    loaded = 0_i32
    open (newunit=u, file=trim(ledger_path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charge_ledger.csv for resume.'
    read (u, '(A)', iostat=ios) header
    if (ios /= 0) error stop 'Failed to read charge_ledger.csv header.'
    do
      read (u, *, iostat=ios) &
        row_batch, species_idx, charge_values, count_values
      if (ios < 0) exit
      if (ios > 0) error stop 'Failed to parse charge_ledger.csv during resume.'
      if (row_batch /= batch_count) error stop 'Resume charge ledger batch count mismatch.'
      if (species_idx < 1_i32 .or. species_idx > nspecies) error stop 'Resume charge ledger species index is invalid.'
      if (seen(species_idx)) error stop 'Resume charge ledger contains duplicate species rows.'
      if (any(.not. ieee_is_finite(charge_values))) error stop 'Resume charge ledger charges must be finite.'
      if (any(count_values < 0_i64)) error stop 'Resume charge ledger counts must be nonnegative.'
      seen(species_idx) = .true.
      ledger%injected_from_remote(species_idx) = charge_values(1)
      ledger%emitted_from_surface(species_idx) = charge_values(2)
      ledger%absorbed_on_surface(species_idx) = charge_values(3)
      ledger%escaped_to_infinity(species_idx) = charge_values(4)
      ledger%discarded_unresolved(species_idx) = charge_values(5)
      ledger%interface_outward_gross(species_idx) = charge_values(6)
      ledger%interface_returned_gross(species_idx) = charge_values(7)
      ledger%injected_count(species_idx) = count_values(1)
      ledger%emitted_count(species_idx) = count_values(2)
      ledger%absorbed_count(species_idx) = count_values(3)
      ledger%escaped_count(species_idx) = count_values(4)
      ledger%discarded_unresolved_count(species_idx) = count_values(5)
      loaded = loaded + 1_i32
    end do
    close (u)
    if (loaded /= nspecies) error stop 'Resume charge ledger species rows are incomplete.'
  end subroutine load_charge_ledger_checkpoint

  !> 保存済み乱数状態を読み戻し、このビルドの RNG 状態へ復元する。
  !! RNG 内部状態の長さが一致しない場合は互換性がないため停止する。
  !! @param[in] path `rng_state.txt` のファイルパス。
  subroutine restore_rng_state(path)
    character(len=*), intent(in) :: path

    integer :: expected_n, file_n, u, ios, i
    integer, allocatable :: seed(:)

    call random_seed(size=expected_n)

    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open rng_state.txt for resume.'

    read (u, *, iostat=ios) file_n
    if (ios /= 0) error stop 'Failed to read rng_state.txt header.'
    if (file_n /= expected_n) then
      error stop 'Resume checkpoint RNG state size does not match this build.'
    end if

    allocate (seed(file_n))
    do i = 1, file_n
      read (u, *, iostat=ios) seed(i)
      if (ios /= 0) error stop 'Failed to parse rng_state.txt.'
    end do
    close (u)

    call random_seed(put=seed)
  end subroutine restore_rng_state

  !> 保存済みマクロ粒子残差を読み戻す。
  !! @param[in] path `macro_residuals.csv` のファイルパス。
  !! @param[inout] state 種別ごとのマクロ粒子残差を書き戻す注入状態。
  subroutine load_macro_residual_file(path, state)
    character(len=*), intent(in) :: path
    type(injection_state), intent(inout) :: state

    integer :: u, ios
    integer(i32) :: species_idx
    real(dp) :: residual
    character(len=512) :: header
    logical, allocatable :: seen(:)

    if (.not. allocated(state%macro_residual)) return

    allocate (seen(size(state%macro_residual)))
    seen = .false.
    state%macro_residual = 0.0d0

    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open macro_residuals.csv for resume.'

    read (u, '(A)', iostat=ios) header
    if (ios /= 0) error stop 'Failed to read macro_residuals.csv header.'

    do
      read (u, *, iostat=ios) species_idx, residual
      if (ios < 0) exit
      if (ios > 0) error stop 'Failed to parse macro_residuals.csv during resume.'
      if (species_idx < 1_i32 .or. species_idx > size(state%macro_residual)) then
        error stop 'Resume checkpoint macro_residuals.csv has an invalid species index.'
      end if
      if (seen(species_idx)) then
        error stop 'Resume checkpoint macro_residuals.csv contains duplicate species rows.'
      end if
      if (.not. ieee_is_finite(residual) .or. residual < 0.0d0 .or. residual >= 1.0d0) then
        error stop 'Resume checkpoint macro_residuals.csv residual values must be finite and in [0, 1).'
      end if
      seen(species_idx) = .true.
      state%macro_residual(species_idx) = residual
    end do
    close (u)
  end subroutine load_macro_residual_file

  !> RNG状態ファイルのパスを返す。MPI複数rank時は rank 接尾辞付きパスへ切り替える。
  function restart_rng_state_path(out_dir, mpi_rank, mpi_size, mpi) result(path)
    character(len=*), intent(in) :: out_dir
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    character(len=1024) :: path
    integer(i32) :: local_rank, world_size

    call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'restart_rng_state_path')
    if (world_size <= 1_i32) then
      path = trim(out_dir)//'/rng_state.txt'
    else
      write (path, '(a,a,i5.5,a)') trim(out_dir), '/rng_state_rank', local_rank, '.txt'
    end if
  end function restart_rng_state_path

  !> MPI rank数によらないglobalマクロ残差ファイルのパスを返す。
  function restart_macro_residual_path(out_dir, mpi_rank, mpi_size, mpi) result(path)
    character(len=*), intent(in) :: out_dir
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    character(len=1024) :: path
    integer(i32) :: local_rank, world_size

    call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'restart_macro_residual_path')
    path = trim(out_dir)//'/macro_residuals.csv'
  end function restart_macro_residual_path

  !> 旧rank別マクロ残差が1個でも存在するかrootで確認し、全rankへ共有する。
  subroutine detect_legacy_ranked_residuals(out_dir, local_rank, world_size, found, mpi)
    character(len=*), intent(in) :: out_dir
    integer(i32), intent(in) :: local_rank, world_size
    logical, intent(out) :: found
    type(mpi_context), intent(in), optional :: mpi

    character(len=1024) :: legacy_path
    integer(i32) :: rank, found_value(1)
    logical :: exists

    found_value = 0_i32
    if (world_size <= 1_i32) then
      found = .false.
      return
    end if

    if (.not. present(mpi) .or. local_rank == 0_i32) then
      do rank = 0_i32, world_size - 1_i32
        write (legacy_path, '(a,a,i5.5,a)') trim(out_dir), '/macro_residuals_rank', rank, '.csv'
        inquire (file=trim(legacy_path), exist=exists)
        if (exists) then
          found_value(1) = 1_i32
          exit
        end if
      end do
    end if
    if (present(mpi)) call mpi_bcast_i32_array(mpi, found_value, 0_i32)
    found = found_value(1) /= 0_i32
  end subroutine detect_legacy_ranked_residuals

  !> summary.txt から復元した統計値が壊れていないことを検証する。
  subroutine validate_summary_stats(stats)
    type(sim_stats), intent(in) :: stats

    if (stats%processed_particles < 0_i64) error stop 'Resume checkpoint processed_particles must be >= 0.'
    if (stats%absorbed < 0_i64) error stop 'Resume checkpoint absorbed must be >= 0.'
    if (stats%escaped < 0_i64) error stop 'Resume checkpoint escaped must be >= 0.'
    if (stats%escaped_boundary < 0_i64) error stop 'Resume checkpoint escaped_boundary must be >= 0.'
    if (stats%survived_max_step < 0_i64) error stop 'Resume checkpoint survived_max_step must be >= 0.'
    if (stats%batches < 0_i32) error stop 'Resume checkpoint batches must be >= 0.'
    if (.not. ieee_is_finite(stats%last_rel_change) .or. stats%last_rel_change < 0.0d0) then
      error stop 'Resume checkpoint last_rel_change must be finite and >= 0.'
    end if
  end subroutine validate_summary_stats

  !> 併存対応のため `mpi_context` と rank/size の両方を受け、最終的なrank/sizeを解決する。
  subroutine resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, caller_name)
    integer(i32), intent(out) :: local_rank, world_size
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    character(len=*), intent(in) :: caller_name

    call mpi_get_rank_size(local_rank, world_size, mpi)
    if (present(mpi_rank)) local_rank = mpi_rank
    if (present(mpi_size)) world_size = mpi_size
    if (world_size <= 0_i32) error stop 'mpi_size must be > 0 in '//trim(caller_name)//'.'
    if (local_rank < 0_i32 .or. local_rank >= world_size) then
      error stop 'mpi_rank out of range in '//trim(caller_name)//'.'
    end if
  end subroutine resolve_parallel_rank_size

end module bem_restart
