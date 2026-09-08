!> チェックポイントファイルの保存/復元を扱う補助モジュール。
module bem_restart
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: i32
  use bem_types, only: sim_stats, mesh_type, injection_state
  use bem_app_config_types, only: app_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_checkpoint_contract, only: checkpoint_schema_is_loadable, inspect_checkpoint_directory
  use bem_mpi, only: mpi_context, mpi_bcast_real_dp_array
  implicit none

  private
  integer(i32), parameter, public :: restart_contract_ok = 0_i32
  integer(i32), parameter, public :: restart_contract_mismatch = 1_i32
  integer(i32), parameter, public :: restart_contract_unsupported_schema = 2_i32
  integer(i32), parameter, public :: restart_contract_malformed = 3_i32
  integer(i32), parameter, public :: restart_contract_configuration_changed = 4_i32
  public :: load_restart_checkpoint
  public :: validate_restart_contract
  public :: write_rng_state_file
  public :: write_macro_residuals_file
  public :: restart_rng_state_path
  public :: restart_macro_residual_path

  interface
    module subroutine load_summary_file(path, expected_nelem, stats, expected_world_size)
      character(len=*), intent(in) :: path
      integer(i32), intent(in) :: expected_nelem
      type(sim_stats), intent(out) :: stats
      integer(i32), intent(in), optional :: expected_world_size
    end subroutine load_summary_file

    module subroutine load_charge_file(path, mesh)
      character(len=*), intent(in) :: path
      type(mesh_type), intent(inout) :: mesh
    end subroutine load_charge_file

    module subroutine load_charge_ledger_checkpoint(summary_path, ledger_path, ledger)
      character(len=*), intent(in) :: summary_path, ledger_path
      type(charge_ledger_type), intent(inout) :: ledger
    end subroutine load_charge_ledger_checkpoint

    module subroutine write_rng_state_file(out_dir, mpi_rank, mpi_size, mpi)
      character(len=*), intent(in) :: out_dir
      integer(i32), intent(in), optional :: mpi_rank, mpi_size
      type(mpi_context), intent(in), optional :: mpi
    end subroutine write_rng_state_file

    module subroutine write_macro_residuals_file(out_dir, state, mpi_rank, mpi_size, mpi)
      character(len=*), intent(in) :: out_dir
      type(injection_state), intent(in) :: state
      integer(i32), intent(in), optional :: mpi_rank, mpi_size
      type(mpi_context), intent(in), optional :: mpi
    end subroutine write_macro_residuals_file

    module subroutine restore_rng_state(path)
      character(len=*), intent(in) :: path
    end subroutine restore_rng_state

    module subroutine load_macro_residual_file(path, state)
      character(len=*), intent(in) :: path
      type(injection_state), intent(inout) :: state
    end subroutine load_macro_residual_file

    module function restart_rng_state_path(out_dir, mpi_rank, mpi_size, mpi) result(path)
      character(len=*), intent(in) :: out_dir
      integer(i32), intent(in), optional :: mpi_rank, mpi_size
      type(mpi_context), intent(in), optional :: mpi
      character(len=1024) :: path
    end function restart_rng_state_path

    module function restart_macro_residual_path(out_dir, mpi_rank, mpi_size, mpi) result(path)
      character(len=*), intent(in) :: out_dir
      integer(i32), intent(in), optional :: mpi_rank, mpi_size
      type(mpi_context), intent(in), optional :: mpi
      character(len=1024) :: path
    end function restart_macro_residual_path

    module subroutine detect_legacy_ranked_residuals(out_dir, local_rank, world_size, found, mpi)
      character(len=*), intent(in) :: out_dir
      integer(i32), intent(in) :: local_rank, world_size
      logical, intent(out) :: found
      type(mpi_context), intent(in), optional :: mpi
    end subroutine detect_legacy_ranked_residuals

    module subroutine validate_restart_contract(path, mesh, app, status, message)
      character(len=*), intent(in) :: path
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: app
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine validate_restart_contract
    module subroutine resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, caller_name)
      integer(i32), intent(out) :: local_rank, world_size
      integer(i32), intent(in), optional :: mpi_rank, mpi_size
      type(mpi_context), intent(in), optional :: mpi
      character(len=*), intent(in) :: caller_name
    end subroutine resolve_parallel_rank_size
  end interface

contains

  !> 既存チェックポイントディレクトリから統計・要素電荷・乱数状態を復元する。
  !! @param[in] out_dir チェックポイントを探索するディレクトリ。
  !! @param[inout] mesh 現在のメッシュ。`q_elem` を復元値で上書きする。
  !! @param[out] stats 復元された統計値。
  !! @param[out] has_restart 復元可能なチェックポイントが存在したか。
  !! @param[inout] state 種別ごとのマクロ粒子残差（指定時のみ復元）。
  subroutine load_restart_checkpoint( &
    out_dir, mesh, stats, has_restart, state, mpi_rank, mpi_size, mpi, require_checkpoint, app, charge_ledger &
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

    character(len=1024) :: summary_path, charges_path, rng_path, residual_path, ledger_path
    character(len=256) :: contract_message
    logical :: has_summary, has_charges, has_rng, has_legacy_residual
    logical :: checkpoint_complete, checkpoint_has_residual, checkpoint_has_ledger
    logical :: must_have_checkpoint
    integer(i32) :: local_rank, world_size, contract_status, residual_species, checkpoint_schema

    stats = sim_stats()
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
    call detect_legacy_ranked_residuals(trim(out_dir), local_rank, world_size, has_legacy_residual, mpi)

    if (has_legacy_residual) then
      error stop 'Resume checkpoint contains legacy rank-local macro residual files; use one global macro_residuals.csv.'
    end if

    if (.not. has_summary .and. .not. has_charges .and. .not. has_rng) then
      if (must_have_checkpoint) error stop 'Resume requested but checkpoint files are missing in checkpoint directory.'
      return
    end if

    call inspect_checkpoint_directory( &
      trim(out_dir), checkpoint_complete, has_macro_residuals=checkpoint_has_residual, &
      has_charge_ledger=checkpoint_has_ledger, schema_version=checkpoint_schema &
      )
    if (has_summary .and. .not. checkpoint_schema_is_loadable(checkpoint_schema)) then
      error stop 'Resume requested but checkpoint schema is unsupported.'
    end if
    if (.not. checkpoint_complete) then
      error stop 'Resume requested but checkpoint files are incomplete in checkpoint directory.'
    end if

    if (present(app)) then
      call validate_restart_contract(trim(summary_path), mesh, app, contract_status, contract_message)
      if (contract_status == restart_contract_configuration_changed) then
        if (local_rank == 0_i32) then
          write (error_unit, '(a)') 'WARNING: resume fingerprint differs: '//trim(contract_message)
          write (error_unit, '(a)') &
            'WARNING: continuing from the saved physical state with the current model/species configuration.'
          flush (error_unit)
        end if
      else if (contract_status /= restart_contract_ok) then
        error stop 'Resume checkpoint contract mismatch: '//trim(contract_message)
      end if
    end if

    call load_summary_file(trim(summary_path), mesh%nelem, stats, expected_world_size=world_size)
    call load_charge_file(trim(charges_path), mesh)
    if (present(charge_ledger)) then
      if (checkpoint_has_ledger) then
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
      if (allocated(state%boundary_macro_residual)) state%boundary_macro_residual = 0.0d0
      if (checkpoint_has_residual .and. allocated(state%macro_residual)) then
        if (.not. present(mpi) .or. local_rank == 0_i32) call load_macro_residual_file(trim(residual_path), state)
        if (present(mpi)) call mpi_bcast_real_dp_array(mpi, state%macro_residual, 0_i32)
        if (present(mpi) .and. allocated(state%boundary_macro_residual)) then
          do residual_species = 1_i32, int(size(state%boundary_macro_residual, 2), i32)
            call mpi_bcast_real_dp_array(mpi, state%boundary_macro_residual(:, residual_species), 0_i32)
          end do
        end if
      end if
    end if
    has_restart = .true.
  end subroutine load_restart_checkpoint

  !> 併存対応のため `mpi_context` と rank/size の両方を受け、最終的なrank/sizeを解決する。

end module bem_restart
