!> accepted batch 境界で二重化した定期チェックポイントを保存・解決する。
module bem_periodic_checkpoint
  use bem_kinds, only: i32
  use bem_types, only: injection_state, mesh_type, sim_stats
  use bem_app_config_types, only: app_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_checkpoint_contract, only: checkpoint_schema_is_loadable, inspect_checkpoint_directory, &
                                     publish_checkpoint_manifest
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_mpi, only: mpi_bcast_i32_array, mpi_context, mpi_is_root, mpi_world_barrier, mpi_world_size
  use bem_output_writer, only: ensure_output_dir, write_checkpoint_state_files
  use bem_restart, only: write_macro_residuals_file, write_rng_state_file
  implicit none
  private

  integer(i32), parameter :: periodic_checkpoint_index_schema = 1_i32
  character(len=*), parameter :: checkpoint_index_name = 'checkpoint_latest.txt'
  character(len=*), parameter :: checkpoint_parent_name = 'checkpoints'

  public :: maybe_write_periodic_checkpoint
  public :: resolve_latest_checkpoint_dir

contains

  !> 設定した stride の accepted batch で、非active slotへ完全な再開状態を書いてからindexを切り替える。
  subroutine maybe_write_periodic_checkpoint(app, mesh, stats, state, mpi, charge_ledger)
    type(app_config), intent(in) :: app
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(injection_state), intent(in), optional :: state
    type(mpi_context), intent(in) :: mpi
    type(charge_ledger_type), intent(in), optional :: charge_ledger

    integer(i32) :: slot_values(1), inactive_slot
    character(len=1024) :: checkpoint_dir
    logical :: has_macro_residuals, has_charge_ledger

    if (.not. app%write_output) return
    if (app%checkpoint_stride <= 0_i32) return
    if (stats%batches <= 0_i32) return
    if (mod(stats%batches, app%checkpoint_stride) /= 0_i32) return

    slot_values = 1_i32
    if (mpi_is_root(mpi)) then
      call read_active_slot(trim(app%output_dir), slot_values(1))
      slot_values(1) = 1_i32 - slot_values(1)
    end if
    call mpi_bcast_i32_array(mpi, slot_values, 0_i32)
    inactive_slot = slot_values(1)
    checkpoint_dir = checkpoint_slot_dir(trim(app%output_dir), inactive_slot)
    has_macro_residuals = .false.
    if (present(state)) has_macro_residuals = allocated(state%macro_residual)
    has_charge_ledger = present(charge_ledger)

    if (mpi_is_root(mpi)) then
      call ensure_output_dir(trim(checkpoint_dir))
      if (present(charge_ledger)) then
        call write_checkpoint_state_files( &
          trim(checkpoint_dir), mesh, stats, app, mpi_world_size=mpi_world_size(mpi), charge_ledger=charge_ledger &
          )
      else
        call write_checkpoint_state_files( &
          trim(checkpoint_dir), mesh, stats, app, mpi_world_size=mpi_world_size(mpi) &
          )
      end if
    end if
    ! 非root rankがslot directoryの作成より先にrank別RNGを書き始めないよう同期する。
    call mpi_world_barrier(mpi)
    call write_rng_state_file(trim(checkpoint_dir), mpi=mpi)
    if (present(state)) call write_macro_residuals_file(trim(checkpoint_dir), state, mpi=mpi)

    ! completion manifestを全rankのファイルが閉じられた後に公開し、続けてadvisory indexを更新する。
    call mpi_world_barrier(mpi)
    if (mpi_is_root(mpi)) then
      call publish_checkpoint_manifest( &
        trim(checkpoint_dir), stats%batches, mpi_world_size(mpi), has_macro_residuals, has_charge_ledger &
        )
      call publish_checkpoint_index(trim(app%output_dir), inactive_slot, stats%batches)
      print '(a,i0,a,a)', 'periodic_checkpoint_batch=', stats%batches, ' dir=', trim(checkpoint_dir)
    end if
    call mpi_world_barrier(mpi)
  end subroutine maybe_write_periodic_checkpoint

  !> base directory直下の最終出力と両定期slotを比較し、batch数が新しい完全checkpointを返す。
  !! indexは同batchのslotを選ぶtie-breakにだけ使う。欠落、破損、古いindexは
  !! complete manifestを持つload可能なslotの回収を妨げない。
  subroutine resolve_latest_checkpoint_dir(base_dir, checkpoint_dir)
    character(len=*), intent(in) :: base_dir
    character(len=*), intent(out) :: checkpoint_dir

    integer(i32) :: root_batch, root_schema, slot_batch, recovered_slot
    logical :: root_complete, slot_found
    character(len=1024) :: slot_dir

    checkpoint_dir = trim(base_dir)
    call inspect_checkpoint_directory( &
      trim(base_dir), root_complete, schema_version=root_schema, batches=root_batch &
      )
    if (.not. root_complete .or. .not. checkpoint_schema_is_loadable(root_schema)) root_batch = -1_i32

    call find_latest_periodic_slot(trim(base_dir), recovered_slot, slot_batch, slot_found)
    if (.not. slot_found) return
    slot_dir = checkpoint_slot_dir(trim(base_dir), recovered_slot)
    if (slot_batch > root_batch) checkpoint_dir = trim(slot_dir)
  end subroutine resolve_latest_checkpoint_dir

  subroutine read_active_slot(base_dir, active_slot)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(out) :: active_slot
    integer(i32) :: recovered_batch
    logical :: recovered

    call find_latest_periodic_slot(base_dir, active_slot, recovered_batch, recovered)
    if (.not. recovered) active_slot = 1_i32
  end subroutine read_active_slot

  !> indexに依存せず両slotを検査し、最新のload可能な完了slotを返す。
  !! batchが同じ場合だけ、summaryと整合するindex指定slotを優先する。
  subroutine find_latest_periodic_slot(base_dir, slot, batch, found)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(out) :: slot, batch
    logical, intent(out) :: found

    integer(i32) :: candidate_slot, candidate_batch, candidate_schema
    integer(i32) :: indexed_slot, indexed_batch
    logical :: candidate_complete, has_index, prefer_candidate
    character(len=1024) :: candidate_dir

    call read_checkpoint_index(base_dir, indexed_slot, indexed_batch, has_index)
    slot = -1_i32
    batch = -1_i32
    found = .false.
    do candidate_slot = 0_i32, 1_i32
      candidate_dir = checkpoint_slot_dir(trim(base_dir), candidate_slot)
      call inspect_checkpoint_directory( &
        trim(candidate_dir), candidate_complete, schema_version=candidate_schema, batches=candidate_batch &
        )
      if (.not. candidate_complete) cycle
      if (.not. checkpoint_schema_is_loadable(candidate_schema)) cycle

      prefer_candidate = .not. found .or. candidate_batch > batch
      if (found .and. candidate_batch == batch .and. has_index) then
        prefer_candidate = candidate_slot == indexed_slot .and. candidate_batch == indexed_batch
      end if
      if (.not. prefer_candidate) cycle
      slot = candidate_slot
      batch = candidate_batch
      found = .true.
    end do
  end subroutine find_latest_periodic_slot

  subroutine read_checkpoint_index(base_dir, active_slot, batch, found)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(out) :: active_slot, batch
    logical, intent(out) :: found

    character(len=1024) :: path
    character(len=256) :: line
    integer :: u, ios, pos
    integer(i32) :: schema
    logical :: file_exists, has_schema, has_slot, has_batch

    path = trim(base_dir)//'/'//checkpoint_index_name
    inquire (file=trim(path), exist=file_exists)
    found = .false.
    active_slot = -1_i32
    batch = -1_i32
    if (.not. file_exists) return

    schema = -1_i32
    has_schema = .false.
    has_slot = .false.
    has_batch = .false.
    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      pos = index(line, '=')
      if (pos <= 1) cycle
      select case (trim(line(:pos - 1)))
      case ('schema_version')
        read (line(pos + 1:), *, iostat=ios) schema
        has_schema = ios == 0
      case ('slot')
        read (line(pos + 1:), *, iostat=ios) active_slot
        has_slot = ios == 0
      case ('batches')
        read (line(pos + 1:), *, iostat=ios) batch
        has_batch = ios == 0
      end select
    end do
    close (u)
    if (.not. has_schema .or. schema /= periodic_checkpoint_index_schema .or. &
        .not. has_slot .or. (active_slot /= 0_i32 .and. active_slot /= 1_i32) .or. &
        .not. has_batch .or. batch < 0_i32) then
      active_slot = -1_i32
      batch = -1_i32
      return
    end if
    found = .true.
  end subroutine read_checkpoint_index

  subroutine publish_checkpoint_index(base_dir, slot, batch)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(in) :: slot, batch

    character(len=1024) :: path, temporary_path
    integer :: u, ios, rename_status

    call ensure_output_dir(trim(base_dir))
    path = trim(base_dir)//'/'//checkpoint_index_name
    temporary_path = trim(path)//'.tmp'
    open (newunit=u, file=trim(temporary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to create temporary periodic checkpoint index.'
    write (u, '(a,i0)') 'schema_version=', periodic_checkpoint_index_schema
    write (u, '(a,i0)') 'slot=', slot
    write (u, '(a,i0)') 'batches=', batch
    close (u)
    call atomic_rename(trim(temporary_path), trim(path), rename_status)
    if (rename_status /= filesystem_success) error stop 'Failed to publish periodic checkpoint index.'
  end subroutine publish_checkpoint_index

  function checkpoint_slot_dir(base_dir, slot) result(path)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(in) :: slot
    character(len=1024) :: path
    character(len=1) :: slot_text

    write (slot_text, '(i1)') slot
    path = trim(base_dir)//'/'//checkpoint_parent_name//'/slot'//slot_text
  end function checkpoint_slot_dir

end module bem_periodic_checkpoint
