!> accepted batch 境界で二重化した定期チェックポイントを保存・解決する。
module bem_periodic_checkpoint
  use bem_kinds, only: i32
  use bem_types, only: injection_state, mesh_type, sim_stats
  use bem_app_config_types, only: app_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_mpi, only: mpi_bcast_i32_array, mpi_context, mpi_is_root, mpi_world_barrier, mpi_world_size
  use bem_output_writer, only: ensure_output_dir, write_checkpoint_state_files
  use bem_restart, only: restart_rng_state_path, write_macro_residuals_file, write_rng_state_file
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

    ! index は全rankのファイルが閉じられた後にだけ公開する。書込み中断時は旧slotがactiveのまま残る。
    call mpi_world_barrier(mpi)
    if (mpi_is_root(mpi)) then
      call publish_checkpoint_index(trim(app%output_dir), inactive_slot, stats%batches)
      print '(a,i0,a,a)', 'periodic_checkpoint_batch=', stats%batches, ' dir=', trim(checkpoint_dir)
    end if
    call mpi_world_barrier(mpi)
  end subroutine maybe_write_periodic_checkpoint

  !> base directory直下の最終出力と定期slotを比較し、batch数が新しい完全checkpointを返す。
  subroutine resolve_latest_checkpoint_dir(base_dir, checkpoint_dir)
    character(len=*), intent(in) :: base_dir
    character(len=*), intent(out) :: checkpoint_dir

    integer(i32) :: active_slot, indexed_batch, root_batch, slot_batch
    logical :: has_index
    character(len=1024) :: slot_dir

    checkpoint_dir = trim(base_dir)
    call read_summary_batch(trim(base_dir)//'/summary.txt', root_batch)
    if (.not. checkpoint_required_files_complete(trim(base_dir))) root_batch = -1_i32
    call read_checkpoint_index(trim(base_dir), active_slot, indexed_batch, has_index)
    if (.not. has_index) return

    slot_dir = checkpoint_slot_dir(trim(base_dir), active_slot)
    call read_summary_batch(trim(slot_dir)//'/summary.txt', slot_batch)
    if (slot_batch /= indexed_batch) then
      error stop 'Periodic checkpoint index does not match its slot summary.'
    end if
    if (.not. checkpoint_required_files_complete(trim(slot_dir))) then
      error stop 'Periodic checkpoint index points to an incomplete slot.'
    end if
    if (slot_batch > root_batch) checkpoint_dir = trim(slot_dir)
  end subroutine resolve_latest_checkpoint_dir

  subroutine read_active_slot(base_dir, active_slot)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(out) :: active_slot
    integer(i32) :: indexed_batch
    logical :: has_index

    call read_checkpoint_index(base_dir, active_slot, indexed_batch, has_index)
    if (.not. has_index) active_slot = 1_i32
  end subroutine read_active_slot

  subroutine read_checkpoint_index(base_dir, active_slot, batch, found)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(out) :: active_slot, batch
    logical, intent(out) :: found

    character(len=1024) :: path
    character(len=256) :: line
    integer :: u, ios, pos
    integer(i32) :: schema
    logical :: has_schema, has_slot, has_batch

    path = trim(base_dir)//'/'//checkpoint_index_name
    inquire (file=trim(path), exist=found)
    active_slot = -1_i32
    batch = -1_i32
    if (.not. found) return

    schema = -1_i32
    has_schema = .false.
    has_slot = .false.
    has_batch = .false.
    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open periodic checkpoint index.'
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
      error stop 'Periodic checkpoint index is malformed or unsupported.'
    end if
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

  subroutine read_summary_batch(path, batch)
    character(len=*), intent(in) :: path
    integer(i32), intent(out) :: batch

    character(len=512) :: line
    integer :: u, ios
    logical :: exists

    batch = -1_i32
    inquire (file=trim(path), exist=exists)
    if (.not. exists) return
    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, 'batches=') /= 1) cycle
      read (line(len('batches=') + 1:), *, iostat=ios) batch
      if (ios /= 0) batch = -1_i32
      exit
    end do
    close (u)
  end subroutine read_summary_batch

  logical function checkpoint_required_files_complete(checkpoint_dir) result(complete)
    character(len=*), intent(in) :: checkpoint_dir

    character(len=1024) :: path
    integer(i32) :: rank, world_size
    logical :: exists

    complete = .false.
    inquire (file=trim(checkpoint_dir)//'/summary.txt', exist=exists)
    if (.not. exists) return
    inquire (file=trim(checkpoint_dir)//'/charges.csv', exist=exists)
    if (.not. exists) return
    call read_summary_world_size(trim(checkpoint_dir)//'/summary.txt', world_size)
    if (world_size <= 0_i32) return
    do rank = 0_i32, world_size - 1_i32
      path = restart_rng_state_path(trim(checkpoint_dir), mpi_rank=rank, mpi_size=world_size)
      inquire (file=trim(path), exist=exists)
      if (.not. exists) return
    end do
    complete = .true.
  end function checkpoint_required_files_complete

  subroutine read_summary_world_size(path, world_size)
    character(len=*), intent(in) :: path
    integer(i32), intent(out) :: world_size

    character(len=512) :: line
    integer :: u, ios

    world_size = -1_i32
    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, 'mpi_world_size=') /= 1) cycle
      read (line(len('mpi_world_size=') + 1:), *, iostat=ios) world_size
      if (ios /= 0) world_size = -1_i32
      exit
    end do
    close (u)
  end subroutine read_summary_world_size

  function checkpoint_slot_dir(base_dir, slot) result(path)
    character(len=*), intent(in) :: base_dir
    integer(i32), intent(in) :: slot
    character(len=1024) :: path
    character(len=1) :: slot_text

    write (slot_text, '(i1)') slot
    path = trim(base_dir)//'/'//checkpoint_parent_name//'/slot'//slot_text
  end function checkpoint_slot_dir

end module bem_periodic_checkpoint
