!> BEACH checkpoint metadata and transactional publication shared by output and restart code.
module bem_checkpoint_contract
  use bem_kinds, only: i32
  use bem_filesystem, only: atomic_rename, filesystem_success
  implicit none
  private

  integer(i32), parameter, public :: checkpoint_schema_version_current = 9_i32
  integer(i32), parameter :: checkpoint_manifest_schema_current = 1_i32
  integer(i32), parameter :: checkpoint_manifest_required_schema = 8_i32
  character(len=*), parameter :: checkpoint_manifest_name = 'checkpoint_complete.txt'

  public :: begin_checkpoint_publish
  public :: checkpoint_schema_is_loadable
  public :: publish_checkpoint_manifest
  public :: inspect_checkpoint_directory

contains

  !> 現行実装が再開候補として扱えるcheckpoint schemaかを返す。
  !! schema keyを持たないlegacy summaryは0として表す。schema v1は正式な
  !! fingerprint contractより前の中間形式であり、現行loaderは対応しない。
  logical function checkpoint_schema_is_loadable(schema_version) result(loadable)
    integer(i32), intent(in) :: schema_version

    loadable = schema_version == 0_i32 .or. &
               (schema_version >= 2_i32 .and. schema_version <= checkpoint_schema_version_current)
  end function checkpoint_schema_is_loadable

  !> schema v9 transactionのcompletion evidenceをrestart-bearing fileの置換前に無効化する。
  subroutine begin_checkpoint_publish(out_dir)
    character(len=*), intent(in) :: out_dir

    character(len=1024) :: path, temporary_path
    integer :: u, ios, rename_status

    path = trim(out_dir)//'/'//checkpoint_manifest_name
    temporary_path = trim(path)//'.tmp'
    open (newunit=u, file=trim(temporary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to create temporary checkpoint manifest.'
    write (u, '(a,i0)') 'schema_version=', checkpoint_manifest_schema_current
    write (u, '(a)') 'state=in_progress'
    close (u, iostat=ios)
    if (ios /= 0) error stop 'Failed to close temporary checkpoint manifest.'
    call atomic_rename(trim(temporary_path), trim(path), rename_status)
    if (rename_status /= filesystem_success) error stop 'Failed to invalidate checkpoint manifest.'
  end subroutine begin_checkpoint_publish

  !> A complete manifest is atomically published only after every restart-bearing file is closed.
  subroutine publish_checkpoint_manifest( &
    out_dir, batches, mpi_world_size, has_macro_residuals, has_charge_ledger &
    )
    character(len=*), intent(in) :: out_dir
    integer(i32), intent(in) :: batches, mpi_world_size
    logical, intent(in) :: has_macro_residuals, has_charge_ledger

    character(len=1024) :: path, temporary_path
    integer :: u, ios, rename_status

    if (batches < 0_i32) error stop 'Checkpoint manifest batches must be >= 0.'
    if (mpi_world_size <= 0_i32) error stop 'Checkpoint manifest mpi_world_size must be > 0.'

    path = trim(out_dir)//'/'//checkpoint_manifest_name
    temporary_path = trim(path)//'.tmp'
    open (newunit=u, file=trim(temporary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to create temporary checkpoint manifest.'
    write (u, '(a,i0)') 'schema_version=', checkpoint_manifest_schema_current
    write (u, '(a)') 'state=complete'
    write (u, '(a,i0)') 'batches=', batches
    write (u, '(a,i0)') 'mpi_world_size=', mpi_world_size
    write (u, '(a,l1)') 'macro_residuals_present=', has_macro_residuals
    write (u, '(a,l1)') 'charge_ledger_present=', has_charge_ledger
    close (u, iostat=ios)
    if (ios /= 0) error stop 'Failed to close temporary checkpoint manifest.'
    call atomic_rename(trim(temporary_path), trim(path), rename_status)
    if (rename_status /= filesystem_success) error stop 'Failed to publish checkpoint manifest.'
  end subroutine publish_checkpoint_manifest

  !> Inspect one directory without mutating it and require every file declared by its checkpoint generation.
  subroutine inspect_checkpoint_directory( &
    checkpoint_dir, complete, schema_version, batches, mpi_world_size, &
    has_macro_residuals, has_charge_ledger &
    )
    character(len=*), intent(in) :: checkpoint_dir
    logical, intent(out) :: complete
    integer(i32), intent(out), optional :: schema_version, batches, mpi_world_size
    logical, intent(out), optional :: has_macro_residuals, has_charge_ledger

    character(len=1024) :: path
    integer(i32) :: saved_schema, saved_batches, saved_world_size, rank
    integer(i32) :: manifest_batches, manifest_world_size
    logical :: summary_ok, summary_has_ledger, manifest_found, manifest_complete
    logical :: manifest_has_residuals, manifest_has_ledger, exists

    complete = .false.
    saved_schema = -1_i32
    saved_batches = -1_i32
    saved_world_size = -1_i32
    manifest_has_residuals = .false.
    manifest_has_ledger = .false.
    summary_has_ledger = .false.

    path = trim(checkpoint_dir)//'/summary.txt'
    inquire (file=trim(path), exist=exists)
    if (.not. exists) then
      call assign_inspection_outputs()
      return
    end if
    call read_summary_checkpoint_metadata( &
      trim(path), saved_schema, saved_batches, saved_world_size, summary_has_ledger, summary_ok &
      )
    if (.not. summary_ok) then
      call assign_inspection_outputs()
      return
    end if
    ! `load_restart_checkpoint` の直接利用でも、future schemaや明示schema v1を
    ! legacy checkpointとして誤読込しない。後続のworld-size分ファイル走査も回避する。
    if (.not. checkpoint_schema_is_loadable(saved_schema)) then
      call assign_inspection_outputs()
      return
    end if

    if (saved_schema >= checkpoint_manifest_required_schema) then
      call read_checkpoint_manifest( &
        trim(checkpoint_dir), manifest_found, manifest_complete, manifest_batches, manifest_world_size, &
        manifest_has_residuals, manifest_has_ledger &
        )
      if (.not. manifest_found .or. .not. manifest_complete) then
        call assign_inspection_outputs()
        return
      end if
      if (manifest_batches /= saved_batches .or. manifest_world_size /= saved_world_size) then
        call assign_inspection_outputs()
        return
      end if
      if (manifest_has_ledger .neqv. summary_has_ledger) then
        call assign_inspection_outputs()
        return
      end if
    else
      inquire (file=trim(checkpoint_dir)//'/macro_residuals.csv', exist=manifest_has_residuals)
      manifest_has_ledger = summary_has_ledger
    end if

    inquire (file=trim(checkpoint_dir)//'/charges.csv', exist=exists)
    if (.not. exists) then
      call assign_inspection_outputs()
      return
    end if
    do rank = 0_i32, saved_world_size - 1_i32
      path = checkpoint_rng_path(trim(checkpoint_dir), rank, saved_world_size)
      inquire (file=trim(path), exist=exists)
      if (.not. exists) then
        call assign_inspection_outputs()
        return
      end if
    end do
    if (manifest_has_residuals) then
      inquire (file=trim(checkpoint_dir)//'/macro_residuals.csv', exist=exists)
      if (.not. exists) then
        call assign_inspection_outputs()
        return
      end if
    end if
    if (summary_has_ledger) then
      inquire (file=trim(checkpoint_dir)//'/charge_ledger.csv', exist=exists)
      if (.not. exists) then
        call assign_inspection_outputs()
        return
      end if
    end if

    complete = .true.
    call assign_inspection_outputs()

  contains

    subroutine assign_inspection_outputs()
      if (present(schema_version)) schema_version = saved_schema
      if (present(batches)) batches = saved_batches
      if (present(mpi_world_size)) mpi_world_size = saved_world_size
      if (present(has_macro_residuals)) has_macro_residuals = manifest_has_residuals
      if (present(has_charge_ledger)) has_charge_ledger = summary_has_ledger
    end subroutine assign_inspection_outputs

  end subroutine inspect_checkpoint_directory

  subroutine read_summary_checkpoint_metadata(path, schema_version, batches, mpi_world_size, has_ledger, valid)
    character(len=*), intent(in) :: path
    integer(i32), intent(out) :: schema_version, batches, mpi_world_size
    logical, intent(out) :: has_ledger, valid

    character(len=512) :: line
    integer :: u, ios, pos
    logical :: found_schema, found_batches, found_world_size

    schema_version = 0_i32
    batches = -1_i32
    mpi_world_size = 1_i32
    has_ledger = .false.
    valid = .false.
    found_schema = .false.
    found_batches = .false.
    found_world_size = .false.

    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      pos = index(line, '=')
      if (pos <= 1) cycle
      select case (trim(adjustl(line(:pos - 1))))
      case ('checkpoint_schema_version')
        read (line(pos + 1:), *, iostat=ios) schema_version
        if (ios /= 0) then
          close (u)
          return
        end if
        found_schema = .true.
      case ('batches')
        read (line(pos + 1:), *, iostat=ios) batches
        if (ios /= 0) then
          close (u)
          return
        end if
        found_batches = .true.
      case ('mpi_world_size')
        read (line(pos + 1:), *, iostat=ios) mpi_world_size
        if (ios /= 0) then
          close (u)
          return
        end if
        found_world_size = .true.
      case default
        if (index(trim(adjustl(line(:pos - 1))), 'charge_ledger_') == 1) has_ledger = .true.
      end select
    end do
    close (u)

    ! 0はschema key自体がないlegacy summaryの内部sentinelとして予約する。
    if (found_schema .and. schema_version == 0_i32) schema_version = -1_i32
    if (.not. found_batches .or. batches < 0_i32) return
    if (schema_version >= 2_i32 .and. .not. found_world_size) return
    if (mpi_world_size <= 0_i32) return
    valid = .true.
  end subroutine read_summary_checkpoint_metadata

  subroutine read_checkpoint_manifest( &
    checkpoint_dir, found, complete, batches, mpi_world_size, has_macro_residuals, has_charge_ledger &
    )
    character(len=*), intent(in) :: checkpoint_dir
    logical, intent(out) :: found, complete, has_macro_residuals, has_charge_ledger
    integer(i32), intent(out) :: batches, mpi_world_size

    character(len=512) :: line
    character(len=64) :: key
    character(len=64) :: state
    integer :: u, ios, pos
    integer(i32) :: schema_version
    logical :: has_schema, has_state, has_batches, has_world_size, has_residual_flag, has_ledger_flag

    inquire (file=trim(checkpoint_dir)//'/'//checkpoint_manifest_name, exist=found)
    complete = .false.
    batches = -1_i32
    mpi_world_size = -1_i32
    has_macro_residuals = .false.
    has_charge_ledger = .false.
    if (.not. found) return

    schema_version = -1_i32
    state = ''
    has_schema = .false.
    has_state = .false.
    has_batches = .false.
    has_world_size = .false.
    has_residual_flag = .false.
    has_ledger_flag = .false.
    open (newunit=u, file=trim(checkpoint_dir)//'/'//checkpoint_manifest_name, status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      pos = index(line, '=')
      if (pos <= 1) cycle
      key = trim(adjustl(line(:pos - 1)))
      select case (trim(key))
      case ('schema_version')
        read (line(pos + 1:), *, iostat=ios) schema_version
        if (ios /= 0) then
          close (u)
          return
        end if
        has_schema = .true.
      case ('state')
        state = trim(adjustl(line(pos + 1:)))
        has_state = .true.
      case ('batches')
        read (line(pos + 1:), *, iostat=ios) batches
        if (ios /= 0) then
          close (u)
          return
        end if
        has_batches = .true.
      case ('mpi_world_size')
        read (line(pos + 1:), *, iostat=ios) mpi_world_size
        if (ios /= 0) then
          close (u)
          return
        end if
        has_world_size = .true.
      case ('macro_residuals_present')
        read (line(pos + 1:), *, iostat=ios) has_macro_residuals
        if (ios /= 0) then
          close (u)
          return
        end if
        has_residual_flag = .true.
      case ('charge_ledger_present')
        read (line(pos + 1:), *, iostat=ios) has_charge_ledger
        if (ios /= 0) then
          close (u)
          return
        end if
        has_ledger_flag = .true.
      end select
    end do
    close (u)

    if (.not. has_schema .or. schema_version /= checkpoint_manifest_schema_current) return
    if (.not. has_state .or. trim(state) /= 'complete') return
    if (.not. has_batches .or. batches < 0_i32) return
    if (.not. has_world_size .or. mpi_world_size <= 0_i32) return
    if (.not. has_residual_flag .or. .not. has_ledger_flag) return
    complete = .true.
  end subroutine read_checkpoint_manifest

  function checkpoint_rng_path(checkpoint_dir, rank, mpi_world_size) result(path)
    character(len=*), intent(in) :: checkpoint_dir
    integer(i32), intent(in) :: rank, mpi_world_size
    character(len=1024) :: path

    if (mpi_world_size <= 1_i32) then
      path = trim(checkpoint_dir)//'/rng_state.txt'
    else
      write (path, '(a,a,i5.5,a)') trim(checkpoint_dir), '/rng_state_rank', rank, '.txt'
    end if
  end function checkpoint_rng_path

end module bem_checkpoint_contract
