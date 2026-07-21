!> Rank-local delayed outer-event queue checkpoint I/O.
module bem_outer_event_queue_io
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_allreduce_sum_i64_array, &
                     mpi_allreduce_sum_real_dp_scalar
  use bem_outer_event_queue, only: outer_event_queue_type, outer_event_record_type, &
                                   outer_event_outcome_return, outer_event_outcome_escape, &
                                   outer_event_queue_global_fingerprint, &
                                   outer_event_queue_fingerprint_is_valid
  implicit none
  private

  integer(i32), parameter, public :: outer_event_queue_checkpoint_schema = 2_i32
  integer(i32), parameter, public :: outer_event_queue_io_ok = 0_i32
  integer(i32), parameter, public :: outer_event_queue_io_error = 1_i32
  integer(i32), parameter, public :: outer_event_queue_unsupported_schema = 2_i32
  integer(i32), parameter, public :: outer_event_queue_contract_mismatch = 3_i32
  integer(i32), parameter, public :: outer_event_queue_malformed = 4_i32

  character(len=*), parameter :: checkpoint_header = &
                                 'event_id,queued_time,due_time,outcome,species_id,origin_rank,origin_batch,'// &
                                 'origin_particle_id,interface_face_index,q,m,w,x1,x2,x3,v1,v2,v3'

  public :: outer_event_queue_checkpoint_path
  public :: write_outer_event_queue_checkpoint
  public :: load_outer_event_queue_checkpoint
  public :: validate_outer_event_queue_inventory

contains

  !> Return the serial or rank-local checkpoint path.
  function outer_event_queue_checkpoint_path(out_dir, mpi_rank, mpi_size, mpi) result(path)
    character(len=*), intent(in) :: out_dir
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    character(len=1024) :: path
    integer(i32) :: local_rank, world_size

    call resolve_parallel_rank_size( &
      local_rank, world_size, mpi_rank, mpi_size, mpi, 'outer_event_queue_checkpoint_path' &
      )
    if (world_size <= 1_i32) then
      path = trim(out_dir)//'/outer_event_queue.csv'
    else
      write (path, '(a,a,i5.5,a)') trim(out_dir), '/outer_event_queue_rank', local_rank, '.csv'
    end if
  end function outer_event_queue_checkpoint_path

  !> Atomically publish one rank's active queue records and ID sequence.
  subroutine write_outer_event_queue_checkpoint( &
    out_dir, queue, completed_batches, status, message, mpi_rank, mpi_size, mpi &
    )
    character(len=*), intent(in) :: out_dir
    type(outer_event_queue_type), intent(in) :: queue
    integer(i32), intent(in) :: completed_batches
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi

    type(outer_event_record_type), allocatable :: events(:)
    integer(i64) :: next_event_id
    integer(i32) :: local_rank, world_size
    integer :: unit, ios, rename_status, idx
    character(len=1024) :: path, temporary_path
    character(len=256) :: io_message, validation_message
    character(len=16) :: local_fingerprint
    logical :: valid

    status = outer_event_queue_io_ok
    message = ''
    if (completed_batches < 0_i32) then
      status = outer_event_queue_malformed
      message = 'completed_batches must be nonnegative'
      return
    end if

    call resolve_parallel_rank_size( &
      local_rank, world_size, mpi_rank, mpi_size, mpi, 'write_outer_event_queue_checkpoint' &
      )
    call queue%snapshot(events, next_event_id)
    call validate_checkpoint_events( &
      events, next_event_id, world_size, valid, validation_message &
      )
    if (.not. valid) then
      status = outer_event_queue_malformed
      message = 'active outer-event queue is not checkpointable: '//trim(validation_message)
      return
    end if
    local_fingerprint = queue%fingerprint(local_rank)

    path = outer_event_queue_checkpoint_path( &
           trim(out_dir), mpi_rank=local_rank, mpi_size=world_size &
           )
    temporary_path = trim(path)//'.tmp'
    open ( &
      newunit=unit, file=trim(temporary_path), status='replace', action='write', &
      form='formatted', iostat=ios, iomsg=io_message &
      )
    if (ios /= 0) then
      status = outer_event_queue_io_error
      message = 'failed to open queue checkpoint temporary file: '//trim(io_message)
      return
    end if

    write (unit, '(a,i0)', iostat=ios, iomsg=io_message) &
      'outer_event_queue_schema=', outer_event_queue_checkpoint_schema
    if (ios == 0) write (unit, '(a,i0)', iostat=ios, iomsg=io_message) 'mpi_world_size=', world_size
    if (ios == 0) write (unit, '(a,i0)', iostat=ios, iomsg=io_message) 'mpi_rank=', local_rank
    if (ios == 0) write (unit, '(a,i0)', iostat=ios, iomsg=io_message) 'completed_batches=', completed_batches
    if (ios == 0) write (unit, '(a,i0)', iostat=ios, iomsg=io_message) 'next_event_id=', next_event_id
    if (ios == 0) write (unit, '(a,i0)', iostat=ios, iomsg=io_message) 'event_count=', size(events)
    if (ios == 0) write (unit, '(a,a)', iostat=ios, iomsg=io_message) &
      'local_fingerprint=', local_fingerprint
    if (ios == 0) write (unit, '(a)', iostat=ios, iomsg=io_message) checkpoint_header
    do idx = 1, size(events)
      if (ios /= 0) exit
      call write_event_record(unit, events(idx), ios, io_message)
    end do
    if (ios /= 0) then
      close (unit)
      call delete_file_if_exists(temporary_path)
      status = outer_event_queue_io_error
      message = 'failed to write queue checkpoint: '//trim(io_message)
      return
    end if

    close (unit, iostat=ios, iomsg=io_message)
    if (ios /= 0) then
      call delete_file_if_exists(temporary_path)
      status = outer_event_queue_io_error
      message = 'failed to close queue checkpoint: '//trim(io_message)
      return
    end if
    call atomic_rename(trim(temporary_path), trim(path), rename_status)
    if (rename_status /= filesystem_success) then
      call delete_file_if_exists(temporary_path)
      status = outer_event_queue_io_error
      message = 'failed to atomically publish queue checkpoint'
    end if
  end subroutine write_outer_event_queue_checkpoint

  !> Load one rank's queue after validating schema and restart identity.
  !! A missing file returns found=.false. and status=outer_event_queue_io_ok.
  !! Any present but invalid file leaves queue empty and returns a nonzero status.
  subroutine load_outer_event_queue_checkpoint( &
    out_dir, expected_completed_batches, queue, found, status, message, mpi_rank, mpi_size, mpi &
    )
    character(len=*), intent(in) :: out_dir
    integer(i32), intent(in) :: expected_completed_batches
    type(outer_event_queue_type), intent(out) :: queue
    logical, intent(out) :: found
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi

    type(outer_event_record_type), allocatable :: events(:)
    integer(i32) :: local_rank, world_size, saved_schema, saved_world_size, saved_rank
    integer(i32) :: saved_completed_batches, event_count
    integer(i64) :: next_event_id, file_size
    integer :: unit, ios, allocation_status, idx
    character(len=1024) :: path
    character(len=512) :: line
    character(len=256) :: io_message, validation_message
    character(len=16) :: saved_local_fingerprint, actual_local_fingerprint
    logical :: exists, valid

    call queue%init()
    found = .false.
    status = outer_event_queue_io_ok
    message = ''
    if (expected_completed_batches < 0_i32) then
      status = outer_event_queue_malformed
      message = 'expected_completed_batches must be nonnegative'
      return
    end if

    call resolve_parallel_rank_size( &
      local_rank, world_size, mpi_rank, mpi_size, mpi, 'load_outer_event_queue_checkpoint' &
      )
    path = outer_event_queue_checkpoint_path( &
           trim(out_dir), mpi_rank=local_rank, mpi_size=world_size &
           )
    inquire (file=trim(path), exist=exists, size=file_size)
    if (.not. exists) return
    found = .true.

    open ( &
      newunit=unit, file=trim(path), status='old', action='read', form='formatted', &
      iostat=ios, iomsg=io_message &
      )
    if (ios /= 0) then
      status = outer_event_queue_io_error
      message = 'failed to open queue checkpoint: '//trim(io_message)
      return
    end if

    call read_i32_metadata(unit, 'outer_event_queue_schema', saved_schema, ios)
    if (ios /= 0) then
      call fail_loaded_file(unit, outer_event_queue_malformed, 'missing or invalid schema metadata', status, message)
      return
    end if
    if (saved_schema /= outer_event_queue_checkpoint_schema) then
      call fail_loaded_file( &
        unit, outer_event_queue_unsupported_schema, 'unsupported outer-event queue checkpoint schema', status, message &
        )
      return
    end if

    call read_i32_metadata(unit, 'mpi_world_size', saved_world_size, ios)
    if (ios == 0) call read_i32_metadata(unit, 'mpi_rank', saved_rank, ios)
    if (ios == 0) call read_i32_metadata(unit, 'completed_batches', saved_completed_batches, ios)
    if (ios == 0) call read_i64_metadata(unit, 'next_event_id', next_event_id, ios)
    if (ios == 0) call read_i32_metadata(unit, 'event_count', event_count, ios)
    if (ios == 0) call read_character_metadata(unit, 'local_fingerprint', saved_local_fingerprint, ios)
    if (ios /= 0) then
      call fail_loaded_file(unit, outer_event_queue_malformed, 'missing or invalid queue metadata', status, message)
      return
    end if
    if (.not. outer_event_queue_fingerprint_is_valid(saved_local_fingerprint)) then
      call fail_loaded_file(unit, outer_event_queue_malformed, 'invalid local queue fingerprint metadata', status, message)
      return
    end if
    if (saved_world_size <= 0_i32 .or. saved_rank < 0_i32 .or. saved_rank >= saved_world_size) then
      call fail_loaded_file(unit, outer_event_queue_malformed, 'invalid saved MPI rank/size', status, message)
      return
    end if
    if (saved_world_size /= world_size) then
      call fail_loaded_file(unit, outer_event_queue_contract_mismatch, 'MPI world size differs', status, message)
      return
    end if
    if (saved_rank /= local_rank) then
      call fail_loaded_file(unit, outer_event_queue_contract_mismatch, 'MPI rank differs', status, message)
      return
    end if
    if (saved_completed_batches /= expected_completed_batches) then
      call fail_loaded_file(unit, outer_event_queue_contract_mismatch, 'completed batch count differs', status, message)
      return
    end if
    if (next_event_id < 1_i64 .or. event_count < 0_i32 .or. int(event_count, i64) > file_size) then
      call fail_loaded_file(unit, outer_event_queue_malformed, 'invalid event count or next event ID', status, message)
      return
    end if

    read (unit, '(A)', iostat=ios, iomsg=io_message) line
    if (ios /= 0) then
      call fail_loaded_file(unit, outer_event_queue_malformed, 'missing or invalid event table header', status, message)
      return
    end if
    if (trim(line) /= checkpoint_header) then
      call fail_loaded_file(unit, outer_event_queue_malformed, 'missing or invalid event table header', status, message)
      return
    end if
    allocate (events(event_count), stat=allocation_status)
    if (allocation_status /= 0) then
      call fail_loaded_file(unit, outer_event_queue_io_error, 'failed to allocate queue checkpoint records', status, message)
      return
    end if
    do idx = 1, event_count
      read (unit, '(A)', iostat=ios, iomsg=io_message) line
      if (ios /= 0) then
        call fail_loaded_file(unit, outer_event_queue_malformed, 'truncated or malformed event record', status, message)
        return
      end if
      if (count_character(line, ',') /= 17_i32) then
        call fail_loaded_file(unit, outer_event_queue_malformed, 'truncated or malformed event record', status, message)
        return
      end if
      call read_event_record(line, events(idx), ios)
      if (ios /= 0) then
        call fail_loaded_file(unit, outer_event_queue_malformed, 'invalid event record value', status, message)
        return
      end if
    end do

    do
      read (unit, '(A)', iostat=ios) line
      if (ios < 0) exit
      if (ios > 0) then
        call fail_loaded_file(unit, outer_event_queue_io_error, 'failed while checking event table end', status, message)
        return
      end if
      if (len_trim(line) > 0) then
        call fail_loaded_file(unit, outer_event_queue_malformed, 'unexpected trailing checkpoint data', status, message)
        return
      end if
    end do
    close (unit)

    call validate_checkpoint_events( &
      events, next_event_id, world_size, valid, validation_message &
      )
    if (.not. valid) then
      status = outer_event_queue_malformed
      message = trim(validation_message)
      return
    end if
    call queue%restore(events, next_event_id)
    actual_local_fingerprint = queue%fingerprint(local_rank)
    if (actual_local_fingerprint /= saved_local_fingerprint) then
      call queue%init()
      status = outer_event_queue_contract_mismatch
      message = 'local outer-event queue fingerprint differs from checkpoint content'
    end if
  end subroutine load_outer_event_queue_checkpoint

  !> Compare the restored global queue inventory with the summary checkpoint.
  !! Event count is exact; signed charge allows only floating-point roundoff from
  !! checkpoint text serialization and the repeated MPI reduction.
  subroutine validate_outer_event_queue_inventory( &
    queue, expected_event_count, expected_signed_charge, expected_fingerprint, status, message, &
    actual_event_count, actual_signed_charge, actual_fingerprint, mpi &
    )
    type(outer_event_queue_type), intent(in) :: queue
    integer(i64), intent(in) :: expected_event_count
    real(dp), intent(in) :: expected_signed_charge
    character(len=*), intent(in) :: expected_fingerprint
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer(i64), intent(out), optional :: actual_event_count
    real(dp), intent(out), optional :: actual_signed_charge
    character(len=*), intent(out), optional :: actual_fingerprint
    type(mpi_context), intent(in), optional :: mpi

    integer(i64) :: event_count_buffer(1), fingerprint_components(2)
    integer(i32) :: local_rank, world_size
    real(dp) :: global_signed_charge, charge_scale, charge_tolerance
    character(len=16) :: global_fingerprint

    status = outer_event_queue_io_ok
    message = ''
    call mpi_get_rank_size(local_rank, world_size, mpi)
    event_count_buffer(1) = int(queue%size(), i64)
    global_signed_charge = queue%signed_charge()
    call queue%fingerprint_components(local_rank, fingerprint_components)
    if (present(mpi)) then
      call mpi_allreduce_sum_i64_array(mpi, event_count_buffer)
      call mpi_allreduce_sum_real_dp_scalar(mpi, global_signed_charge)
      call mpi_allreduce_sum_i64_array(mpi, fingerprint_components)
    end if
    global_fingerprint = outer_event_queue_global_fingerprint(fingerprint_components, world_size)
    if (present(actual_event_count)) actual_event_count = event_count_buffer(1)
    if (present(actual_signed_charge)) actual_signed_charge = global_signed_charge
    if (present(actual_fingerprint)) actual_fingerprint = global_fingerprint

    if (expected_event_count < 0_i64 .or. .not. ieee_is_finite(expected_signed_charge) .or. &
        .not. outer_event_queue_fingerprint_is_valid(expected_fingerprint)) then
      status = outer_event_queue_malformed
      message = 'summary outer-event queue inventory is invalid'
      return
    end if
    if (.not. ieee_is_finite(global_signed_charge)) then
      status = outer_event_queue_malformed
      message = 'restored outer-event queue signed charge is not finite'
      return
    end if
    if (event_count_buffer(1) /= expected_event_count) then
      status = outer_event_queue_contract_mismatch
      write (message, '(a,i0,a,i0)') &
        'global outer-event queue count differs: checkpoint=', event_count_buffer(1), &
        ' summary=', expected_event_count
      return
    end if

    charge_scale = max(tiny(1.0_dp), abs(global_signed_charge), abs(expected_signed_charge))
    charge_tolerance = 512.0_dp*epsilon(1.0_dp)*charge_scale
    if (abs(global_signed_charge - expected_signed_charge) > charge_tolerance) then
      status = outer_event_queue_contract_mismatch
      write (message, '(a,es24.16,a,es24.16,a,es24.16)') &
        'global outer-event queue signed charge differs: checkpoint=', global_signed_charge, &
        ' summary=', expected_signed_charge, ' tolerance=', charge_tolerance
      return
    end if
    if (global_fingerprint /= expected_fingerprint) then
      status = outer_event_queue_contract_mismatch
      write (message, '(a,a,a,a)') &
        'global outer-event queue fingerprint differs: checkpoint=', global_fingerprint, &
        ' summary=', expected_fingerprint
    end if
  end subroutine validate_outer_event_queue_inventory

  subroutine write_event_record(unit, event, ios, io_message)
    integer, intent(in) :: unit
    type(outer_event_record_type), intent(in) :: event
    integer, intent(out) :: ios
    character(len=*), intent(out) :: io_message

    write ( &
      unit, &
      '(i0,2(",",es26.17e3),6(",",i0),9(",",es26.17e3))', &
      iostat=ios, iomsg=io_message &
      ) event%event_id, event%queued_time, event%due_time, event%outcome, event%species_id, &
      event%origin_rank, event%origin_batch, event%origin_particle_id, event%interface_face_index, &
      event%q, event%m, event%w, event%x, event%v
  end subroutine write_event_record

  subroutine read_event_record(line, event, ios)
    character(len=*), intent(in) :: line
    type(outer_event_record_type), intent(out) :: event
    integer, intent(out) :: ios

    event = outer_event_record_type()
    read (line, *, iostat=ios) event%event_id, event%queued_time, event%due_time, event%outcome, &
      event%species_id, event%origin_rank, event%origin_batch, event%origin_particle_id, &
      event%interface_face_index, event%q, event%m, event%w, event%x, event%v
  end subroutine read_event_record

  subroutine read_i32_metadata(unit, expected_key, value, ios)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: expected_key
    integer(i32), intent(out) :: value
    integer, intent(out) :: ios
    character(len=512) :: line
    integer :: separator

    value = 0_i32
    read (unit, '(A)', iostat=ios) line
    if (ios /= 0) return
    separator = index(line, '=')
    if (separator <= 1) then
      ios = 1
      return
    end if
    if (trim(line(:separator - 1)) /= expected_key) then
      ios = 1
      return
    end if
    read (line(separator + 1:), *, iostat=ios) value
  end subroutine read_i32_metadata

  subroutine read_i64_metadata(unit, expected_key, value, ios)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: expected_key
    integer(i64), intent(out) :: value
    integer, intent(out) :: ios
    character(len=512) :: line
    integer :: separator

    value = 0_i64
    read (unit, '(A)', iostat=ios) line
    if (ios /= 0) return
    separator = index(line, '=')
    if (separator <= 1) then
      ios = 1
      return
    end if
    if (trim(line(:separator - 1)) /= expected_key) then
      ios = 1
      return
    end if
    read (line(separator + 1:), *, iostat=ios) value
  end subroutine read_i64_metadata

  subroutine read_character_metadata(unit, expected_key, value, ios)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: expected_key
    character(len=*), intent(out) :: value
    integer, intent(out) :: ios
    character(len=512) :: line
    integer :: separator

    value = ''
    read (unit, '(A)', iostat=ios) line
    if (ios /= 0) return
    separator = index(line, '=')
    if (separator <= 1) then
      ios = 1
      return
    end if
    if (trim(line(:separator - 1)) /= expected_key) then
      ios = 1
      return
    end if
    value = trim(adjustl(line(separator + 1:)))
    if (len_trim(line(separator + 1:)) > len(value)) ios = 1
  end subroutine read_character_metadata

  subroutine validate_checkpoint_events(events, next_event_id, world_size, valid, message)
    type(outer_event_record_type), intent(in) :: events(:)
    integer(i64), intent(in) :: next_event_id
    integer(i32), intent(in) :: world_size
    logical, intent(out) :: valid
    character(len=*), intent(out) :: message
    integer(i64), allocatable :: ids(:)
    integer(i64) :: maximum_id
    integer :: idx

    valid = .false.
    message = ''
    if (next_event_id < 1_i64) then
      message = 'next event ID must be positive'
      return
    end if
    maximum_id = 0_i64
    allocate (ids(size(events)))
    do idx = 1, size(events)
      if (.not. event_payload_is_valid(events(idx), world_size)) then
        message = 'invalid active event payload'
        return
      end if
      if (idx > 1) then
        if (event_precedes(events(idx), events(idx - 1))) then
          message = 'active events are not ordered by due time and event ID'
          return
        end if
      end if
      ids(idx) = events(idx)%event_id
      maximum_id = max(maximum_id, events(idx)%event_id)
    end do
    call sort_i64(ids)
    do idx = 2, size(ids)
      if (ids(idx) == ids(idx - 1)) then
        message = 'active event IDs are not unique'
        return
      end if
    end do
    if (next_event_id <= maximum_id) then
      message = 'next event ID does not exceed all active IDs'
      return
    end if
    valid = .true.
  end subroutine validate_checkpoint_events

  logical function event_payload_is_valid(event, world_size) result(valid)
    type(outer_event_record_type), intent(in) :: event
    integer(i32), intent(in) :: world_size

    valid = event%event_id > 0_i64
    valid = valid .and. ieee_is_finite(event%queued_time) .and. ieee_is_finite(event%due_time)
    valid = valid .and. event%queued_time >= 0.0_dp .and. event%due_time > event%queued_time
    valid = valid .and. ( &
            event%outcome == outer_event_outcome_return .or. event%outcome == outer_event_outcome_escape &
            )
    valid = valid .and. event%species_id > 0_i32
    valid = valid .and. event%origin_rank >= 0_i32 .and. event%origin_rank < world_size
    valid = valid .and. event%origin_batch > 0_i32 .and. event%origin_particle_id > 0_i64
    valid = valid .and. ieee_is_finite(event%q)
    valid = valid .and. ieee_is_finite(event%m) .and. event%m > 0.0_dp
    valid = valid .and. ieee_is_finite(event%w) .and. event%w > 0.0_dp
    valid = valid .and. all(ieee_is_finite(event%x)) .and. all(ieee_is_finite(event%v))
  end function event_payload_is_valid

  pure logical function event_precedes(lhs, rhs) result(precedes)
    type(outer_event_record_type), intent(in) :: lhs, rhs

    precedes = lhs%due_time < rhs%due_time .or. &
               (lhs%due_time == rhs%due_time .and. lhs%event_id < rhs%event_id)
  end function event_precedes

  subroutine sort_i64(values)
    integer(i64), intent(inout) :: values(:)
    integer(i64), allocatable :: work(:)
    integer :: width, left, middle, right, src_left, src_right, dst, count

    count = size(values)
    if (count < 2) return
    allocate (work(count))
    width = 1
    do while (width < count)
      left = 1
      do while (left <= count)
        middle = min(left + width, count + 1)
        right = min(left + 2*width - 1, count)
        src_left = left
        src_right = middle
        do dst = left, right
          if (src_left >= middle) then
            work(dst) = values(src_right)
            src_right = src_right + 1
          else if (src_right > right) then
            work(dst) = values(src_left)
            src_left = src_left + 1
          else if (values(src_right) < values(src_left)) then
            work(dst) = values(src_right)
            src_right = src_right + 1
          else
            work(dst) = values(src_left)
            src_left = src_left + 1
          end if
        end do
        left = right + 1
      end do
      values = work
      if (width > count/2) exit
      width = 2*width
    end do
  end subroutine sort_i64

  integer(i32) function count_character(text, character_to_count) result(character_count)
    character(len=*), intent(in) :: text
    character(len=1), intent(in) :: character_to_count
    integer :: idx

    character_count = 0_i32
    do idx = 1, len_trim(text)
      if (text(idx:idx) == character_to_count) character_count = character_count + 1_i32
    end do
  end function count_character

  subroutine fail_loaded_file(unit, failure_status, failure_message, status, message)
    integer, intent(in) :: unit
    integer(i32), intent(in) :: failure_status
    character(len=*), intent(in) :: failure_message
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    close (unit)
    status = failure_status
    message = failure_message
  end subroutine fail_loaded_file

  subroutine delete_file_if_exists(path)
    character(len=*), intent(in) :: path
    logical :: exists
    integer :: unit, ios

    inquire (file=trim(path), exist=exists)
    if (.not. exists) return
    open (newunit=unit, file=trim(path), status='old', action='readwrite', iostat=ios)
    if (ios == 0) close (unit, status='delete')
  end subroutine delete_file_if_exists

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

end module bem_outer_event_queue_io
