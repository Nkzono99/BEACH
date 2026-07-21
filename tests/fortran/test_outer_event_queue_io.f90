program test_outer_event_queue_io
  use bem_kinds, only: dp, i32, i64
  use bem_outer_event_queue, only: outer_event_queue_type, outer_event_record_type, &
                                   outer_event_outcome_return, outer_event_outcome_escape, &
                                   outer_event_queue_global_fingerprint
  use bem_outer_event_queue_io, only: &
    outer_event_queue_checkpoint_schema, outer_event_queue_io_ok, &
    outer_event_queue_unsupported_schema, outer_event_queue_contract_mismatch, &
    outer_event_queue_malformed, outer_event_queue_checkpoint_path, &
    write_outer_event_queue_checkpoint, load_outer_event_queue_checkpoint, &
    validate_outer_event_queue_inventory
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, &
                          assert_allclose_1d, delete_file_if_exists, ensure_directory, &
                          remove_empty_directory
  implicit none

  character(len=*), parameter :: out_dir = 'test_outer_event_queue_io_tmp'
  character(len=*), parameter :: table_header = &
                                 'event_id,queued_time,due_time,outcome,species_id,origin_rank,origin_batch,'// &
                                 'origin_particle_id,interface_face_index,q,m,w,x1,x2,x3,v1,v2,v3'
  type(outer_event_queue_type) :: queue, restored, stale_queue
  type(outer_event_record_type) :: event
  type(outer_event_record_type), allocatable :: expected(:), actual(:), popped(:)
  integer(i64) :: expected_next_id, actual_next_id, actual_event_count
  real(dp) :: actual_signed_charge
  integer(i32) :: status
  character(len=256) :: message
  character(len=512) :: fixture_row
  character(len=16) :: expected_fingerprint, actual_fingerprint, stale_fingerprint
  character(len=1024) :: serial_path, rank1_path, rank2_path
  logical :: found, exists
  integer :: idx

  serial_path = outer_event_queue_checkpoint_path(out_dir)
  rank1_path = outer_event_queue_checkpoint_path(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  rank2_path = outer_event_queue_checkpoint_path(out_dir, mpi_rank=2_i32, mpi_size=4_i32)
  call cleanup_files()
  call ensure_directory(out_dir)

  call test_init(14)

  call test_begin('serial_and_ranked_paths_and_missing_file')
  call assert_true( &
    trim(serial_path) == out_dir//'/outer_event_queue.csv', 'serial queue checkpoint path mismatch' &
    )
  call assert_true( &
    trim(rank1_path) == out_dir//'/outer_event_queue_rank00001.csv', 'ranked queue checkpoint path mismatch' &
    )
  call load_outer_event_queue_checkpoint( &
    out_dir, 0_i32, restored, found, status, message &
    )
  call assert_true(.not. found, 'missing queue checkpoint must report found=false')
  call assert_equal_i32(status, outer_event_queue_io_ok, 'missing queue checkpoint should not be an I/O error')
  call assert_equal_i32(restored%size(), 0_i32, 'missing queue checkpoint must leave an empty queue')
  call test_end()

  call test_begin('roundtrip_preserves_all_fields_sequence_and_atomic_publish')
  call queue%init(first_event_id=41_i64)
  event = make_event(3.0_dp, outer_event_outcome_return, 2_i32, 1_i32, 7_i32, 7001_i64, 13_i32, 1.0_dp)
  call queue%enqueue(event)
  event = make_event(2.0_dp, outer_event_outcome_escape, 3_i32, 1_i32, 6_i32, 6002_i64, 17_i32, 2.0_dp)
  call queue%enqueue(event)
  event = make_event(1.0_dp, outer_event_outcome_escape, 4_i32, 1_i32, 7_i32, 7003_i64, 19_i32, 3.0_dp)
  call queue%enqueue(event)
  call queue%pop_due(1.5_dp, popped)
  call assert_equal_i32(int(size(popped), i32), 1_i32, 'roundtrip fixture should pop one earlier event')
  call queue%snapshot(expected, expected_next_id)
  call write_outer_event_queue_checkpoint( &
    out_dir, queue, 7_i32, status, message, mpi_rank=1_i32, mpi_size=4_i32 &
    )
  call assert_equal_i32(status, outer_event_queue_io_ok, 'queue checkpoint write failed')
  inquire (file=trim(rank1_path), exist=exists)
  call assert_true(exists, 'final queue checkpoint was not published')
  inquire (file=trim(rank1_path)//'.tmp', exist=exists)
  call assert_true(.not. exists, 'atomic queue checkpoint temporary file remains')
  call load_outer_event_queue_checkpoint( &
    out_dir, 7_i32, restored, found, status, message, mpi_rank=1_i32, mpi_size=4_i32 &
    )
  call assert_true(found, 'published queue checkpoint was not found')
  call assert_equal_i32(status, outer_event_queue_io_ok, 'valid queue checkpoint did not load')
  call restored%snapshot(actual, actual_next_id)
  call assert_equal_i32(int(size(actual), i32), int(size(expected), i32), 'restored event count mismatch')
  call assert_equal_i64(actual_next_id, expected_next_id, 'next event ID was not preserved')
  do idx = 1, size(expected)
    call assert_record_equal(actual(idx), expected(idx))
  end do
  call test_end()

  call test_begin('global_inventory_match_is_accepted')
  expected_fingerprint = queue_global_fingerprint(restored, 0_i32, 1_i32)
  call validate_outer_event_queue_inventory( &
    restored, int(queue%size(), i64), queue%signed_charge(), expected_fingerprint, status, message, &
    actual_event_count=actual_event_count, actual_signed_charge=actual_signed_charge, &
    actual_fingerprint=actual_fingerprint &
    )
  call assert_equal_i32(status, outer_event_queue_io_ok, 'matching global queue inventory was rejected')
  call assert_equal_i64(actual_event_count, int(queue%size(), i64), 'measured global queue count mismatch')
  call assert_close_dp(actual_signed_charge, queue%signed_charge(), 0.0_dp, &
                       'measured global queue signed charge mismatch')
  call assert_true(actual_fingerprint == expected_fingerprint, 'measured global queue fingerprint mismatch')
  call test_end()

  call test_begin('global_inventory_charge_roundoff_is_accepted')
  call validate_outer_event_queue_inventory( &
    restored, int(queue%size(), i64), &
    queue%signed_charge()*(1.0_dp + 64.0_dp*epsilon(1.0_dp)), expected_fingerprint, status, message &
    )
  call assert_equal_i32(status, outer_event_queue_io_ok, &
                        'roundoff-scale global queue signed-charge difference was rejected')
  call test_end()

  call test_begin('stale_checkpoint_count_is_rejected')
  call validate_outer_event_queue_inventory( &
    restored, int(queue%size(), i64) + 1_i64, queue%signed_charge(), expected_fingerprint, status, message &
    )
  call assert_equal_i32(status, outer_event_queue_contract_mismatch, &
                        'global queue count mismatch was not rejected')
  call assert_true(len_trim(message) > 0, 'global queue count mismatch should explain the rejection')
  call test_end()

  call test_begin('stale_checkpoint_signed_charge_is_rejected')
  call validate_outer_event_queue_inventory( &
    restored, int(queue%size(), i64), queue%signed_charge() + 1.0e-9_dp, expected_fingerprint, status, message &
    )
  call assert_equal_i32(status, outer_event_queue_contract_mismatch, &
                        'global queue signed-charge mismatch was not rejected')
  call assert_true(len_trim(message) > 0, 'global queue signed-charge mismatch should explain the rejection')
  call test_end()

  call test_begin('stale_checkpoint_fingerprint_is_rejected_with_matching_inventory')
  call restored%snapshot(actual, actual_next_id)
  actual(1)%due_time = actual(1)%due_time + 0.125_dp
  call stale_queue%restore(actual, actual_next_id)
  stale_fingerprint = queue_global_fingerprint(stale_queue, 0_i32, 1_i32)
  call validate_outer_event_queue_inventory( &
    restored, int(restored%size(), i64), restored%signed_charge(), stale_fingerprint, status, message &
    )
  call assert_equal_i32(status, outer_event_queue_contract_mismatch, &
                        'global queue fingerprint mismatch was not rejected')
  call assert_true(len_trim(message) > 0, 'global queue fingerprint mismatch should explain the rejection')
  call test_end()

  call test_begin('empty_queue_roundtrip_preserves_next_id')
  call queue%init(first_event_id=987654321_i64)
  call write_outer_event_queue_checkpoint(out_dir, queue, 0_i32, status, message)
  call assert_equal_i32(status, outer_event_queue_io_ok, 'empty queue checkpoint write failed')
  call load_outer_event_queue_checkpoint(out_dir, 0_i32, restored, found, status, message)
  call assert_true(found, 'empty queue checkpoint was not found')
  call assert_equal_i32(status, outer_event_queue_io_ok, 'empty queue checkpoint did not load')
  call assert_equal_i32(restored%size(), 0_i32, 'empty queue checkpoint restored an event')
  call assert_equal_i64(restored%next_id(), 987654321_i64, 'empty queue next event ID mismatch')
  call test_end()

  call test_begin('completed_batch_mismatch_fails_closed')
  call load_outer_event_queue_checkpoint( &
    out_dir, 8_i32, restored, found, status, message, mpi_rank=1_i32, mpi_size=4_i32 &
    )
  call assert_true(found, 'batch-mismatch fixture was not found')
  call assert_equal_i32(status, outer_event_queue_contract_mismatch, 'batch mismatch was not rejected')
  call assert_equal_i32(restored%size(), 0_i32, 'batch mismatch exposed checkpoint records')
  call test_end()

  call test_begin('world_size_mismatch_fails_closed')
  call load_outer_event_queue_checkpoint( &
    out_dir, 7_i32, restored, found, status, message, mpi_rank=1_i32, mpi_size=3_i32 &
    )
  call assert_true(found, 'world-size-mismatch fixture was not found')
  call assert_equal_i32(status, outer_event_queue_contract_mismatch, 'world size mismatch was not rejected')
  call assert_equal_i32(restored%size(), 0_i32, 'world size mismatch exposed checkpoint records')
  call test_end()

  call test_begin('rank_mismatch_fails_closed')
  call write_fixture(rank2_path, outer_event_queue_checkpoint_schema, 4_i32, 1_i32, 7_i32, 44_i64, 0_i32)
  call load_outer_event_queue_checkpoint( &
    out_dir, 7_i32, restored, found, status, message, mpi_rank=2_i32, mpi_size=4_i32 &
    )
  call assert_true(found, 'rank-mismatch fixture was not found')
  call assert_equal_i32(status, outer_event_queue_contract_mismatch, 'rank mismatch was not rejected')
  call assert_equal_i32(restored%size(), 0_i32, 'rank mismatch exposed checkpoint records')
  call test_end()

  call test_begin('unsupported_schema_fails_closed')
  call write_fixture(serial_path, 99_i32, 1_i32, 0_i32, 0_i32, 1_i64, 0_i32)
  call load_outer_event_queue_checkpoint(out_dir, 0_i32, restored, found, status, message)
  call assert_true(found, 'unsupported-schema fixture was not found')
  call assert_equal_i32(status, outer_event_queue_unsupported_schema, 'unsupported schema was not rejected')
  call assert_equal_i32(restored%size(), 0_i32, 'unsupported schema exposed checkpoint records')
  call test_end()

  call test_begin('local_fingerprint_mismatch_fails_closed')
  call queue%init()
  event = make_event(2.0_dp, outer_event_outcome_return, 1_i32, 0_i32, 1_i32, 1001_i64, 6_i32, 1.0_dp)
  call queue%enqueue(event)
  call queue%snapshot(expected, expected_next_id)
  call format_event_record(expected(1), fixture_row)
  stale_fingerprint = queue%fingerprint(0_i32)
  if (stale_fingerprint(1:1) == '0') then
    stale_fingerprint(1:1) = '1'
  else
    stale_fingerprint(1:1) = '0'
  end if
  call write_fixture( &
    serial_path, outer_event_queue_checkpoint_schema, 1_i32, 0_i32, 0_i32, expected_next_id, 1_i32, &
    trim(fixture_row), stale_fingerprint &
    )
  call load_outer_event_queue_checkpoint(out_dir, 0_i32, restored, found, status, message)
  call assert_true(found, 'local-fingerprint-mismatch fixture was not found')
  call assert_equal_i32(status, outer_event_queue_contract_mismatch, &
                        'local fingerprint mismatch was not rejected')
  call assert_equal_i32(restored%size(), 0_i32, 'local fingerprint mismatch exposed checkpoint records')
  call test_end()

  call test_begin('malformed_record_fails_closed')
  call write_fixture( &
    serial_path, outer_event_queue_checkpoint_schema, 1_i32, 0_i32, 0_i32, 2_i64, 1_i32, '1,0.0' &
    )
  call load_outer_event_queue_checkpoint(out_dir, 0_i32, restored, found, status, message)
  call assert_true(found, 'malformed-record fixture was not found')
  call assert_equal_i32(status, outer_event_queue_malformed, 'malformed event record was not rejected')
  call assert_equal_i32(restored%size(), 0_i32, 'malformed record exposed checkpoint records')
  call test_end()

  call cleanup_files()
  call test_summary()

contains

  function make_event( &
    due_time, outcome, species_id, origin_rank, origin_batch, particle_id, face_index, scale &
    ) result(record)
    real(dp), intent(in) :: due_time, scale
    integer(i32), intent(in) :: outcome, species_id, origin_rank, origin_batch, face_index
    integer(i64), intent(in) :: particle_id
    type(outer_event_record_type) :: record

    record%queued_time = 0.012345678901234567_dp*scale
    record%due_time = due_time
    record%outcome = outcome
    record%species_id = species_id
    record%origin_rank = origin_rank
    record%origin_batch = origin_batch
    record%origin_particle_id = particle_id
    record%interface_face_index = face_index
    record%q = -1.602176634e-19_dp*scale
    record%m = 9.1093837139e-31_dp*scale
    record%w = 1.2345678901234567e7_dp*scale
    record%x = [0.125_dp, -2.5e-9_dp, 3.75e2_dp]*scale
    record%v = [-4.5e5_dp, 6.25_dp, -7.875e-4_dp]*scale
  end function make_event

  subroutine assert_record_equal(actual_record, expected_record)
    type(outer_event_record_type), intent(in) :: actual_record, expected_record

    call assert_equal_i64(actual_record%event_id, expected_record%event_id, 'event ID mismatch')
    call assert_close_dp(actual_record%queued_time, expected_record%queued_time, 0.0_dp, 'queued time mismatch')
    call assert_close_dp(actual_record%due_time, expected_record%due_time, 0.0_dp, 'due time mismatch')
    call assert_equal_i32(actual_record%outcome, expected_record%outcome, 'outcome mismatch')
    call assert_equal_i32(actual_record%species_id, expected_record%species_id, 'species ID mismatch')
    call assert_equal_i32(actual_record%origin_rank, expected_record%origin_rank, 'origin rank mismatch')
    call assert_equal_i32(actual_record%origin_batch, expected_record%origin_batch, 'origin batch mismatch')
    call assert_equal_i64( &
      actual_record%origin_particle_id, expected_record%origin_particle_id, 'origin particle ID mismatch' &
      )
    call assert_equal_i32( &
      actual_record%interface_face_index, expected_record%interface_face_index, 'interface face mismatch' &
      )
    call assert_close_dp(actual_record%q, expected_record%q, 0.0_dp, 'particle charge mismatch')
    call assert_close_dp(actual_record%m, expected_record%m, 0.0_dp, 'particle mass mismatch')
    call assert_close_dp(actual_record%w, expected_record%w, 0.0_dp, 'macro weight mismatch')
    call assert_allclose_1d(actual_record%x, expected_record%x, 0.0_dp, 'position mismatch')
    call assert_allclose_1d(actual_record%v, expected_record%v, 0.0_dp, 'velocity mismatch')
  end subroutine assert_record_equal

  subroutine write_fixture(path, schema, world_size, rank, batches, next_id, event_count, row, local_fingerprint)
    character(len=*), intent(in) :: path
    integer(i32), intent(in) :: schema, world_size, rank, batches, event_count
    integer(i64), intent(in) :: next_id
    character(len=*), intent(in), optional :: row
    character(len=*), intent(in), optional :: local_fingerprint
    integer :: unit, ios
    character(len=16) :: fixture_fingerprint

    fixture_fingerprint = '0000000000000000'
    if (present(local_fingerprint)) fixture_fingerprint = local_fingerprint
    open (newunit=unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create queue checkpoint fixture'
    write (unit, '(a,i0)') 'outer_event_queue_schema=', schema
    write (unit, '(a,i0)') 'mpi_world_size=', world_size
    write (unit, '(a,i0)') 'mpi_rank=', rank
    write (unit, '(a,i0)') 'completed_batches=', batches
    write (unit, '(a,i0)') 'next_event_id=', next_id
    write (unit, '(a,i0)') 'event_count=', event_count
    write (unit, '(a,a)') 'local_fingerprint=', fixture_fingerprint
    write (unit, '(a)') table_header
    if (present(row)) write (unit, '(a)') trim(row)
    close (unit)
  end subroutine write_fixture

  subroutine format_event_record(record, row)
    type(outer_event_record_type), intent(in) :: record
    character(len=*), intent(out) :: row

    write (row, '(i0,2(",",es26.17e3),6(",",i0),9(",",es26.17e3))') &
      record%event_id, record%queued_time, record%due_time, record%outcome, record%species_id, &
      record%origin_rank, record%origin_batch, record%origin_particle_id, record%interface_face_index, &
      record%q, record%m, record%w, record%x, record%v
  end subroutine format_event_record

  function queue_global_fingerprint(active_queue, rank, world_size) result(fingerprint)
    type(outer_event_queue_type), intent(in) :: active_queue
    integer(i32), intent(in) :: rank, world_size
    character(len=16) :: fingerprint
    integer(i64) :: components(2)

    call active_queue%fingerprint_components(rank, components)
    fingerprint = outer_event_queue_global_fingerprint(components, world_size)
  end function queue_global_fingerprint

  subroutine cleanup_files()
    call delete_file_if_exists(trim(serial_path))
    call delete_file_if_exists(trim(serial_path)//'.tmp')
    call delete_file_if_exists(trim(rank1_path))
    call delete_file_if_exists(trim(rank1_path)//'.tmp')
    call delete_file_if_exists(trim(rank2_path))
    call delete_file_if_exists(trim(rank2_path)//'.tmp')
    call remove_empty_directory(out_dir)
  end subroutine cleanup_files

end program test_outer_event_queue_io
