program test_outer_event_queue
  use bem_kinds, only: dp, i32, i64
  use bem_outer_event_queue, only: outer_event_queue_type, outer_event_record_type, &
                                   outer_event_outcome_return, outer_event_outcome_escape, &
                                   outer_event_queue_fingerprint_is_valid
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, &
                          assert_close_dp, assert_equal_i32, assert_equal_i64
  implicit none

  type(outer_event_queue_type) :: queue, restored
  type(outer_event_record_type) :: event
  type(outer_event_record_type), allocatable :: events(:), snapshot(:)
  integer(i64) :: event_id, next_event_id
  character(len=16) :: fingerprint, changed_fingerprint

  call test_init(6)

  call test_begin('enqueue_orders_by_due_time_then_stable_id')
  call queue%init(1_i32)
  event = make_event(3.0_dp, outer_event_outcome_return, 1_i32, -1.0_dp, 2.0_dp)
  event%origin_particle_id = 101_i64
  call queue%enqueue(event, event_id)
  call assert_equal_i64(event_id, 1_i64, 'first queue ID mismatch')
  event = make_event(1.0_dp, outer_event_outcome_return, 1_i32, -1.0_dp, 3.0_dp)
  event%origin_particle_id = 102_i64
  call queue%enqueue(event)
  event = make_event(1.0_dp, outer_event_outcome_escape, 2_i32, 1.0_dp, 4.0_dp)
  event%origin_particle_id = 103_i64
  call queue%enqueue(event)
  event = make_event(2.0_dp, outer_event_outcome_return, 1_i32, -1.0_dp, 5.0_dp)
  event%origin_particle_id = 104_i64
  call queue%enqueue(event)
  call queue%snapshot(snapshot, next_event_id)
  call assert_equal_i32(queue%size(), 4_i32, 'queue size mismatch after growth')
  call assert_equal_i64(next_event_id, 5_i64, 'next queue ID mismatch')
  call assert_equal_i64(snapshot(1)%event_id, 2_i64, 'earliest event ID mismatch')
  call assert_equal_i64(snapshot(2)%event_id, 3_i64, 'equal-time enqueue order was not stable')
  call assert_equal_i64(snapshot(3)%event_id, 4_i64, 'middle event ID mismatch')
  call assert_equal_i64(snapshot(4)%event_id, 1_i64, 'latest event ID mismatch')
  call assert_equal_i64(snapshot(2)%origin_particle_id, 103_i64, 'event metadata changed during ordering')
  call test_end()

  call test_begin('pop_due_includes_time_boundary_and_preserves_order')
  call queue%pop_due(2.0_dp, events)
  call assert_equal_i32(int(size(events), i32), 3_i32, 'due event count mismatch')
  call assert_equal_i64(events(1)%event_id, 2_i64, 'first popped event mismatch')
  call assert_equal_i64(events(2)%event_id, 3_i64, 'second popped event mismatch')
  call assert_equal_i64(events(3)%event_id, 4_i64, 'boundary-time event was not popped')
  call assert_equal_i32(queue%size(), 1_i32, 'remaining queue size mismatch')
  call queue%snapshot(snapshot)
  call assert_equal_i64(snapshot(1)%event_id, 1_i64, 'wrong event remained after pop')
  call test_end()

  call test_begin('inventory_uses_physical_charge_times_macro_weight')
  call queue%clear(.true.)
  event = make_event(1.0_dp, outer_event_outcome_return, 7_i32, -2.0_dp, 3.0_dp)
  call queue%enqueue(event)
  event = make_event(2.0_dp, outer_event_outcome_escape, 7_i32, -2.0_dp, 4.0_dp)
  call queue%enqueue(event)
  event = make_event(3.0_dp, outer_event_outcome_return, 8_i32, 5.0_dp, 2.0_dp)
  call queue%enqueue(event)
  call assert_close_dp(queue%signed_charge(), -4.0_dp, 0.0_dp, 'total signed-charge inventory mismatch')
  call assert_close_dp(queue%signed_charge(7_i32), -14.0_dp, 0.0_dp, 'species charge inventory mismatch')
  call assert_close_dp(queue%photoelectron_number(7_i32), 7.0_dp, 0.0_dp, 'photoelectron inventory mismatch')
  call assert_close_dp(queue%species_particle_number(8_i32), 2.0_dp, 0.0_dp, 'species particle inventory mismatch')
  call test_end()

  call test_begin('snapshot_restore_preserves_sequence_and_sorts')
  call queue%snapshot(snapshot, next_event_id)
  snapshot = snapshot([3, 1, 2])
  call restored%restore(snapshot, next_event_id, 8_i32)
  call restored%snapshot(events)
  call assert_equal_i64(events(1)%event_id, 1_i64, 'restore did not sort first due event')
  call assert_equal_i64(events(2)%event_id, 2_i64, 'restore did not sort second due event')
  call assert_equal_i64(events(3)%event_id, 3_i64, 'restore did not sort last due event')
  event = make_event(4.0_dp, outer_event_outcome_return, 7_i32, -2.0_dp, 1.0_dp)
  call restored%enqueue(event, event_id)
  call assert_equal_i64(event_id, next_event_id, 'restore did not preserve the next event ID')
  call test_end()

  call test_begin('clear_retains_or_resets_id_sequence_explicitly')
  call restored%clear()
  call assert_equal_i32(restored%size(), 0_i32, 'clear did not empty the queue')
  call assert_equal_i64(restored%next_id(), next_event_id + 1_i64, 'clear unexpectedly reset the ID sequence')
  call restored%clear(.true.)
  call assert_equal_i64(restored%next_id(), 1_i64, 'explicit sequence reset failed')
  call restored%pop_due(10.0_dp, events)
  call assert_true(allocated(events), 'empty pop must return an allocated zero-length array')
  call assert_equal_i32(int(size(events), i32), 0_i32, 'empty pop returned an event')
  call test_end()

  call test_begin('fingerprint_covers_rank_sequence_and_active_payload')
  call queue%init(first_event_id=17_i64)
  event = make_event(3.0_dp, outer_event_outcome_return, 2_i32, -2.0_dp, 4.0_dp)
  call queue%enqueue(event)
  fingerprint = queue%fingerprint(0_i32)
  call assert_true(outer_event_queue_fingerprint_is_valid(fingerprint), 'queue fingerprint format is invalid')
  call assert_true(queue%fingerprint(0_i32) == fingerprint, 'queue fingerprint must be deterministic')
  call assert_true(queue%fingerprint(1_i32) /= fingerprint, 'queue fingerprint must include rank ownership')
  call queue%snapshot(snapshot, next_event_id)
  call restored%restore(snapshot, next_event_id + 1_i64)
  call assert_true(restored%fingerprint(0_i32) /= fingerprint, 'queue fingerprint must include next event ID')
  snapshot(1)%v(3) = snapshot(1)%v(3) + 0.5_dp
  call restored%restore(snapshot, next_event_id)
  changed_fingerprint = restored%fingerprint(0_i32)
  call assert_true(changed_fingerprint /= fingerprint, 'queue fingerprint must include active event payload')
  call test_end()

  call test_summary()

contains

  function make_event(due_time, outcome, species_id, charge, weight) result(record)
    real(dp), intent(in) :: due_time, charge, weight
    integer(i32), intent(in) :: outcome, species_id
    type(outer_event_record_type) :: record

    record%queued_time = 0.0_dp
    record%due_time = due_time
    record%outcome = outcome
    record%species_id = species_id
    record%origin_batch = 1_i32
    record%origin_particle_id = 1_i64
    record%interface_face_index = 6_i32
    record%q = charge
    record%m = 10.0_dp
    record%w = weight
    record%x = [1.0_dp, 2.0_dp, 3.0_dp]
    record%v = [-1.0_dp, -2.0_dp, -3.0_dp]
  end function make_event

end program test_outer_event_queue
