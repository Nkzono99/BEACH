program test_interface_particle_buffer
  use bem_kinds, only: dp, i32
  use bem_interface_particle_buffer, only: interface_particle_buffer_type, interface_particle_record_type
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32
  implicit none

  type(interface_particle_buffer_type) :: buffer
  type(interface_particle_record_type) :: record

  call test_init(2)
  call test_begin('append_grows_and_preserves_records')
  call buffer%init(1_i32)
  record%charge = -1.0_dp
  record%species_id = 2_i32
  call buffer%append(record)
  record%charge = -2.0_dp
  call buffer%append(record)
  call assert_equal_i32(buffer%count, 2_i32, 'buffer count mismatch')
  call assert_close_dp(buffer%records(1)%charge, -1.0_dp, 0.0_dp, 'first record changed during growth')
  call assert_close_dp(buffer%records(2)%charge, -2.0_dp, 0.0_dp, 'second record mismatch')
  call test_end()

  call test_begin('clear_keeps_capacity')
  call buffer%clear()
  call assert_equal_i32(buffer%count, 0_i32, 'clear count mismatch')
  call assert_equal_i32(int(size(buffer%records), i32), 2_i32, 'clear must retain capacity')
  call test_end()
  call test_summary()
end program test_interface_particle_buffer
