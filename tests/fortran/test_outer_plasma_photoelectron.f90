program test_outer_plasma_photoelectron
  use bem_kinds, only: dp, i32, i64
  use bem_outer_plasma_photoelectron, only: photoelectron_histogram_type, photoelectron_coupling_state_type, &
                                            validate_photoelectron_linear_applicability, photoelectron_closure_ok, &
                                            photoelectron_closure_not_applicable
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32, &
                          assert_equal_i64
  implicit none

  type(photoelectron_histogram_type) :: histogram, batch
  type(photoelectron_coupling_state_type) :: state
  integer(i32) :: status

  call test_init(3)
  call test_begin('histogram_charge_energy_moments')
  call histogram%init(4_i32, 4.0_dp)
  call histogram%add(-1.0_dp, 2.0_dp, 3.0_dp, [1.0_dp, -2.0_dp, 1.0_dp])
  call assert_close_dp(histogram%total_signed_charge(), -3.0_dp, 0.0_dp, 'histogram charge mismatch')
  call assert_close_dp(histogram%total_kinetic_energy(), 18.0_dp, 0.0_dp, 'histogram energy mismatch')
  call assert_equal_i64(histogram%total_count(), 1_i64, 'histogram count mismatch')
  call assert_close_dp(sum(histogram%tangential_momentum_x), 6.0_dp, 0.0_dp, 'x momentum mismatch')
  call assert_close_dp(sum(histogram%tangential_momentum_y), -12.0_dp, 0.0_dp, 'y momentum mismatch')
  call test_end()

  call test_begin('explicit_previous_batch_ownership')
  call state%init(4_i32, 4.0_dp)
  call state%begin_batch(batch)
  call batch%add(-1.0_dp, 1.0_dp, 2.0_dp, [0.0_dp, 0.0_dp, 1.0_dp])
  call state%commit_batch(1_i32, batch)
  call assert_equal_i32(state%last_completed_batch, 1_i32, 'completed batch mismatch')
  call assert_close_dp(state%previous_batch%total_signed_charge(), -2.0_dp, 0.0_dp, 'previous batch mismatch')
  call assert_close_dp(state%cumulative%total_signed_charge(), -2.0_dp, 0.0_dp, 'cumulative mismatch')
  call test_end()

  call test_begin('strong_photoelectron_guard')
  call validate_photoelectron_linear_applicability(0.1_dp, 1.0_dp, 0.2_dp, status)
  call assert_equal_i32(status, photoelectron_closure_ok, 'weak photoelectron closure should be applicable')
  call validate_photoelectron_linear_applicability(0.3_dp, 1.0_dp, 0.2_dp, status)
  call assert_equal_i32(status, photoelectron_closure_not_applicable, 'strong photoelectron closure must be rejected')
  call test_end()
  call test_summary()
end program test_outer_plasma_photoelectron
