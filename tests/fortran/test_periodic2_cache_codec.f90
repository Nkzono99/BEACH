program test_periodic2_cache_codec
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_periodic_cache, only: write_periodic_operator_cache, read_periodic_operator_cache, &
                                            periodic_cache_ok, periodic_cache_miss, periodic_cache_invalid
  use test_support, only: test_begin, test_end, test_summary, assert_true, assert_equal_i32, assert_allclose_1d
  implicit none

  character(len=*), parameter :: cache_path = 'test_periodic2_cache_codec.bin'
  integer(i32) :: status
  integer(i32) :: target_nodes(2)
  integer(i32), allocatable :: loaded_nodes(:)
  real(dp) :: operator(3, 3, 2)
  real(dp), allocatable :: loaded_operator(:, :, :)
  integer :: unit, ios

  call delete_cache()
  target_nodes = [2_i32, 5_i32]
  operator = reshape([(real(status, dp), status=1, size(operator))], shape(operator))

  call test_begin('missing_cache')
  call read_periodic_operator_cache(cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status)
  call assert_equal_i32(status, periodic_cache_miss, 'missing cache status')
  call test_end()

  call test_begin('roundtrip_and_fingerprint')
  call write_periodic_operator_cache(cache_path, 'fingerprint-a', target_nodes, operator, status)
  call assert_equal_i32(status, periodic_cache_ok, 'cache write status')
  call read_periodic_operator_cache(cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status)
  call assert_equal_i32(status, periodic_cache_ok, 'cache read status')
  call assert_true(all(loaded_nodes == target_nodes), 'cached target nodes mismatch')
  call assert_allclose_1d(reshape(loaded_operator, [size(loaded_operator)]), reshape(operator, [size(operator)]), &
                          0.0_dp, 'cached operator mismatch')
  deallocate (loaded_nodes, loaded_operator)
  call read_periodic_operator_cache( &
    cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status, expected_ncoef=4_i32, expected_ntarget=2_i32 &
    )
  call assert_equal_i32(status, periodic_cache_invalid, 'unexpected operator shape must invalidate cache')
  call read_periodic_operator_cache(cache_path, 'fingerprint-b', loaded_nodes, loaded_operator, status)
  call assert_equal_i32(status, periodic_cache_invalid, 'fingerprint mismatch must invalidate cache')
  call test_end()

  call test_begin('corruption_is_rejected')
  open (newunit=unit, file=cache_path, access='stream', form='unformatted', status='old', action='write', iostat=ios)
  call assert_equal_i32(int(ios, i32), 0_i32, 'open cache for corruption')
  write (unit, pos=1, iostat=ios) 'BROKEN-CACHE!!!'
  close (unit)
  call assert_equal_i32(int(ios, i32), 0_i32, 'corrupt cache write')
  call read_periodic_operator_cache(cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status)
  call assert_equal_i32(status, periodic_cache_invalid, 'corrupt cache must be rejected')
  call test_end()

  call delete_cache()
  call test_summary()

contains

  subroutine delete_cache()
    logical :: exists
    integer :: delete_unit
    inquire (file=cache_path, exist=exists)
    if (.not. exists) return
    open (newunit=delete_unit, file=cache_path, status='old')
    close (delete_unit, status='delete')
  end subroutine delete_cache
end program test_periodic2_cache_codec
