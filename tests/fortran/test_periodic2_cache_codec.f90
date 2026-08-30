program test_periodic2_cache_codec
  use, intrinsic :: iso_fortran_env, only: int8
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_periodic_cache, only: write_periodic_operator_cache, read_periodic_operator_cache, &
                                            periodic_cache_ok, periodic_cache_miss, periodic_cache_invalid
  use test_support, only: test_begin, test_end, test_init, test_summary, assert_true, assert_equal_i32, &
                          assert_allclose_1d, delete_file_if_exists
  implicit none

  character(len=*), parameter :: cache_path = 'test_periodic2_cache_codec.bin'
  integer(i32) :: status
  integer(i32) :: target_nodes(2)
  integer(i32), allocatable :: loaded_nodes(:)
  real(dp) :: operator(3, 3, 2)
  real(dp), allocatable :: loaded_operator(:, :, :)
  integer(int8) :: payload_byte
  integer :: unit, ios, file_size

  call delete_file_if_exists(cache_path)
  target_nodes = [2_i32, 5_i32]
  operator = reshape([(real(status, dp), status=1, size(operator))], shape(operator))
  call test_init(3)

  call test_begin('missing_cache')
  call read_periodic_operator_cache(cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status)
  call assert_equal_i32(status, periodic_cache_miss, 'missing cache status')
  call test_end()

  call test_begin('roundtrip_and_identity_validation')
  call write_periodic_operator_cache(cache_path, 'fingerprint-a', target_nodes, operator, status)
  call assert_equal_i32(status, periodic_cache_ok, 'cache write status')
  call read_periodic_operator_cache( &
    cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status, expected_ncoef=3_i32, expected_ntarget=2_i32 &
    )
  call assert_equal_i32(status, periodic_cache_ok, 'cache read status')
  call assert_true(all(loaded_nodes == target_nodes), 'cached target nodes mismatch')
  call assert_allclose_1d(reshape(loaded_operator, [size(loaded_operator)]), reshape(operator, [size(operator)]), &
                          0.0_dp, 'cached operator mismatch')
  deallocate (loaded_nodes, loaded_operator)
  call read_periodic_operator_cache( &
    cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status, expected_ncoef=4_i32, expected_ntarget=2_i32 &
    )
  call assert_equal_i32(status, periodic_cache_invalid, 'unexpected operator shape must invalidate cache')
  call read_periodic_operator_cache( &
    cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status, expected_ncoef=3_i32, expected_ntarget=1_i32 &
    )
  call assert_equal_i32(status, periodic_cache_invalid, 'unexpected target count must invalidate cache')
  call read_periodic_operator_cache(cache_path, 'fingerprint-b', loaded_nodes, loaded_operator, status)
  call assert_equal_i32(status, periodic_cache_invalid, 'fingerprint mismatch must invalidate cache')
  call test_end()

  call test_begin('payload_corruption_is_rejected')
  inquire (file=cache_path, size=file_size, iostat=ios)
  call assert_equal_i32(int(ios, i32), 0_i32, 'inspect cache size for corruption')
  open (newunit=unit, file=cache_path, access='stream', form='unformatted', status='old', action='readwrite', iostat=ios)
  call assert_equal_i32(int(ios, i32), 0_i32, 'open cache for corruption')
  read (unit, pos=file_size, iostat=ios) payload_byte
  call assert_equal_i32(int(ios, i32), 0_i32, 'read cached payload byte')
  write (unit, pos=file_size, iostat=ios) ieor(payload_byte, 1_int8)
  close (unit)
  call assert_equal_i32(int(ios, i32), 0_i32, 'corrupt cached payload byte')
  call read_periodic_operator_cache(cache_path, 'fingerprint-a', loaded_nodes, loaded_operator, status)
  call assert_equal_i32(status, periodic_cache_invalid, 'checksum mismatch must invalidate cache')
  call test_end()

  call delete_file_if_exists(cache_path)
  call test_summary()
end program test_periodic2_cache_codec
