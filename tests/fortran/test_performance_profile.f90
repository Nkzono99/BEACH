!> 計測プロファイル CSV の生成と基本ヘッダを検証する。
program test_performance_profile
  use bem_kinds, only: i32
  use bem_mpi, only: mpi_context
  use bem_performance_profile, only: perf_reset, perf_configure, perf_set_output_context, perf_add_elapsed, &
                                     perf_write_outputs, perf_region_program_total, perf_region_simulation_total
  use test_support, only: assert_true, delete_file_if_exists, ensure_directory, remove_empty_directory
  implicit none

  type(mpi_context) :: mpi
  character(len=*), parameter :: out_dir = 'test_performance_profile_tmp'
  character(len=*), parameter :: profile_path = 'test_performance_profile_tmp/performance_profile.csv'
  character(len=512) :: line
  logical :: exists, saw_header, saw_program_total, saw_simulation_total
  integer :: u, ios

  call ensure_directory(out_dir)
  call delete_file_if_exists(profile_path)

  call perf_configure(.true.)
  call perf_set_output_context(out_dir, .true.)
  call perf_add_elapsed(perf_region_program_total, 1.5d0, 1_i32)
  call perf_add_elapsed(perf_region_simulation_total, 1.0d0, 2_i32)
  call perf_write_outputs(mpi)

  inquire (file=profile_path, exist=exists)
  call assert_true(exists, 'performance_profile.csv should be created')

  saw_header = .false.
  saw_program_total = .false.
  saw_simulation_total = .false.

  open (newunit=u, file=profile_path, status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open performance_profile.csv'
  do
    read (u, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (index(line, 'region,calls_sum,calls_mean,rank_min_s,rank_mean_s,rank_max_s,imbalance_ratio') == 1) then
      saw_header = .true.
    end if
    if (index(line, 'program_total,') == 1) saw_program_total = .true.
    if (index(line, 'simulation_total,') == 1) saw_simulation_total = .true.
  end do
  close (u)

  call assert_true(saw_header, 'performance_profile.csv header row is missing')
  call assert_true(saw_program_total, 'program_total row is missing from performance_profile.csv')
  call assert_true(saw_simulation_total, 'simulation_total row is missing from performance_profile.csv')

  call perf_reset()
  call delete_file_if_exists(profile_path)
  call remove_empty_directory(out_dir)
end program test_performance_profile
