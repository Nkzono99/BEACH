program benchmark_periodic2_runtime
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_plan, update_state, &
                                  eval_point, destroy_plan, destroy_state
  implicit none

  integer(i32), parameter :: repeat_count = 200000_i32
  type(fmm_options_type) :: options
  type(fmm_plan_type) :: cached_plan, shell_plan
  type(fmm_state_type) :: cached_state, shell_state
  real(dp) :: src_pos(3, 4), q(4), target(3), field(3), accumulator
  real(dp) :: cached_start, cached_end, shell_start, shell_end
  integer(i32) :: repeat_idx
  integer :: clock_count
  character(len=512) :: cache_dir, cache_path

  src_pos(:, 1) = [0.15_dp, 0.20_dp, -0.10_dp]
  src_pos(:, 2) = [0.75_dp, 0.25_dp, 0.08_dp]
  src_pos(:, 3) = [0.30_dp, 0.80_dp, 0.12_dp]
  src_pos(:, 4) = [0.82_dp, 0.72_dp, -0.06_dp]
  q = [1.0_dp, -0.8_dp, 0.6_dp, -0.6_dp]
  target = [0.55_dp, 0.48_dp, 0.18_dp]

  options%theta = 0.4_dp
  options%leaf_max = 2_i32
  options%order = 4_i32
  options%use_periodic2 = .true.
  options%periodic_axes = [1_i32, 2_i32]
  options%periodic_len = [1.0_dp, 1.0_dp]
  options%periodic_image_layers = 2_i32
  options%periodic_far_correction = 'cached_kneq0'
  options%periodic_ewald_layers = 3_i32
  call system_clock(count=clock_count)
  write (cache_dir, '(a,i0)') 'benchmark_periodic2_runtime_tmp_', clock_count
  options%periodic_cache_dir = trim(cache_dir)
  options%periodic_generation_tolerance = 1.0e-8_dp
  options%target_box_min = [0.0_dp, 0.0_dp, -0.5_dp]
  options%target_box_max = [1.0_dp, 1.0_dp, 0.5_dp]
  call build_plan(cached_plan, src_pos, options)
  cache_path = cached_plan%periodic_cache_path
  call update_state(cached_plan, cached_state, q)

  options%periodic_far_correction = 'none'
  call build_plan(shell_plan, src_pos, options)
  call update_state(shell_plan, shell_state, q)

  accumulator = 0.0_dp
  call cpu_time(cached_start)
  do repeat_idx = 1_i32, repeat_count
    target(1) = modulo(0.55_dp + 1.0e-5_dp*real(repeat_idx, dp), 1.0_dp)
    call eval_point(cached_plan, cached_state, target, field)
    accumulator = accumulator + field(1)
  end do
  call cpu_time(cached_end)

  call cpu_time(shell_start)
  do repeat_idx = 1_i32, repeat_count
    target(1) = modulo(0.55_dp + 1.0e-5_dp*real(repeat_idx, dp), 1.0_dp)
    call eval_point(shell_plan, shell_state, target, field)
    accumulator = accumulator + field(1)
  end do
  call cpu_time(shell_end)

  if (cached_end <= cached_start .or. shell_end <= shell_start) error stop 'benchmark timing is not measurable'
  if (abs(accumulator) >= huge(1.0_dp)) error stop 'benchmark accumulator is not finite'
  write (*, '(a,i0)') 'periodic benchmark repeats=', repeat_count
  write (*, '(a,3(es12.5,1x))') 'periodic eval timing(cached,shell,ratio)=', &
    cached_end - cached_start, shell_end - shell_start, &
    (cached_end - cached_start)/max(shell_end - shell_start, tiny(1.0_dp))

  call destroy_state(cached_state)
  call destroy_state(shell_state)
  call destroy_plan(cached_plan)
  call destroy_plan(shell_plan)
  call delete_file_if_exists(cache_path)
  call delete_file_if_exists(trim(cache_path)//'.lock')
  call remove_empty_directory(cache_dir)

contains

  subroutine delete_file_if_exists(path)
    character(len=*), intent(in) :: path
    logical :: exists
    integer :: unit, ios

    inquire (file=trim(path), exist=exists)
    if (.not. exists) return
    open (newunit=unit, file=trim(path), status='old', action='readwrite', iostat=ios)
    if (ios == 0) close (unit, status='delete')
  end subroutine delete_file_if_exists

  subroutine remove_empty_directory(path)
    character(len=*), intent(in) :: path
    character(len=1024) :: command
    integer :: exit_status

    command = 'rmdir "'//trim(path)//'"'
    call execute_command_line(trim(command), wait=.true., exitstat=exit_status)
    if (exit_status /= 0) error stop 'failed to remove benchmark cache directory'
  end subroutine remove_empty_directory
end program benchmark_periodic2_runtime
