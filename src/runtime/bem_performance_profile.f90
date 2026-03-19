!> 実行フェーズごとの壁時計計測と MPI 集約出力を担う軽量プロファイラ。
module bem_performance_profile
!$ use omp_lib
  use, intrinsic :: iso_fortran_env, only: int64, output_unit
  use bem_kinds, only: dp, i32
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_is_root, &
                     mpi_allreduce_sum_real_dp_array, mpi_allreduce_sum_i32_array, &
                     mpi_allreduce_min_real_dp_array, mpi_allreduce_max_real_dp_array
  implicit none
  private

  integer(i32), parameter, public :: perf_region_program_total = 1_i32
  integer(i32), parameter, public :: perf_region_load_or_init = 2_i32
  integer(i32), parameter, public :: perf_region_history_open = 3_i32
  integer(i32), parameter, public :: perf_region_simulation_total = 4_i32
  integer(i32), parameter, public :: perf_region_field_solver_init = 5_i32
  integer(i32), parameter, public :: perf_region_batch_total = 6_i32
  integer(i32), parameter, public :: perf_region_prepare_batch = 7_i32
  integer(i32), parameter, public :: perf_region_field_refresh = 8_i32
  integer(i32), parameter, public :: perf_region_particle_batch = 9_i32
  integer(i32), parameter, public :: perf_region_commit_charge = 10_i32
  integer(i32), parameter, public :: perf_region_count_outcomes = 11_i32
  integer(i32), parameter, public :: perf_region_mpi_reduce = 12_i32
  integer(i32), parameter, public :: perf_region_stats_update = 13_i32
  integer(i32), parameter, public :: perf_region_history_write = 14_i32
  integer(i32), parameter, public :: perf_region_write_results = 15_i32
  integer(i32), parameter, public :: perf_region_write_checkpoint = 16_i32
  integer(i32), parameter, public :: perf_region_particle_field_eval = 17_i32
  integer(i32), parameter, public :: perf_region_particle_push = 18_i32
  integer(i32), parameter, public :: perf_region_particle_collision = 19_i32
  integer(i32), parameter, public :: perf_region_count = 19_i32

  type :: perf_region_stats
    real(dp) :: total_s = 0.0d0
    integer(i32) :: call_count = 0_i32
  end type perf_region_stats

  type :: perf_state_type
    logical :: enabled = .false.
    logical :: detail_enabled = .false.
    logical :: write_files = .false.
    integer(i32) :: omp_max_threads = 1_i32
    character(len=1024) :: output_dir = ''
    type(perf_region_stats) :: regions(perf_region_count)
  end type perf_state_type

  type(perf_state_type), save :: perf_state = perf_state_type()

  character(len=32), parameter :: perf_region_names(perf_region_count) = [character(len=32) :: &
                                                                          'program_total', &
                                                                          'load_or_init', &
                                                                          'history_open', &
                                                                          'simulation_total', &
                                                                          'field_solver_init', &
                                                                          'batch_total', &
                                                                          'prepare_batch', &
                                                                          'field_refresh', &
                                                                          'particle_batch', &
                                                                          'commit_charge', &
                                                                          'count_outcomes', &
                                                                          'mpi_reduce', &
                                                                          'stats_update', &
                                                                          'history_write', &
                                                                          'write_results', &
                                                                          'write_checkpoint', &
                                                                          'particle_field_eval', &
                                                                          'particle_push', &
                                                                          'particle_collision' &
                                                                          ]

  public :: perf_reset
  public :: perf_configure
  public :: perf_configure_from_env
  public :: perf_set_output_context
  public :: perf_is_enabled
  public :: perf_is_detail_enabled
  public :: perf_wall_time_seconds
  public :: perf_region_begin
  public :: perf_region_end
  public :: perf_add_elapsed
  public :: perf_write_outputs

contains

  !> プロファイラ状態を既定値へ戻す。
  subroutine perf_reset()
    perf_state = perf_state_type()
  end subroutine perf_reset

  !> 明示フラグで計測状態を設定する。
  subroutine perf_configure(enabled, detail_enabled)
    logical, intent(in) :: enabled
    logical, intent(in) :: detail_enabled

    call perf_reset()
    perf_state%enabled = enabled .or. detail_enabled
    perf_state%detail_enabled = detail_enabled
    perf_state%omp_max_threads = 1_i32
!$  perf_state%omp_max_threads = int(max(1, omp_get_max_threads()), i32)
  end subroutine perf_configure

  !> 環境変数 `BEACH_PROFILE` / `BEACH_PROFILE_DETAIL` から計測状態を初期化する。
  subroutine perf_configure_from_env()
    logical :: enabled, detail_enabled

    enabled = env_flag_enabled('BEACH_PROFILE')
    detail_enabled = env_flag_enabled('BEACH_PROFILE_DETAIL')
    call perf_configure(enabled, detail_enabled)
  end subroutine perf_configure_from_env

  !> 出力先ディレクトリとファイル書き出し可否を登録する。
  subroutine perf_set_output_context(output_dir, write_files)
    character(len=*), intent(in) :: output_dir
    logical, intent(in) :: write_files

    perf_state%output_dir = ''
    perf_state%output_dir = trim(output_dir)
    perf_state%write_files = write_files
  end subroutine perf_set_output_context

  !> 粗粒度プロファイルが有効かを返す。
  logical function perf_is_enabled()
    perf_is_enabled = perf_state%enabled
  end function perf_is_enabled

  !> 詳細計測が有効かを返す。
  logical function perf_is_detail_enabled()
    perf_is_detail_enabled = perf_state%detail_enabled
  end function perf_is_detail_enabled

  !> OpenMP有効時は `omp_get_wtime`、それ以外は `system_clock` を使う壁時計。
  real(dp) function perf_wall_time_seconds()
    integer(int64) :: count, rate

!$  perf_wall_time_seconds = omp_get_wtime()
!$  return

    call system_clock(count=count, count_rate=rate)
    if (rate > 0_int64) then
      perf_wall_time_seconds = real(count, dp)/real(rate, dp)
    else
      call cpu_time(perf_wall_time_seconds)
    end if
  end function perf_wall_time_seconds

  !> フェーズ開始時刻を取得する。
  subroutine perf_region_begin(region_id, t0)
    integer(i32), intent(in) :: region_id
    real(dp), intent(out) :: t0

    if (.not. perf_state%enabled) then
      t0 = 0.0d0
      return
    end if
    if (.not. valid_region(region_id)) then
      t0 = 0.0d0
      return
    end if
    t0 = perf_wall_time_seconds()
  end subroutine perf_region_begin

  !> フェーズ終了時刻との差分を累積する。
  subroutine perf_region_end(region_id, t0)
    integer(i32), intent(in) :: region_id
    real(dp), intent(in) :: t0

    if (.not. perf_state%enabled) return
    if (.not. valid_region(region_id)) return
    call perf_add_elapsed(region_id, perf_wall_time_seconds() - t0)
  end subroutine perf_region_end

  !> 経過時間と呼び出し回数を累積する。
  subroutine perf_add_elapsed(region_id, elapsed_s, call_count)
    integer(i32), intent(in) :: region_id
    real(dp), intent(in) :: elapsed_s
    integer(i32), intent(in), optional :: call_count
    integer(i32) :: calls

    if (.not. perf_state%enabled) return
    if (.not. valid_region(region_id)) return

    calls = 1_i32
    if (present(call_count)) calls = max(0_i32, call_count)
    perf_state%regions(region_id)%total_s = perf_state%regions(region_id)%total_s + max(0.0d0, elapsed_s)
    perf_state%regions(region_id)%call_count = perf_state%regions(region_id)%call_count + calls
  end subroutine perf_add_elapsed

  !> 集計済みプロファイルを標準出力および CSV へ書き出す。
  subroutine perf_write_outputs(mpi)
    type(mpi_context), intent(in) :: mpi

    integer(i32) :: rank, world_size, region_idx
    integer(i32) :: call_sum(perf_region_count)
    real(dp) :: total_sum(perf_region_count), total_min(perf_region_count), total_max(perf_region_count)
    real(dp) :: total_mean, imbalance
    integer :: u, ios
    character(len=1024) :: path

    if (.not. perf_state%enabled) return

    call mpi_get_rank_size(rank, world_size, mpi)

    do region_idx = 1, perf_region_count
      total_sum(region_idx) = perf_state%regions(region_idx)%total_s
      total_min(region_idx) = perf_state%regions(region_idx)%total_s
      total_max(region_idx) = perf_state%regions(region_idx)%total_s
      call_sum(region_idx) = perf_state%regions(region_idx)%call_count
    end do

    call mpi_allreduce_sum_real_dp_array(mpi, total_sum)
    call mpi_allreduce_min_real_dp_array(mpi, total_min)
    call mpi_allreduce_max_real_dp_array(mpi, total_max)
    call mpi_allreduce_sum_i32_array(mpi, call_sum)

    if (.not. mpi_is_root(mpi)) return

    write (output_unit, '(a,es12.4,a,es12.4,a,i0,a,i0,a,l1)') &
      'performance: program_total(rank_max)=', total_max(perf_region_program_total), &
      ' simulation_total(rank_max)=', total_max(perf_region_simulation_total), &
      ' mpi=', world_size, ' omp=', perf_state%omp_max_threads, ' detail=', perf_state%detail_enabled

    if (.not. perf_state%write_files .or. len_trim(perf_state%output_dir) == 0) then
      flush (output_unit)
      return
    end if

    path = trim(perf_state%output_dir)//'/performance_profile.csv'
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open performance_profile.csv.'

    write (u, '(a)') '# BEACH performance profile'
    write (u, '(a,i0)') '# mpi_world_size=', world_size
    write (u, '(a,i0)') '# omp_max_threads=', perf_state%omp_max_threads
    write (u, '(a,l1)') '# detail_enabled=', perf_state%detail_enabled
    write (u, '(a)') '# use rank_max_s of simulation_total for scaling comparisons'
    write (u, '(a)') 'region,calls_sum,calls_mean,rank_min_s,rank_mean_s,rank_max_s,imbalance_ratio'

    do region_idx = 1, perf_region_count
      total_mean = total_sum(region_idx)/real(max(1_i32, world_size), dp)
      if (total_mean > tiny(1.0d0)) then
        imbalance = total_max(region_idx)/total_mean
      else
        imbalance = 0.0d0
      end if
      write (u, '(a,a,i0,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
        trim(perf_region_names(region_idx)), ',', call_sum(region_idx), ',', &
        real(call_sum(region_idx), dp)/real(max(1_i32, world_size), dp), ',', &
        total_min(region_idx), ',', total_mean, ',', total_max(region_idx), ',', imbalance
    end do

    close (u)
    write (output_unit, '(a,a)') 'performance profile written to ', trim(path)
    flush (output_unit)
  end subroutine perf_write_outputs

  !> 計測対象 region 番号かどうかを返す。
  logical function valid_region(region_id)
    integer(i32), intent(in) :: region_id

    valid_region = (region_id >= 1_i32 .and. region_id <= perf_region_count)
  end function valid_region

  !> 真偽値環境変数を読み取る。
  logical function env_flag_enabled(name)
    character(len=*), intent(in) :: name
    character(len=32) :: value
    integer :: status, value_len

    env_flag_enabled = .false.
    value = ''
    call get_environment_variable(name, value, length=value_len, status=status)
    if (status /= 0 .or. value_len <= 0) return

    select case (trim(lower_ascii(value(1:value_len))))
    case ('1', 'true', 'yes', 'on')
      env_flag_enabled = .true.
    end select
  end function env_flag_enabled

  !> ASCII英字だけを小文字化する。
  pure function lower_ascii(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i, code

    out = s
    do i = 1, len(s)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
    end do
  end function lower_ascii

end module bem_performance_profile
