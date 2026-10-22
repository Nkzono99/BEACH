!> 逐次に引く PCG64 の乱数列（Fortran 組み込み `random_number` の置き換え）。
!>
!> `seed_random_stream` で seed を決めると、各 OpenMP thread は (seed, thread 番号) から作った自分の PCG64 stream を
!> 初回の draw で開く。同じ seed・同じ thread で引けば同じ列になる。バッチ粒子の生成と粒子追跡は
!> `bem_counter_rng` の counter-based 乱数を使い、この列は単体の sampling API と test だけが使う。
module bem_random_stream
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int64_t
  use bem_kinds, only: dp, i32, i64
  use bem_counter_rng, only: counter_key
  use m_pcg64_bind, only: pcg64_create, pcg64_destroy, pcg64_random_double53
!$ use omp_lib, only: omp_get_thread_num
  implicit none
  private

  public :: seed_random_stream
  public :: random_uniform

  interface random_uniform
    module procedure random_uniform_scalar
    module procedure random_uniform_1d
    module procedure random_uniform_2d
  end interface random_uniform

  !> 共有の seed と、その世代番号（seed を変えるたびに増え、各 thread は古い stream を開き直す）。
  integer(i64), save :: stream_seed_value = 12345_i64
  integer(i64), save :: stream_generation = 1_i64
  type(c_ptr), save :: thread_handle = c_null_ptr
  integer(i64), save :: thread_generation = 0_i64
!$omp threadprivate(thread_handle, thread_generation)

contains

  !> 乱数列の seed を設定する。複数の値は順に混ぜて 1 つの seed にする。
  subroutine seed_random_stream(seed_values)
    integer(i32), intent(in) :: seed_values(:)
    integer(i64) :: state_seed, stream

    if (size(seed_values) == 0) error stop 'seed_random_stream requires at least one seed value.'
    call counter_key(seed_values(1), int(seed_values, i64), state_seed, stream)
    !$omp critical (beach_random_stream_seed)
    stream_seed_value = ieor(state_seed, stream)
    stream_generation = stream_generation + 1_i64
    !$omp end critical (beach_random_stream_seed)
  end subroutine seed_random_stream

  !> [0, 1) の一様乱数を 1 つ引く。
  subroutine random_uniform_scalar(x)
    real(dp), intent(out) :: x

    call ensure_thread_stream()
    x = pcg64_random_double53(thread_handle)
  end subroutine random_uniform_scalar

  !> 配列を要素順に [0, 1) の一様乱数で埋める。
  subroutine random_uniform_1d(x)
    real(dp), intent(out) :: x(:)
    integer :: i

    call ensure_thread_stream()
    do i = 1, size(x)
      x(i) = pcg64_random_double53(thread_handle)
    end do
  end subroutine random_uniform_1d

  !> 配列を列優先の要素順に [0, 1) の一様乱数で埋める。
  subroutine random_uniform_2d(x)
    real(dp), intent(out) :: x(:, :)
    integer :: i, j

    call ensure_thread_stream()
    do j = 1, size(x, 2)
      do i = 1, size(x, 1)
        x(i, j) = pcg64_random_double53(thread_handle)
      end do
    end do
  end subroutine random_uniform_2d

  subroutine ensure_thread_stream()
    integer(i64) :: state_seed, stream, generation, seed_value
    integer(i32) :: thread_id

    !$omp critical (beach_random_stream_seed)
    generation = stream_generation
    seed_value = stream_seed_value
    !$omp end critical (beach_random_stream_seed)
    if (c_associated(thread_handle) .and. thread_generation == generation) return
    thread_id = 0_i32
!$  thread_id = int(omp_get_thread_num(), i32)
    call counter_key(int(iand(seed_value, 2147483647_i64), i32), [seed_value, int(thread_id, i64), 211_i64], &
                     state_seed, stream)
    if (c_associated(thread_handle)) call pcg64_destroy(thread_handle)
    thread_handle = pcg64_create(int(state_seed, c_int64_t), int(stream, c_int64_t))
    if (.not. c_associated(thread_handle)) error stop 'failed to create a PCG64 generator.'
    thread_generation = generation
  end subroutine ensure_thread_stream

end module bem_random_stream
