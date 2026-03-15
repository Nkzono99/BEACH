!> MPIの初期化・集約を抽象化し、非MPIビルドでは単一ランク動作へフォールバックする。
module bem_mpi
  use bem_kinds, only: dp, i32
  implicit none
  private

  type :: mpi_context
    integer(i32) :: rank = 0_i32
    integer(i32) :: size = 1_i32
    logical :: enabled = .false.
#ifdef USE_MPI
    logical :: initialized_here = .false.
#endif
  end type mpi_context

  public :: mpi_context
  public :: mpi_initialize
  public :: mpi_shutdown
  public :: mpi_is_root
  public :: mpi_world_size
  public :: mpi_get_rank_size
  public :: mpi_split_count
  public :: mpi_allreduce_sum_real_dp_array
  public :: mpi_allreduce_sum_real_dp_scalar
  public :: mpi_allreduce_min_real_dp_array
  public :: mpi_allreduce_max_real_dp_array
  public :: mpi_allreduce_sum_i32_array
  public :: mpi_allreduce_sum_i32_scalar
  public :: mpi_world_barrier

contains

  !> MPIを初期化して rank / size を取得する。非MPIビルドでは単一ランクを返す。
  subroutine mpi_initialize(ctx)
    type(mpi_context), intent(out) :: ctx
#ifdef USE_MPI
    include 'mpif.h'
    logical :: is_initialized
    integer :: ierr
    integer :: rank_int, size_int

    ctx = mpi_context()

    call MPI_Initialized(is_initialized, ierr)
    if (.not. is_initialized) then
      call MPI_Init(ierr)
      ctx%initialized_here = .true.
    end if

    call MPI_Comm_rank(MPI_COMM_WORLD, rank_int, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size_int, ierr)
    ctx%rank = int(rank_int, i32)
    ctx%size = int(size_int, i32)
    ctx%enabled = (ctx%size > 1_i32)
#else
    ctx = mpi_context()
#endif
  end subroutine mpi_initialize

  !> `mpi_initialize` が実際に初期化した場合のみ `MPI_Finalize` を呼ぶ。
  subroutine mpi_shutdown(ctx)
    type(mpi_context), intent(inout) :: ctx
#ifdef USE_MPI
    include 'mpif.h'
    logical :: is_finalized
    integer :: ierr

    if (ctx%initialized_here) then
      call MPI_Finalized(is_finalized, ierr)
      if (.not. is_finalized) call MPI_Finalize(ierr)
      ctx%initialized_here = .false.
    end if
#endif
    ctx%rank = 0_i32
    ctx%size = 1_i32
    ctx%enabled = .false.
  end subroutine mpi_shutdown

  !> root rank (rank=0) かどうかを返す。
  logical function mpi_is_root(ctx)
    type(mpi_context), intent(in) :: ctx

    mpi_is_root = (ctx%rank == 0_i32)
  end function mpi_is_root

  !> MPI world size を返す（size<=0 は 1 へ補正）。
  integer(i32) function mpi_world_size(ctx)
    type(mpi_context), intent(in), optional :: ctx

    mpi_world_size = 1_i32
    if (present(ctx)) mpi_world_size = max(1_i32, ctx%size)
  end function mpi_world_size

  !> `mpi_context` から rank/size を取得する。未指定時は単一rank(0/1)。
  subroutine mpi_get_rank_size(rank, size, ctx)
    integer(i32), intent(out) :: rank, size
    type(mpi_context), intent(in), optional :: ctx

    rank = 0_i32
    size = 1_i32
    if (present(ctx)) then
      rank = ctx%rank
      size = max(1_i32, ctx%size)
    end if
    if (rank < 0_i32 .or. rank >= size) then
      error stop 'mpi_get_rank_size detected an invalid rank/size pair.'
    end if
  end subroutine mpi_get_rank_size

  !> 総数 `total_count` をrankへ均等分割したときの局所個数を返す。
  integer(i32) function mpi_split_count(total_count, rank, size) result(local_count)
    integer(i32), intent(in) :: total_count, rank, size
    integer(i32) :: base_count, n_remainder

    if (total_count < 0_i32) error stop 'mpi_split_count requires total_count >= 0.'
    if (size <= 0_i32) error stop 'mpi_split_count requires size > 0.'
    if (rank < 0_i32 .or. rank >= size) error stop 'mpi_split_count rank out of range.'

    base_count = total_count/size
    n_remainder = modulo(total_count, size)
    local_count = base_count
    if (rank < n_remainder) local_count = local_count + 1_i32
  end function mpi_split_count

  !> 倍精度配列の総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_real_dp_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: values(:)
#ifdef USE_MPI
    include 'mpif.h'
    real(dp), allocatable :: recvbuf(:)
    integer :: ierr

    if (.not. ctx%enabled) return
    allocate(recvbuf(size(values)))
    call MPI_Allreduce(values, recvbuf, size(values), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    values = recvbuf
#endif
  end subroutine mpi_allreduce_sum_real_dp_array

  !> 倍精度スカラの総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_real_dp_scalar(ctx, value)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: value
#ifdef USE_MPI
    include 'mpif.h'
    real(dp) :: recvval
    integer :: ierr

    if (.not. ctx%enabled) return
    call MPI_Allreduce(value, recvval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    value = recvval
#endif
  end subroutine mpi_allreduce_sum_real_dp_scalar

  !> 倍精度配列の最小値Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_min_real_dp_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: values(:)
#ifdef USE_MPI
    include 'mpif.h'
    real(dp), allocatable :: recvbuf(:)
    integer :: ierr

    if (.not. ctx%enabled) return
    allocate(recvbuf(size(values)))
    call MPI_Allreduce(values, recvbuf, size(values), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    values = recvbuf
#endif
  end subroutine mpi_allreduce_min_real_dp_array

  !> 倍精度配列の最大値Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_max_real_dp_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: values(:)
#ifdef USE_MPI
    include 'mpif.h'
    real(dp), allocatable :: recvbuf(:)
    integer :: ierr

    if (.not. ctx%enabled) return
    allocate(recvbuf(size(values)))
    call MPI_Allreduce(values, recvbuf, size(values), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    values = recvbuf
#endif
  end subroutine mpi_allreduce_max_real_dp_array

  !> 32bit整数配列の総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_i32_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: values(:)
#ifdef USE_MPI
    include 'mpif.h'
    integer, allocatable :: sendbuf(:), recvbuf(:)
    integer :: ierr

    if (.not. ctx%enabled) return
    allocate(sendbuf(size(values)), recvbuf(size(values)))
    sendbuf = int(values, kind=kind(0))
    call MPI_Allreduce(sendbuf, recvbuf, size(values), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    values = int(recvbuf, kind=i32)
#endif
  end subroutine mpi_allreduce_sum_i32_array

  !> 32bit整数スカラの総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_i32_scalar(ctx, value)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: value
#ifdef USE_MPI
    include 'mpif.h'
    integer :: sendval, recvval, ierr

    if (.not. ctx%enabled) return
    sendval = int(value, kind=kind(0))
    call MPI_Allreduce(sendval, recvval, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    value = int(recvval, kind=i32)
#endif
  end subroutine mpi_allreduce_sum_i32_scalar

  !> 全rankの同期ポイント。
  subroutine mpi_world_barrier(ctx)
    type(mpi_context), intent(in) :: ctx
#ifdef USE_MPI
    include 'mpif.h'
    integer :: ierr

    if (.not. ctx%enabled) return
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_world_barrier

end module bem_mpi
