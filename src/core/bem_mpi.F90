!> MPIの初期化・集約を抽象化し、非MPIビルドでは単一ランク動作へフォールバックする。
module bem_mpi
  use bem_kinds, only: dp, i32, i64
#ifdef USE_MPI
  use mpi
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
#endif
  implicit none
  private

  type :: mpi_context
    integer(i32) :: rank = 0_i32
    integer(i32) :: size = 1_i32
    logical :: enabled = .false.
    !> `MPI_Init_thread` が返した thread support level（非 MPI build では 0）。
    integer(i32) :: thread_level = 0_i32
#ifdef USE_MPI
    logical :: initialized_here = .false.
#endif
  end type mpi_context

  !> rank 0 が保持する 64 bit カウンタを全 rank が fetch-and-add で取り合う、バッチ粒子 slot pool の claim 窓口。
  !! `claim` は OpenMP の critical 区間から呼ぶ前提で、MPI は `MPI_THREAD_SERIALIZED` 以上を要求する。
  !! 非 MPI build や 1 rank では rank 内のカウンタに落ちる。
  type, public :: mpi_shared_counter_type
    logical :: active = .false.
    logical :: distributed = .false.
    !> この rank が窓を持つ（rank 0）。remote の fetch-and-add は窓を持つ rank の MPI 進行で処理される。
    logical :: hosts_window = .false.
    !> 次の chunk の fetch-and-add を先に発行してある（結果は `prefetch_result`）。
    logical :: prefetch_pending = .false.
    integer(i64) :: local_value = 0_i64
#ifdef USE_MPI
    integer :: win = MPI_WIN_NULL
#endif
  contains
    procedure :: init => shared_counter_init
    procedure :: reset => shared_counter_reset
    procedure :: claim => shared_counter_claim
    procedure :: progress => shared_counter_progress
    procedure :: complete_prefetch => shared_counter_complete_prefetch
    procedure :: destroy => shared_counter_destroy
  end type mpi_shared_counter_type

  !> 先行発行した fetch-and-add の結果バッファ。MPI が非同期に書くので volatile にし、run 中は 1 pool だけ使う。
  integer(i64), save, volatile, target :: prefetch_result = 0_i64
  integer(i64), save, target :: prefetch_count = 0_i64

  public :: mpi_context
  public :: mpi_initialize
  public :: mpi_shutdown
  public :: mpi_is_root
  public :: mpi_world_size
  public :: mpi_get_rank_size
  public :: mpi_split_count
  public :: mpi_select_lowest_rank_i32_values
  public :: mpi_allreduce_sum_real_dp_array
  public :: mpi_allreduce_sum_real_dp_scalar
  public :: mpi_allreduce_min_real_dp_array
  public :: mpi_allreduce_max_real_dp_array
  public :: mpi_allreduce_sum_i32_array
  public :: mpi_allreduce_sum_i64_array
  public :: mpi_allreduce_sum_i32_scalar
  public :: mpi_allreduce_min_i32_scalar
  public :: mpi_allreduce_max_i32_scalar
  public :: mpi_bcast_i32_array
  public :: mpi_bcast_real_dp_array
  public :: mpi_gatherv_real_dp_array
  public :: mpi_allgatherv_inplace_real_dp
  public :: mpi_allgatherv_inplace_i32
  public :: mpi_world_barrier

contains

  !> MPIを初期化して rank / size を取得する。非MPIビルドでは単一ランクを返す。
  subroutine mpi_initialize(ctx)
    type(mpi_context), intent(out) :: ctx
#ifdef USE_MPI
    logical :: is_initialized
    integer :: ierr
    integer :: rank_int, size_int, provided
#endif
    ctx = mpi_context()
#ifdef USE_MPI
    call MPI_Initialized(is_initialized, ierr)
    if (.not. is_initialized) then
      ! shared particle pool は OpenMP critical 区間から RMA を呼ぶので SERIALIZED を求める。
      call MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
      ctx%initialized_here = .true.
    else
      call MPI_Query_thread(provided, ierr)
    end if
    ctx%thread_level = int(provided, i32)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank_int, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size_int, ierr)
    ctx%rank = int(rank_int, i32)
    ctx%size = int(size_int, i32)
    ctx%enabled = (ctx%size > 1_i32)
#endif
    ! MPI ビルド時に launcher 環境変数から rank/size を補完して
    ! root 専用ログの重複を避ける。非 MPI build では collective が no-op
    ! になるため、複数 task launcher を検出したら明示的に停止する。
#ifdef USE_MPI
    if (ctx%size <= 1_i32) call infer_launcher_rank_size(ctx)
#else
    call infer_launcher_rank_size(ctx)
    if (ctx%size > 1_i32) then
      error stop 'BEACH was built without MPI but launcher reports multiple ranks. Rebuild with USE_MPI or run one task.'
    end if
    ctx = mpi_context()
#endif
  end subroutine mpi_initialize

  !> `mpi_initialize` が実際に初期化した場合のみ `MPI_Finalize` を呼ぶ。
  subroutine mpi_shutdown(ctx)
    type(mpi_context), intent(inout) :: ctx
#ifdef USE_MPI
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

  !> 非MPIビルド時に launcher 環境変数から rank / size を補完する。
  subroutine infer_launcher_rank_size(ctx)
    type(mpi_context), intent(inout) :: ctx
    logical :: found

    call try_launcher_env_pair(ctx, 'OMPI_COMM_WORLD_RANK', 'OMPI_COMM_WORLD_SIZE', found)
    if (found) return
    call try_launcher_env_pair(ctx, 'PMI_RANK', 'PMI_SIZE', found)
    if (found) return
    call try_launcher_env_pair(ctx, 'MV2_COMM_WORLD_RANK', 'MV2_COMM_WORLD_SIZE', found)
  end subroutine infer_launcher_rank_size

  !> rank/size の環境変数ペアを解釈できたときだけ `ctx` を更新する。
  subroutine try_launcher_env_pair(ctx, rank_name, size_name, found)
    type(mpi_context), intent(inout) :: ctx
    character(len=*), intent(in) :: rank_name, size_name
    logical, intent(out) :: found
    integer(i32) :: rank_value, size_value
    logical :: has_rank, has_size

    found = .false.
    call read_env_i32(rank_name, rank_value, has_rank)
    call read_env_i32(size_name, size_value, has_size)
    if (.not. has_rank .or. .not. has_size) return
    if (size_value <= 0_i32) return
    if (rank_value < 0_i32 .or. rank_value >= size_value) return

    ctx%rank = rank_value
    ctx%size = size_value
    ctx%enabled = (ctx%size > 1_i32)
    found = .true.
  end subroutine try_launcher_env_pair

  !> 整数環境変数を読み取る。未設定や parse 失敗時は `found=.false.`。
  subroutine read_env_i32(name, value, found)
    character(len=*), intent(in) :: name
    integer(i32), intent(out) :: value
    logical, intent(out) :: found
    integer :: status, length, ios
    character(len=64) :: raw

    value = 0_i32
    found = .false.
    raw = ''
    call get_environment_variable(name, raw, length=length, status=status)
    if (status /= 0 .or. length <= 0 .or. length > len(raw)) return

    read (raw(:length), *, iostat=ios) value
    if (ios /= 0) return
    found = .true.
  end subroutine read_env_i32

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

  !> Select values supplied by the lowest MPI rank whose local flag is true.
  subroutine mpi_select_lowest_rank_i32_values(ctx, local_present, local_values, selected_rank, selected_values)
    type(mpi_context), intent(in) :: ctx
    logical, intent(in) :: local_present
    integer(i32), intent(in) :: local_values(:)
    integer(i32), intent(out) :: selected_rank
    integer(i32), intent(out) :: selected_values(:)

    integer :: base, nvalues, stride
    integer(i32) :: candidate_rank, local_rank, world_size
    integer(i32), allocatable :: packed(:)

    if (size(selected_values) /= size(local_values)) then
      error stop 'mpi_select_lowest_rank_i32_values requires matching value array sizes.'
    end if

    call mpi_get_rank_size(local_rank, world_size, ctx)
    nvalues = size(local_values)
    stride = nvalues + 1
    allocate (packed(int(world_size)*stride))
    packed = 0_i32
    if (local_present) then
      base = int(local_rank)*stride
      packed(base + 1) = 1_i32
      if (nvalues > 0) packed(base + 2:base + 1 + nvalues) = local_values
    end if

    call mpi_allreduce_sum_i32_array(ctx, packed)

    selected_rank = -1_i32
    selected_values = 0_i32
    do candidate_rank = 0_i32, world_size - 1_i32
      base = int(candidate_rank)*stride
      if (packed(base + 1) <= 0_i32) cycle
      selected_rank = candidate_rank
      if (nvalues > 0) selected_values = packed(base + 2:base + 1 + nvalues)
      exit
    end do
  end subroutine mpi_select_lowest_rank_i32_values

  !> 倍精度配列の総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_real_dp_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: values(:)
#ifdef USE_MPI
    integer :: ierr

    if (.not. ctx%enabled) return
    call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_allreduce_sum_real_dp_array

  !> 倍精度スカラの総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_real_dp_scalar(ctx, value)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: value
#ifdef USE_MPI
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
    integer :: ierr

    if (.not. ctx%enabled) return
    call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_allreduce_min_real_dp_array

  !> 倍精度配列の最大値Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_max_real_dp_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: values(:)
#ifdef USE_MPI
    integer :: ierr

    if (.not. ctx%enabled) return
    call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_allreduce_max_real_dp_array

  !> 32bit整数配列の総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_i32_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: values(:)
#ifdef USE_MPI
    integer :: ierr

    if (.not. ctx%enabled) return
    call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_allreduce_sum_i32_array

  !> 64bit整数配列の総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_i64_array(ctx, values)
    type(mpi_context), intent(in) :: ctx
    integer(i64), intent(inout) :: values(:)
#ifdef USE_MPI
    integer :: ierr

    if (.not. ctx%enabled) return
    call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_allreduce_sum_i64_array

  !> 32bit整数スカラの総和Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_sum_i32_scalar(ctx, value)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: value
#ifdef USE_MPI
    integer :: sendval, recvval, ierr

    if (.not. ctx%enabled) return
    sendval = int(value, kind=kind(0))
    call MPI_Allreduce(sendval, recvval, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    value = int(recvval, kind=i32)
#endif
  end subroutine mpi_allreduce_sum_i32_scalar

  !> 32bit整数スカラの最小値Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_min_i32_scalar(ctx, value)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: value
#ifdef USE_MPI
    integer :: sendval, recvval, ierr

    if (.not. ctx%enabled) return
    sendval = int(value, kind=kind(0))
    call MPI_Allreduce(sendval, recvval, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    value = int(recvval, kind=i32)
#endif
  end subroutine mpi_allreduce_min_i32_scalar

  !> 32bit整数スカラの最大値Allreduceをin-placeで実行する。
  subroutine mpi_allreduce_max_i32_scalar(ctx, value)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: value
#ifdef USE_MPI
    integer :: sendval, recvval, ierr

    if (.not. ctx%enabled) return
    sendval = int(value, kind=kind(0))
    call MPI_Allreduce(sendval, recvval, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    value = int(recvval, kind=i32)
#endif
  end subroutine mpi_allreduce_max_i32_scalar

  !> root rankの32bit整数配列を全rankへbroadcastする。
  subroutine mpi_bcast_i32_array(ctx, values, root)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: values(:)
    integer(i32), intent(in) :: root
#ifdef USE_MPI
    integer, allocatable :: buffer(:)
    integer :: ierr
#endif

    if (root < 0_i32 .or. root >= max(1_i32, ctx%size)) then
      error stop 'mpi_bcast_i32_array root out of range.'
    end if
#ifdef USE_MPI
    if (.not. ctx%enabled) return
    allocate (buffer(size(values)))
    buffer = int(values, kind=kind(0))
    call MPI_Bcast(buffer, size(buffer), MPI_INTEGER, int(root), MPI_COMM_WORLD, ierr)
    values = int(buffer, kind=i32)
#endif
  end subroutine mpi_bcast_i32_array

  !> root rankの倍精度配列を全rankへbroadcastする。
  subroutine mpi_bcast_real_dp_array(ctx, values, root)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: values(:)
    integer(i32), intent(in) :: root
#ifdef USE_MPI
    integer :: ierr
#endif

    if (root < 0_i32 .or. root >= max(1_i32, ctx%size)) then
      error stop 'mpi_bcast_real_dp_array root out of range.'
    end if
#ifdef USE_MPI
    if (.not. ctx%enabled) return
    call MPI_Bcast(values, size(values), MPI_DOUBLE_PRECISION, int(root), MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_bcast_real_dp_array

  !> 各rankの可変長倍精度配列をrootへrank順に連結する。
  !!
  !! 非rootでは`global_values`を長さ0で返す。非MPIまたは単一rankでは
  !! local_valuesをそのまま複製する。
  subroutine mpi_gatherv_real_dp_array(ctx, local_values, global_values, root)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(in) :: local_values(:)
    real(dp), allocatable, intent(out) :: global_values(:)
    integer(i32), intent(in) :: root
#ifdef USE_MPI
    integer, allocatable :: counts(:), displacements(:)
    integer :: ierr, local_count, rank_index, total_count
#endif

    if (root < 0_i32 .or. root >= max(1_i32, ctx%size)) then
      error stop 'mpi_gatherv_real_dp_array root out of range.'
    end if
#ifdef USE_MPI
    if (ctx%enabled) then
      allocate (counts(int(ctx%size)), displacements(int(ctx%size)))
      counts = 0
      displacements = 0
      local_count = size(local_values)
      call MPI_Gather( &
        local_count, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, int(root), MPI_COMM_WORLD, ierr &
        )
      if (ctx%rank == root) then
        total_count = 0
        do rank_index = 1, int(ctx%size)
          displacements(rank_index) = total_count
          total_count = total_count + counts(rank_index)
        end do
        allocate (global_values(total_count))
      else
        allocate (global_values(0))
      end if
      call MPI_Gatherv( &
        local_values, local_count, MPI_DOUBLE_PRECISION, global_values, counts, displacements, &
        MPI_DOUBLE_PRECISION, int(root), MPI_COMM_WORLD, ierr &
        )
      return
    end if
#endif
    if (ctx%rank == root) then
      allocate (global_values(size(local_values)))
      global_values = local_values
    else
      allocate (global_values(0))
    end if
  end subroutine mpi_gatherv_real_dp_array

  !> 各rankが `values` の自分の区間 (`displacements(rank+1)+1` から `counts(rank+1)` 要素) を持ち寄り、
  !! 全rankで `values` の全区間を揃える。非MPIまたは単一rankでは何もしない。
  subroutine mpi_allgatherv_inplace_real_dp(ctx, values, counts, displacements)
    type(mpi_context), intent(in) :: ctx
    real(dp), intent(inout) :: values(*)
    integer(i32), intent(in) :: counts(:), displacements(:)
#ifdef USE_MPI
    integer :: ierr

    if (.not. ctx%enabled .or. ctx%size <= 1_i32) return
    call MPI_Allgatherv( &
      MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, values, int(counts), int(displacements), MPI_DOUBLE_PRECISION, &
      MPI_COMM_WORLD, ierr &
      )
#endif
  end subroutine mpi_allgatherv_inplace_real_dp

  !> `mpi_allgatherv_inplace_real_dp` の 32 bit 整数版。
  subroutine mpi_allgatherv_inplace_i32(ctx, values, counts, displacements)
    type(mpi_context), intent(in) :: ctx
    integer(i32), intent(inout) :: values(*)
    integer(i32), intent(in) :: counts(:), displacements(:)
#ifdef USE_MPI
    integer :: ierr

    if (.not. ctx%enabled .or. ctx%size <= 1_i32) return
    call MPI_Allgatherv( &
      MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, values, int(counts), int(displacements), MPI_INTEGER4, &
      MPI_COMM_WORLD, ierr &
      )
#endif
  end subroutine mpi_allgatherv_inplace_i32

  !> rank 0 に 1 要素の RMA 窓を作り、全 rank で passive-target lock を張る（collective）。
  subroutine shared_counter_init(self, ctx)
    class(mpi_shared_counter_type), intent(inout) :: self
    type(mpi_context), intent(in) :: ctx
#ifdef USE_MPI
    integer(MPI_ADDRESS_KIND) :: window_bytes
    type(c_ptr) :: base_ptr
    integer(i64), pointer :: base_value
    integer :: ierr
#endif

    call self%destroy()
    self%active = .true.
    self%distributed = .false.
    self%local_value = 0_i64
#ifdef USE_MPI
    if (.not. ctx%enabled) return
    if (ctx%thread_level < MPI_THREAD_SERIALIZED) then
      error stop 'shared particle pool needs MPI_THREAD_SERIALIZED, but MPI_Init_thread provided less.'
    end if
    window_bytes = 0_MPI_ADDRESS_KIND
    if (ctx%rank == 0_i32) window_bytes = 8_MPI_ADDRESS_KIND
    call MPI_Win_allocate(window_bytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, base_ptr, self%win, ierr)
    if (ierr /= MPI_SUCCESS) error stop 'MPI_Win_allocate failed for the shared particle pool counter.'
    if (ctx%rank == 0_i32) then
      call c_f_pointer(base_ptr, base_value)
      base_value = 0_i64
    end if
    call MPI_Win_lock_all(MPI_MODE_NOCHECK, self%win, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    self%distributed = .true.
    self%hosts_window = ctx%rank == 0_i32
#endif
  end subroutine shared_counter_init

  !> カウンタを 0 に戻す（collective）。前バッチの claim が全て終わってから呼ぶ。
  subroutine shared_counter_reset(self, ctx)
    class(mpi_shared_counter_type), intent(inout) :: self
    type(mpi_context), intent(in) :: ctx
#ifdef USE_MPI
    integer(i64) :: zero_value, previous
    integer :: ierr
#endif

    if (.not. self%active) error stop 'shared particle pool counter is not initialized.'
    self%local_value = 0_i64
#ifdef USE_MPI
    if (.not. self%distributed) return
    call self%complete_prefetch()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (ctx%rank == 0_i32) then
      zero_value = 0_i64
      call MPI_Fetch_and_op(zero_value, previous, MPI_INTEGER8, 0, 0_MPI_ADDRESS_KIND, MPI_REPLACE, self%win, ierr)
      call MPI_Win_flush(0, self%win, ierr)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  end subroutine shared_counter_reset

  !> `count` 個を確保し、確保前のカウンタ値（0 始まりの先頭 slot）を返す。
  subroutine shared_counter_claim(self, count, first)
    class(mpi_shared_counter_type), intent(inout) :: self
    integer(i64), intent(in) :: count
    integer(i64), intent(out) :: first
#ifdef USE_MPI
    integer :: ierr
#endif

#ifdef USE_MPI
    if (self%distributed) then
      ! The claim for the next chunk is issued right away and completed at the following call, so the
      ! round trip to the window owner overlaps with the tracking of the current chunk.
      if (self%prefetch_pending .and. prefetch_count == count) then
        call MPI_Win_flush(0, self%win, ierr)
        first = prefetch_result
      else
        if (self%prefetch_pending) call MPI_Win_flush(0, self%win, ierr)
        call MPI_Fetch_and_op(count, first, MPI_INTEGER8, 0, 0_MPI_ADDRESS_KIND, MPI_SUM, self%win, ierr)
        call MPI_Win_flush(0, self%win, ierr)
      end if
      prefetch_count = count
      call MPI_Fetch_and_op( &
        prefetch_count, prefetch_result, MPI_INTEGER8, 0, 0_MPI_ADDRESS_KIND, MPI_SUM, self%win, ierr &
        )
      self%prefetch_pending = .true.
      return
    end if
#endif
    first = self%local_value
    self%local_value = self%local_value + count
  end subroutine shared_counter_claim

  !> 窓を持つ rank が MPI の進行エンジンを回し、他 rank から届いた fetch-and-add を処理する。
  !! passive-target RMA は target が MPI を呼ばないと進まないので、rank 0 の粒子ループから定期的に呼ぶ。
  subroutine shared_counter_progress(self)
    class(mpi_shared_counter_type), intent(inout) :: self
#ifdef USE_MPI
    integer :: ierr
    logical :: flag
    integer :: status(MPI_STATUS_SIZE)

    if (.not. self%distributed .or. .not. self%hosts_window) return
    call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, status, ierr)
    call MPI_Win_flush(0, self%win, ierr)
#endif
  end subroutine shared_counter_progress

  !> 先行発行した fetch-and-add があれば完了させて捨てる（barrier の前に呼ぶ）。
  subroutine shared_counter_complete_prefetch(self)
    class(mpi_shared_counter_type), intent(inout) :: self
#ifdef USE_MPI
    integer :: ierr

    if (self%distributed .and. self%prefetch_pending) call MPI_Win_flush(0, self%win, ierr)
#endif
    self%prefetch_pending = .false.
  end subroutine shared_counter_complete_prefetch

  !> RMA 窓を解放する（collective。`init` した全 rank が呼ぶ）。
  subroutine shared_counter_destroy(self)
    class(mpi_shared_counter_type), intent(inout) :: self
#ifdef USE_MPI
    integer :: ierr

    if (self%distributed) then
      call self%complete_prefetch()
      call MPI_Win_unlock_all(self%win, ierr)
      call MPI_Win_free(self%win, ierr)
      self%win = MPI_WIN_NULL
    end if
#endif
    self%distributed = .false.
    self%hosts_window = .false.
    self%active = .false.
    self%local_value = 0_i64
  end subroutine shared_counter_destroy

  !> 全rankの同期ポイント。
  subroutine mpi_world_barrier(ctx)
    type(mpi_context), intent(in) :: ctx
#ifdef USE_MPI
    integer :: ierr
    if (.not. ctx%enabled) return
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  end subroutine mpi_world_barrier

end module bem_mpi
