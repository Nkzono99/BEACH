!> periodic2 root operator のセル単位 cache の読み書きと排他。前計算の数値部分（bem_coulomb_fmm_periodic_root_ops）
!! からファイル I/O とロックを分けたもの。呼ぶのは MPI root だけで、collective は使わない。
module bem_coulomb_fmm_periodic_cell_cache
  use bem_kinds, only: dp, i32
  use bem_filesystem, only: create_directories, filesystem_success, acquire_file_lock, release_file_lock
  use bem_coulomb_fmm_types, only: fmm_plan_type
  use bem_coulomb_fmm_periodic_cache, only: read_periodic_operator_cache, write_periodic_operator_cache, periodic_cache_ok
  implicit none
  private

  !> 欠けたセルを作って publish するまで保持する cache ロック。
  type, public :: periodic_cell_cache_lock_type
    logical :: held = .false.
    integer :: fd = -1
  end type periodic_cell_cache_lock_type

  public :: lookup_cell_cache
  public :: publish_cell_cache_and_unlock

contains

  !> required のセルを cache から読む。足りなければ cache ディレクトリを作ってロックを取り、待つ間に他の job が
  !! 足したセルを読み直す。それでも足りなければロックを保持したまま返し、呼び出し側が作ったセルを
  !! publish_cell_cache_and_unlock で書く。have は cache から得たセル。
  subroutine lookup_cell_cache(plan, nanchor, required, cached_cells, cached_blocks, have, lock)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: nanchor
    logical, intent(in) :: required(:)
    integer(i32), allocatable, intent(out) :: cached_cells(:)
    real(dp), allocatable, intent(out) :: cached_blocks(:, :, :)
    logical, intent(out) :: have(:)
    type(periodic_cell_cache_lock_type), intent(out) :: lock
    integer(i32) :: ncell
    integer :: directory_status, lock_status

    ncell = int(size(required), i32)
    call read_cell_cache(plan, nanchor, ncell, cached_cells, cached_blocks)
    have = .false.
    have(cached_cells) = .true.
    if (.not. any(required .and. .not. have)) return
    call create_directories(trim(plan%options%periodic_cache_dir), directory_status)
    if (directory_status /= filesystem_success) error stop 'Failed to create periodic operator cache directory.'
    call acquire_file_lock(trim(plan%periodic_cache_path)//'.lock', lock%fd, lock_status)
    if (lock_status /= filesystem_success) error stop 'Failed to lock periodic operator cache.'
    lock%held = .true.
    ! Another job may have added cells while this one waited for the lock.
    call read_cell_cache(plan, nanchor, ncell, cached_cells, cached_blocks)
    have = .false.
    have(cached_cells) = .true.
    if (.not. any(required .and. .not. have)) call release_cell_cache_lock(lock)
  end subroutine lookup_cell_cache

  !> cache にあったセルと新しく作ったセルを合わせて書き、ロックを外す。
  subroutine publish_cell_cache_and_unlock(plan, nanchor, cached_cells, cached_blocks, new_cells, new_blocks, lock)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: nanchor
    integer(i32), allocatable, intent(in) :: cached_cells(:), new_cells(:)
    real(dp), allocatable, intent(in) :: cached_blocks(:, :, :), new_blocks(:, :, :)
    type(periodic_cell_cache_lock_type), intent(inout) :: lock

    call publish_cell_cache(plan, nanchor, cached_cells, cached_blocks, new_cells, new_blocks)
    call release_cell_cache_lock(lock)
  end subroutine publish_cell_cache_and_unlock

  subroutine release_cell_cache_lock(lock)
    type(periodic_cell_cache_lock_type), intent(inout) :: lock
    integer :: lock_status

    call release_file_lock(lock%fd, lock_status)
    if (lock_status /= filesystem_success) error stop 'Failed to release periodic operator cache lock.'
    lock%held = .false.
  end subroutine release_cell_cache_lock

  !> Read the cell blocks stored for this fingerprint. Invalid or foreign files count as empty.
  !! The integer payload holds the source cell of every pair; each cell owns `nanchor` consecutive pairs.
  subroutine read_cell_cache(plan, nanchor, ncell, cells, blocks)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: nanchor, ncell
    integer(i32), allocatable, intent(out) :: cells(:)
    real(dp), allocatable, intent(out) :: blocks(:, :, :)
    integer(i32), allocatable :: pair_cells(:)
    logical, allocatable :: seen(:)
    integer(i32) :: cache_status, nblock, b, c

    allocate (cells(0), blocks(plan%ncoef, plan%ncoef, 0))
    call read_periodic_operator_cache( &
      trim(plan%periodic_cache_path), trim(plan%periodic_cache_fingerprint), pair_cells, blocks, cache_status, &
      expected_ncoef=plan%ncoef &
      )
    if (cache_status /= periodic_cache_ok) then
      if (allocated(blocks)) deallocate (blocks)
      allocate (blocks(plan%ncoef, plan%ncoef, 0))
      return
    end if
    if (mod(size(pair_cells), nanchor) /= 0) go to 900
    nblock = int(size(pair_cells)/nanchor, i32)
    allocate (seen(ncell))
    seen = .false.
    deallocate (cells)
    allocate (cells(nblock))
    do b = 1_i32, nblock
      c = pair_cells((b - 1_i32)*nanchor + 1_i32)
      if (c < 1_i32 .or. c > ncell) go to 900
      if (seen(c)) go to 900
      if (any(pair_cells((b - 1_i32)*nanchor + 1_i32:b*nanchor) /= c)) go to 900
      seen(c) = .true.
      cells(b) = c
    end do
    return
900 continue
    if (allocated(cells)) deallocate (cells)
    allocate (cells(0))
    deallocate (blocks)
    allocate (blocks(plan%ncoef, plan%ncoef, 0))
  end subroutine read_cell_cache

  !> Write the union of the cached and newly built cells. The caller holds the cache lock.
  subroutine publish_cell_cache(plan, nanchor, cached_cells, cached_blocks, new_cells, new_blocks)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: nanchor
    integer(i32), allocatable, intent(in) :: cached_cells(:), new_cells(:)
    real(dp), allocatable, intent(in) :: cached_blocks(:, :, :), new_blocks(:, :, :)
    integer(i32), allocatable :: pair_cells(:), cells(:)
    real(dp), allocatable :: blocks(:, :, :)
    integer(i32) :: ncached, nnew, b, cache_status

    ncached = 0_i32
    if (allocated(cached_cells)) ncached = int(size(cached_cells), i32)
    nnew = int(size(new_cells), i32)
    allocate (cells(ncached + nnew), blocks(plan%ncoef, plan%ncoef, nanchor*(ncached + nnew)))
    if (ncached > 0_i32) then
      cells(1:ncached) = cached_cells
      blocks(:, :, 1:nanchor*ncached) = cached_blocks
    end if
    cells(ncached + 1_i32:) = new_cells
    blocks(:, :, nanchor*ncached + 1_i32:) = new_blocks
    allocate (pair_cells(nanchor*(ncached + nnew)))
    do b = 1_i32, ncached + nnew
      pair_cells((b - 1_i32)*nanchor + 1_i32:b*nanchor) = cells(b)
    end do
    call write_periodic_operator_cache( &
      trim(plan%periodic_cache_path), trim(plan%periodic_cache_fingerprint), pair_cells, blocks, cache_status &
      )
    if (cache_status /= periodic_cache_ok) error stop 'Failed to publish periodic operator cache.'
  end subroutine publish_cell_cache

end module bem_coulomb_fmm_periodic_cell_cache
