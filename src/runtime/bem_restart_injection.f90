!> 再開に必要な乱数状態・マクロ粒子端数の保存と復元。
submodule(bem_restart) bem_restart_injection
  use bem_kinds, only: dp
  use bem_mpi, only: mpi_get_rank_size, mpi_bcast_i32_array
  implicit none
contains

  !> 現在の Fortran 乱数状態をファイルへ保存する。
  !! checkpoint transactionのownerは、全rankがこの手続きへ入る前に
  !! completion manifestを無効化し、書き終わり後に全rankを同期する。
  !! @param[in] out_dir 出力ディレクトリ。
  module procedure write_rng_state_file

  character(len=1024) :: path
  integer :: n, u, ios, i
  integer, allocatable :: seed(:)
  integer(i32) :: local_rank, world_size

  call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'write_rng_state_file')
  call random_seed(size=n)
  allocate (seed(n))
  call random_seed(get=seed)

  path = restart_rng_state_path(trim(out_dir), mpi_rank=local_rank, mpi_size=world_size)
  open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'Failed to open rng_state.txt.'

  write (u, '(i0)') n
  do i = 1, n
    write (u, '(i0)') seed(i)
  end do
  close (u)
  end procedure write_rng_state_file

  !> マクロ粒子残差を `macro_residuals.csv` として保存する。
  !! completion manifestの無効化と最終公開はcheckpoint transactionのownerが行う。
  !! @param[in] out_dir 出力ディレクトリ。
  !! @param[in] state 種別ごとのマクロ粒子残差を保持した注入状態。
  module procedure write_macro_residuals_file

  character(len=1024) :: path
  integer :: u, ios, i, face
  integer(i32) :: local_rank, world_size

  if (.not. allocated(state%macro_residual)) return
  if (allocated(state%boundary_macro_residual)) then
    if (size(state%boundary_macro_residual, 1) /= 6 .or. &
        size(state%boundary_macro_residual, 2) /= size(state%macro_residual)) then
      error stop 'injection_state boundary residual shape must be (6, nspecies).'
    end if
  end if

  call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'write_macro_residuals_file')
  if (local_rank /= 0_i32) return
  path = restart_macro_residual_path(trim(out_dir), mpi_rank=local_rank, mpi_size=world_size)
  open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'Failed to open macro_residuals.csv.'

  write (u, '(a)') 'species_idx,face,residual'
  do i = 1, size(state%macro_residual)
    write (u, '(i0,a,i0,a,es24.16)') i, ',', 0, ',', state%macro_residual(i)
    if (.not. allocated(state%boundary_macro_residual)) cycle
    do face = 1, 6
      write (u, '(i0,a,i0,a,es24.16)') i, ',', face, ',', state%boundary_macro_residual(face, i)
    end do
  end do
  close (u)
  end procedure write_macro_residuals_file

  !> 保存済み乱数状態を読み戻し、このビルドの RNG 状態へ復元する。
  !! RNG 内部状態の長さが一致しない場合は互換性がないため停止する。
  !! @param[in] path `rng_state.txt` のファイルパス。
  module procedure restore_rng_state

  integer :: expected_n, file_n, u, ios, i
  integer, allocatable :: seed(:)

  call random_seed(size=expected_n)

  open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'Failed to open rng_state.txt for resume.'

  read (u, *, iostat=ios) file_n
  if (ios /= 0) error stop 'Failed to read rng_state.txt header.'
  if (file_n /= expected_n) then
    error stop 'Resume checkpoint RNG state size does not match this build.'
  end if

  allocate (seed(file_n))
  do i = 1, file_n
    read (u, *, iostat=ios) seed(i)
    if (ios /= 0) error stop 'Failed to parse rng_state.txt.'
  end do
  close (u)

  call random_seed(put=seed)
  end procedure restore_rng_state

  !> 保存済みマクロ粒子残差を読み戻す。
  !! @param[in] path `macro_residuals.csv` のファイルパス。
  !! @param[inout] state 種別ごとのマクロ粒子残差を書き戻す注入状態。
  module procedure load_macro_residual_file

  integer :: u, ios
  integer(i32) :: species_idx, face
  real(dp) :: residual
  character(len=512) :: header
  logical :: extended_format
  logical, allocatable :: seen(:), boundary_seen(:, :)

  if (.not. allocated(state%macro_residual)) return
  if (allocated(state%boundary_macro_residual)) then
    if (size(state%boundary_macro_residual, 1) /= 6 .or. &
        size(state%boundary_macro_residual, 2) /= size(state%macro_residual)) then
      error stop 'injection_state boundary residual shape must be (6, nspecies).'
    end if
  end if

  allocate (seen(size(state%macro_residual)))
  seen = .false.
  state%macro_residual = 0.0d0
  if (allocated(state%boundary_macro_residual)) then
    allocate (boundary_seen(size(state%boundary_macro_residual, 1), size(state%boundary_macro_residual, 2)))
    boundary_seen = .false.
    state%boundary_macro_residual = 0.0d0
  end if

  open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'Failed to open macro_residuals.csv for resume.'

  read (u, '(A)', iostat=ios) header
  if (ios /= 0) error stop 'Failed to read macro_residuals.csv header.'
  select case (trim(header))
  case ('species_idx,face,residual')
    extended_format = .true.
  case ('species_idx,residual')
    extended_format = .false.
  case default
    error stop 'Resume checkpoint macro_residuals.csv has an unsupported header.'
  end select

  do
    if (extended_format) then
      read (u, *, iostat=ios) species_idx, face, residual
    else
      read (u, *, iostat=ios) species_idx, residual
      face = 0_i32
    end if
    if (ios < 0) exit
    if (ios > 0) error stop 'Failed to parse macro_residuals.csv during resume.'
    if (species_idx < 1_i32 .or. species_idx > size(state%macro_residual)) then
      error stop 'Resume checkpoint macro_residuals.csv has an invalid species index.'
    end if
    if (.not. ieee_is_finite(residual) .or. residual < 0.0d0 .or. residual >= 1.0d0) then
      error stop 'Resume checkpoint macro_residuals.csv residual values must be finite and in [0, 1).'
    end if
    if (face == 0_i32) then
      if (seen(species_idx)) error stop 'Resume checkpoint macro_residuals.csv contains duplicate source rows.'
      seen(species_idx) = .true.
      state%macro_residual(species_idx) = residual
    else
      if (.not. extended_format .or. face < 1_i32 .or. face > 6_i32) then
        error stop 'Resume checkpoint macro_residuals.csv has an invalid boundary face index.'
      end if
      if (.not. allocated(state%boundary_macro_residual)) cycle
      if (boundary_seen(face, species_idx)) then
        error stop 'Resume checkpoint macro_residuals.csv contains duplicate boundary rows.'
      end if
      boundary_seen(face, species_idx) = .true.
      state%boundary_macro_residual(face, species_idx) = residual
    end if
  end do
  close (u)
  if (.not. all(seen)) then
    error stop 'Resume checkpoint macro_residuals.csv is missing source rows.'
  end if
  if (extended_format .and. allocated(state%boundary_macro_residual)) then
    if (.not. all(boundary_seen)) then
      error stop 'Resume checkpoint macro_residuals.csv is missing boundary rows.'
    end if
  end if
  end procedure load_macro_residual_file

  !> RNG状態ファイルのパスを返す。MPI複数rank時は rank 接尾辞付きパスへ切り替える。
  module procedure restart_rng_state_path
  integer(i32) :: local_rank, world_size

  call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'restart_rng_state_path')
  if (world_size <= 1_i32) then
    path = trim(out_dir)//'/rng_state.txt'
  else
    write (path, '(a,a,i5.5,a)') trim(out_dir), '/rng_state_rank', local_rank, '.txt'
  end if
  end procedure restart_rng_state_path

  !> MPI rank数によらないglobalマクロ残差ファイルのパスを返す。
  module procedure restart_macro_residual_path
  integer(i32) :: local_rank, world_size

  call resolve_parallel_rank_size(local_rank, world_size, mpi_rank, mpi_size, mpi, 'restart_macro_residual_path')
  path = trim(out_dir)//'/macro_residuals.csv'
  end procedure restart_macro_residual_path

  !> 旧rank別マクロ残差が1個でも存在するかrootで確認し、全rankへ共有する。
  module procedure detect_legacy_ranked_residuals

  character(len=1024) :: legacy_path
  integer(i32) :: rank, found_value(1)
  logical :: exists

  found_value = 0_i32
  if (world_size <= 1_i32) then
    found = .false.
    return
  end if

  if (.not. present(mpi) .or. local_rank == 0_i32) then
    do rank = 0_i32, world_size - 1_i32
      write (legacy_path, '(a,a,i5.5,a)') trim(out_dir), '/macro_residuals_rank', rank, '.csv'
      inquire (file=trim(legacy_path), exist=exists)
      if (exists) then
        found_value(1) = 1_i32
        exit
      end if
    end do
  end if
  if (present(mpi)) call mpi_bcast_i32_array(mpi, found_value, 0_i32)
  found = found_value(1) /= 0_i32
  end procedure detect_legacy_ranked_residuals

  module procedure resolve_parallel_rank_size

  call mpi_get_rank_size(local_rank, world_size, mpi)
  if (present(mpi_rank)) local_rank = mpi_rank
  if (present(mpi_size)) world_size = mpi_size
  if (world_size <= 0_i32) error stop 'mpi_size must be > 0 in '//trim(caller_name)//'.'
  if (local_rank < 0_i32 .or. local_rank >= world_size) then
    error stop 'mpi_rank out of range in '//trim(caller_name)//'.'
  end if
  end procedure resolve_parallel_rank_size

end submodule bem_restart_injection
