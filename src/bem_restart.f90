!> 出力ディレクトリに保存したチェックポイントの保存/復元を扱う補助モジュール。
module bem_restart
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_stats, mesh_type, injection_state
  implicit none

  private
  public :: load_restart_checkpoint
  public :: write_rng_state_file
  public :: write_macro_residuals_file

contains

  !> 既存出力ディレクトリから統計・要素電荷・乱数状態を復元する。
  !! @param[in] out_dir チェックポイントを読む出力ディレクトリ。
  !! @param[inout] mesh 現在のメッシュ。`q_elem` を復元値で上書きする。
  !! @param[out] stats 復元された統計値。
  !! @param[out] has_restart 復元可能なチェックポイントが存在したか。
  subroutine load_restart_checkpoint(out_dir, mesh, stats, has_restart, state)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(inout) :: mesh
    type(sim_stats), intent(out) :: stats
    logical, intent(out) :: has_restart
    type(injection_state), intent(inout), optional :: state

    character(len=1024) :: summary_path, charges_path, rng_path, residual_path
    logical :: has_summary, has_charges, has_rng, has_residual

    stats = sim_stats()
    has_restart = .false.

    summary_path = trim(out_dir) // '/summary.txt'
    charges_path = trim(out_dir) // '/charges.csv'
    rng_path = trim(out_dir) // '/rng_state.txt'
    residual_path = trim(out_dir) // '/macro_residuals.csv'

    inquire(file=trim(summary_path), exist=has_summary)
    inquire(file=trim(charges_path), exist=has_charges)
    inquire(file=trim(rng_path), exist=has_rng)
    inquire(file=trim(residual_path), exist=has_residual)

    if (.not. has_summary .and. .not. has_charges .and. .not. has_rng) return

    if (.not. (has_summary .and. has_charges .and. has_rng)) then
      error stop 'Resume requested but checkpoint files are incomplete in output directory.'
    end if

    call load_summary_file(trim(summary_path), mesh%nelem, stats)
    call load_charge_file(trim(charges_path), mesh)
    call restore_rng_state(trim(rng_path))
    if (present(state)) then
      if (allocated(state%macro_residual)) state%macro_residual = 0.0d0
      if (has_residual) call load_macro_residual_file(trim(residual_path), state)
    end if
    has_restart = .true.
  end subroutine load_restart_checkpoint

  !> 現在の Fortran 乱数状態をファイルへ保存する。
  !! @param[in] out_dir 出力ディレクトリ。
  subroutine write_rng_state_file(out_dir)
    character(len=*), intent(in) :: out_dir

    character(len=1024) :: path
    integer :: n, u, ios, i
    integer, allocatable :: seed(:)

    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)

    path = trim(out_dir) // '/rng_state.txt'
    open(newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open rng_state.txt.'

    write(u, '(i0)') n
    do i = 1, n
      write(u, '(i0)') seed(i)
    end do
    close(u)
  end subroutine write_rng_state_file

  !> マクロ粒子残差を `macro_residuals.csv` として保存する。
  subroutine write_macro_residuals_file(out_dir, state)
    character(len=*), intent(in) :: out_dir
    type(injection_state), intent(in) :: state

    character(len=1024) :: path
    integer :: u, ios, i

    if (.not. allocated(state%macro_residual)) return

    path = trim(out_dir) // '/macro_residuals.csv'
    open(newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open macro_residuals.csv.'

    write(u, '(a)') 'species_idx,residual'
    do i = 1, size(state%macro_residual)
      write(u, '(i0,a,es24.16)') i, ',', state%macro_residual(i)
    end do
    close(u)
  end subroutine write_macro_residuals_file

  !> `summary.txt` を読み込み、必須キーの存在と要素数整合を検証する。
  !! 欠落キーやメッシュ要素数不一致は再開不能として停止する。
  subroutine load_summary_file(path, expected_nelem, stats)
    character(len=*), intent(in) :: path
    integer(i32), intent(in) :: expected_nelem
    type(sim_stats), intent(out) :: stats

    integer :: u, ios, pos
    integer(i32) :: mesh_nelem
    character(len=512) :: line
    character(len=64) :: key
    character(len=256) :: value
    logical :: found_mesh, found_processed, found_absorbed, found_escaped
    logical :: found_batches, found_rel

    stats = sim_stats()
    mesh_nelem = -1_i32
    found_mesh = .false.
    found_processed = .false.
    found_absorbed = .false.
    found_escaped = .false.
    found_batches = .false.
    found_rel = .false.

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary.txt for resume.'

    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = trim(adjustl(line))
      if (len_trim(line) == 0) cycle
      pos = index(line, '=')
      if (pos <= 0) cycle

      key = trim(adjustl(line(:pos - 1)))
      value = trim(adjustl(line(pos + 1:)))

      select case (trim(key))
      case ('mesh_nelem')
        read(value, *) mesh_nelem
        found_mesh = .true.
      case ('processed_particles')
        read(value, *) stats%processed_particles
        found_processed = .true.
      case ('absorbed')
        read(value, *) stats%absorbed
        found_absorbed = .true.
      case ('escaped')
        read(value, *) stats%escaped
        found_escaped = .true.
      case ('batches')
        read(value, *) stats%batches
        found_batches = .true.
      case ('escaped_boundary')
        read(value, *) stats%escaped_boundary
      case ('survived_max_step')
        read(value, *) stats%survived_max_step
      case ('last_rel_change')
        read(value, *) stats%last_rel_change
        found_rel = .true.
      end select
    end do
    close(u)

    if (.not. (found_mesh .and. found_processed .and. found_absorbed .and. found_escaped .and. found_batches .and. found_rel)) then
      error stop 'Resume checkpoint summary is missing required keys.'
    end if
    if (mesh_nelem /= expected_nelem) then
      error stop 'Resume checkpoint mesh element count does not match current mesh.'
    end if
  end subroutine load_summary_file

  !> `charges.csv` を読み込み、各要素の電荷をメッシュへ復元する。
  !! 行重複や要素数不足を検出した場合は停止する。
  subroutine load_charge_file(path, mesh)
    character(len=*), intent(in) :: path
    type(mesh_type), intent(inout) :: mesh

    integer :: u, ios
    integer(i32) :: elem_idx, n_loaded
    real(dp) :: charge
    character(len=512) :: header
    logical, allocatable :: seen(:)

    if (.not. allocated(mesh%q_elem)) error stop 'Mesh charges are not allocated.'

    allocate(seen(mesh%nelem))
    seen = .false.
    mesh%q_elem = 0.0d0
    n_loaded = 0_i32

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charges.csv for resume.'

    read(u, '(A)', iostat=ios) header
    if (ios /= 0) error stop 'Failed to read charges.csv header.'

    do
      read(u, *, iostat=ios) elem_idx, charge
      if (ios < 0) exit
      if (ios > 0) error stop 'Failed to parse charges.csv during resume.'
      if (elem_idx < 1_i32 .or. elem_idx > mesh%nelem) then
        error stop 'Resume checkpoint charges.csv has an invalid element index.'
      end if
      if (seen(elem_idx)) then
        error stop 'Resume checkpoint charges.csv contains duplicate element rows.'
      end if
      seen(elem_idx) = .true.
      mesh%q_elem(elem_idx) = charge
      n_loaded = n_loaded + 1_i32
    end do
    close(u)

    if (n_loaded /= mesh%nelem) then
      error stop 'Resume checkpoint charges.csv does not match the current mesh.'
    end if
  end subroutine load_charge_file

  !> 保存済み乱数状態を読み戻し、このビルドの RNG 状態へ復元する。
  !! RNG 内部状態の長さが一致しない場合は互換性がないため停止する。
  subroutine restore_rng_state(path)
    character(len=*), intent(in) :: path

    integer :: expected_n, file_n, u, ios, i
    integer, allocatable :: seed(:)

    call random_seed(size=expected_n)

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open rng_state.txt for resume.'

    read(u, *, iostat=ios) file_n
    if (ios /= 0) error stop 'Failed to read rng_state.txt header.'
    if (file_n /= expected_n) then
      error stop 'Resume checkpoint RNG state size does not match this build.'
    end if

    allocate(seed(file_n))
    do i = 1, file_n
      read(u, *, iostat=ios) seed(i)
      if (ios /= 0) error stop 'Failed to parse rng_state.txt.'
    end do
    close(u)

    call random_seed(put=seed)
  end subroutine restore_rng_state

  !> 保存済みマクロ粒子残差を読み戻す。
  subroutine load_macro_residual_file(path, state)
    character(len=*), intent(in) :: path
    type(injection_state), intent(inout) :: state

    integer :: u, ios
    integer(i32) :: species_idx
    real(dp) :: residual
    character(len=512) :: header
    logical, allocatable :: seen(:)

    if (.not. allocated(state%macro_residual)) return

    allocate(seen(size(state%macro_residual)))
    seen = .false.
    state%macro_residual = 0.0d0

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Failed to open macro_residuals.csv for resume.'

    read(u, '(A)', iostat=ios) header
    if (ios /= 0) error stop 'Failed to read macro_residuals.csv header.'

    do
      read(u, *, iostat=ios) species_idx, residual
      if (ios < 0) exit
      if (ios > 0) error stop 'Failed to parse macro_residuals.csv during resume.'
      if (species_idx < 1_i32 .or. species_idx > size(state%macro_residual)) then
        error stop 'Resume checkpoint macro_residuals.csv has an invalid species index.'
      end if
      if (seen(species_idx)) then
        error stop 'Resume checkpoint macro_residuals.csv contains duplicate species rows.'
      end if
      seen(species_idx) = .true.
      state%macro_residual(species_idx) = residual
    end do
    close(u)
  end subroutine load_macro_residual_file

end module bem_restart
