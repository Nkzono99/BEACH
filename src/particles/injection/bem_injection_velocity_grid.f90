!> CSV速度grid reservoirの読込、補間、sampling。
module bem_injection_velocity_grid
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp
  use bem_string_utils, only: lower_ascii
  use bem_injection_geometry, only: resolve_face_axes, resolve_face_geometry
  implicit none
  private

  type :: velocity_grid_snapshot_entry
    character(len=:), allocatable :: path
    real(dp), allocatable :: grid_v(:, :)
    real(dp), allocatable :: f(:)
  end type velocity_grid_snapshot_entry

  type(velocity_grid_snapshot_entry), allocatable, save :: velocity_grid_snapshots(:)

  public :: sample_reservoir_velocity_grid_particles
  public :: get_velocity_grid_snapshot
  public :: reset_velocity_grid_snapshot_cache

contains

  !> 速度グリッド分布から reservoir_face 粒子をサンプルする。
  !! `velocity_grid_pdf_kind="phase_space"` では `max(v_n,0) f(v)` で流入粒子を選び、
  !! `"flux_weighted"` では入力 `f` を流入粒子分布として扱う。
  subroutine sample_reservoir_velocity_grid_particles( &
    box_min, box_max, inject_face, pos_low, pos_high, velocity_grid_path, velocity_grid_pdf_kind, batch_duration, x, v, &
    barrier_normal_energy, vmin_normal, position_jitter_dt, apply_barrier_energy_shift, velocity_grid_sampling &
    )
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    character(len=*), intent(in) :: velocity_grid_path, velocity_grid_pdf_kind
    real(dp), intent(in) :: batch_duration
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(in), optional :: barrier_normal_energy
    real(dp), intent(in), optional :: vmin_normal
    real(dp), intent(in), optional :: position_jitter_dt
    logical, intent(in), optional :: apply_barrier_energy_shift
    character(len=*), intent(in), optional :: velocity_grid_sampling

    integer :: i, axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), vn_floor, barrier, jitter_dt, vn_inf, vn_out
    real(dp), allocatable :: u(:, :), tau(:)
    logical :: apply_energy_shift
    character(len=16) :: grid_sampling

    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop "reservoir particle arrays must have first dimension 3"
    if (size(x, 2) /= size(v, 2)) error stop "reservoir x/v size mismatch"
    if (batch_duration < 0.0_dp) error stop "batch_duration must be >= 0"
    if (present(position_jitter_dt)) then
      if (position_jitter_dt < 0.0_dp) error stop "position_jitter_dt must be >= 0"
      jitter_dt = position_jitter_dt
    else
      jitter_dt = 0.0_dp
    end if
    apply_energy_shift = .true.
    if (present(apply_barrier_energy_shift)) apply_energy_shift = apply_barrier_energy_shift
    grid_sampling = 'auto'
    if (present(velocity_grid_sampling)) grid_sampling = lower_ascii(trim(velocity_grid_sampling))
    if (size(x, 2) == 0) return

    call resolve_face_geometry(box_min, box_max, inject_face, axis_n, boundary_value, inward_normal)
    call resolve_face_axes(inject_face, axis_t1, axis_t2)

    barrier = 0.0_dp
    if (present(barrier_normal_energy)) barrier = barrier_normal_energy
    vn_floor = 0.0_dp
    if (barrier > 0.0_dp) vn_floor = sqrt(barrier)
    if (present(vmin_normal)) vn_floor = max(vn_floor, max(0.0_dp, vmin_normal))

    call sample_velocity_grid_distribution(velocity_grid_path, velocity_grid_pdf_kind, grid_sampling, inward_normal, vn_floor, v)
    if (apply_energy_shift) then
      do i = 1, size(v, 2)
        vn_inf = dot_product(v(:, i), inward_normal)
        vn_out = sqrt(max(0.0_dp, vn_inf*vn_inf - barrier))
        v(:, i) = v(:, i) - inward_normal*vn_inf + inward_normal*vn_out
      end do
    end if

    allocate (u(2, size(x, 2)))
    call random_number(u)
    if (jitter_dt > 0.0_dp) then
      allocate (tau(size(x, 2)))
      ! Keep the historical draw count so MPI ranks and restarted runs retain the same RNG stream.
      call random_number(tau)
    end if

    do i = 1, size(x, 2)
      x(:, i) = 0.0_dp
      ! One representable step inward avoids a zero-time boundary event without an untracked flight segment.
      x(axis_n, i) = nearest(boundary_value, inward_normal(axis_n))
      x(axis_t1, i) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1))*u(1, i)
      x(axis_t2, i) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2))*u(2, i)
    end do
  end subroutine sample_reservoir_velocity_grid_particles

  !> CSV 速度グリッドから、流入条件で重み付けした速度をサンプルする。
  subroutine sample_velocity_grid_distribution(path, pdf_kind, grid_sampling, inward_normal, vmin_normal, v)
    character(len=*), intent(in) :: path, pdf_kind, grid_sampling
    real(dp), intent(in) :: inward_normal(3), vmin_normal
    real(dp), intent(out) :: v(:, :)

    real(dp), allocatable :: grid_v(:, :), f(:), weights(:), cdf(:)
    real(dp) :: f_sum, w_sum, vn, draw
    integer :: i, j, ngrid, idx
    character(len=16) :: kind

    if (size(v, 1) /= 3) error stop "v first dimension must be 3"
    if (size(v, 2) == 0) return
    call get_velocity_grid_snapshot(path, grid_v, f)
    ngrid = size(f)
    if (ngrid <= 0) error stop "velocity grid must contain at least one row"

    f_sum = sum(f)
    if (.not. ieee_is_finite(f_sum) .or. f_sum <= 0.0_dp) error stop "velocity grid f sum must be > 0"
    f = f/f_sum
    kind = lower_ascii(trim(pdf_kind))
    select case (trim(lower_ascii(grid_sampling)))
    case ('auto')
      if (try_sample_velocity_grid_interpolated(grid_v, f, kind, inward_normal, vmin_normal, v)) return
    case ('rectilinear')
      if (.not. try_sample_velocity_grid_interpolated(grid_v, f, kind, inward_normal, vmin_normal, v)) then
        error stop 'velocity_grid_sampling="rectilinear" requires a complete rectilinear grid with inward support.'
      end if
      return
    case ('discrete')
      continue
    case default
      error stop 'velocity_grid_sampling must be "auto", "rectilinear", or "discrete"'
    end select

    allocate (weights(ngrid), cdf(ngrid))
    weights = 0.0_dp
    do i = 1, ngrid
      vn = dot_product(grid_v(:, i), inward_normal)
      select case (trim(kind))
      case ('phase_space')
        if (vn >= vmin_normal .and. vn > 0.0_dp) weights(i) = vn*f(i)
      case ('flux_weighted')
        if (vn >= vmin_normal .and. vn > 0.0_dp) weights(i) = f(i)
      case default
        error stop 'velocity_grid_pdf_kind must be "phase_space" or "flux_weighted"'
      end select
    end do
    w_sum = sum(weights)
    if (.not. ieee_is_finite(w_sum) .or. w_sum <= 0.0_dp) then
      error stop "velocity grid has no inward entries after vmin_normal filtering"
    end if

    cdf(1) = weights(1)
    do i = 2, ngrid
      cdf(i) = cdf(i - 1) + weights(i)
    end do
    do j = 1, size(v, 2)
      call random_number(draw)
      draw = draw*w_sum
      idx = find_cdf_index(cdf, ngrid, draw)
      v(:, j) = grid_v(:, idx)
    end do
  end subroutine sample_velocity_grid_distribution

  !> 完全な直交速度グリッドなら、三線形補間したセル内分布からサンプルする。
  logical function try_sample_velocity_grid_interpolated( &
    grid_v, f, pdf_kind, inward_normal, vmin_normal, v &
    ) result(sampled)
    real(dp), intent(in) :: grid_v(:, :), f(:), inward_normal(3), vmin_normal
    character(len=*), intent(in) :: pdf_kind
    real(dp), intent(out) :: v(:, :)

    real(dp), allocatable :: vx(:), vy(:), vz(:), fgrid(:, :, :)
    integer, allocatable :: cell_ix(:), cell_iy(:), cell_iz(:)
    real(dp), allocatable :: cell_upper(:), cdf(:)
    integer :: nx, ny, nz, ncx, ncy, ncz, ncells, nvalid
    integer :: ix, iy, iz, j, idx, attempt
    real(dp) :: upper, measure, w_sum, draw, accept_draw, density
    real(dp) :: u(3), tx, ty, tz

    sampled = .false.
    if (.not. build_rectilinear_velocity_grid(grid_v, f, vx, vy, vz, fgrid)) return
    nx = size(vx)
    ny = size(vy)
    nz = size(vz)
    if (nx == 1 .and. ny == 1 .and. nz == 1) return

    ncx = max(1, nx - 1)
    ncy = max(1, ny - 1)
    ncz = max(1, nz - 1)
    ncells = ncx*ncy*ncz
    allocate (cell_ix(ncells), cell_iy(ncells), cell_iz(ncells), cell_upper(ncells), cdf(ncells))

    nvalid = 0
    do iz = 1, ncz
      do iy = 1, ncy
        do ix = 1, ncx
          upper = cell_density_upper(vx, vy, vz, fgrid, ix, iy, iz, pdf_kind, inward_normal, vmin_normal)
          if (upper <= 0.0_dp) cycle
          measure = axis_cell_width(vx, ix)*axis_cell_width(vy, iy)*axis_cell_width(vz, iz)
          if (measure <= 0.0_dp) cycle
          nvalid = nvalid + 1
          cell_ix(nvalid) = ix
          cell_iy(nvalid) = iy
          cell_iz(nvalid) = iz
          cell_upper(nvalid) = upper
          cdf(nvalid) = upper*measure
        end do
      end do
    end do
    if (nvalid <= 0) return

    do idx = 2, nvalid
      cdf(idx) = cdf(idx) + cdf(idx - 1)
    end do
    w_sum = cdf(nvalid)
    if (.not. ieee_is_finite(w_sum) .or. w_sum <= 0.0_dp) return

    do j = 1, size(v, 2)
      do attempt = 1, 10000
        call random_number(draw)
        idx = find_cdf_index(cdf, nvalid, draw*w_sum)
        call random_number(u)
        call sample_velocity_cell_point( &
          vx, vy, vz, cell_ix(idx), cell_iy(idx), cell_iz(idx), u, v(:, j), tx, ty, tz &
          )
        density = interpolated_velocity_density( &
                  vx, vy, vz, fgrid, cell_ix(idx), cell_iy(idx), cell_iz(idx), tx, ty, tz, &
                  pdf_kind, inward_normal, vmin_normal &
                  )
        if (density <= 0.0_dp) cycle
        call random_number(accept_draw)
        if (accept_draw*cell_upper(idx) <= density) exit
      end do
      if (attempt > 10000) error stop "velocity grid interpolation rejection sampler did not converge"
    end do
    sampled = .true.
  end function try_sample_velocity_grid_interpolated

  !> CSV 行集合が完全な直交格子なら `f(ix,iy,iz)` に詰め替える。
  logical function build_rectilinear_velocity_grid(grid_v, f, vx, vy, vz, fgrid) result(ok)
    real(dp), intent(in) :: grid_v(:, :), f(:)
    real(dp), allocatable, intent(out) :: vx(:), vy(:), vz(:), fgrid(:, :, :)

    logical, allocatable :: filled(:, :, :)
    integer :: i, ix, iy, iz, nx, ny, nz

    ok = .false.
    call unique_sorted_values(grid_v(1, :), vx)
    call unique_sorted_values(grid_v(2, :), vy)
    call unique_sorted_values(grid_v(3, :), vz)
    nx = size(vx)
    ny = size(vy)
    nz = size(vz)
    if (nx*ny*nz /= size(f)) return

    allocate (fgrid(nx, ny, nz), filled(nx, ny, nz))
    fgrid = 0.0_dp
    filled = .false.
    do i = 1, size(f)
      ix = axis_index(vx, grid_v(1, i))
      iy = axis_index(vy, grid_v(2, i))
      iz = axis_index(vz, grid_v(3, i))
      if (ix <= 0 .or. iy <= 0 .or. iz <= 0) return
      if (filled(ix, iy, iz)) error stop "velocity grid CSV contains duplicate grid point"
      fgrid(ix, iy, iz) = f(i)
      filled(ix, iy, iz) = .true.
    end do
    if (.not. all(filled)) return
    ok = .true.
  end function build_rectilinear_velocity_grid

  !> 速度セルの上界密度を返す。rejection sampling の envelope に使う。
  real(dp) function cell_density_upper( &
    vx, vy, vz, fgrid, ix, iy, iz, pdf_kind, inward_normal, vmin_normal &
    ) result(upper)
    real(dp), intent(in) :: vx(:), vy(:), vz(:), fgrid(:, :, :)
    integer, intent(in) :: ix, iy, iz
    character(len=*), intent(in) :: pdf_kind
    real(dp), intent(in) :: inward_normal(3), vmin_normal

    integer :: ox, oy, oz, ixc, iyc, izc
    real(dp) :: max_f, max_vn, vn, vel(3)

    max_f = 0.0_dp
    max_vn = -huge(1.0_dp)
    do oz = 0, 1
      izc = min(iz + oz, size(vz))
      do oy = 0, 1
        iyc = min(iy + oy, size(vy))
        do ox = 0, 1
          ixc = min(ix + ox, size(vx))
          max_f = max(max_f, fgrid(ixc, iyc, izc))
          vel = [vx(ixc), vy(iyc), vz(izc)]
          vn = dot_product(vel, inward_normal)
          max_vn = max(max_vn, vn)
        end do
      end do
    end do

    upper = 0.0_dp
    if (max_f <= 0.0_dp) return
    if (max_vn < vmin_normal .or. max_vn <= 0.0_dp) return
    select case (trim(pdf_kind))
    case ('phase_space')
      upper = max_f*max_vn
    case ('flux_weighted')
      upper = max_f
    case default
      error stop 'velocity_grid_pdf_kind must be "phase_space" or "flux_weighted"'
    end select
  end function cell_density_upper

  !> 速度セル内で一様な候補点を作り、三線形補間の局所座標も返す。
  subroutine sample_velocity_cell_point(vx, vy, vz, ix, iy, iz, u, vel, tx, ty, tz)
    real(dp), intent(in) :: vx(:), vy(:), vz(:), u(3)
    integer, intent(in) :: ix, iy, iz
    real(dp), intent(out) :: vel(3), tx, ty, tz

    call sample_axis_cell_value(vx, ix, u(1), vel(1), tx)
    call sample_axis_cell_value(vy, iy, u(2), vel(2), ty)
    call sample_axis_cell_value(vz, iz, u(3), vel(3), tz)
  end subroutine sample_velocity_cell_point

  !> セル内候補点における補間済み流入密度を返す。
  real(dp) function interpolated_velocity_density( &
    vx, vy, vz, fgrid, ix, iy, iz, tx, ty, tz, pdf_kind, inward_normal, vmin_normal &
    ) result(density)
    real(dp), intent(in) :: vx(:), vy(:), vz(:), fgrid(:, :, :)
    integer, intent(in) :: ix, iy, iz
    real(dp), intent(in) :: tx, ty, tz, inward_normal(3), vmin_normal
    character(len=*), intent(in) :: pdf_kind

    real(dp) :: f_interp, vel(3), vn

    vel(1) = interpolated_axis_value(vx, ix, tx)
    vel(2) = interpolated_axis_value(vy, iy, ty)
    vel(3) = interpolated_axis_value(vz, iz, tz)
    vn = dot_product(vel, inward_normal)
    if (vn < vmin_normal .or. vn <= 0.0_dp) then
      density = 0.0_dp
      return
    end if

    f_interp = max(0.0_dp, trilinear_f(fgrid, ix, iy, iz, tx, ty, tz))
    select case (trim(pdf_kind))
    case ('phase_space')
      density = f_interp*vn
    case ('flux_weighted')
      density = f_interp
    case default
      error stop 'velocity_grid_pdf_kind must be "phase_space" or "flux_weighted"'
    end select
  end function interpolated_velocity_density

  !> 三線形補間で `f` を評価する。1点だけの軸は固定軸として扱う。
  real(dp) function trilinear_f(fgrid, ix, iy, iz, tx, ty, tz) result(value)
    real(dp), intent(in) :: fgrid(:, :, :)
    integer, intent(in) :: ix, iy, iz
    real(dp), intent(in) :: tx, ty, tz

    integer :: ix1, iy1, iz1
    real(dp) :: c00, c10, c01, c11, c0, c1

    ix1 = min(ix + 1, size(fgrid, 1))
    iy1 = min(iy + 1, size(fgrid, 2))
    iz1 = min(iz + 1, size(fgrid, 3))
    c00 = (1.0_dp - tx)*fgrid(ix, iy, iz) + tx*fgrid(ix1, iy, iz)
    c10 = (1.0_dp - tx)*fgrid(ix, iy1, iz) + tx*fgrid(ix1, iy1, iz)
    c01 = (1.0_dp - tx)*fgrid(ix, iy, iz1) + tx*fgrid(ix1, iy, iz1)
    c11 = (1.0_dp - tx)*fgrid(ix, iy1, iz1) + tx*fgrid(ix1, iy1, iz1)
    c0 = (1.0_dp - ty)*c00 + ty*c10
    c1 = (1.0_dp - ty)*c01 + ty*c11
    value = (1.0_dp - tz)*c0 + tz*c1
  end function trilinear_f

  subroutine sample_axis_cell_value(axis, idx, u, value, t)
    real(dp), intent(in) :: axis(:), u
    integer, intent(in) :: idx
    real(dp), intent(out) :: value, t

    if (size(axis) == 1) then
      value = axis(1)
      t = 0.0_dp
    else
      t = u
      value = axis(idx) + (axis(idx + 1) - axis(idx))*t
    end if
  end subroutine sample_axis_cell_value

  real(dp) function interpolated_axis_value(axis, idx, t) result(value)
    real(dp), intent(in) :: axis(:), t
    integer, intent(in) :: idx

    if (size(axis) == 1) then
      value = axis(1)
    else
      value = axis(idx) + (axis(idx + 1) - axis(idx))*t
    end if
  end function interpolated_axis_value

  real(dp) function axis_cell_width(axis, idx) result(width)
    real(dp), intent(in) :: axis(:)
    integer, intent(in) :: idx

    if (size(axis) == 1) then
      width = 1.0_dp
    else
      width = axis(idx + 1) - axis(idx)
    end if
  end function axis_cell_width

  !> CDF 配列から二分探索でサンプル位置を返す。
  integer function find_cdf_index(cdf, n, draw) result(idx)
    real(dp), intent(in) :: cdf(:), draw
    integer, intent(in) :: n
    integer :: lo, hi, mid

    lo = 1
    hi = n
    do while (lo < hi)
      mid = (lo + hi)/2
      if (draw <= cdf(mid)) then
        hi = mid
      else
        lo = mid + 1
      end if
    end do
    idx = lo
  end function find_cdf_index

  !> 実数配列から昇順 unique 値を作る。
  subroutine unique_sorted_values(values, unique)
    real(dp), intent(in) :: values(:)
    real(dp), allocatable, intent(out) :: unique(:)

    real(dp), allocatable :: tmp(:)
    real(dp) :: key
    integer :: i, j, n_unique
    logical :: found

    allocate (tmp(size(values)))
    n_unique = 0
    do i = 1, size(values)
      found = .false.
      do j = 1, n_unique
        if (tmp(j) == values(i)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        n_unique = n_unique + 1
        tmp(n_unique) = values(i)
      end if
    end do

    do i = 2, n_unique
      key = tmp(i)
      j = i - 1
      do while (j >= 1)
        if (tmp(j) <= key) exit
        tmp(j + 1) = tmp(j)
        j = j - 1
      end do
      tmp(j + 1) = key
    end do
    allocate (unique(n_unique))
    unique = tmp(1:n_unique)
  end subroutine unique_sorted_values

  integer function axis_index(axis, value) result(idx)
    real(dp), intent(in) :: axis(:), value
    integer :: i

    idx = 0
    do i = 1, size(axis)
      if (axis(i) == value) then
        idx = i
        return
      end if
    end do
  end function axis_index

  !> 最初の読込み結果をpathごとに固定し、同一runのsamplingとfingerprintで共有する。
  subroutine get_velocity_grid_snapshot(path, grid_v, f)
    character(len=*), intent(in) :: path
    real(dp), allocatable, intent(out) :: grid_v(:, :)
    real(dp), allocatable, intent(out) :: f(:)

    type(velocity_grid_snapshot_entry), allocatable :: grown(:)
    real(dp), allocatable :: loaded_grid_v(:, :), loaded_f(:)
    character(len=:), allocatable :: normalized_path
    integer :: entry, old_size

    normalized_path = trim(path)
    if (len(normalized_path) == 0) error stop 'velocity_grid_path must not be empty.'
    if (allocated(velocity_grid_snapshots)) then
      do entry = 1, size(velocity_grid_snapshots)
        if (velocity_grid_snapshots(entry)%path == normalized_path) then
          allocate (grid_v, source=velocity_grid_snapshots(entry)%grid_v)
          allocate (f, source=velocity_grid_snapshots(entry)%f)
          return
        end if
      end do
      old_size = size(velocity_grid_snapshots)
    else
      old_size = 0
    end if

    call read_velocity_grid_csv(normalized_path, loaded_grid_v, loaded_f)
    allocate (grown(old_size + 1))
    if (old_size > 0) grown(:old_size) = velocity_grid_snapshots
    grown(old_size + 1)%path = normalized_path
    grown(old_size + 1)%grid_v = loaded_grid_v
    grown(old_size + 1)%f = loaded_f
    call move_alloc(grown, velocity_grid_snapshots)
    allocate (grid_v, source=loaded_grid_v)
    allocate (f, source=loaded_f)
  end subroutine get_velocity_grid_snapshot

  !> 独立runを同一processで開始するtest/support用途にsnapshot cacheを解放する。
  subroutine reset_velocity_grid_snapshot_cache()
    if (allocated(velocity_grid_snapshots)) deallocate (velocity_grid_snapshots)
  end subroutine reset_velocity_grid_snapshot_cache

  !> `vx,vy,vz,f` CSV を一度だけ走査する。先頭の非数値行は header とみなして無視する。
  subroutine read_velocity_grid_csv(path, grid_v, f)
    character(len=*), intent(in) :: path
    real(dp), allocatable, intent(out) :: grid_v(:, :)
    real(dp), allocatable, intent(out) :: f(:)

    integer :: u, ios, parse_ios, row, capacity
    character(len=512) :: line
    real(dp) :: vx, vy, vz, weight
    logical :: skipped_header
    real(dp), allocatable :: grid_work(:, :), f_work(:)

    capacity = 64
    allocate (grid_work(3, capacity), f_work(capacity))
    row = 0
    skipped_header = .false.
    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop "could not open velocity_grid_path"
    do
      read (u, '(A)', iostat=ios) line
      if (ios < 0) exit
      if (ios > 0) then
        close (u)
        error stop 'failed to read velocity grid CSV'
      end if
      if (is_blank_or_comment(line)) cycle
      read (line, *, iostat=parse_ios) vx, vy, vz, weight
      if (parse_ios /= 0) then
        if (.not. skipped_header .and. row == 0) then
          skipped_header = .true.
          cycle
        end if
        close (u)
        error stop "invalid velocity grid CSV row"
      end if
      if (.not. all(ieee_is_finite([vx, vy, vz, weight]))) then
        close (u)
        error stop "velocity grid values must be finite"
      end if
      if (weight < 0.0_dp) then
        close (u)
        error stop "velocity grid f values must be >= 0"
      end if
      if (row >= capacity) call grow_velocity_grid_buffers(grid_work, f_work, capacity)
      row = row + 1
      grid_work(:, row) = [vx, vy, vz]
      f_work(row) = weight
    end do
    close (u)
    if (row <= 0) error stop "velocity grid CSV contains no numeric rows"
    allocate (grid_v(3, row), f(row))
    grid_v = grid_work(:, :row)
    f = f_work(:row)
  end subroutine read_velocity_grid_csv

  subroutine grow_velocity_grid_buffers(grid_v, f, capacity)
    real(dp), allocatable, intent(inout) :: grid_v(:, :), f(:)
    integer, intent(inout) :: capacity
    real(dp), allocatable :: grown_grid(:, :), grown_f(:)
    integer :: new_capacity

    new_capacity = 2*capacity
    if (new_capacity <= capacity) error stop 'velocity grid CSV row capacity overflow'
    allocate (grown_grid(3, new_capacity), grown_f(new_capacity))
    grown_grid(:, :capacity) = grid_v
    grown_f(:capacity) = f
    call move_alloc(grown_grid, grid_v)
    call move_alloc(grown_f, f)
    capacity = new_capacity
  end subroutine grow_velocity_grid_buffers

  !> 空行または `#` コメント行かを判定する。
  pure logical function is_blank_or_comment(line) result(is_skip)
    character(len=*), intent(in) :: line
    character(len=:), allocatable :: trimmed

    trimmed = adjustl(trim(line))
    is_skip = len_trim(trimmed) == 0
    if (.not. is_skip) is_skip = trimmed(1:1) == '#'
  end function is_blank_or_comment

end module bem_injection_velocity_grid
