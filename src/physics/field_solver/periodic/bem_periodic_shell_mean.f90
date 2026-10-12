!> periodic2 の有限画像殻の面平均を、無限面電荷の平均場から差し引いた滑らかな残り。
!!
!! 有限画像殻（中心セルと周囲 N 層の画像）の電位をセル `C` 上で面平均すると、無限面電荷の平均場
!! `k0(h) = -2*pi*|h|/A` と、滑らかな残り `s(h)` の和になる。`s` は 1/r を画像ブロック
!! `[c-(N+1/2)L, c+(N+1/2)L]^2` 上で積分した閉じた式から `k0` を引いたもので、`h=0` でも解析的である。
!! cached_kneq0 の遠方補正は面平均 0 の残差だけを当てはめ、`s` は実行時に高さの表から引く。
!! 単位は FMM core と同じ raw（1/r、Coulomb 定数なし）。
module bem_periodic_shell_mean
  use bem_kinds, only: dp, i32
  implicit none
  private

  real(dp), parameter :: pi = acos(-1.0_dp)
  ! Dunavant の 5 次の三角形積分則（重心座標の 1 成分目、重み）。
  real(dp), parameter :: dunavant_a(7) = [ &
                         1.0_dp/3.0_dp, 0.059715871789770_dp, 0.470142064105115_dp, 0.470142064105115_dp, &
                         0.797426985353087_dp, 0.101286507323456_dp, 0.101286507323456_dp]
  real(dp), parameter :: dunavant_b(7) = [ &
                         1.0_dp/3.0_dp, 0.470142064105115_dp, 0.059715871789770_dp, 0.470142064105115_dp, &
                         0.101286507323456_dp, 0.797426985353087_dp, 0.101286507323456_dp]
  real(dp), parameter :: dunavant_w(7) = [ &
                         0.225_dp, 0.132394152788506_dp, 0.132394152788506_dp, 0.132394152788506_dp, &
                         0.125939180544827_dp, 0.125939180544827_dp, 0.125939180544827_dp]

  !> 線源ごとの `s` の高さ表。高さの格子は等間隔。
  type, public :: periodic_shell_mean_plan_type
    logical :: ready = .false.
    integer(i32) :: nz = 0_i32
    integer(i32) :: periodic_axes(2) = 0_i32
    integer(i32) :: free_axis = 0_i32
    real(dp) :: z0 = 0.0_dp
    real(dp) :: dz = 0.0_dp
    !> 平均をとるセルの下端と周期長、画像の層数、softening の二乗。
    real(dp) :: cell_lo(2) = 0.0_dp
    real(dp) :: cell_len(2) = 0.0_dp
    integer(i32) :: image_layers = 0_i32
    real(dp) :: soft2 = 0.0_dp
    !> `(nz, nsrc)`。値と高さ微分。
    real(dp), allocatable :: value(:, :)
    real(dp), allocatable :: slope(:, :)
  end type periodic_shell_mean_plan_type

  !> 電荷を重ねた `S(z) = sum_i q_i s_i(z)` と `dS/dz` の格子値。
  type, public :: periodic_shell_mean_state_type
    real(dp), allocatable :: value(:)
    real(dp), allocatable :: slope(:)
  end type periodic_shell_mean_state_type

  public :: periodic_shell_mean_excess
  public :: periodic_shell_mean_point_excess
  public :: build_periodic_shell_mean_plan
  public :: refresh_periodic_shell_mean_state
  public :: eval_periodic_shell_mean
  public :: reset_periodic_shell_mean_plan
  public :: reset_periodic_shell_mean_state

contains

  !> 単位電荷 1 個の `s(h)` と `ds/dh`。
  !! @param[in] h 評価点の高さ − 線源の高さ。
  !! @param[in] x1,x2,y1,y2 線源から見た画像ブロックの範囲（x1<0<x2、y1<0<y2）。
  !! @param[in] area セルの面積。
  !! @param[in] soft2 softening の二乗（P0 panel では 0）。
  pure subroutine periodic_shell_mean_excess(h, x1, x2, y1, y2, area, soft2, s, dsdh)
    real(dp), intent(in) :: h, x1, x2, y1, y2, area, soft2
    real(dp), intent(out) :: s, dsdh
    real(dp) :: h_eff, log_part, atan_part, smooth

    h_eff = sqrt(h*h + soft2)
    log_part = corner_log(x2, y2, h_eff) - corner_log(x1, y2, h_eff) - corner_log(x2, y1, h_eff) + &
               corner_log(x1, y1, h_eff)
    atan_part = corner_atan(x2, y2, h_eff) + corner_atan(x1, y2, h_eff) + corner_atan(x2, y1, h_eff) + &
                corner_atan(x1, y1, h_eff)
    ! smooth = mu(h_eff) + 2*pi*h_eff/A。四隅の atan(xy/(h r)) の和を 2*pi - sum atan(h r/|xy|) に直し、
    ! |h| の折れ目を解析的に打ち消している。
    smooth = (log_part + h_eff*atan_part)/area
    s = smooth + 2.0_dp*pi*(abs(h) - h_eff)/area
    if (h_eff > 0.0_dp) then
      dsdh = (h/h_eff)*atan_part/area + 2.0_dp*pi*(sign(1.0_dp, h) - h/h_eff)/area
    else
      dsdh = 0.0_dp
    end if
  end subroutine periodic_shell_mean_excess

  !> 点線源 1 個の `s` と `ds/dh`。平均をとるセルと画像の層数は plan の値を使う。
  pure subroutine periodic_shell_mean_point_excess(plan, source, target_z, s, dsdh)
    type(periodic_shell_mean_plan_type), intent(in) :: plan
    real(dp), intent(in) :: source(3), target_z
    real(dp), intent(out) :: s, dsdh
    real(dp) :: lo(2), hi(2), xs, ys

    lo = plan%cell_lo - real(plan%image_layers, dp)*plan%cell_len
    hi = plan%cell_lo + real(plan%image_layers + 1_i32, dp)*plan%cell_len
    xs = source(plan%periodic_axes(1))
    ys = source(plan%periodic_axes(2))
    call periodic_shell_mean_excess( &
      target_z - source(plan%free_axis), lo(1) - xs, hi(1) - xs, lo(2) - ys, hi(2) - ys, &
      plan%cell_len(1)*plan%cell_len(2), plan%soft2, s, dsdh &
      )
  end subroutine periodic_shell_mean_point_excess

  !> 線源ごとの `s` を高さの格子に表として持つ。
  !! @param[in] periodic_axes 周期軸、free_axis 非周期軸。
  !! @param[in] cell_center 平均をとるセルの周期軸方向の中心。線源から画像ブロックの端まで約 L 離すため線源の中心に置く。
  !! @param[in] cell_len 周期長、image_layers 画像の層数、softening P0 panel では 0。
  !! @param[in] z_min,z_max 評価点がとる高さの範囲、spacing 格子間隔の上限。
  !! @param[in] vertices `(3,3,nsrc)` 三角形の頂点。panel でなければ 3 頂点とも同じ点を渡す。
  subroutine build_periodic_shell_mean_plan( &
    periodic_axes, free_axis, cell_center, cell_len, image_layers, softening, z_min, z_max, spacing, vertices, plan &
    )
    integer(i32), intent(in) :: periodic_axes(2), free_axis, image_layers
    real(dp), intent(in) :: cell_center(2), cell_len(2), softening, z_min, z_max, spacing
    real(dp), intent(in) :: vertices(:, :, :)
    type(periodic_shell_mean_plan_type), intent(inout) :: plan
    integer(i32) :: nsrc, k, i, q, nq
    real(dp) :: point(3), weight, s, dsdh, z
    logical :: is_point

    call reset_periodic_shell_mean_plan(plan)
    if (z_max < z_min .or. spacing <= 0.0_dp) error stop 'periodic shell-mean plan has an invalid height range.'
    nsrc = int(size(vertices, 3), i32)
    plan%periodic_axes = periodic_axes
    plan%free_axis = free_axis
    plan%cell_lo = cell_center - 0.5_dp*cell_len
    plan%cell_len = cell_len
    plan%image_layers = image_layers
    plan%soft2 = softening*softening
    plan%nz = max(2_i32, int(ceiling((z_max - z_min)/spacing), i32) + 1_i32)
    plan%z0 = z_min
    plan%dz = max(z_max - z_min, tiny(1.0_dp))/real(plan%nz - 1_i32, dp)
    allocate (plan%value(plan%nz, nsrc), plan%slope(plan%nz, nsrc))
    plan%value = 0.0_dp
    plan%slope = 0.0_dp
    !$omp parallel do default(shared) private(i, k, q, nq, point, weight, s, dsdh, z, is_point) schedule(static)
    do i = 1_i32, nsrc
      is_point = all(vertices(:, 2, i) == vertices(:, 1, i)) .and. all(vertices(:, 3, i) == vertices(:, 1, i))
      nq = merge(1_i32, 7_i32, is_point)
      do q = 1_i32, nq
        if (is_point) then
          point = vertices(:, 1, i)
          weight = 1.0_dp
        else
          point = dunavant_a(q)*vertices(:, 1, i) + dunavant_b(q)*vertices(:, 2, i) + &
                  (1.0_dp - dunavant_a(q) - dunavant_b(q))*vertices(:, 3, i)
          weight = dunavant_w(q)
        end if
        do k = 1_i32, plan%nz
          z = plan%z0 + real(k - 1_i32, dp)*plan%dz
          call periodic_shell_mean_point_excess(plan, point, z, s, dsdh)
          plan%value(k, i) = plan%value(k, i) + weight*s
          plan%slope(k, i) = plan%slope(k, i) + weight*dsdh
        end do
      end do
    end do
    !$omp end parallel do
    plan%ready = .true.
  end subroutine build_periodic_shell_mean_plan

  !> 線源電荷から `S` と `dS/dz` の格子値を作る。
  subroutine refresh_periodic_shell_mean_state(plan, charge, state)
    type(periodic_shell_mean_plan_type), intent(in) :: plan
    real(dp), intent(in) :: charge(:)
    type(periodic_shell_mean_state_type), intent(inout) :: state

    if (.not. plan%ready) return
    if (size(charge) /= size(plan%value, 2)) error stop 'periodic shell-mean charge size mismatch.'
    if (.not. allocated(state%value)) allocate (state%value(plan%nz), state%slope(plan%nz))
    if (size(state%value) /= plan%nz) then
      deallocate (state%value, state%slope)
      allocate (state%value(plan%nz), state%slope(plan%nz))
    end if
    state%value = matmul(plan%value, charge)
    state%slope = matmul(plan%slope, charge)
  end subroutine refresh_periodic_shell_mean_state

  !> 3 次 Hermite 補間した `S(z)` と、同じ多項式の微分 `dS/dz`。電場は電位の勾配と一致する。
  pure subroutine eval_periodic_shell_mean(plan, state, z, value, slope)
    type(periodic_shell_mean_plan_type), intent(in) :: plan
    type(periodic_shell_mean_state_type), intent(in) :: state
    real(dp), intent(in) :: z
    real(dp), intent(out) :: value, slope
    integer(i32) :: k
    real(dp) :: t, t2, t3, h00, h10, h01, h11, d00, d10, d01, d11

    value = 0.0_dp
    slope = 0.0_dp
    if (.not. plan%ready .or. .not. allocated(state%value)) return
    t = (z - plan%z0)/plan%dz
    k = min(max(int(floor(t), i32), 0_i32), plan%nz - 2_i32)
    t = min(max(t - real(k, dp), 0.0_dp), 1.0_dp)
    k = k + 1_i32
    t2 = t*t
    t3 = t2*t
    h00 = 2.0_dp*t3 - 3.0_dp*t2 + 1.0_dp
    h10 = t3 - 2.0_dp*t2 + t
    h01 = -2.0_dp*t3 + 3.0_dp*t2
    h11 = t3 - t2
    d00 = 6.0_dp*t2 - 6.0_dp*t
    d10 = 3.0_dp*t2 - 4.0_dp*t + 1.0_dp
    d01 = -6.0_dp*t2 + 6.0_dp*t
    d11 = 3.0_dp*t2 - 2.0_dp*t
    value = h00*state%value(k) + h10*plan%dz*state%slope(k) + h01*state%value(k + 1_i32) + &
            h11*plan%dz*state%slope(k + 1_i32)
    slope = (d00*state%value(k) + d01*state%value(k + 1_i32))/plan%dz + d10*state%slope(k) + &
            d11*state%slope(k + 1_i32)
  end subroutine eval_periodic_shell_mean

  subroutine reset_periodic_shell_mean_plan(plan)
    type(periodic_shell_mean_plan_type), intent(inout) :: plan

    if (allocated(plan%value)) deallocate (plan%value)
    if (allocated(plan%slope)) deallocate (plan%slope)
    plan%ready = .false.
    plan%nz = 0_i32
    plan%z0 = 0.0_dp
    plan%dz = 0.0_dp
  end subroutine reset_periodic_shell_mean_plan

  subroutine reset_periodic_shell_mean_state(state)
    type(periodic_shell_mean_state_type), intent(inout) :: state

    if (allocated(state%value)) deallocate (state%value)
    if (allocated(state%slope)) deallocate (state%slope)
  end subroutine reset_periodic_shell_mean_state

  !> 四隅の対数項 x*log(y+r) + y*log(x+r)。負の引数では a+r = (r^2-a^2)/(r-a) で桁落ちを避ける。
  pure real(dp) function corner_log(x, y, h) result(value)
    real(dp), intent(in) :: x, y, h
    real(dp) :: r

    r = sqrt(x*x + y*y + h*h)
    value = x*log_a_plus_r(y, r, x*x + h*h) + y*log_a_plus_r(x, r, y*y + h*h)
  end function corner_log

  pure real(dp) function log_a_plus_r(a, r, rest) result(value)
    real(dp), intent(in) :: a, r, rest

    if (a >= 0.0_dp) then
      value = log(a + r)
    else
      value = log(rest) - log(r - a)
    end if
  end function log_a_plus_r

  !> 四隅の atan(h r/|x y|)。h -> 0 で 0 になり、|h| の折れ目を消した形。
  pure real(dp) function corner_atan(x, y, h) result(value)
    real(dp), intent(in) :: x, y, h

    value = atan(h*sqrt(x*x + y*y + h*h)/abs(x*y))
  end function corner_atan

end module bem_periodic_shell_mean
