!> 粒子SoAデータ構造の初期化を提供するモジュール。
module bem_particles
  use bem_kinds, only: dp, i32
  use bem_types, only: particles_soa
  implicit none
contains

  !> 位置・速度・電荷・質量(と任意重み)配列から `particles_soa` を検証付きで構築する。
  !! @param[out] pcls 検証済み配列を内部に保持した `particles_soa` 構造体。
  !! @param[in] x 粒子位置配列 `x(3,n)` [m]。
  !! @param[in] v 粒子速度配列 `v(3,n)` [m/s]。
  !! @param[in] q 粒子電荷配列 `q(n)` [C]。
  !! @param[in] m 粒子質量配列 `m(n)` [kg]。
  !! @param[in] w マクロ粒子重み配列 `w(n)`（省略時は1）。
  subroutine init_particles(pcls, x, v, q, m, w)
    type(particles_soa), intent(out) :: pcls
    real(dp), intent(in) :: x(:, :), v(:, :), q(:), m(:)
    real(dp), intent(in), optional :: w(:)
    integer(i32) :: n

    n = size(q)
    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) then
      error stop "particle input first dimension must be 3"
    end if
    if (size(x, 2) /= n .or. size(v, 2) /= n .or. size(m) /= n) then
      error stop "particle input size mismatch"
    end if

    pcls%n = n
    allocate (pcls%x(3, n), pcls%v(3, n), pcls%q(n), pcls%m(n), pcls%w(n), pcls%alive(n))
    pcls%x = x
    pcls%v = v
    pcls%q = q
    pcls%m = m
    if (present(w)) then
      if (size(w) /= n) error stop "w size mismatch"
      pcls%w = w
    else
      pcls%w = 1.0d0
    end if
    pcls%alive = .true.
  end subroutine init_particles

end module bem_particles
