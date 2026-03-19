!> 荷電粒子の時間発展にBoris法を適用する運動方程式ソルバ。
module bem_pusher
  use bem_kinds, only: dp
  implicit none
contains

  !> 電場半ステップ加速と磁場回転を組み合わせ、1タイムステップ後の位置・速度を計算する。
  !! @param[in] x 現在時刻の粒子位置ベクトル `(x,y,z)` [m]。
  !! @param[in] v 現在時刻の粒子速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] q 粒子1個あたりの電荷 [C]。
  !! @param[in] m 粒子1個あたりの質量 [kg]。
  !! @param[in] dt 時間刻み幅 [s]。
  !! @param[in] e 粒子位置での電場ベクトル [V/m]。
  !! @param[in] b 粒子位置での磁束密度ベクトル [T]。
  !! @param[out] x_new 1ステップ後の更新粒子位置 [m]。
  !! @param[out] v_new 1ステップ後の更新粒子速度 [m/s]。
  subroutine boris_push(x, v, q, m, dt, e, b, x_new, v_new)
    real(dp), intent(in) :: x(3), v(3), q, m, dt, e(3), b(3)
    real(dp), intent(out) :: x_new(3), v_new(3)
    real(dp) :: qm, v_minus(3), t(3), s(3), v_prime(3), v_plus(3), t2

    qm = q/m
    v_minus = v + qm*e*(0.5d0*dt)
    t = qm*b*(0.5d0*dt)
    t2 = sum(t*t)
    s = 2.0d0*t/(1.0d0 + t2)
    v_prime = v_minus + cross(v_minus, t)
    v_plus = v_minus + cross(v_prime, s)
    v_new = v_plus + qm*e*(0.5d0*dt)
    x_new = x + v_new*dt
  end subroutine boris_push

  !> 3次元ベクトルの外積を返す基本演算。
  !! @param[in] a 左オペランドの3次元ベクトル。
  !! @param[in] b 右オペランドの3次元ベクトル。
  !! @return c 関数の戻り値。
  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

end module bem_pusher
