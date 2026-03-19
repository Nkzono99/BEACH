!> 境界要素に蓄積した電荷から観測点の電場を評価する場計算モジュール。
module bem_field
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type
  implicit none
contains

  !> 全要素電荷を点電荷近似で総和し、softening付きで観測点 `r` の電場ベクトルを返す。
  !! @param[in] mesh 要素重心 `centers` と要素電荷 `q_elem` を保持したメッシュ情報。
  !! @param[in] r 電場を評価する観測点座標 `(x,y,z)` [m]。
  !! @param[in] softening 特異点回避のために距離2乗へ加える softening 長さ [m]。
  !! @param[out] e 観測点 `r` における電場ベクトル `(Ex,Ey,Ez)` [V/m]。
  subroutine electric_field_at(mesh, r, softening, e)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: softening
    real(dp), intent(out) :: e(3)

    integer(i32) :: i
    real(dp) :: soft2, r2, inv_r3, ex, ey, ez
    real(dp) :: rx, ry, rz, dx, dy, dz, qi

    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    soft2 = softening*softening
    rx = r(1)
    ry = r(2)
    rz = r(3)

    !$omp simd reduction(+:ex,ey,ez) private(dx,dy,dz,r2,inv_r3,qi)
    do i = 1, mesh%nelem
      dx = rx - mesh%center_x(i)
      dy = ry - mesh%center_y(i)
      dz = rz - mesh%center_z(i)
      r2 = dx*dx + dy*dy + dz*dz + soft2
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      qi = mesh%q_elem(i)
      ex = ex + qi*inv_r3*dx
      ey = ey + qi*inv_r3*dy
      ez = ez + qi*inv_r3*dz
    end do

    e(1) = k_coulomb*ex
    e(2) = k_coulomb*ey
    e(3) = k_coulomb*ez
  end subroutine electric_field_at

  !> 全要素電荷を点電荷近似で総和し、softening付きで観測点 `r` の電位を返す。
  !! @param[in] mesh 要素重心 `centers` と要素電荷 `q_elem` を保持したメッシュ情報。
  !! @param[in] r 電位を評価する観測点座標 `(x,y,z)` [m]。
  !! @param[in] softening 特異点回避のために距離2乗へ加える softening 長さ [m]。
  !! @param[out] phi 観測点 `r` における電位 [V]。
  subroutine electric_potential_at(mesh, r, softening, phi)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: softening
    real(dp), intent(out) :: phi

    integer(i32) :: i
    real(dp) :: soft2, r2, inv_r, phi_sum
    real(dp) :: rx, ry, rz, dx, dy, dz

    phi_sum = 0.0d0
    soft2 = softening*softening
    rx = r(1)
    ry = r(2)
    rz = r(3)

    !$omp simd reduction(+:phi_sum) private(dx,dy,dz,r2,inv_r)
    do i = 1, mesh%nelem
      dx = rx - mesh%center_x(i)
      dy = ry - mesh%center_y(i)
      dz = rz - mesh%center_z(i)
      r2 = dx*dx + dy*dy + dz*dz + soft2
      inv_r = 1.0d0/sqrt(r2)
      phi_sum = phi_sum + mesh%q_elem(i)*inv_r
    end do

    phi = k_coulomb*phi_sum
  end subroutine electric_potential_at

end module bem_field
