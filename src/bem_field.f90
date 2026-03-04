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
    use omp_lib
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: softening
    real(dp), intent(out) :: e(3)

    integer(i32) :: i
    real(dp) :: dr(3), r2, inv_r3, ex, ey, ez

    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0

    !$omp parallel do default(none) private(i,dr,r2,inv_r3) reduction(+:ex,ey,ez) shared(mesh,r,softening)
    do i = 1, mesh%nelem
      dr = r - mesh%centers(:, i)
      r2 = sum(dr * dr) + softening * softening
      inv_r3 = 1.0d0 / (sqrt(r2) * r2)
      ex = ex + mesh%q_elem(i) * inv_r3 * dr(1)
      ey = ey + mesh%q_elem(i) * inv_r3 * dr(2)
      ez = ez + mesh%q_elem(i) * inv_r3 * dr(3)
    end do
    !$omp end parallel do

    e(1) = k_coulomb * ex
    e(2) = k_coulomb * ey
    e(3) = k_coulomb * ez
  end subroutine electric_field_at

end module bem_field
