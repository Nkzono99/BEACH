!> シミュレーションで使用する物理定数を定義する。
module bem_constants
  use bem_kinds, only: dp
  implicit none
  real(dp), parameter :: k_coulomb = 8.9875517923d9
  real(dp), parameter :: k_boltzmann = 1.380649d-23
  real(dp), parameter :: pi = 3.1415926535897932384626433832795d0
  real(dp), parameter :: eps0 = 8.8541878128d-12
  real(dp), parameter :: qe = 1.602176634d-19
end module bem_constants
