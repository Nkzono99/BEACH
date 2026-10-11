!> z-high を外向きに通過した光電子の法線運動エネルギー [eV] を数える固定 bin。
!! 外部シースへ渡す bin 光電子源と、checkpoint に保存する外部源の両方がこの bin を使う。
!! 最初の bin は [0, 1 meV)、以降は 1 meV から 1 keV までの対数 bin で、1 keV 以上は最後の bin に数える。
module bem_outflow_spectrum
  use bem_kinds, only: dp, i32
  implicit none
  private

  integer(i32), parameter, public :: outflow_spectrum_bins_per_decade = 24_i32
  real(dp), parameter :: lowest_edge_ev = 1.0e-3_dp
  integer(i32), parameter :: log_decades = 6_i32
  integer(i32), parameter, public :: outflow_spectrum_bins = 1_i32 + log_decades*outflow_spectrum_bins_per_decade

  public :: outflow_spectrum_edges_ev, outflow_spectrum_bin

contains

  !> bin 境界 [eV]。要素数は bin 数 + 1。
  pure function outflow_spectrum_edges_ev() result(edges)
    real(dp) :: edges(outflow_spectrum_bins + 1)
    integer(i32) :: j

    edges(1) = 0.0_dp
    do j = 0_i32, outflow_spectrum_bins - 1_i32
      edges(j + 2_i32) = lowest_edge_ev*10.0_dp**(real(j, dp)/real(outflow_spectrum_bins_per_decade, dp))
    end do
  end function outflow_spectrum_edges_ev

  !> 法線運動エネルギーが属する bin。範囲外は端の bin に数える。
  pure integer(i32) function outflow_spectrum_bin(energy_ev) result(bin)
    real(dp), intent(in) :: energy_ev

    if (.not. energy_ev >= lowest_edge_ev) then
      bin = 1_i32
      return
    end if
    bin = 2_i32 + int(floor(log10(energy_ev/lowest_edge_ev)*real(outflow_spectrum_bins_per_decade, dp)), i32)
    bin = min(max(bin, 2_i32), outflow_spectrum_bins)
  end function outflow_spectrum_bin

end module bem_outflow_spectrum
