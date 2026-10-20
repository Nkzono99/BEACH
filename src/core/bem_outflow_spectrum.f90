!> z-high を外向きに通過した電子の法線運動エネルギー [eV] を数える固定 bin。
!! sheath-model の既定の対数格子（[0, 1 meV) と 1 meV から 1 keV までの 24 bin/decade、計 145 bin）をそのまま使い、
!! 外部シースへ渡す bin 放出源と、checkpoint に保存する外部源の両方がこの bin を使う。1 keV 以上は最後の bin に数える。
module bem_outflow_spectrum
  use bem_kinds, only: dp, i32
  use sheath_model, only: log_energy_edges, log_energy_bin, LOG_ENERGY_DEFAULT_BINS, LOG_ENERGY_DEFAULT_BINS_PER_DECADE
  implicit none
  private

  integer(i32), parameter, public :: outflow_spectrum_bins_per_decade = int(LOG_ENERGY_DEFAULT_BINS_PER_DECADE, i32)
  integer(i32), parameter, public :: outflow_spectrum_bins = int(LOG_ENERGY_DEFAULT_BINS, i32)

  public :: outflow_spectrum_edges_ev, outflow_spectrum_bin

contains

  !> bin 境界 [eV]。要素数は bin 数 + 1。
  pure function outflow_spectrum_edges_ev() result(edges)
    real(dp) :: edges(outflow_spectrum_bins + 1)

    edges = log_energy_edges()
  end function outflow_spectrum_edges_ev

  !> 法線運動エネルギーが属する bin。範囲外と NaN は端の bin に数える。
  pure integer(i32) function outflow_spectrum_bin(energy_ev) result(bin)
    real(dp), intent(in) :: energy_ev

    bin = int(log_energy_bin(energy_ev), i32)
  end function outflow_spectrum_bin

end module bem_outflow_spectrum
