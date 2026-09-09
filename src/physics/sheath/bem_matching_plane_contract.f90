!> Matching-plane 応答の 5 入力・6 出力の並び。モデル、CSV、MPI の実装には依存しない。
module bem_matching_plane_contract
  use bem_kinds, only: i32
  implicit none
  private

  integer(i32), parameter, public :: matching_plane_response_input_count = 5_i32
  integer(i32), parameter, public :: matching_plane_response_output_count = 6_i32

  integer(i32), parameter, public :: matching_plane_input_displacement = 1_i32
  integer(i32), parameter, public :: matching_plane_input_photoelectron_outward_flux = 2_i32
  integer(i32), parameter, public :: matching_plane_input_photoelectron_mean_normal_energy = 3_i32
  integer(i32), parameter, public :: matching_plane_input_electron_outward_flux = 4_i32
  integer(i32), parameter, public :: matching_plane_input_ion_outward_flux = 5_i32

  integer(i32), parameter, public :: matching_plane_output_matching_potential = 1_i32
  integer(i32), parameter, public :: matching_plane_output_electron_inward_flux = 2_i32
  integer(i32), parameter, public :: matching_plane_output_ion_inward_flux = 3_i32
  integer(i32), parameter, public :: matching_plane_output_electron_access_potential = 4_i32
  integer(i32), parameter, public :: matching_plane_output_ion_access_potential = 5_i32
  integer(i32), parameter, public :: matching_plane_output_photoelectron_barrier_potential = 6_i32

end module bem_matching_plane_contract
