!> root の応答テーブルを配信し、全 rank で同じ補間軸と値を使う。
submodule(bem_matching_plane_response) bem_matching_plane_response_mpi
  use bem_mpi, only: mpi_is_root, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  implicit none
contains

  module procedure load_matching_plane_response_mpi
  integer(i32) :: header(7)
  integer(i32), allocatable :: path_codes(:)
  real(dp) :: height(1)
  integer :: axis, component, character_index

  table = matching_plane_response_table_type()
  header = 0_i32
  call accept(status, message)
  if (mpi_is_root(mpi)) then
    call get_matching_plane_response_snapshot(path, table, status, message)
    header(1) = status
    if (status == matching_plane_response_ok) then
      header(2:6) = table%axis_sizes
      header(7) = len(table%source_path)
    end if
  end if
  call mpi_bcast_i32_array(mpi, header, 0_i32)
  status = header(1)
  if (status /= matching_plane_response_ok) then
    if (.not. mpi_is_root(mpi)) message = 'matching-plane response failed to load on MPI root.'
    return
  end if

  if (.not. mpi_is_root(mpi)) then
    table%axis_sizes = header(2:6)
    allocate (table%axes(maxval(table%axis_sizes), matching_plane_response_input_count))
    allocate (table%response_values(matching_plane_response_output_count, product(table%axis_sizes)))
    allocate (character(len=header(7)) :: table%source_path)
  end if
  height = [table%matching_plane_z_m]
  call mpi_bcast_real_dp_array(mpi, height, 0_i32)
  table%matching_plane_z_m = height(1)
  do axis = 1, matching_plane_response_input_count
    call mpi_bcast_real_dp_array(mpi, table%axes(:, axis), 0_i32)
  end do
  do component = 1, matching_plane_response_output_count
    call mpi_bcast_real_dp_array(mpi, table%response_values(component, :), 0_i32)
  end do
  allocate (path_codes(header(7)), source=0_i32)
  if (mpi_is_root(mpi)) then
    do character_index = 1, size(path_codes)
      path_codes(character_index) = iachar(table%source_path(character_index:character_index))
    end do
  end if
  call mpi_bcast_i32_array(mpi, path_codes, 0_i32)
  if (.not. mpi_is_root(mpi)) then
    do character_index = 1, size(path_codes)
      table%source_path(character_index:character_index) = achar(path_codes(character_index))
    end do
  end if
  table%loaded = .true.
  end procedure load_matching_plane_response_mpi

end submodule bem_matching_plane_response_mpi
