!> Zhao A/B/C branchのsolvabilityを独立評価するoffline atlas生成器。
module bem_matching_plane_zhao_atlas
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use bem_kinds, only: dp, i32
  use bem_app_config_types, only: app_config
  use bem_matching_plane_response, only: matching_plane_response_input_count, &
                                         matching_plane_response_output_count
  use bem_matching_plane_response_provider, only: matching_plane_response_provider_type, &
                                                  matching_plane_provider_ok
  use bem_matching_plane_zhao, only: matching_plane_zhao_diagnostics_type, &
                                     matching_plane_zhao_ok, &
                                     matching_plane_zhao_invalid_argument, &
                                     matching_plane_zhao_no_physical_solution, &
                                     matching_plane_zhao_numerical_failure, &
                                     matching_plane_zhao_ambiguous_solution
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_mpi, only: mpi_context
  use bem_string_utils, only: lower_ascii
  use bem_matching_plane_query_io, only: read_matching_plane_query_csv, &
                                         matching_plane_query_ok, matching_plane_query_io_error, &
                                         matching_plane_query_invalid_grid
  implicit none
  private

  integer, parameter :: atlas_input_count = 3
  integer, parameter :: branch_count = 3
  character(len=1), parameter :: branch_names(branch_count) = ['a', 'b', 'c']

  integer(i32), parameter, public :: matching_plane_atlas_ok = matching_plane_query_ok
  integer(i32), parameter, public :: matching_plane_atlas_invalid_argument = 1_i32
  integer(i32), parameter, public :: matching_plane_atlas_io_error = matching_plane_query_io_error
  integer(i32), parameter, public :: matching_plane_atlas_invalid_grid = matching_plane_query_invalid_grid
  integer(i32), parameter, public :: matching_plane_atlas_initialization_failure = 4_i32

  character(len=*), parameter, public :: matching_plane_zhao_atlas_query_csv_header = &
                                         'displacement_c_m2,photoelectron_outward_number_flux_m2_s,'// &
                                         'photoelectron_outward_mean_normal_energy_ev'
  character(len=*), parameter, public :: matching_plane_zhao_atlas_csv_header = &
                                         'displacement_c_m2,photoelectron_outward_number_flux_m2_s,'// &
                                         'photoelectron_outward_mean_normal_energy_ev,branch,status,matching_potential_v,'// &
                                         'electron_inward_number_flux_m2_s,ion_inward_number_flux_m2_s,'// &
                                         'electron_access_potential_v,ion_access_potential_v,'// &
                                         'photoelectron_barrier_potential_v,residual_norm,minimum_field_squared_hat,'// &
                                         'nonlinear_iterations'

  public :: generate_matching_plane_zhao_atlas

contains

  subroutine generate_matching_plane_zhao_atlas(cfg, query_path, output_path, status, message)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: query_path, output_path
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(app_config) :: branch_cfg
    type(matching_plane_response_provider_type) :: providers(branch_count)
    type(matching_plane_zhao_diagnostics_type) :: diagnostics
    type(mpi_context) :: serial_mpi
    real(dp), allocatable :: queries(:, :)
    real(dp) :: input(matching_plane_response_input_count)
    real(dp) :: output(matching_plane_response_output_count)
    integer(i32) :: provider_status, zhao_status
    integer :: branch, row, output_unit, ios, close_ios, rename_status
    character(len=512) :: provider_message, zhao_message
    character(len=:), allocatable :: temporary_output_path

    call accept_atlas(status, message)
    if (trim(lower_ascii(cfg%surface_current%model)) /= 'matching_plane_quasistatic' .or. &
        trim(lower_ascii(cfg%surface_current%response_backend)) /= 'zhao_online') then
      call reject_atlas( &
        matching_plane_atlas_invalid_argument, &
        'Zhao atlas requires matching_plane_quasistatic with response_backend="zhao_online".', &
        status, message &
        )
      return
    end if
    if (len_trim(query_path) == 0 .or. len_trim(output_path) == 0 .or. &
        trim(query_path) == trim(output_path)) then
      call reject_atlas( &
        matching_plane_atlas_invalid_argument, &
        'atlas query and output paths must be nonempty and different.', status, message &
        )
      return
    end if
    temporary_output_path = trim(output_path)//'.beach-zhao-atlas.tmp'
    if (trim(query_path) == temporary_output_path) then
      call reject_atlas( &
        matching_plane_atlas_invalid_argument, &
        'atlas query path conflicts with the temporary output path.', status, message &
        )
      return
    end if

    serial_mpi = mpi_context()
    do branch = 1, branch_count
      branch_cfg = cfg
      branch_cfg%surface_current%zhao_branch = branch_names(branch)
      branch_cfg%surface_current%zhao_root_selection = 'require_unique'
      call providers(branch)%initialize( &
        branch_cfg, serial_mpi, provider_status, provider_message &
        )
      if (provider_status /= matching_plane_provider_ok) then
        call reject_atlas( &
          matching_plane_atlas_initialization_failure, &
          'Zhao-'//branch_names(branch)//' initialization failed: '//trim(provider_message), &
          status, message &
          )
        return
      end if
    end do

    call read_matching_plane_query_csv( &
      query_path, matching_plane_zhao_atlas_query_csv_header, queries, status, message &
      )
    if (status /= matching_plane_atlas_ok) return

    open (newunit=output_unit, file=temporary_output_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      call reject_atlas( &
        matching_plane_atlas_io_error, 'could not open Zhao atlas output: '//trim(output_path), &
        status, message &
        )
      return
    end if
    write (output_unit, '(a)', iostat=ios) matching_plane_zhao_atlas_csv_header
    do row = 1, size(queries, 2)
      if (ios /= 0) exit
      input = [queries(:, row), 0.0_dp, 0.0_dp]
      do branch = 1, branch_count
        call providers(branch)%evaluate_zhao_local( &
          input, output, zhao_status, zhao_message, diagnostics &
          )
        call write_atlas_row( &
          output_unit, queries(:, row), branch_names(branch), zhao_status, output, &
          diagnostics, ios &
          )
        if (ios /= 0) exit
      end do
    end do
    close (output_unit, iostat=close_ios)
    if (ios == 0) ios = close_ios
    if (ios /= 0) then
      call reject_atlas( &
        matching_plane_atlas_io_error, 'failed while writing Zhao atlas output.', status, message &
        )
      return
    end if
    call atomic_rename(temporary_output_path, trim(output_path), rename_status)
    if (rename_status /= filesystem_success) then
      call reject_atlas( &
        matching_plane_atlas_io_error, 'failed to atomically publish Zhao atlas output.', &
        status, message &
        )
    end if
  end subroutine generate_matching_plane_zhao_atlas

  subroutine write_atlas_row(unit_id, query, branch, status, output, diagnostics, ios)
    integer, intent(in) :: unit_id
    real(dp), intent(in) :: query(atlas_input_count)
    character(len=1), intent(in) :: branch
    integer(i32), intent(in) :: status
    real(dp), intent(in) :: output(matching_plane_response_output_count)
    type(matching_plane_zhao_diagnostics_type), intent(in) :: diagnostics
    integer, intent(out) :: ios

    real(dp) :: values(matching_plane_response_output_count + 2)
    integer(i32) :: iterations
    integer :: column

    if (status == matching_plane_zhao_ok) then
      values = [output, diagnostics%residual_norm, diagnostics%minimum_field_squared_hat]
      iterations = diagnostics%nonlinear_iterations
    else
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      iterations = 0_i32
    end if
    write (unit_id, '(es24.16,",",es24.16,",",es24.16,",",a,",",a)', &
           advance='no', iostat=ios) query, upper_ascii(branch), trim(atlas_status_text(status))
    do column = 1, size(values)
      if (ios /= 0) return
      write (unit_id, '(",",es24.16)', advance='no', iostat=ios) values(column)
    end do
    if (ios == 0) write (unit_id, '(",",i0)', iostat=ios) iterations
  end subroutine write_atlas_row

  pure function atlas_status_text(status) result(text)
    integer(i32), intent(in) :: status
    character(len=32) :: text

    select case (status)
    case (matching_plane_zhao_ok)
      text = 'ok'
    case (matching_plane_zhao_invalid_argument)
      text = 'invalid_input'
    case (matching_plane_zhao_no_physical_solution)
      text = 'no_physical_solution'
    case (matching_plane_zhao_numerical_failure)
      text = 'numerical_failure'
    case (matching_plane_zhao_ambiguous_solution)
      text = 'ambiguous_within_branch'
    case default
      text = 'numerical_failure'
    end select
  end function atlas_status_text

  pure function upper_ascii(value) result(upper)
    character(len=1), intent(in) :: value
    character(len=1) :: upper

    upper = value
    if (upper >= 'a' .and. upper <= 'z') upper = achar(iachar(upper) - 32)
  end function upper_ascii

  subroutine accept_atlas(status, message)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = matching_plane_atlas_ok
    message = ''
  end subroutine accept_atlas

  subroutine reject_atlas(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = code
    message = text
  end subroutine reject_atlas

end module bem_matching_plane_zhao_atlas
