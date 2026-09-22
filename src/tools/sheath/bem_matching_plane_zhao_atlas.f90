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
  use bem_source_kinetic_sheath, only: source_kinetic_root, source_kinetic_excluded, source_kinetic_numerical_failure
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

    if (cfg%surface_current%density_model == 'source_kinetic') then
      call generate_source_atlas(cfg, query_path, output_path, temporary_output_path, status, message)
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

  !> Source-connected mode retains all detected roots before applying a selector.
  subroutine generate_source_atlas(cfg, query_path, output_path, temporary_path, status, message)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: query_path, output_path, temporary_path
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(app_config) :: query_cfg
    type(matching_plane_response_provider_type) :: provider
    type(matching_plane_zhao_diagnostics_type) :: diagnostics
    type(source_kinetic_root) :: root
    type(mpi_context) :: serial_mpi
    real(dp), allocatable :: queries(:, :)
    real(dp) :: output(matching_plane_response_output_count)
    integer :: unit_id, ios, close_ios, row, i, count, rejected_count, rename_status
    integer(i32) :: query_status
    character(len=512) :: query_message
    character(len=24) :: classification

    query_cfg = cfg
    query_cfg%surface_current%zhao_branch = 'auto'
    query_cfg%surface_current%zhao_root_selection = 'require_unique'
    serial_mpi = mpi_context()
    call provider%initialize(query_cfg, serial_mpi, status, message)
    if (status /= matching_plane_provider_ok) return
    call read_matching_plane_query_csv(query_path, matching_plane_zhao_atlas_query_csv_header, queries, status, message)
    if (status /= matching_plane_atlas_ok) return
    open (newunit=unit_id, file=temporary_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      call reject_atlas(matching_plane_atlas_io_error, 'Could not open source_kinetic atlas.', status, message)
      return
    end if
    write (unit_id, '(a)', iostat=ios) &
      'displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,'// &
      'classification,root_count,root_index,accepted,branch,phi_h_te,phi_min_te,electron_amplitude,'// &
      'escaping_flux_hat,electron_inward_flux_hat,ion_inward_flux_hat,current_hat,neutrality_residual,'// &
      'field_squared_residual,outer_residual,minimum_field_squared,minimum_curvature,edge_coefficient,'// &
      'ion_turning_margin,barrier_relative_roundoff,depth_min,depth_max,search_points,search_complete,'// &
      'boundary_turning,deep_root,barrier_resolution_limited,numerical_failures,rejection'
    do row = 1, size(queries, 2)
      if (ios /= 0) exit
      call provider%evaluate_zhao_local([queries(:, row), 0._dp, 0._dp], output, query_status, query_message, diagnostics)
      count = 0
      rejected_count = 0
      if (allocated(diagnostics%kinetic%roots)) count = size(diagnostics%kinetic%roots)
      if (allocated(diagnostics%kinetic%rejected)) rejected_count = size(diagnostics%kinetic%rejected)
      classification = 'search_unresolved'
      if (count > 0) classification = 'solutions'
      if (diagnostics%kinetic%status == source_kinetic_excluded) classification = 'analytically_excluded'
      if (diagnostics%kinetic%status == source_kinetic_numerical_failure) classification = 'numerical_exception'
      if (query_status == matching_plane_zhao_invalid_argument) classification = 'invalid_input'
      do i = 1, max(1, count + rejected_count)
        root = source_kinetic_root()
        if (i <= count) then
          root = diagnostics%kinetic%roots(i)
        else if (i <= count + rejected_count) then
          root = diagnostics%kinetic%rejected(i - count)
        else
          root%phi_h = ieee_value(0._dp, ieee_quiet_nan)
          root%phi_min = root%phi_h
          root%rejection = trim(query_message)
        end if
        call write_source_row(unit_id, queries(:, row), classification, count, i, root, diagnostics, ios)
        if (ios /= 0) exit
      end do
    end do
    close (unit_id, iostat=close_ios)
    if (ios == 0) ios = close_ios
    if (ios /= 0) then
      call reject_atlas(matching_plane_atlas_io_error, 'Failed writing source_kinetic atlas.', status, message)
      return
    end if
    call atomic_rename(temporary_path, output_path, rename_status)
    if (rename_status /= filesystem_success) then
      call reject_atlas(matching_plane_atlas_io_error, 'Failed publishing source_kinetic atlas.', status, message)
    end if
  end subroutine

  subroutine write_source_row(unit_id, query, classification, count, index, root, diagnostics, ios)
    integer, intent(in) :: unit_id, count, index
    real(dp), intent(in) :: query(3)
    character(len=*), intent(in) :: classification
    type(source_kinetic_root), intent(in) :: root
    type(matching_plane_zhao_diagnostics_type), intent(in) :: diagnostics
    integer, intent(out) :: ios
    write (unit_id, '(3(es24.16,","),a,2(",",i0),",",l1,",",a)', advance='no', iostat=ios) &
      query, trim(classification), count, index, root%accepted, trim(root%branch)
    if (ios /= 0) return
    write (unit_id, '(17(",",es24.16))', advance='no', iostat=ios) &
      root%phi_h, root%phi_min, root%amplitude, root%escaping_flux, root%electron_flux, root%ion_flux, root%current, &
      root%neutrality_residual, root%field_squared_residual, root%outer_residual, root%minimum_field_squared, &
      root%minimum_curvature, root%edge_coefficient, root%ion_turning_margin, root%barrier_roundoff, &
      diagnostics%kinetic%depth_min, diagnostics%kinetic%depth_max
    if (ios /= 0) return
    write (unit_id, '(",",i0,4(",",l1),",",i0,",",a)', iostat=ios) &
      diagnostics%kinetic%search_points, diagnostics%kinetic%search_complete, root%boundary_turning, &
      -min(root%phi_h, root%phi_min) > 1.e4_dp, root%barrier_roundoff > 1.e-7_dp, &
      diagnostics%kinetic%numerical_failures, '"'//trim(root%rejection)//'"'
  end subroutine

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
