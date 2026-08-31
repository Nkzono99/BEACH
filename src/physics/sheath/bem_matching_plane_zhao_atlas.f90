!> Zhao A/B/C branchのsolvabilityを独立評価するoffline atlas生成器。
module bem_matching_plane_zhao_atlas
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
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
  use bem_string_utils, only: lower_ascii, is_decimal_real_token
  implicit none
  private

  integer, parameter :: atlas_input_count = 3
  integer, parameter :: branch_count = 3
  integer, parameter :: line_length = 4096
  integer, parameter :: initial_capacity = 64
  character(len=1), parameter :: branch_names(branch_count) = ['a', 'b', 'c']

  integer(i32), parameter, public :: matching_plane_atlas_ok = 0_i32
  integer(i32), parameter, public :: matching_plane_atlas_invalid_argument = 1_i32
  integer(i32), parameter, public :: matching_plane_atlas_io_error = 2_i32
  integer(i32), parameter, public :: matching_plane_atlas_invalid_grid = 3_i32
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

    call read_atlas_query_grid(query_path, queries, status, message)
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

  subroutine read_atlas_query_grid(path, queries, status, message)
    character(len=*), intent(in) :: path
    real(dp), allocatable, intent(out) :: queries(:, :)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp), allocatable :: buffer(:, :), grown(:, :)
    real(dp) :: query(atlas_input_count)
    character(len=line_length) :: line
    character(len=:), allocatable :: record
    integer :: unit_id, ios, line_number, row_count, capacity
    logical :: header_found

    call accept_atlas(status, message)
    capacity = initial_capacity
    allocate (buffer(atlas_input_count, capacity))
    row_count = 0
    line_number = 0
    header_found = .false.
    open (newunit=unit_id, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      allocate (queries(atlas_input_count, 0))
      call reject_atlas( &
        matching_plane_atlas_io_error, 'could not open Zhao atlas query CSV: '//trim(path), &
        status, message &
        )
      return
    end if

    do
      read (unit_id, '(a)', iostat=ios) line
      if (ios < 0) exit
      line_number = line_number + 1
      if (ios > 0) then
        close (unit_id)
        allocate (queries(atlas_input_count, 0))
        call reject_atlas( &
          matching_plane_atlas_io_error, &
          'failed to read Zhao atlas query CSV line '//trim(integer_text(line_number))//'.', &
          status, message &
          )
        return
      end if
      record = trim(adjustl(line))
      if (len(record) == 0) cycle
      if (record(1:1) == '#') cycle
      if (.not. header_found) then
        if (record /= matching_plane_zhao_atlas_query_csv_header) then
          close (unit_id)
          allocate (queries(atlas_input_count, 0))
          call reject_atlas( &
            matching_plane_atlas_invalid_grid, &
            'Zhao atlas query CSV header does not match the required three-column contract.', &
            status, message &
            )
          return
        end if
        header_found = .true.
        cycle
      end if
      if (count_character(record, ',') /= atlas_input_count - 1) then
        close (unit_id)
        allocate (queries(atlas_input_count, 0))
        call reject_atlas( &
          matching_plane_atlas_invalid_grid, &
          'Zhao atlas query row '//trim(integer_text(line_number))//' must contain exactly three values.', &
          status, message &
          )
        return
      end if
      call parse_query_record(record, query, ios)
      if (ios /= 0 .or. any(.not. ieee_is_finite(query))) then
        close (unit_id)
        allocate (queries(atlas_input_count, 0))
        call reject_atlas( &
          matching_plane_atlas_invalid_grid, &
          'Zhao atlas query row '//trim(integer_text(line_number))//' contains invalid values.', &
          status, message &
          )
        return
      end if
      if (row_count == capacity) then
        allocate (grown(atlas_input_count, 2*capacity))
        grown(:, :capacity) = buffer
        call move_alloc(grown, buffer)
        capacity = 2*capacity
      end if
      row_count = row_count + 1
      buffer(:, row_count) = query
    end do
    close (unit_id)

    if (.not. header_found .or. row_count == 0) then
      allocate (queries(atlas_input_count, 0))
      call reject_atlas( &
        matching_plane_atlas_invalid_grid, &
        'Zhao atlas query CSV must contain the required header and at least one row.', &
        status, message &
        )
      return
    end if
    allocate (queries(atlas_input_count, row_count), source=buffer(:, :row_count))
  end subroutine read_atlas_query_grid

  subroutine parse_query_record(record, query, ios)
    character(len=*), intent(in) :: record
    real(dp), intent(out) :: query(atlas_input_count)
    integer, intent(out) :: ios

    character(len=:), allocatable :: token
    integer :: column, comma, first, last

    query = 0.0_dp
    ios = 0
    first = 1
    do column = 1, atlas_input_count
      if (column < atlas_input_count) then
        comma = index(record(first:), ',')
        if (comma <= 0) then
          ios = 1
          return
        end if
        last = first + comma - 2
      else
        if (index(record(first:), ',') /= 0) then
          ios = 1
          return
        end if
        last = len(record)
      end if
      if (last < first) then
        ios = 1
        return
      end if
      token = trim(adjustl(record(first:last)))
      if (len(token) == 0 .or. scan(token, ' '//achar(9)) > 0 .or. &
          .not. is_decimal_real_token(token)) then
        ios = 1
        return
      end if
      read (token, *, iostat=ios) query(column)
      if (ios /= 0) return
      first = last + 2
    end do
  end subroutine parse_query_record

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

  pure integer function count_character(text, needle) result(count)
    character(len=*), intent(in) :: text
    character(len=1), intent(in) :: needle
    integer :: position

    count = 0
    do position = 1, len(text)
      if (text(position:position) == needle) count = count + 1
    end do
  end function count_character

  pure function integer_text(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write (text, '(i0)') value
  end function integer_text

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
