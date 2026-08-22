!> Matching-plane outer-sheath response table の immutable snapshot と補間。
module bem_matching_plane_response
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_allreduce_min_i32_scalar, mpi_allreduce_max_i32_scalar, &
                     mpi_bcast_i32_array
  use bem_string_utils, only: is_decimal_real_token
  implicit none
  private

  integer(i32), parameter, public :: matching_plane_response_input_count = 5_i32
  integer(i32), parameter, public :: matching_plane_response_output_count = 6_i32

  integer(i32), parameter, public :: matching_plane_response_ok = 0_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_argument = 1_i32
  integer(i32), parameter, public :: matching_plane_response_io_error = 2_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_header = 3_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_metadata = 4_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_row = 5_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_grid = 6_i32
  integer(i32), parameter, public :: matching_plane_response_out_of_range = 7_i32

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

  character(len=*), parameter, public :: matching_plane_response_csv_header = &
                                         'displacement_c_m2,'// &
                                         'photoelectron_outward_number_flux_m2_s,'// &
                                         'photoelectron_outward_mean_normal_energy_ev,'// &
                                         'electron_outward_number_flux_m2_s,'// &
                                         'ion_outward_number_flux_m2_s,'// &
                                         'matching_potential_v,'// &
                                         'electron_inward_number_flux_m2_s,'// &
                                         'ion_inward_number_flux_m2_s,'// &
                                         'electron_access_potential_v,'// &
                                         'ion_access_potential_v,'// &
                                         'photoelectron_barrier_potential_v'

  integer, parameter :: response_column_count = 11
  integer, parameter :: response_line_length = 4096
  integer, parameter :: initial_row_capacity = 64
  integer(i64), parameter :: content_hash_modulus = 2147483647_i64
  integer(i64), parameter :: content_hash_multiplier_a = 65599_i64
  integer(i64), parameter :: content_hash_multiplier_b = 131071_i64
  character(len=*), parameter :: matching_plane_z_prefix = '# matching_plane_z_m='
  character(len=52), parameter :: input_names(matching_plane_response_input_count) = [character(len=52) :: &
                                                                                      'displacement_c_m2', &
                                                                                      'photoelectron_outward_number_flux_m2_s', &
                                                                                    'photoelectron_outward_mean_normal_energy_ev', &
                                                                                      'electron_outward_number_flux_m2_s', &
                                                                                      'ion_outward_number_flux_m2_s']

  type :: matching_plane_axis_type
    real(dp), allocatable :: values(:)
  end type matching_plane_axis_type

  !> 検証済み response table。components は外部から変更できない。
  type, public :: matching_plane_response_table_type
    private
    logical :: loaded = .false.
    character(len=:), allocatable :: source_path
    real(dp) :: matching_plane_z_m = 0.0_dp
    integer(i32) :: axis_sizes(matching_plane_response_input_count) = 0_i32
    real(dp), allocatable :: axes(:, :)
    real(dp), allocatable :: response_values(:, :)
    character(len=16) :: content_fingerprint = ''
  contains
    procedure, public :: evaluate => evaluate_matching_plane_response
    procedure, public :: get_fingerprint_data => get_matching_plane_fingerprint_data
    procedure, public :: get_matching_plane_z => get_matching_plane_z
    procedure, public :: get_source_path => get_matching_plane_source_path
    procedure, public :: is_loaded => matching_plane_response_is_loaded
  end type matching_plane_response_table_type

  type :: matching_plane_response_snapshot_entry
    character(len=:), allocatable :: path
    type(matching_plane_response_table_type) :: table
  end type matching_plane_response_snapshot_entry

  type :: matching_plane_content_hash_state_type
    integer(i64) :: a = 146959810_i64
    integer(i64) :: b = 109951162_i64
  end type matching_plane_content_hash_state_type

  type(matching_plane_response_snapshot_entry), allocatable, save :: response_snapshots(:)

  public :: get_matching_plane_response_snapshot
  public :: get_matching_plane_response_content_fingerprint
  public :: preflight_matching_plane_response_mpi
  public :: reset_matching_plane_response_snapshot_cache

contains

  !> path ごとに最初の検証済み内容を固定し、同一 process 内で共有する。
  subroutine get_matching_plane_response_snapshot(path, table, status, message)
    character(len=*), intent(in) :: path
    type(matching_plane_response_table_type), intent(out) :: table
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer :: entry

    table = matching_plane_response_table_type()
    call ensure_matching_plane_response_snapshot(path, entry, status, message)
    if (status /= matching_plane_response_ok) return
    table = response_snapshots(entry)%table
  end subroutine get_matching_plane_response_snapshot

  !> canonical table内容のcache済みfingerprintを配列copyなしで返す。
  subroutine get_matching_plane_response_content_fingerprint(path, fingerprint, status, message)
    character(len=*), intent(in) :: path
    character(len=16), intent(out) :: fingerprint
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer :: entry

    fingerprint = ''
    call ensure_matching_plane_response_snapshot(path, entry, status, message)
    if (status /= matching_plane_response_ok) return
    fingerprint = response_snapshots(entry)%table%content_fingerprint
  end subroutine get_matching_plane_response_content_fingerprint

  !> active flag・load成否・canonical contentを全MPI rankで照合する。
  subroutine preflight_matching_plane_response_mpi(active, path, mpi, fingerprint, status, message, table)
    logical, intent(in) :: active
    character(len=*), intent(in) :: path
    type(mpi_context), intent(in) :: mpi
    character(len=16), intent(out) :: fingerprint
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(matching_plane_response_table_type), intent(out), optional :: table

    integer :: entry, character_index
    integer(i32) :: active_min, active_max, load_ok, mismatch
    integer(i32) :: local_codes(16), root_codes(16)

    fingerprint = ''
    if (present(table)) table = matching_plane_response_table_type()
    call accept(status, message)

    active_min = merge(1_i32, 0_i32, active)
    active_max = active_min
    call mpi_allreduce_min_i32_scalar(mpi, active_min)
    call mpi_allreduce_max_i32_scalar(mpi, active_max)
    if (active_min /= active_max) then
      call reject( &
        matching_plane_response_invalid_argument, &
        'matching-plane activation differs across MPI ranks.', status, message &
        )
      return
    end if
    if (active_max == 0_i32) return

    call ensure_matching_plane_response_snapshot(path, entry, status, message)
    load_ok = merge(1_i32, 0_i32, status == matching_plane_response_ok)
    call mpi_allreduce_min_i32_scalar(mpi, load_ok)
    if (load_ok == 0_i32) then
      if (status == matching_plane_response_ok) then
        call reject( &
          matching_plane_response_io_error, &
          'matching-plane response failed to load on another MPI rank.', status, message &
          )
      end if
      return
    end if

    fingerprint = response_snapshots(entry)%table%content_fingerprint
    do character_index = 1, len(fingerprint)
      local_codes(character_index) = int(iachar(fingerprint(character_index:character_index)), i32)
    end do
    root_codes = 0_i32
    if (mpi_is_root(mpi)) root_codes = local_codes
    call mpi_bcast_i32_array(mpi, root_codes, 0_i32)
    mismatch = merge(1_i32, 0_i32, any(local_codes /= root_codes))
    call mpi_allreduce_max_i32_scalar(mpi, mismatch)
    if (mismatch /= 0_i32) then
      fingerprint = ''
      call reject( &
        matching_plane_response_invalid_grid, &
        'matching-plane response content differs across MPI ranks.', status, message &
        )
      return
    end if
    if (present(table)) table = response_snapshots(entry)%table
  end subroutine preflight_matching_plane_response_mpi

  subroutine ensure_matching_plane_response_snapshot(path, entry_index, status, message)
    character(len=*), intent(in) :: path
    integer, intent(out) :: entry_index
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(matching_plane_response_snapshot_entry), allocatable :: grown(:)
    type(matching_plane_response_table_type) :: loaded_table
    character(len=:), allocatable :: normalized_path
    integer :: entry, old_size

    entry_index = 0
    call accept(status, message)
    normalized_path = trim(path)
    if (len(normalized_path) == 0) then
      call reject( &
        matching_plane_response_invalid_argument, 'matching-plane response path must not be empty.', status, message &
        )
      return
    end if

    if (allocated(response_snapshots)) then
      do entry = 1, size(response_snapshots)
        if (response_snapshots(entry)%path == normalized_path) then
          entry_index = entry
          return
        end if
      end do
      old_size = size(response_snapshots)
    else
      old_size = 0
    end if

    call read_matching_plane_response_csv(normalized_path, loaded_table, status, message)
    if (status /= matching_plane_response_ok) return

    allocate (grown(old_size + 1))
    do entry = 1, old_size
      call move_alloc(response_snapshots(entry)%path, grown(entry)%path)
      call move_matching_plane_response_table(response_snapshots(entry)%table, grown(entry)%table)
    end do
    grown(old_size + 1)%path = normalized_path
    call move_matching_plane_response_table(loaded_table, grown(old_size + 1)%table)
    call move_alloc(grown, response_snapshots)
    entry_index = old_size + 1
  end subroutine ensure_matching_plane_response_snapshot

  !> allocatable componentsをcopyせずcache entryへ移す。
  subroutine move_matching_plane_response_table(source, destination)
    type(matching_plane_response_table_type), intent(inout) :: source
    type(matching_plane_response_table_type), intent(out) :: destination

    destination = matching_plane_response_table_type()
    destination%loaded = source%loaded
    destination%matching_plane_z_m = source%matching_plane_z_m
    destination%axis_sizes = source%axis_sizes
    destination%content_fingerprint = source%content_fingerprint
    if (allocated(source%source_path)) call move_alloc(source%source_path, destination%source_path)
    if (allocated(source%axes)) call move_alloc(source%axes, destination%axes)
    if (allocated(source%response_values)) call move_alloc(source%response_values, destination%response_values)
    source = matching_plane_response_table_type()
  end subroutine move_matching_plane_response_table

  !> 独立 run や unit test の開始時に immutable snapshot cache を解放する。
  subroutine reset_matching_plane_response_snapshot_cache()
    if (allocated(response_snapshots)) deallocate (response_snapshots)
  end subroutine reset_matching_plane_response_snapshot_cache

  !> 5入力から6出力を多重線形補間する。singleton axis は inactive dimension として扱う。
  subroutine evaluate_matching_plane_response(self, input, output, status, message)
    class(matching_plane_response_table_type), intent(in) :: self
    real(dp), intent(in) :: input(matching_plane_response_input_count)
    real(dp), intent(out) :: output(matching_plane_response_output_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer(i32) :: lower(matching_plane_response_input_count)
    integer(i32) :: upper(matching_plane_response_input_count)
    integer(i32) :: index_at_corner(matching_plane_response_input_count)
    real(dp) :: upper_weight(matching_plane_response_input_count)
    real(dp) :: query_value, range_tolerance
    real(dp) :: corner_weight
    integer(i64) :: linear_index, stride
    integer :: axis, corner

    output = 0.0_dp
    call accept(status, message)
    if (.not. self%loaded) then
      call reject( &
        matching_plane_response_invalid_argument, 'matching-plane response table is not loaded.', status, message &
        )
      return
    end if
    if (any(.not. ieee_is_finite(input))) then
      call reject( &
        matching_plane_response_invalid_argument, 'matching-plane response query must be finite.', status, message &
        )
      return
    end if

    do axis = 1, matching_plane_response_input_count
      if (self%axis_sizes(axis) == 1_i32) then
        ! A singleton explicitly removes this dimension from the response.  Its
        ! coordinate documents the generation point but does not restrict queries.
        lower(axis) = 1_i32
        upper(axis) = 1_i32
        upper_weight(axis) = 0.0_dp
        cycle
      end if
      query_value = input(axis)
      range_tolerance = 64.0_dp*epsilon(1.0_dp)*max( &
                        1.0_dp, abs(self%axes(1, axis)), &
                        abs(self%axes(self%axis_sizes(axis), axis)) &
                        )
      if (query_value < self%axes(1, axis) - range_tolerance .or. &
          query_value > self%axes(self%axis_sizes(axis), axis) + range_tolerance) then
        call reject_out_of_range( &
          axis, query_value, self%axes(1, axis), self%axes(self%axis_sizes(axis), axis), status, message &
          )
        return
      end if
      ! 丸め幅だけ端点を越えた値は端点そのものとして補間する。この幅を
      ! 越える外挿は上で拒否し、物理的な範囲外入力をclampしない。
      if (query_value < self%axes(1, axis)) query_value = self%axes(1, axis)
      if (query_value > self%axes(self%axis_sizes(axis), axis)) then
        query_value = self%axes(self%axis_sizes(axis), axis)
      end if
      call locate_axis_interval( &
        self%axes(:self%axis_sizes(axis), axis), query_value, lower(axis), upper(axis), upper_weight(axis) &
        )
    end do

    do corner = 0, 2**matching_plane_response_input_count - 1
      corner_weight = 1.0_dp
      do axis = 1, matching_plane_response_input_count
        if (btest(corner, axis - 1)) then
          if (self%axis_sizes(axis) == 1_i32) then
            corner_weight = 0.0_dp
            exit
          end if
          index_at_corner(axis) = upper(axis)
          corner_weight = corner_weight*upper_weight(axis)
        else
          index_at_corner(axis) = lower(axis)
          corner_weight = corner_weight*(1.0_dp - upper_weight(axis))
        end if
      end do
      if (corner_weight == 0.0_dp) cycle

      linear_index = 1_i64
      stride = 1_i64
      do axis = 1, matching_plane_response_input_count
        linear_index = linear_index + int(index_at_corner(axis) - 1_i32, i64)*stride
        stride = stride*int(self%axis_sizes(axis), i64)
      end do
      output = output + corner_weight*self%response_values(:, int(linear_index))
    end do

    if (any(.not. ieee_is_finite(output))) then
      output = 0.0_dp
      call reject( &
        matching_plane_response_invalid_grid, &
        'matching-plane response interpolation produced a non-finite output.', status, message &
        )
    end if
  end subroutine evaluate_matching_plane_response

  !> Fingerprint 用にcanonical axesとaxis-1-fastest出力配列のcopyを返す。
  subroutine get_matching_plane_fingerprint_data( &
    self, axis_sizes, axis_values, response_values, matching_plane_z_m, status, message &
    )
    class(matching_plane_response_table_type), intent(in) :: self
    integer(i32), allocatable, intent(out) :: axis_sizes(:)
    real(dp), allocatable, intent(out) :: axis_values(:)
    real(dp), allocatable, intent(out), optional :: response_values(:, :)
    real(dp), intent(out) :: matching_plane_z_m
    integer(i32), intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    integer :: axis, first, last, total_axis_values

    matching_plane_z_m = 0.0_dp
    if (.not. self%loaded) then
      allocate (axis_sizes(0), axis_values(0))
      if (present(response_values)) allocate (response_values(0, 0))
      call assign_optional_status( &
        matching_plane_response_invalid_argument, 'matching-plane response table is not loaded.', status, message &
        )
      return
    end if

    allocate (axis_sizes(matching_plane_response_input_count), source=self%axis_sizes)
    total_axis_values = sum(self%axis_sizes)
    allocate (axis_values(total_axis_values))
    first = 1
    do axis = 1, matching_plane_response_input_count
      last = first + self%axis_sizes(axis) - 1
      axis_values(first:last) = self%axes(:self%axis_sizes(axis), axis)
      first = last + 1
    end do
    if (present(response_values)) allocate (response_values, source=self%response_values)
    matching_plane_z_m = self%matching_plane_z_m
    call assign_optional_status(matching_plane_response_ok, '', status, message)
  end subroutine get_matching_plane_fingerprint_data

  !> CSV metadata のmatching-plane高さを返す。
  subroutine get_matching_plane_z(self, matching_plane_z_m, status, message)
    class(matching_plane_response_table_type), intent(in) :: self
    real(dp), intent(out) :: matching_plane_z_m
    integer(i32), intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    matching_plane_z_m = 0.0_dp
    if (.not. self%loaded) then
      call assign_optional_status( &
        matching_plane_response_invalid_argument, 'matching-plane response table is not loaded.', status, message &
        )
      return
    end if
    matching_plane_z_m = self%matching_plane_z_m
    call assign_optional_status(matching_plane_response_ok, '', status, message)
  end subroutine get_matching_plane_z

  !> snapshot のcache keyとなったpathを返す。
  subroutine get_matching_plane_source_path(self, path, status, message)
    class(matching_plane_response_table_type), intent(in) :: self
    character(len=:), allocatable, intent(out) :: path
    integer(i32), intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    if (.not. self%loaded) then
      path = ''
      call assign_optional_status( &
        matching_plane_response_invalid_argument, 'matching-plane response table is not loaded.', status, message &
        )
      return
    end if
    path = self%source_path
    call assign_optional_status(matching_plane_response_ok, '', status, message)
  end subroutine get_matching_plane_source_path

  pure logical function matching_plane_response_is_loaded(self) result(loaded)
    class(matching_plane_response_table_type), intent(in) :: self

    loaded = self%loaded
  end function matching_plane_response_is_loaded

  subroutine read_matching_plane_response_csv(path, table, status, message)
    character(len=*), intent(in) :: path
    type(matching_plane_response_table_type), intent(out) :: table
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp), allocatable :: rows(:, :)
    real(dp) :: parsed_row(response_column_count), matching_plane_z_m
    character(len=response_line_length) :: line
    character(len=:), allocatable :: left_adjusted
    integer :: unit_id, ios, line_number, row_count, capacity
    logical :: header_found, matching_plane_z_found

    table = matching_plane_response_table_type()
    call accept(status, message)
    capacity = initial_row_capacity
    allocate (rows(response_column_count, capacity))
    row_count = 0
    line_number = 0
    header_found = .false.
    matching_plane_z_found = .false.
    matching_plane_z_m = 0.0_dp

    open (newunit=unit_id, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      call reject( &
        matching_plane_response_io_error, 'could not open matching-plane response CSV: '//trim(path), status, message &
        )
      return
    end if

    do
      read (unit_id, '(a)', iostat=ios) line
      if (ios < 0) exit
      line_number = line_number + 1
      if (ios > 0) then
        close (unit_id)
        call reject_line( &
          matching_plane_response_io_error, line_number, 'failed to read matching-plane response CSV.', status, message &
          )
        return
      end if
      if (len_trim(line) == 0) cycle
      left_adjusted = adjustl(trim(line))

      if (left_adjusted(1:1) == '#') then
        if (index(left_adjusted, '# matching_plane_z_m') == 1) then
          if (header_found) then
            close (unit_id)
            call reject_line( &
              matching_plane_response_invalid_metadata, line_number, &
              'matching_plane_z_m metadata must precede the CSV header.', status, message &
              )
            return
          end if
          if (matching_plane_z_found) then
            close (unit_id)
            call reject_line( &
              matching_plane_response_invalid_metadata, line_number, &
              'matching_plane_z_m metadata must appear exactly once.', status, message &
              )
            return
          end if
          call parse_matching_plane_z(left_adjusted, matching_plane_z_m, status, message)
          if (status /= matching_plane_response_ok) then
            close (unit_id)
            call prefix_line_number(line_number, message)
            return
          end if
          matching_plane_z_found = .true.
        end if
        cycle
      end if

      if (.not. header_found) then
        if (trim(line) /= matching_plane_response_csv_header) then
          close (unit_id)
          call reject_line( &
            matching_plane_response_invalid_header, line_number, &
            'matching-plane response CSV header does not match the exact v1 contract.', status, message &
            )
          return
        end if
        if (.not. matching_plane_z_found) then
          close (unit_id)
          call reject_line( &
            matching_plane_response_invalid_metadata, line_number, &
            'matching_plane_z_m metadata is required before the CSV header.', status, message &
            )
          return
        end if
        header_found = .true.
        cycle
      end if

      call parse_numeric_csv_row(line, parsed_row, status, message)
      if (status /= matching_plane_response_ok) then
        close (unit_id)
        call prefix_line_number(line_number, message)
        return
      end if
      if (row_count >= capacity) call grow_row_buffer(rows, capacity)
      row_count = row_count + 1
      rows(:, row_count) = parsed_row
    end do
    close (unit_id)

    if (.not. header_found) then
      call reject( &
        matching_plane_response_invalid_header, 'matching-plane response CSV contains no v1 header.', status, message &
        )
      return
    end if
    if (row_count <= 0) then
      call reject( &
        matching_plane_response_invalid_grid, 'matching-plane response CSV contains no numeric rows.', status, message &
        )
      return
    end if

    call build_matching_plane_response_table( &
      path, matching_plane_z_m, rows(:, :row_count), table, status, message &
      )
  end subroutine read_matching_plane_response_csv

  subroutine build_matching_plane_response_table(path, matching_plane_z_m, rows, table, status, message)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: matching_plane_z_m
    real(dp), intent(in) :: rows(:, :)
    type(matching_plane_response_table_type), intent(out) :: table
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(matching_plane_axis_type) :: axis_data(matching_plane_response_input_count)
    logical, allocatable :: filled(:)
    integer(i32) :: coordinate_index(matching_plane_response_input_count)
    integer(i64) :: expected_points, linear_index, stride
    integer :: axis, row, max_axis_size, point_count

    table = matching_plane_response_table_type()
    call accept(status, message)
    if (size(rows, 1) /= response_column_count .or. size(rows, 2) <= 0) then
      call reject( &
        matching_plane_response_invalid_grid, 'matching-plane response rows have an invalid shape.', status, message &
        )
      return
    end if
    if (.not. ieee_is_finite(matching_plane_z_m)) then
      call reject( &
        matching_plane_response_invalid_metadata, 'matching_plane_z_m metadata must be finite.', status, message &
        )
      return
    end if
    if (any(.not. ieee_is_finite(rows))) then
      call reject( &
        matching_plane_response_invalid_row, 'matching-plane response values must all be finite.', status, message &
        )
      return
    end if
    if (any(rows(2:5, :) < 0.0_dp) .or. any(rows(7:8, :) < 0.0_dp)) then
      call reject( &
        matching_plane_response_invalid_row, &
        'matching-plane response fluxes and normal energies must be non-negative.', status, message &
        )
      return
    end if

    expected_points = 1_i64
    do axis = 1, matching_plane_response_input_count
      call unique_sorted_values(rows(axis, :), axis_data(axis)%values)
      table%axis_sizes(axis) = int(size(axis_data(axis)%values), i32)
      if (axis >= matching_plane_input_photoelectron_outward_flux .and. &
          .not. any(axis_data(axis)%values == 0.0_dp)) then
        call reject( &
          matching_plane_response_invalid_grid, &
          'every matching-plane feedback axis must include zero.', status, message &
          )
        return
      end if
      if (expected_points > huge(expected_points)/int(table%axis_sizes(axis), i64)) then
        call reject( &
          matching_plane_response_invalid_grid, 'matching-plane Cartesian grid size overflow.', status, message &
          )
        return
      end if
      expected_points = expected_points*int(table%axis_sizes(axis), i64)
    end do
    if (expected_points > int(huge(0), i64)) then
      call reject( &
        matching_plane_response_invalid_grid, 'matching-plane Cartesian grid is too large.', status, message &
        )
      return
    end if
    if (expected_points /= int(size(rows, 2), i64)) then
      call reject( &
        matching_plane_response_invalid_grid, &
        'matching-plane response rows do not form one complete Cartesian product.', status, message &
        )
      return
    end if

    max_axis_size = maxval(table%axis_sizes)
    point_count = int(expected_points)
    allocate (table%axes(max_axis_size, matching_plane_response_input_count))
    table%axes = 0.0_dp
    do axis = 1, matching_plane_response_input_count
      table%axes(:table%axis_sizes(axis), axis) = axis_data(axis)%values
    end do
    allocate (table%response_values(matching_plane_response_output_count, point_count))
    allocate (filled(point_count), source=.false.)
    table%response_values = 0.0_dp

    do row = 1, size(rows, 2)
      do axis = 1, matching_plane_response_input_count
        coordinate_index(axis) = exact_axis_index(axis_data(axis)%values, rows(axis, row))
        if (coordinate_index(axis) <= 0_i32) then
          call reject( &
            matching_plane_response_invalid_grid, &
            'matching-plane row could not be mapped to its canonical axis.', status, message &
            )
          return
        end if
      end do
      linear_index = 1_i64
      stride = 1_i64
      do axis = 1, matching_plane_response_input_count
        linear_index = linear_index + int(coordinate_index(axis) - 1_i32, i64)*stride
        stride = stride*int(table%axis_sizes(axis), i64)
      end do
      if (filled(int(linear_index))) then
        call reject( &
          matching_plane_response_invalid_grid, &
          'matching-plane response grid contains a duplicate input coordinate.', status, message &
          )
        return
      end if
      filled(int(linear_index)) = .true.
      table%response_values(:, int(linear_index)) = rows(6:11, row)
    end do
    if (size(rows, 2) /= point_count .or. .not. all(filled)) then
      call reject( &
        matching_plane_response_invalid_grid, &
        'matching-plane response rows must form one complete Cartesian product.', status, message &
        )
      return
    end if

    table%source_path = trim(path)
    table%matching_plane_z_m = matching_plane_z_m
    table%content_fingerprint = compute_matching_plane_content_fingerprint(table)
    table%loaded = .true.
  end subroutine build_matching_plane_response_table

  subroutine parse_matching_plane_z(line, matching_plane_z_m, status, message)
    character(len=*), intent(in) :: line
    real(dp), intent(out) :: matching_plane_z_m
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=:), allocatable :: value_text
    integer :: ios

    matching_plane_z_m = 0.0_dp
    call accept(status, message)
    if (index(line, matching_plane_z_prefix) /= 1) then
      call reject( &
        matching_plane_response_invalid_metadata, &
        'matching_plane_z_m metadata must use "# matching_plane_z_m=<finite>".', status, message &
        )
      return
    end if
    value_text = trim(adjustl(line(len(matching_plane_z_prefix) + 1:)))
    if (len(value_text) == 0 .or. scan(value_text, ' ,'//achar(9)) > 0 .or. &
        .not. is_decimal_real_token(value_text)) then
      call reject( &
        matching_plane_response_invalid_metadata, 'matching_plane_z_m metadata has an invalid value.', status, message &
        )
      return
    end if
    read (value_text, *, iostat=ios) matching_plane_z_m
    if (ios /= 0 .or. .not. ieee_is_finite(matching_plane_z_m)) then
      matching_plane_z_m = 0.0_dp
      call reject( &
        matching_plane_response_invalid_metadata, 'matching_plane_z_m metadata must be finite.', status, message &
        )
    end if
  end subroutine parse_matching_plane_z

  subroutine parse_numeric_csv_row(line, values, status, message)
    character(len=*), intent(in) :: line
    real(dp), intent(out) :: values(response_column_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=:), allocatable :: record, token
    integer :: column, comma, first, last, ios

    values = 0.0_dp
    call accept(status, message)
    record = trim(line)
    first = 1
    do column = 1, response_column_count
      if (column < response_column_count) then
        comma = index(record(first:), ',')
        if (comma == 0) then
          call reject( &
            matching_plane_response_invalid_row, &
            'matching-plane response row must contain exactly 11 comma-separated values.', status, message &
            )
          return
        end if
        last = first + comma - 2
      else
        if (index(record(first:), ',') /= 0) then
          call reject( &
            matching_plane_response_invalid_row, &
            'matching-plane response row must contain exactly 11 comma-separated values.', status, message &
            )
          return
        end if
        last = len(record)
      end if
      if (last < first) then
        call reject( &
          matching_plane_response_invalid_row, 'matching-plane response CSV values must not be empty.', status, message &
          )
        return
      end if
      token = trim(adjustl(record(first:last)))
      if (len(token) == 0 .or. scan(token, ' '//achar(9)) > 0 .or. &
          .not. is_decimal_real_token(token)) then
        call reject( &
          matching_plane_response_invalid_row, 'matching-plane response CSV contains an invalid numeric token.', &
          status, message &
          )
        return
      end if
      read (token, *, iostat=ios) values(column)
      if (ios /= 0) then
        call reject( &
          matching_plane_response_invalid_row, 'matching-plane response CSV contains an invalid numeric token.', &
          status, message &
          )
        return
      end if
      first = last + 2
    end do
    if (any(.not. ieee_is_finite(values))) then
      call reject( &
        matching_plane_response_invalid_row, 'matching-plane response values must all be finite.', status, message &
        )
    end if
  end subroutine parse_numeric_csv_row

  subroutine grow_row_buffer(rows, capacity)
    real(dp), allocatable, intent(inout) :: rows(:, :)
    integer, intent(inout) :: capacity

    real(dp), allocatable :: grown(:, :)
    integer :: new_capacity

    if (capacity > huge(capacity)/2) error stop 'matching-plane response row capacity overflow.'
    new_capacity = 2*capacity
    allocate (grown(response_column_count, new_capacity))
    grown(:, :capacity) = rows
    call move_alloc(grown, rows)
    capacity = new_capacity
  end subroutine grow_row_buffer

  subroutine unique_sorted_values(values, unique)
    real(dp), intent(in) :: values(:)
    real(dp), allocatable, intent(out) :: unique(:)

    real(dp), allocatable :: sorted(:), compact(:)
    integer :: i, count

    allocate (sorted, source=values)
    call merge_sort_real(sorted)
    allocate (compact(size(sorted)))
    count = 1
    compact(1) = sorted(1)
    do i = 2, size(sorted)
      if (sorted(i) /= compact(count)) then
        count = count + 1
        compact(count) = sorted(i)
      end if
    end do
    allocate (unique(count), source=compact(:count))
  end subroutine unique_sorted_values

  subroutine merge_sort_real(values)
    real(dp), intent(inout) :: values(:)

    real(dp), allocatable :: work(:)
    integer :: width, left, middle, right, i, j, k, n

    n = size(values)
    if (n <= 1) return
    allocate (work(n))
    width = 1
    do while (width < n)
      left = 1
      do while (left <= n)
        middle = min(left + width, n + 1)
        right = min(left + 2*width - 1, n)
        i = left
        j = middle
        do k = left, right
          if (i >= middle) then
            work(k) = values(j)
            j = j + 1
          else if (j > right) then
            work(k) = values(i)
            i = i + 1
          else if (values(i) <= values(j)) then
            work(k) = values(i)
            i = i + 1
          else
            work(k) = values(j)
            j = j + 1
          end if
        end do
        left = left + 2*width
      end do
      values = work
      if (width > n/2) exit
      width = 2*width
    end do
  end subroutine merge_sort_real

  pure integer(i32) function exact_axis_index(axis_values, value) result(index_value)
    real(dp), intent(in) :: axis_values(:), value

    integer :: low, high, middle

    index_value = 0_i32
    low = 1
    high = size(axis_values)
    do while (low <= high)
      middle = (low + high)/2
      if (value == axis_values(middle)) then
        index_value = int(middle, i32)
        return
      else if (value < axis_values(middle)) then
        high = middle - 1
      else
        low = middle + 1
      end if
    end do
  end function exact_axis_index

  pure subroutine locate_axis_interval(axis_values, value, lower, upper, upper_weight)
    real(dp), intent(in) :: axis_values(:), value
    integer(i32), intent(out) :: lower, upper
    real(dp), intent(out) :: upper_weight

    integer :: low_index, high_index, middle

    if (value == axis_values(size(axis_values))) then
      lower = int(size(axis_values) - 1, i32)
      upper = int(size(axis_values), i32)
      upper_weight = 1.0_dp
      return
    end if

    low_index = 1
    high_index = size(axis_values) - 1
    do while (low_index <= high_index)
      middle = (low_index + high_index)/2
      if (value < axis_values(middle)) then
        high_index = middle - 1
      else if (value >= axis_values(middle + 1)) then
        low_index = middle + 1
      else
        lower = int(middle, i32)
        upper = int(middle + 1, i32)
        upper_weight = (value - axis_values(middle))/(axis_values(middle + 1) - axis_values(middle))
        return
      end if
    end do

    ! evaluate側のclosed-range検証後なので、ここへ到達するのは丸め境界だけ。
    lower = 1_i32
    upper = 2_i32
    upper_weight = 0.0_dp
  end subroutine locate_axis_interval

  subroutine reject_out_of_range(axis, value, lower, upper, status, message)
    integer, intent(in) :: axis
    real(dp), intent(in) :: value, lower, upper
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=512) :: text

    write (text, '(a,a,a,es16.8,a,es16.8,a,es16.8,a)') &
      'matching-plane query ', trim(input_names(axis)), '=', value, ' is outside [', lower, ',', upper, '].'
    call reject(matching_plane_response_out_of_range, trim(text), status, message)
  end subroutine reject_out_of_range

  subroutine reject_line(code, line_number, text, status, message)
    integer(i32), intent(in) :: code
    integer, intent(in) :: line_number
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=512) :: contextual

    write (contextual, '(a,i0,a,a)') 'line ', line_number, ': ', trim(text)
    call reject(code, trim(contextual), status, message)
  end subroutine reject_line

  subroutine prefix_line_number(line_number, message)
    integer, intent(in) :: line_number
    character(len=*), intent(inout) :: message

    character(len=512) :: contextual

    write (contextual, '(a,i0,a,a)') 'line ', line_number, ': ', trim(message)
    message = trim(contextual)
  end subroutine prefix_line_number

  function compute_matching_plane_content_fingerprint(table) result(fingerprint)
    type(matching_plane_response_table_type), intent(in) :: table
    character(len=16) :: fingerprint
    type(matching_plane_content_hash_state_type) :: hash
    integer :: axis, point

    call content_hash_feed_string(hash, 'matching_plane_response_csv_v1')
    call content_hash_feed_real(hash, table%matching_plane_z_m)
    call content_hash_feed_integer(hash, matching_plane_response_input_count)
    do axis = 1, matching_plane_response_input_count
      call content_hash_feed_integer(hash, table%axis_sizes(axis))
      call content_hash_feed_real_vector(hash, table%axes(:table%axis_sizes(axis), axis))
    end do
    call content_hash_feed_integer(hash, matching_plane_response_output_count)
    call content_hash_feed_integer(hash, int(size(table%response_values, 2), i32))
    do point = 1, size(table%response_values, 2)
      call content_hash_feed_real_vector(hash, table%response_values(:, point))
    end do
    write (fingerprint, '(z8.8,z8.8)') hash%a, hash%b
  end function compute_matching_plane_content_fingerprint

  subroutine content_hash_feed_string(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    character(len=*), intent(in) :: value
    integer :: index

    call content_hash_feed_byte(hash, len_trim(value))
    do index = 1, len_trim(value)
      call content_hash_feed_byte(hash, iachar(value(index:index)))
    end do
  end subroutine content_hash_feed_string

  subroutine content_hash_feed_integer(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    integer(i32), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(i0)') value
    call content_hash_feed_string(hash, trim(encoded))
  end subroutine content_hash_feed_integer

  subroutine content_hash_feed_real(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    real(dp), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(es24.16e3)') value
    call content_hash_feed_string(hash, trim(adjustl(encoded)))
  end subroutine content_hash_feed_real

  subroutine content_hash_feed_real_vector(hash, values)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    real(dp), intent(in) :: values(:)
    integer :: index

    call content_hash_feed_integer(hash, int(size(values), i32))
    do index = 1, size(values)
      call content_hash_feed_real(hash, values(index))
    end do
  end subroutine content_hash_feed_real_vector

  subroutine content_hash_feed_byte(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    integer, intent(in) :: value

    hash%a = modulo(hash%a*content_hash_multiplier_a + int(value, i64) + 1_i64, content_hash_modulus)
    hash%b = modulo(hash%b*content_hash_multiplier_b + int(value, i64) + 1_i64, content_hash_modulus)
  end subroutine content_hash_feed_byte

  pure subroutine accept(status, message)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = matching_plane_response_ok
    message = ''
  end subroutine accept

  pure subroutine reject(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = code
    message = text
  end subroutine reject

  pure subroutine assign_optional_status(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    if (present(status)) status = code
    if (present(message)) message = text
  end subroutine assign_optional_status

end module bem_matching_plane_response
