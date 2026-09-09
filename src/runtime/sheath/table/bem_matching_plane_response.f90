!> Matching-plane outer-sheath response table の immutable snapshot と補間。
module bem_matching_plane_response
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_mpi, only: mpi_context
  use bem_matching_plane_contract, only: &
    matching_plane_response_input_count, &
    matching_plane_response_output_count, &
    matching_plane_input_displacement, &
    matching_plane_input_photoelectron_outward_flux, &
    matching_plane_input_photoelectron_mean_normal_energy, &
    matching_plane_input_electron_outward_flux, &
    matching_plane_input_ion_outward_flux, &
    matching_plane_output_matching_potential, &
    matching_plane_output_electron_inward_flux, &
    matching_plane_output_ion_inward_flux, &
    matching_plane_output_electron_access_potential, &
    matching_plane_output_ion_access_potential, &
    matching_plane_output_photoelectron_barrier_potential
  implicit none
  private

  ! The table and analytic model share the same response coordinates.
  public :: matching_plane_response_input_count
  public :: matching_plane_response_output_count
  public :: matching_plane_input_displacement
  public :: matching_plane_input_photoelectron_outward_flux
  public :: matching_plane_input_photoelectron_mean_normal_energy
  public :: matching_plane_input_electron_outward_flux
  public :: matching_plane_input_ion_outward_flux
  public :: matching_plane_output_matching_potential
  public :: matching_plane_output_electron_inward_flux
  public :: matching_plane_output_ion_inward_flux
  public :: matching_plane_output_electron_access_potential
  public :: matching_plane_output_ion_access_potential
  public :: matching_plane_output_photoelectron_barrier_potential

  integer(i32), parameter, public :: matching_plane_response_ok = 0_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_argument = 1_i32
  integer(i32), parameter, public :: matching_plane_response_io_error = 2_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_header = 3_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_metadata = 4_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_row = 5_i32
  integer(i32), parameter, public :: matching_plane_response_invalid_grid = 6_i32
  integer(i32), parameter, public :: matching_plane_response_out_of_range = 7_i32

  character(len=*), parameter, public :: matching_plane_response_query_csv_header = &
                                         'displacement_c_m2,'// &
                                         'photoelectron_outward_number_flux_m2_s,'// &
                                         'photoelectron_outward_mean_normal_energy_ev,'// &
                                         'electron_outward_number_flux_m2_s,'// &
                                         'ion_outward_number_flux_m2_s'

  character(len=*), parameter, public :: matching_plane_response_csv_header = &
                                         matching_plane_response_query_csv_header//','// &
                                         'matching_potential_v,'// &
                                         'electron_inward_number_flux_m2_s,'// &
                                         'ion_inward_number_flux_m2_s,'// &
                                         'electron_access_potential_v,'// &
                                         'ion_access_potential_v,'// &
                                         'photoelectron_barrier_potential_v'

  character(len=52), parameter :: input_names(matching_plane_response_input_count) = &
                                  [character(len=52) :: &
                                   'displacement_c_m2', &
                                   'photoelectron_outward_number_flux_m2_s', &
                                   'photoelectron_outward_mean_normal_energy_ev', &
                                   'electron_outward_number_flux_m2_s', &
                                   'ion_outward_number_flux_m2_s']

  !> 検証済み response table。components は外部から変更できない。
  type, public :: matching_plane_response_table_type
    private
    logical :: loaded = .false.
    character(len=:), allocatable :: source_path
    real(dp) :: matching_plane_z_m = 0.0_dp
    integer(i32) :: axis_sizes(matching_plane_response_input_count) = 0_i32
    real(dp), allocatable :: axes(:, :)
    real(dp), allocatable :: response_values(:, :)
  contains
    procedure, public :: evaluate => evaluate_matching_plane_response
    procedure, public :: get_axis_data => get_matching_plane_axis_data
    procedure, public :: get_matching_plane_z => get_matching_plane_z
    procedure, public :: get_source_path => get_matching_plane_source_path
    procedure, public :: is_loaded => matching_plane_response_is_loaded
  end type matching_plane_response_table_type

  type :: matching_plane_response_snapshot_entry
    character(len=:), allocatable :: path
    type(matching_plane_response_table_type) :: table
  end type matching_plane_response_snapshot_entry

  type(matching_plane_response_snapshot_entry), allocatable, save :: response_snapshots(:)

  public :: get_matching_plane_response_snapshot
  public :: load_matching_plane_response_mpi
  public :: reset_matching_plane_response_snapshot_cache

  interface
    !> root が読み込んだテーブルを全 rank に配信する。
    module subroutine load_matching_plane_response_mpi(path, mpi, table, status, message)
      character(len=*), intent(in) :: path
      type(mpi_context), intent(in) :: mpi
      type(matching_plane_response_table_type), intent(out) :: table
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine load_matching_plane_response_mpi

    module subroutine read_matching_plane_response_csv(path, table, status, message)
      character(len=*), intent(in) :: path
      type(matching_plane_response_table_type), intent(out) :: table
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine read_matching_plane_response_csv

    pure module subroutine accept(status, message)
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine accept

    pure module subroutine reject(code, text, status, message)
      integer(i32), intent(in) :: code
      character(len=*), intent(in) :: text
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine reject
  end interface

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

  !> フィードバックの範囲と尺度を決めるために補間軸を返す。
  subroutine get_matching_plane_axis_data( &
    self, axis_sizes, axis_values, matching_plane_z_m, status, message &
    )
    class(matching_plane_response_table_type), intent(in) :: self
    integer(i32), allocatable, intent(out) :: axis_sizes(:)
    real(dp), allocatable, intent(out) :: axis_values(:)
    real(dp), intent(out) :: matching_plane_z_m
    integer(i32), intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    integer :: axis, first, last, total_axis_values

    matching_plane_z_m = 0.0_dp
    if (.not. self%loaded) then
      allocate (axis_sizes(0), axis_values(0))
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
    matching_plane_z_m = self%matching_plane_z_m
    call assign_optional_status(matching_plane_response_ok, '', status, message)
  end subroutine get_matching_plane_axis_data

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

  pure subroutine assign_optional_status(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    if (present(status)) status = code
    if (present(message)) message = text
  end subroutine assign_optional_status

end module bem_matching_plane_response
