!> summary.txt の外部シース（matching-plane）receipt のキー表。
!! summary の writer と checkpoint の loader が同じ表からキー名・順序・必須条件を取る。
module bem_matching_plane_summary
  use bem_kinds, only: dp, i32
  use bem_types, only: matching_plane_stats
  use bem_string_utils, only: is_decimal_real_token, is_decimal_integer_token, is_logical_token
  implicit none
  private

  public :: matching_plane_summary_key_count
  public :: write_matching_plane_summary
  public :: read_matching_plane_summary_key
  public :: require_matching_plane_summary_keys
  public :: read_matching_summary_real, read_matching_summary_integer, read_matching_summary_logical

  integer, parameter :: value_real = 1, value_integer = 2, value_logical = 3, value_spectrum = 4

  !> summary の 1 キー。unique は重複を拒否し、required は schema 9 以降の checkpoint で必須、
  !! restored は再開時に読み戻す（top_cell_mean は snapshot から再計算するので読まない）。
  type :: receipt_key
    character(len=56) :: name = ''
    integer :: value_kind = value_real
    logical :: unique = .false.
    logical :: required = .false.
    logical :: restored = .true.
  end type receipt_key

  ! summary.txt に書く順。順序と名前は出力形式の一部なので変えない。
  type(receipt_key), parameter :: receipt_keys(*) = [ &
                                  receipt_key('matching_plane_secondary_share', value_real, .false., .false., .true.), &
                                  receipt_key('secondary_budget_residual_a_m2', value_real, .false., .false., .true.), &
                                  receipt_key('matching_plane_state_valid', value_logical, .true., .true., .true.), &
                                  receipt_key('matching_plane_displacement_C_m2', value_real, .true., .true., .true.), &
                                  receipt_key('matching_plane_phi_V', value_real, .true., .true., .true.), &
                                  receipt_key('matching_plane_top_cell_mean_V_m', value_real, .false., .false., .false.), &
                                  receipt_key('matching_plane_electron_inward_flux_m2_s', value_real, .true., .true., .true.), &
                                  receipt_key('matching_plane_ion_inward_flux_m2_s', value_real, .true., .true., .true.), &
                                  receipt_key('matching_plane_electron_access_potential_V', value_real, .true., .true., .true.), &
                                  receipt_key('matching_plane_ion_access_potential_V', value_real, .true., .true., .true.), &
                                  receipt_key( &
                                  'matching_plane_photoelectron_barrier_potential_V', value_real, .true., .true., .true.), &
                                  receipt_key( &
                                  'matching_plane_photoelectron_outward_flux_m2_s', value_real, .true., .true., .true.), &
                                  receipt_key( &
                                  'matching_plane_photoelectron_mean_normal_energy_eV', value_real, .true., .true., .true.), &
                                  receipt_key( &
                                  'matching_plane_photoelectron_spectrum_m2_s', value_spectrum, .true., .false., .true.), &
                                  receipt_key('matching_plane_electron_outward_flux_m2_s', value_real, .true., .true., .true.), &
                                  receipt_key('matching_plane_ion_outward_flux_m2_s', value_real, .true., .true., .true.), &
                                  receipt_key( &
                                  'matching_plane_photoelectron_return_flux_m2_s', value_real, .true., .true., .true.), &
                                  receipt_key( &
                                  'matching_plane_photoelectron_escape_flux_m2_s', value_real, .true., .true., .true.), &
                                  receipt_key('matching_plane_iterations', value_integer, .true., .true., .true.), &
                                  receipt_key('matching_plane_residual', value_real, .true., .true., .true.) &
                                  ]
  integer, parameter :: matching_plane_summary_key_count = size(receipt_keys)

contains

  !> receipt を表の順に summary.txt へ書く。
  subroutine write_matching_plane_summary(u, receipt)
    integer, intent(in) :: u
    type(matching_plane_stats), intent(in) :: receipt
    integer :: k

    do k = 1, size(receipt_keys)
      select case (receipt_keys(k)%value_kind)
      case (value_logical)
        write (u, '(a,l1)') trim(receipt_keys(k)%name)//'=', receipt%state_valid
      case (value_integer)
        write (u, '(a,i0)') trim(receipt_keys(k)%name)//'=', receipt%iterations
      case (value_spectrum)
        write (u, '(a,*(1x,es24.16))') trim(receipt_keys(k)%name)//'=', receipt%photoelectron_spectrum_m2_s
      case default
        write (u, '(a,es24.16)') trim(receipt_keys(k)%name)//'=', receipt_real(receipt, k)
      end select
    end do
  end subroutine write_matching_plane_summary

  !> key が receipt のキーなら value を読み戻して handled=.true.。found は表と同じ添字で、重複の検出に使う。
  subroutine read_matching_plane_summary_key(key, value, receipt, found, handled)
    character(len=*), intent(in) :: key, value
    type(matching_plane_stats), intent(inout) :: receipt
    logical, intent(inout) :: found(matching_plane_summary_key_count)
    logical, intent(out) :: handled
    integer :: k, ios
    real(dp) :: x

    handled = .false.
    k = key_index(key)
    if (k == 0) return
    if (.not. receipt_keys(k)%restored) return
    handled = .true.
    if (receipt_keys(k)%unique) then
      if (found(k)) error stop 'Resume checkpoint summary contains duplicate matching-plane key: '//trim(key)
    end if
    found(k) = .true.
    select case (receipt_keys(k)%value_kind)
    case (value_logical)
      call read_matching_summary_logical(value, key, receipt%state_valid)
    case (value_integer)
      call read_matching_summary_integer(value, key, receipt%iterations)
    case (value_spectrum)
      read (value, *, iostat=ios) receipt%photoelectron_spectrum_m2_s
      if (ios /= 0) error stop 'Resume checkpoint matching-plane photoelectron spectrum could not be parsed.'
    case default
      call read_matching_summary_real(value, key, x)
      call set_receipt_real(receipt, k, x)
    end select
  end subroutine read_matching_plane_summary_key

  !> schema 9 以降の checkpoint では required のキーが全て揃っていなければ再開しない。
  subroutine require_matching_plane_summary_keys(found)
    logical, intent(in) :: found(matching_plane_summary_key_count)

    if (.not. all(found .or. .not. receipt_keys%required)) then
      error stop 'Resume checkpoint schema-v9 summary is missing required matching-plane keys.'
    end if
  end subroutine require_matching_plane_summary_keys

  integer function key_index(key) result(k)
    character(len=*), intent(in) :: key

    do k = 1, size(receipt_keys)
      if (trim(receipt_keys(k)%name) == trim(key)) return
    end do
    k = 0
  end function key_index

  !> 実数キー k の receipt の値。
  real(dp) function receipt_real(receipt, k) result(x)
    type(matching_plane_stats), intent(in) :: receipt
    integer, intent(in) :: k

    select case (trim(receipt_keys(k)%name))
    case ('matching_plane_secondary_share')
      x = receipt%secondary_share
    case ('secondary_budget_residual_a_m2')
      x = receipt%secondary_budget_residual_a_m2
    case ('matching_plane_displacement_C_m2')
      x = receipt%displacement_c_m2
    case ('matching_plane_phi_V')
      x = receipt%phi_v
    case ('matching_plane_top_cell_mean_V_m')
      x = receipt%top_cell_mean_v_m
    case ('matching_plane_electron_inward_flux_m2_s')
      x = receipt%response(2)
    case ('matching_plane_ion_inward_flux_m2_s')
      x = receipt%response(3)
    case ('matching_plane_electron_access_potential_V')
      x = receipt%response(4)
    case ('matching_plane_ion_access_potential_V')
      x = receipt%response(5)
    case ('matching_plane_photoelectron_barrier_potential_V')
      x = receipt%response(6)
    case ('matching_plane_photoelectron_outward_flux_m2_s')
      x = receipt%feedback(1)
    case ('matching_plane_photoelectron_mean_normal_energy_eV')
      x = receipt%feedback(2)
    case ('matching_plane_electron_outward_flux_m2_s')
      x = receipt%feedback(3)
    case ('matching_plane_ion_outward_flux_m2_s')
      x = receipt%feedback(4)
    case ('matching_plane_photoelectron_return_flux_m2_s')
      x = receipt%photoelectron_return_flux_m2_s
    case ('matching_plane_photoelectron_escape_flux_m2_s')
      x = receipt%photoelectron_escape_flux_m2_s
    case ('matching_plane_residual')
      x = receipt%residual
    case default
      error stop 'matching-plane summary key has no real value.'
    end select
  end function receipt_real

  !> 実数キー k へ値を戻す。phi_V は応答の電位（response(1)）も兼ねる。
  subroutine set_receipt_real(receipt, k, x)
    type(matching_plane_stats), intent(inout) :: receipt
    integer, intent(in) :: k
    real(dp), intent(in) :: x

    select case (trim(receipt_keys(k)%name))
    case ('matching_plane_secondary_share')
      receipt%secondary_share = x
    case ('secondary_budget_residual_a_m2')
      receipt%secondary_budget_residual_a_m2 = x
    case ('matching_plane_displacement_C_m2')
      receipt%displacement_c_m2 = x
    case ('matching_plane_phi_V')
      receipt%phi_v = x
      receipt%response(1) = x
    case ('matching_plane_electron_inward_flux_m2_s')
      receipt%response(2) = x
    case ('matching_plane_ion_inward_flux_m2_s')
      receipt%response(3) = x
    case ('matching_plane_electron_access_potential_V')
      receipt%response(4) = x
    case ('matching_plane_ion_access_potential_V')
      receipt%response(5) = x
    case ('matching_plane_photoelectron_barrier_potential_V')
      receipt%response(6) = x
    case ('matching_plane_photoelectron_outward_flux_m2_s')
      receipt%feedback(1) = x
    case ('matching_plane_photoelectron_mean_normal_energy_eV')
      receipt%feedback(2) = x
    case ('matching_plane_electron_outward_flux_m2_s')
      receipt%feedback(3) = x
    case ('matching_plane_ion_outward_flux_m2_s')
      receipt%feedback(4) = x
    case ('matching_plane_photoelectron_return_flux_m2_s')
      receipt%photoelectron_return_flux_m2_s = x
    case ('matching_plane_photoelectron_escape_flux_m2_s')
      receipt%photoelectron_escape_flux_m2_s = x
    case ('matching_plane_residual')
      receipt%residual = x
    case default
      error stop 'matching-plane summary key has no real value.'
    end select
  end subroutine set_receipt_real

  subroutine read_matching_summary_real(value, key, result_value)
    character(len=*), intent(in) :: value, key
    real(dp), intent(out) :: result_value
    integer :: ios

    if (.not. is_decimal_real_token(trim(value))) then
      error stop 'Resume checkpoint matching-plane key has an invalid real token: '//trim(key)
    end if
    read (value, *, iostat=ios) result_value
    if (ios /= 0) error stop 'Resume checkpoint matching-plane real value could not be parsed: '//trim(key)
  end subroutine read_matching_summary_real

  subroutine read_matching_summary_integer(value, key, result_value)
    character(len=*), intent(in) :: value, key
    integer(i32), intent(out) :: result_value
    integer :: ios

    if (.not. is_decimal_integer_token(trim(value))) then
      error stop 'Resume checkpoint matching-plane key has an invalid integer token: '//trim(key)
    end if
    read (value, *, iostat=ios) result_value
    if (ios /= 0) error stop 'Resume checkpoint matching-plane integer value could not be parsed: '//trim(key)
  end subroutine read_matching_summary_integer

  subroutine read_matching_summary_logical(value, key, result_value)
    character(len=*), intent(in) :: value, key
    logical, intent(out) :: result_value
    integer :: ios

    if (.not. is_logical_token(trim(value))) then
      error stop 'Resume checkpoint matching-plane key has an invalid logical token: '//trim(key)
    end if
    read (value, *, iostat=ios) result_value
    if (ios /= 0) error stop 'Resume checkpoint matching-plane logical value could not be parsed: '//trim(key)
  end subroutine read_matching_summary_logical

end module bem_matching_plane_summary
