!> Fortranテスト用の簡易アサートと一時ファイル補助。
module test_support
  use iso_fortran_env, only: error_unit
  use bem_kinds, only: dp, i32
  implicit none
  private

  public :: fail_test
  public :: assert_true
  public :: assert_equal_i32
  public :: assert_close_dp
  public :: assert_allclose_1d
  public :: delete_file_if_exists
  public :: ensure_directory
  public :: remove_empty_directory

contains

  !> 失敗理由を表示してテストを停止する。
  subroutine fail_test(message)
    character(len=*), intent(in) :: message

    write (error_unit, '(a)') 'TEST FAILURE: ' // trim(message)
    error stop 1
  end subroutine fail_test

  !> 条件が偽ならテストを失敗させる。
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail_test(message)
  end subroutine assert_true

  !> 32bit整数の一致を確認する。
  subroutine assert_equal_i32(actual, expected, message)
    integer(i32), intent(in) :: actual
    integer(i32), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (actual /= expected) call fail_test(message)
  end subroutine assert_equal_i32

  !> 倍精度実数の近接一致を確認する。
  subroutine assert_close_dp(actual, expected, tol, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tol) call fail_test(message)
  end subroutine assert_close_dp

  !> 1次元実数配列の要素ごとの近接一致を確認する。
  subroutine assert_allclose_1d(actual, expected, tol, message)
    real(dp), intent(in) :: actual(:)
    real(dp), intent(in) :: expected(:)
    real(dp), intent(in) :: tol
    character(len=*), intent(in) :: message
    integer :: i

    if (size(actual) /= size(expected)) call fail_test(message)
    do i = 1, size(actual)
      if (abs(actual(i) - expected(i)) > tol) call fail_test(message)
    end do
  end subroutine assert_allclose_1d

  !> 既存ファイルがあれば削除する。
  subroutine delete_file_if_exists(path)
    character(len=*), intent(in) :: path
    logical :: exists
    integer :: u, ios

    inquire (file=trim(path), exist=exists)
    if (.not. exists) return

    open (newunit=u, file=trim(path), status='old', action='readwrite', iostat=ios)
    if (ios /= 0) return
    close (u, status='delete')
  end subroutine delete_file_if_exists

  !> ディレクトリを作成する。
  subroutine ensure_directory(path)
    character(len=*), intent(in) :: path
    character(len=1024) :: cmd
    integer :: ios

    cmd = 'mkdir -p "' // trim(path) // '"'
    call execute_command_line(trim(cmd), wait=.true., exitstat=ios)
    if (ios /= 0) call fail_test('failed to create directory: ' // trim(path))
  end subroutine ensure_directory

  !> 空ディレクトリを削除する。
  subroutine remove_empty_directory(path)
    character(len=*), intent(in) :: path
    character(len=1024) :: cmd
    integer :: ios

    cmd = 'if [ -d "' // trim(path) // '" ]; then rmdir "' // trim(path) // '"; fi'
    call execute_command_line(trim(cmd), wait=.true., exitstat=ios)
  end subroutine remove_empty_directory

end module test_support
