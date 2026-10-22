!> 倍精度実数と32bit整数のkind定義を集約する基盤モジュール。
module bem_kinds
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  implicit none
  integer, parameter :: dp = real64
  !> Single precision, used only by the optional near-field panel tables.
  integer, parameter :: sp = real32
  integer, parameter :: i32 = int32
  integer, parameter :: i64 = int64
end module bem_kinds
